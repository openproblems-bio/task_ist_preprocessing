include { checkItemAllowed } from "${meta.resources_dir}/helper.nf"


workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:

  /****************************************
   *      INITIALIZE DATA STRUCTURES      *
   ****************************************/
  init_ch = input_ch
    | map { id, state ->
      // Load method_parameters from YAML file if provided in state
      def method_parameters = null
      if (state.containsKey("method_parameters_yaml") && state.method_parameters_yaml) {
        def yaml_blob = readYaml(state.method_parameters_yaml)
        if (yaml_blob instanceof Map && yaml_blob.containsKey('parameters')) {
          method_parameters = yaml_blob.parameters
        } else if (yaml_blob instanceof Map) {
          // assume the file itself is the mapping
          method_parameters = yaml_blob
        }
      }

      // Initialize steps
      def steps = [
        [type: "dataset", dataset_id: id]
      ]

      // Create new state
      def new_state = state + [
        orig_id: id,
        steps: steps,
        method_parameters: method_parameters
      ]
      [id, new_state]
    }

    // extract the dataset metadata
    | extract_uns_metadata.run(
      fromState: [input: "input_sc"],
      toState: { id, output, state ->
        state + [
          dataset_uns: readYaml(output.output).uns
        ]
      }
    )

  /****************************************
   *        CONTROL METHODS               *
   ****************************************/
  control_methods = [
    identity,
    permute_celltype_annotations
  ]
  control_ch = init_ch
    | runEach(
      components: control_methods,
      id: { id, state, comp ->
        id + "/control_" + comp.name
      },
      fromState: [
        input_scrnaseq_reference: "input_sc"
      ],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "control",
            component_id: comp.config.name,
            component_variant: comp.name,
            run_id: id
          ]],
          output_correction: out_dict.output,
          output_qc_filter: out_dict.output_qc_col,
          output_assignment: out_dict.output_transcript_assignments
        ]
      }
    )

  /****************************************
   *       RUN SEGMENTATION METHODS       *
   ****************************************/
  // base segmentation methods
  segm_methods = [
    custom_segmentation.run(
      args: [labels_key: "cell_labels"]
    ),
    cellpose,
    binning,
    stardist,
    watershed
  ]
  
  segm_ch = init_ch
    | expandChannelWithParameterSets(segm_methods, "segm", "segmentation_methods")
    | runEach(
      components: segm_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [input: state.input_sp] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
          steps: state.steps + [[
            type: "segmentation",
            component_id: state.current_method_id,
            component_variant: state.current_method_variant,
            run_id: id
          ]],
          output_segmentation: out_dict.output
        ]
      }
    )

  /****************************************
   *   RUN ASSIGNMENT AFTER SEGMENTATION  *
   ****************************************/
  segm_ass_methods = [
    basic_transcript_assignment.run(
      args: [
        transcripts_key: "transcripts",
        coordinate_system: "global"
      ]
    ),
    baysor,
    clustermap,
    pciseq,
    comseg,
    proseg
  ]
  
  segm_ass_ch = segm_ch
    | expandChannelWithParameterSets(segm_ass_methods, "ass", "transcript_assignment_methods")
    | runEach(
      components: segm_ass_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [
          input_ist: state.input_sp,
          input_scrnaseq: state.input_sc,
          input_segmentation: state.output_segmentation
        ] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
          steps: state.steps + [[
            type: "assignment",
            component_id: state.current_method_id,
            component_variant: state.current_method_variant,
            run_id: id
          ]],
          output_assignment: out_dict.output
        ]
      }
    )

  /****************************************
   *         RUN DIRECT ASSIGNMENT        *
   ****************************************/

  // TODO: implement this when direct assignment methods are added
  direct_ass_methods = []
  direct_ass_ch = Channel.empty()
  // direct_ass_ch = init_ch
  //   | runEach(
  //     components: segm_ass_methods,
  //     id: { id, state, comp ->
  //       id + "/ass_" + comp.name
  //     },
  //     fromState: [
  //       input_ist: "input_sp",
  //       input_scrnaseq: "input_sc"
  //     ],
  //     toState: { id, out_dict, state, comp ->
  //       state + [
  //         steps: state.steps + [[
  //           type: "assignment",
  //           run_id: id,
  //           component_id: comp.name,
  //           input_state: state,
  //           output_dict: out_dict
  //         ]],
  //         output_assignment: out_dict.output
  //       ]
  //     }
  //   )

  /****************************************
   *          COMBINE ASSIGNMENT          *
   ****************************************/
  assignment_ch = segm_ass_ch.mix(direct_ass_ch)

  /****************************************
   *          COUNT AGGREGATION           *
   ****************************************/
  count_aggr_methods = [
    basic_count_aggregation
  ]
  
  count_aggr_ch = assignment_ch
    | expandChannelWithParameterSets(count_aggr_methods, "aggr", "count_aggregation_methods")
    | runEach(
      components: count_aggr_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [input: state.output_assignment] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
          steps: state.steps + [[
            type: "count_aggregation",
            component_id: state.current_method_id,
            component_variant: state.current_method_variant,
            run_id: id
          ]],
          output_count_aggregation: out_dict.output
        ]
      }
    )


  /************************************
   *          QC FILTERING            *
   ************************************/
  qc_filter_methods = [
    basic_qc_filter
  ]
  
  qc_filter_ch = count_aggr_ch
    | expandChannelWithParameterSets(qc_filter_methods, "qc_filter", "qc_filtering_methods")
    | runEach(
      components: qc_filter_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [input: state.output_count_aggregation] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
          steps: state.steps + [[
            type: "qc_filter",
            component_id: state.current_method_id,
            component_variant: state.current_method_variant,
            run_id: id
          ]],
          output_qc_filter: out_dict.output
        ]
      }
    )


  /****************************************
   *          VOLUME CALCULATION          *
   ****************************************/
  cell_vol_methods = [
    alpha_shapes
  ]
  
  cell_vol_ch = qc_filter_ch
    | expandChannelWithParameterSets(cell_vol_methods, "cell_vol", "volume_calculation_methods")
    | runEach(
      components: cell_vol_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [input: state.output_assignment] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
          steps: state.steps + [[
            type: "calculate_cell_volume",
            component_id: state.current_method_id,
            component_variant: state.current_method_variant,
            run_id: id
          ]],
          output_cell_volume: out_dict.output
        ]
      }
    )

  /****************************************
   *        NORMALIZATION BY VOLUME       *
   ****************************************/
  vol_norm_methods = [
    normalize_by_volume
  ]
  
  vol_norm_ch = cell_vol_ch
    | expandChannelWithParameterSets(vol_norm_methods, "norm", "normalization_methods")
    | runEach(
      components: vol_norm_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [
          input_spatial_aggregated_counts: state.output_count_aggregation,
          input_cell_volumes: state.output_cell_volume
        ] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
          steps: state.steps + [[
            type: "normalization",
            component_id: state.current_method_id,
            component_variant: state.current_method_variant,
            run_id: id
          ]],
          output_normalization: out_dict.output
        ]
      }
    )


  /****************************************
   *         DIRECT NORMALIZATION         *
   ****************************************/

   // TODO: implement this when direct normalization methods are added
  direct_norm_methods = [
    normalize_by_counts,
    spanorm
  ]
  
  direct_norm_ch = qc_filter_ch
    | expandChannelWithParameterSets(direct_norm_methods, "norm", "normalization_methods")
    | runEach(
      components: direct_norm_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [input_spatial_aggregated_counts: state.output_count_aggregation] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
          steps: state.steps + [[
            type: "normalization",
            component_id: state.current_method_id,
            component_variant: state.current_method_variant,
            run_id: id
          ]],
          output_normalization: out_dict.output
        ]
      }
    )

  /****************************************
   *          COMBINE NORMALIZATION        *
   ****************************************/
  normalization_ch = vol_norm_ch.mix(direct_norm_ch)


  /****************************************
   *         CELL TYPE ANNOTATION         *
   ****************************************/
  cta_methods = [
    ssam,
    tacco,
    moscot,
    mapmycells,
    tangram
  ]
  
  cta_ch = normalization_ch
    | expandChannelWithParameterSets(cta_methods, "cta", "celltype_annotation_methods")
    | runEach(
      components: cta_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [
          input_spatial_normalized_counts: state.output_normalization,
          input_transcript_assignments: state.output_assignment,
          input_scrnaseq_reference: state.input_sc
        ] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
          steps: state.steps + [[
            type: "cell_type_assignment",
            component_id: state.current_method_id,
            component_variant: state.current_method_variant,
            run_id: id
          ]],
          output_cta: out_dict.output
        ]
      }
    )

  /****************************************
   *         EXPRESSION CORRECTION        *
   ****************************************/
  expr_corr_methods = [
    no_correction,
    gene_efficiency_correction,
    resolvi_correction
  ]
  
  expr_corr_ch = cta_ch
    | expandChannelWithParameterSets(expr_corr_methods, "corr", "expression_correction_methods")
    | runEach(
      components: expr_corr_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [
          input_spatial_with_cell_types: state.output_cta,
          input_scrnaseq_reference: state.input_sc
        ] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
          steps: state.steps + [[
            type: "expression_correction",
            component_id: state.current_method_id,
            component_variant: state.current_method_variant,
            run_id: id
          ]],
          output_correction: out_dict.output
        ]
      }
    )
    // aggregate spatial data for quality metrics
    | aggregate_spatial_data.run(
      fromState: [
        input_raw_sp: "input_sp",
        input_transcript_assignments: "output_assignment",
        input_qc_col: "output_qc_filter",
        input_spatial_corrected_counts: "output_correction"
      ],
      toState: [output_agg_spatial_data: "output"],
      auto: params.save_spatial_data ? [publish: true] : [:],
      directives: params.save_spatial_data ? [
        publishDir: [
          path: "${params.publish_dir}/spatial_data/",
          mode: "copy"
        ]
      ] : [:]
    )



  /****************************************
   *          COMBINE WITH CONTROL        *
   ****************************************/

  expr_corr_and_control_ch = expr_corr_ch.mix(control_ch)


  /****************************************
   *                METRICS               *
   ****************************************/
  metrics = [
    similarity,
    quality
  ]

  metric_ch = expr_corr_and_control_ch
    | runEach(
      components: metrics,
      filter: { id, state, comp ->
        // only run quality metrics if they have been computed
        if (comp.config.name == "quality" && !state.containsKey("output_agg_spatial_data")) {
          return false
        }
        return true
      },
      id: { id, state, comp ->
        id + "/metric_" + comp.name
      },
      fromState: { id, state, comp ->
        if (comp.config.name == "quality") {
          return [
            input: state.output_agg_spatial_data,
            input_qc_col: state.output_qc_filter
          ]
        } else {
          return [
            input: state.output_correction,
            input_qc_col: state.output_qc_filter,
            input_sc: state.input_sc,
            input_transcript_assignments: state.output_assignment
          ]
        }
      },
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "metric",
            component_id: comp.name,
            run_id: id
          ]],
          output_metric: out_dict.output
        ]
      }
    )

    // extract the scores
    | extract_uns_metadata.run(
      key: "extract_uns_scores",
      fromState: [input: "output_metric"],
      args: [
        uns_length_cutoff: 100
      ],
      toState: { id, output, state ->
        state + [
          score_uns: readYaml(output.output).uns
        ]
      }
    )

    | joinStates { ids, states ->
      // TODO: determine what to store in the score_uns file

      // store the scores in a file
      def score_uns = states.collect{state ->
        def method_ids =
          state.steps.collectMany{step ->
            step.type in ["dataset", "metric"] ? [] : [step.component_id]
          }
        return [
          dataset_id: state.orig_id,
          // dataset_sc_id: ..., // todo: extract this from the dataset
          // dataset_sp_id: ..., // todo: extract this from the dataset
          method_ids: method_ids,
          steps: state.steps,
          metric_ids: state.score_uns['metric_ids'],
          metrics_values: state.score_uns['metric_values']
        ]
      }
      def score_uns_yaml_blob = toYamlBlob(score_uns)
      def score_uns_file = tempFile("score_uns.yaml")
      score_uns_file.write(score_uns_yaml_blob)

      ["output", [output_scores: score_uns_file]]
    }


  /******************************
   * GENERATE OUTPUT YAML FILES *
   ******************************/
  // extract the dataset metadata
  meta_ch = init_ch

    // store join id
    | map{ id, state ->
      [id, state + ["_meta": [join_id: id]]]
    }

    // only keep one of the normalization methods
    | joinStates { ids, states ->
      // TODO: determine what to store in the dataset_uns file

      // store the dataset metadata in a file
      def dataset_uns = states.collect{it.dataset_uns}
      def dataset_uns_yaml_blob = toYamlBlob(dataset_uns)
      def dataset_uns_file = tempFile("dataset_uns.yaml")
      dataset_uns_file.write(dataset_uns_yaml_blob)

      // store the method configs in a file
      def methods =
        segm_methods + segm_ass_methods + direct_ass_methods + count_aggr_methods +
        qc_filter_methods + cell_vol_methods + vol_norm_methods + direct_norm_methods +
        cta_methods + expr_corr_methods
      def method_configs = methods.collect{it.config}
      def method_configs_yaml_blob = toYamlBlob(method_configs)
      def method_configs_file = tempFile("method_configs.yaml")
      method_configs_file.write(method_configs_yaml_blob)

      // store the metric configs in a file
      def metric_configs = metrics.collect{it.config}
      def metric_configs_yaml_blob = toYamlBlob(metric_configs)
      def metric_configs_file = tempFile("metric_configs.yaml")
      metric_configs_file.write(metric_configs_yaml_blob)

      // retrieve task info
      def viash_file = meta.resources_dir.resolve("_viash.yaml")

      // create state
      def new_state = [
        output_dataset_info: dataset_uns_file,
        output_method_configs: method_configs_file,
        output_metric_configs: metric_configs_file,
        output_task_info: viash_file,
        _meta: states[0]._meta
      ]

      ["output", new_state]
    }

  output_ch = metric_ch
    | mix(meta_ch)
    | joinStates{ ids, states ->
      def mergedStates = states.inject([:]) { acc, m -> acc + m }
      [ids[0], mergedStates]
    }

  emit:
  output_ch
}

/**
 * Remove specified keys from a map.
 * 
 * @param map The map to remove keys from
 * @param keys List of keys to remove
 * @return New map with specified keys removed
 */
def removeKeys(map, keys) {
  def result = map.clone()
  keys.each { key -> result.remove(key) }
  return result
}

/**
 * Expand a channel with parameter sets for each method.
 * 
 * This function returns a workflow that expands a channel by creating multiple variants
 * for each method that has parameter sets defined in state.method_parameters.
 * Each variant gets a unique ID and stores the method configuration in the state.
 * 
 * @param methods List of component objects to potentially expand
 * @param method_type String identifying the method type (e.g., "segmentation")
 * @param selected_methods_key Key in state containing the list of selected methods
 * @return Workflow that takes a channel and emits expanded [id, state] tuples
 */
def expandChannelWithParameterSets(methods, method_type, selected_methods_key) {
  workflow expand_ch_wf {
    take: input_ch
    
    main:
    output_ch = input_ch
      | flatMap { id, state ->
        def selected_methods = state[selected_methods_key]
        def default_methods = state.default_methods ?: []
        def results = []
        
        methods.each { comp ->
          def comp_name = comp.config.name
          
          // Check if method is selected
          if (!checkItemAllowed(comp_name, selected_methods, null, selected_methods_key, "NA")) {
            return
          }
          
          // Check constraint: at most one non-default method or parameter set variant
          // Count non-default items in previous steps (both non-default methods and variants)
          def has_non_default = state.steps.any { step ->
            step.component_id != null && 
            (!(step.component_id in default_methods) || 
             (step.component_variant != null && step.component_id != step.component_variant))
          }
          
          // Get parameter sets for this method from state
          def method_parameters = state.method_parameters
          def spec = method_parameters ? method_parameters[comp_name] : null
          
          if (!spec) {
            // No parameter sets defined, check if this is a non-default method
            def is_non_default_method = !(comp_name in default_methods)
            
            // Skip if we already have a non-default item and this would be another
            if (has_non_default && is_non_default_method) {
              return
            }
            
            // Use default
            def new_id = "${id}/${method_type}_${comp_name}"
            def new_state = state + [
              current_method_id: comp_name,
              current_method_variant: comp_name,
              current_method_args: [:]
            ]
            results << [new_id, new_state]
          } else {
            def default_args = spec['default'] ?: [:]
            def is_non_default_method = !(comp_name in default_methods)
            
            // Skip if we already have a non-default item and this is a non-default method
            if (has_non_default && is_non_default_method) {
              return
            }
            
            // Add default variant
            def new_id = "${id}/${method_type}_${comp_name}"
            def new_state = state + [
              current_method_id: comp_name,
              current_method_variant: comp_name,
              current_method_args: default_args
            ]
            results << [new_id, new_state]
            
            // Add parameter sweep variants only if no non-default item exists yet
            def sweep = spec['sweep']
            if (sweep instanceof Map && !has_non_default) {
              sweep.each { parName, parVals ->
                if (!(parVals instanceof Collection)) parVals = [parVals]
                parVals.each { val ->
                  def args = default_args + [(parName): val]
                  def variant_name = "${comp_name}_${parName}_${val}"
                  def variant_id = "${id}/${method_type}_${variant_name}"
                  def variant_state = state + [
                    current_method_id: comp_name,
                    current_method_variant: variant_name,
                    current_method_args: args
                  ]
                  results << [variant_id, variant_state]
                }
              }
            }
          }
        }
        
        return results
      }
    
    emit:
    output_ch
  }
  return expand_ch_wf.cloneWithName("expand_channel_for_" + method_type)
}
