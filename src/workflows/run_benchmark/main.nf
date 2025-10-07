include { checkItemAllowed } from "${meta.resources_dir}/helper.nf"

/**
 * Check if a component should be run based on method selection and constraint rules.
 * 
 * This function enforces three constraints for the benchmark workflow:
 * 1. Method selection: The component must be in the selected_methods list (if provided)
 * 2. Parameter set constraint: At most one non-default parameter set can be used across the entire pipeline
 * 3. Default methods constraint: At most one non-default method can be used in the pipeline
 * 
 * A non-default parameter set is detected when comp.config.name != comp.getName(), which
 * occurs when expandMethodsWithParameterSets creates variants with modified parameters.
 * These variants are named "${comp_name}_${parName}_${val}" while the default variant
 * keeps the original component name.
 * 
 * Default methods are those specified in the default_methods list. This constraint allows
 * running all default methods together, or one non-default method with all defaults.
 * 
 * @param comp The component object to check, must have config.name and getName() methods
 * @param selected_methods List of method names that are allowed to run (null means all allowed)
 * @param steps List of step objects representing methods already in the pipeline, each with
 *              component_id (original name) and component_variant (name after expansion)
 * @param default_methods List of method names that are considered "default" methods
 * @return true if the component should run, false if it violates any constraints
 */
Boolean checkMethodConstraints(comp, List selected_methods, List steps, List default_methods) {
  // Constraint 1: Check if this method is in the selected methods list
  if (!checkItemAllowed(comp.config.name, selected_methods, null, "selected_methods", "NA")) {
    return false
  }

  // Constraint 2: Check if this component uses a non-default parameter set
  def is_non_default_parameterset = comp.config.name != comp.getName()
  
  if (is_non_default_parameterset) {
    // Count how many non-default parameter sets are already in the steps
    // A step uses a non-default parameter set if component_id != component_variant
    def non_default_parameterset_count = steps.findAll { step ->
      step.component_id != null && step.component_variant != null && 
      step.component_id != step.component_variant
    }.size()
    
    // Only allow this method if there are no other non-default parameter sets
    if (non_default_parameterset_count > 0) {
      return false
    }
  }
  
  // Constraint 3: Check the default methods constraint (at most one non-default method)
  if (default_methods && default_methods.size() > 0 && steps && steps.size() > 0) {
    // Extract component_ids from steps and count non-default methods
    def component_ids = steps.collect{it.component_id}.findAll{it != null}
    def non_default_count = component_ids.findAll { !(it in default_methods) }.size()
    
    if (non_default_count > 1) {
      return false
    }
  }
  
  return true
}

// Expand a list of components into variants according to params.method_parameters
List expandMethodsWithParameterSets(List methods, Map parameters_map) {
  if (!parameters_map) {
    return methods
  }

  def expanded = []

  methods.each { comp ->
    def comp_name = comp.config.name
    def spec = parameters_map[comp_name]

    if (!spec) {
      expanded << comp
    } else {
      def default_args = spec['default'] ?: [:]

      // include default variant
      expanded << comp.run(key: comp_name, args: default_args)

      def sweep = spec['sweep']
      if (sweep instanceof Map) {
        // vary one parameter at a time
        sweep.each { parName, parVals ->
          if (!(parVals instanceof Collection)) parVals = [parVals]
          parVals.each { val ->
            def args = default_args + [(parName): val]
            def k = "${comp_name}_${parName}_${val}"
            expanded << comp.run(key: k, args: args)
          }
        }
      }
    }
  }

  return expanded
}

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
      def new_state = state + [
        orig_id: id,
        steps: [
          [type: "dataset", dataset_id: id]
        ]
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

  // If a method_parameters_yaml is provided, load it and set params.method_parameters
  if (params.containsKey("method_parameters_yaml") && params.method_parameters_yaml) {
    def yaml_blob = readYaml(file(params.method_parameters_yaml))
    if (yaml_blob instanceof Map && yaml_blob.containsKey('parameters')) {
      params.method_parameters = yaml_blob.parameters
    } else if (yaml_blob instanceof Map) {
      // assume the file itself is the mapping
      params.method_parameters = yaml_blob
    }
  }

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
  _segm_methods_base = [
    custom_segmentation.run(
      args: [labels_key: "cell_labels"]
    ),
    cellpose,
    binning,
    stardist,
    watershed
  ]
  segm_methods = expandMethodsWithParameterSets(_segm_methods_base, params.method_parameters)
  segm_ch = init_ch
    | runEach(
      components: segm_methods,
      filter: { id, state, comp ->
        checkMethodConstraints(
          comp,
          state.segmentation_methods,
          state.steps,
          state.default_methods
        )
      },
      id: { id, state, comp ->
        id + "/segm_" + comp.name
      },
      fromState: ["input": "input_sp"],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "segmentation",
            component_id: comp.config.name,
            component_variant: comp.name,
            run_id: id
          ]],
          output_segmentation: out_dict.output
        ]
      }
    )

  /****************************************
   *   RUN ASSIGNMENT AFTER SEGMENTATION  *
   ****************************************/
  _segm_ass_methods_base = [
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
  segm_ass_methods = expandMethodsWithParameterSets(_segm_ass_methods_base, params.method_parameters)
  segm_ass_ch = segm_ch
    | runEach(
      components: segm_ass_methods,
      filter: { id, state, comp ->
        checkMethodConstraints(
          comp,
          state.transcript_assignment_methods,
          state.steps,
          state.default_methods
        )
      },
      id: { id, state, comp ->
        id + "/ass_" + comp.name
      },
      fromState: [
        input_ist: "input_sp",
        input_scrnaseq: "input_sc",
        input_segmentation: "output_segmentation"
      ],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "assignment",
            component_id: comp.config.name,
            component_variant: comp.name,
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
  _count_aggr_methods_base = [
    basic_count_aggregation
  ]
  count_aggr_methods = expandMethodsWithParameterSets(_count_aggr_methods_base, params.method_parameters)
  count_aggr_ch = assignment_ch
    | runEach(
      components: count_aggr_methods,
      filter: { id, state, comp ->
        checkMethodConstraints(
          comp,
          state.count_aggregation_methods,
          state.steps,
          state.default_methods
        )
      },
      id: { id, state, comp ->
        id + "/aggr_" + comp.name
      },
      fromState: [
        input: "output_assignment"
      ],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "count_aggregation",
            component_id: comp.config.name,
            component_variant: comp.name,
            run_id: id
          ]],
          output_count_aggregation: out_dict.output
        ]
      }
    )


  /************************************
   *          QC FILTERING            *
   ************************************/
  _qc_filter_methods_base = [
    basic_qc_filter
  ]
  qc_filter_methods = expandMethodsWithParameterSets(_qc_filter_methods_base, params.method_parameters)
  qc_filter_ch = count_aggr_ch
    | runEach(
      components: qc_filter_methods,
      filter: { id, state, comp ->
        checkMethodConstraints(
          comp,
          state.qc_filtering_methods,
          state.steps,
          state.default_methods
        )
      },
      id: { id, state, comp ->
        id + "/qc_filter_" + comp.name
      },
      fromState: [
        input: "output_count_aggregation"
      ],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "qc_filter",
            component_id: comp.config.name,
            component_variant: comp.name,
            run_id: id
          ]],
          output_qc_filter: out_dict.output
        ]
      }
    )


  /****************************************
   *          VOLUME CALCULATION          *
   ****************************************/
  _cell_vol_methods_base = [
    alpha_shapes
  ]
  cell_vol_methods = expandMethodsWithParameterSets(_cell_vol_methods_base, params.method_parameters)
  cell_vol_ch = qc_filter_ch
    | runEach(
      components: cell_vol_methods,
      filter: { id, state, comp ->
        checkMethodConstraints(
          comp,
          state.volume_calculation_methods,
          state.steps,
          state.default_methods
        )
      },
      id: { id, state, comp ->
        id + "/cell_vol_" + comp.name
      },
      fromState: [
        input: "output_assignment"
      ],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "calculate_cell_volume",
            component_id: comp.config.name,
            component_variant: comp.name,
            run_id: id
          ]],
          output_cell_volume: out_dict.output
        ]
      }
    )

  /****************************************
   *        NORMALIZATION BY VOLUME       *
   ****************************************/
  _vol_norm_methods_base = [
    normalize_by_volume
  ]
  vol_norm_methods = expandMethodsWithParameterSets(_vol_norm_methods_base, params.method_parameters)
  vol_norm_ch = cell_vol_ch
    | runEach(
      components: vol_norm_methods,
      filter: { id, state, comp ->
        checkMethodConstraints(
          comp,
          state.normalization_methods,
          state.steps,
          state.default_methods
        )
      },
      id: { id, state, comp ->
        id + "/norm_" + comp.name
      },
      fromState: [
        input_spatial_aggregated_counts: "output_count_aggregation",
        input_cell_volumes: "output_cell_volume"
      ],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "normalization",
            component_id: comp.config.name,
            component_variant: comp.name,
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
  _direct_norm_methods_base = [
    normalize_by_counts,
    spanorm
  ]
  direct_norm_methods = expandMethodsWithParameterSets(_direct_norm_methods_base, params.method_parameters)
  direct_norm_ch = qc_filter_ch
    | runEach(
      components: direct_norm_methods,
      filter: { id, state, comp ->
        checkMethodConstraints(
          comp,
          state.normalization_methods,
          state.steps,
          state.default_methods
        )
      },
      id: { id, state, comp ->
        id + "/norm_" + comp.name
      },
      fromState: [
        input_spatial_aggregated_counts: "output_count_aggregation",
      ],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "normalization",
            component_id: comp.config.name,
            component_variant: comp.name,
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
  _cta_methods_base = [
    ssam,
    tacco,
    moscot
  ]
  cta_methods = expandMethodsWithParameterSets(_cta_methods_base, params.method_parameters)
  cta_ch = normalization_ch
    | runEach(
      components: cta_methods,
      filter: { id, state, comp ->
        checkMethodConstraints(
          comp,
          state.celltype_annotation_methods,
          state.steps,
          state.default_methods
        )
      },
      id: { id, state, comp ->
        id + "/cta_" + comp.name
      },
      fromState: [
        input_spatial_normalized_counts: "output_normalization",
        input_transcript_assignments: "output_assignment",
        input_scrnaseq_reference: "input_sc"
      ],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "cell_type_assignment",
            component_id: comp.config.name,
            component_variant: comp.name,
            run_id: id
          ]],
          output_cta: out_dict.output
        ]
      }
    )

  /****************************************
   *         EXPRESSION CORRECTION        *
   ****************************************/
  _expr_corr_methods_base = [
    no_correction,
    gene_efficiency_correction,
    resolvi_correction
  ]
  expr_corr_methods = expandMethodsWithParameterSets(_expr_corr_methods_base, params.method_parameters)
  expr_corr_ch = cta_ch
    | runEach(
      components: expr_corr_methods,
      filter: { id, state, comp ->
        checkMethodConstraints(
          comp,
          state.expression_correction_methods,
          state.steps,
          state.default_methods
        )
      },
      id: { id, state, comp ->
        id + "/corr_" + comp.name
      },
      fromState: [
        input_spatial_with_cell_types: "output_cta",
        input_scrnaseq_reference: "input_sc"
      ],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "expression_correction",
            component_id: comp.config.name,
            component_variant: comp.name,
            run_id: id
          ]],
          output_correction: out_dict.output
        ]
      }
    )

  /****************************************
   *          COMBINE WITH CONTROL        *
   ****************************************/

  expr_corr_and_control_ch = expr_corr_ch.mix(control_ch)


  /****************************************
   *                METRICS               *
   ****************************************/
  _metrics_base = [
    similarity
  ]
  metrics = expandMethodsWithParameterSets(_metrics_base, params.method_parameters)
  metric_ch = expr_corr_and_control_ch
    | runEach(
      components: metrics,
      id: { id, state, comp ->
        id + "/metric_" + comp.name
      },
      fromState: [
        input: "output_correction",
        input_qc_col: "output_qc_filter",
        input_sc: "input_sc",
        input_transcript_assignments: "output_assignment"
      ],
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
