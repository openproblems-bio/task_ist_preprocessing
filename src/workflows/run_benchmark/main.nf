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
            component_id: comp.name,
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
    | runEach(
      components: segm_methods,
      id: { id, state, comp ->
        id + "/segm_" + comp.name
      },
      fromState: ["input": "input_sp"],
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "segmentation",
            component_id: comp.name,
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
    pciseq
  ]
  segm_ass_ch = segm_ch
    | runEach(
      components: segm_ass_methods,
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
            component_id: comp.name,
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
    | runEach(
      components: count_aggr_methods,
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
            component_id: comp.name,
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
    | runEach(
      components: qc_filter_methods,
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
            component_id: comp.name,
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
    | runEach(
      components: cell_vol_methods,
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
            component_id: comp.name,
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
    | runEach(
      components: vol_norm_methods,
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
            component_id: comp.name,
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
    spanorm
  ]
  //direct_norm_ch = Channel.empty()
   direct_norm_ch = qc_filter_ch
     | runEach(
       components: direct_norm_methods,
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
            component_id: comp.name,
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
    moscot
  ]
  cta_ch = normalization_ch
    | runEach(
      components: cta_methods,
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
            component_id: comp.name,
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
    gene_efficiency_correction,
    resolvi_correction
  ]
  expr_corr_ch = cta_ch
    | runEach(
      components: expr_corr_methods,
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
            component_id: comp.name,
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
  metrics = [
    similarity
  ]
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
