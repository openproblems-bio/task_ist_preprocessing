include { findArgumentSchema } from "${meta.resources_dir}/helper.nf"

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
  output_ch = input_ch

    // copy id to the state
    | map{ id, state ->
      def new_state = state + [dataset_id: id]
      [id, new_state]
    }

    | lu_human_liver_cancer_sc.run(
      fromState: [
        "dataset_id",
        "dataset_name",
        "dataset_url",
        "dataset_reference",
        "dataset_summary",
        "dataset_description",
        "dataset_organism",
      ],
      toState: [
        "output_raw": "output"
      ]
    )

    | log_cp.run(
      key: "log_cp10k",
      fromState: [
        "input": "output_raw"
      ],
      args: [
        "normalization_id": "log_cp10k",
        "n_cp": 10000
      ],
      toState: [
        "output_normalized": "output"
      ]
    )
    | hvg.run(
      fromState: ["input": "output_normalized"],
      toState: ["output_hvg": "output"]
    )

    | pca.run(
      fromState: ["input": "output_hvg"],
      toState: ["output_pca": "output" ]
    )

    | knn.run(
      fromState: ["input": "output_pca"],
      toState: ["output_knn": "output"]
    )
    // add synonym
    | map{ id, state ->
      [id, state + [output_dataset: state.output_knn]]
    }

    | extract_uns_metadata.run(
      fromState: { id, state ->
        def schema = findArgumentSchema(meta.config, "output_dataset")
        // workaround: convert GString to String
        schema = iterateMap(schema, { it instanceof GString ? it.toString() : it })
        def schemaYaml = tempFile("schema.yaml")
        writeYaml(schema, schemaYaml)
        [
          "input": state.output_dataset,
          "schema": schemaYaml
        ]
      },
      toState: ["output_meta": "output"]
    )

    | setState([
      "output_dataset": "output_dataset",
      "output_meta": "output_meta"
    ])

  emit:
  output_ch
}

