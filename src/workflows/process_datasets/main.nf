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

    | process_dataset_comp.run(
      fromState: [
        "input_sc",
        "input_sp",
        "dataset_id",
        "dataset_name",
        "dataset_url",
        "dataset_reference",
        "dataset_summary",
        "dataset_description",
        "dataset_organism"
      ],
      toState: [
        "output_sc",
        "output_sp"
      ]
    )

    | setState(["output_sp", "output_sc"])

  emit:
  output_ch
}

