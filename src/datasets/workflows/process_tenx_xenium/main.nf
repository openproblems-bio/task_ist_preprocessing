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

    | tenx_xenium.run(
      fromState: [
        "input",
        "replicate_id",
        "segmentation_id",
        "dataset_id",
        "dataset_name",
        "dataset_url",
        "dataset_reference",
        "dataset_summary",
        "dataset_description",
        "dataset_organism",
      ],
      toState: ["output"]
    )

    | setState(["output"])

  emit:
  output_ch
}

