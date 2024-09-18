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

    | crop_region.run(
      runIf: { id, state -> state.crop_region_min_x },
      fromState: [
        "input": "output",
        "replicate_id": "replicate_id",
        "crop_region_min_x": "crop_region_min_x",
        "crop_region_min_y": "crop_region_min_y",
        "crop_region_max_x": "crop_region_max_x",
        "crop_region_max_y": "crop_region_max_y"
      ],
      toState: ["output"]
    )

    | setState(["output"])

  emit:
  output_ch
}

