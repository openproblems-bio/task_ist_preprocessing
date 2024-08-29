include { findArgumentSchema } from "${meta.resources_dir}/helper.nf"

workflow auto {
  findStatesTemp(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // example of channel event: 
    //   ["my_id", ["input_sc": file("..."), "input_sp": file("...")]]

    | process_datasets.run(
      fromState: ["input_sc", "input_sp"],
      toState: ["output_sc", "output_sp"]
    )

    // example of channel event at this point:
    //   ["my_id", ["input_sc": ..., "input_sp": ..., 
    //              "output_sc": file("..."), "output_sp": file("...")]]

    | setState(["output_sp", "output_sc"])

    // example of channel event at this point:
    //   ["my_id", ["output_sc": file("..."), "output_sp": file("...")]]

  emit:
  output_ch
}

