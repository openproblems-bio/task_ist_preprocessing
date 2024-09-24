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

    // example of channel event: 
    //   ["mouse_brain_combined", [
    //      "input_sc": file("resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad"),
    //      "input_sp": file("resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr")]]

    | process_dataset_comp.run(
      fromState: [
        input_ist: "input_sp",
        input_scrnaseq: "input_sc"
      ],
      toState: [
        output_sc: "output_scrnaseq",
        output_sp: "output_ist"
      ]
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

