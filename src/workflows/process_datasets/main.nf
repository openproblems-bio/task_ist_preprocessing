meta.resources_dir = 'target/executable/data_processors/process_dataset',
meta.config = 'target/executable/data_processors/process_dataset/.config.vsh.yaml'

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

    | process_dataset_comp.run(
      fromState: [
        input_sp: "input_sp",
        input_sc: "input_sc"
      ],
      toState: [
        output_sc: "output_sc",
        output_sp: "output_sp"
      ]
    )

    | setState(["output_sp", "output_sc"])

  emit:
  output_ch
}

