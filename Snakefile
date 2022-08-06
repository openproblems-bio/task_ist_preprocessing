from TxsimConfig import *

configfile: 'configs/config.yaml'
parsed = ParsedConfig(config)
final_files = parsed.gen_file_names()

#Ensures dataset is a name not a file path, and id_code is an int
wildcard_constraints:
    dataset="[^\/]+",
    id_code="\d+"

def get_params(
    method: str, 
    id_code: int,
    param_name: str
) -> dict:
    """Get parameters for a given method

    Parameters
    ----------
    method : str
        The name of the method
    id_code : int
        The specific ID code
    param_name : str
        The parameter to get. Leave as 'p' for method-specific hyperparameters    

    Returns
    -------
    dict
        Dictionary of parameters
    """
    if parsed.get_method_params(method, id_code) is not None:
        return parsed.get_method_params(method, id_code).get(param_name)
    return None

#Main rule (will always be run)
rule all:
    input:
        final_files

#Rules corresponding to each method
rule custom_segmentation:
    conda:
        "envs/base-env.yaml"
    input:
        img = lambda w: parsed.get_data_file(w.dataset, 'segmented_image')
    output:
        '{results}/{dataset}/segments_custom-{id_code}.tif',
        '{results}/{dataset}/areas_custom-{id_code}.csv'
    params:
        hyp = lambda w: get_params('custom', int(w.id_code), 'p'),
        exp = lambda w: get_params('custom', int(w.id_code), 'expand')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyp}\" "
        "-o {wildcards.results}/{wildcards.dataset} "
        "-s custom "
        "-id {wildcards.id_code} "
        "-e {params.exp} "
        "-b "

rule watershed:
    conda:
        "envs/base-env.yaml"
    input:
        img = lambda w: parsed.get_data_file(w.dataset, 'image')
    output:
        '{results}/{dataset}/segments_watershed-{id_code}.tif',
        '{results}/{dataset}/areas_watershed-{id_code}.csv'
    params:
        hyp = lambda w: get_params('watershed', int(w.id_code), 'p'),
        exp = lambda w: get_params('watershed', int(w.id_code), 'expand')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyp}\" "
        "-o {wildcards.results}/{wildcards.dataset} "
        "-s watershed "
        "-id {wildcards.id_code} "
        "-e {params.exp} "

rule cellpose:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt
    conda: 
        "envs/cellpose-env.yaml"
    input:
        img = lambda w: parsed.get_data_file(w.dataset, 'image')
    output:
        '{results}/{dataset}/segments_cellpose-{id_code}.tif',
        '{results}/{dataset}/areas_cellpose-{id_code}.csv'
    params:
        hyp = lambda w: get_params('cellpose', int(w.id_code), 'p'),
        exp = lambda w: get_params('cellpose', int(w.id_code), 'expand')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyp}\" "
        "-o {wildcards.results}/{wildcards.dataset} "
        "-s cellpose "
        "-id {wildcards.id_code} "
        "-e {params.exp} "

rule pciSeq:
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt
    conda:
        "envs/pciSeq-env.yaml"
    input:
        '{results}/{dataset}/segments_{seg}.tif',
        mol = lambda w: parsed.get_data_file(w.dataset, 'molecules'),
        scd = lambda w: parsed.get_data_file(w.dataset, 'sc_data')
    params:
        hyp = lambda w: get_params('pciSeq', int(w.id_code), 'p')
    output:
        '{results}/{dataset}/assignments_{seg}_pciSeq-{id_code}.csv',
        '{results}/{dataset}/celltypes_{seg}_pciSeq-{id_code}.csv'
    shell:
        "python3 scripts/run_pciseq.py "
        "-m {input.mol} "
        "-p \"{params.hyp}\" "
        "-sc {input.scd} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-s {wildcards.seg} "
        "-id {wildcards.id_code} "

rule basic_assign:
    conda:
        "envs/base-env.yaml"
    input:
        '{results}/{dataset}/segments_{seg}.tif',
        mol = lambda w: parsed.get_data_file(w.dataset, 'molecules')
    params:
        hyp = lambda w: get_params('basic', int(w.id_code), 'p')
    output:
        '{results}/{dataset}/assignments_{seg}_basic-{id_code}.csv'
    shell:
        "python3 scripts/basic_assignment.py "
        "-m {input.mol} "
        "-p \"{params.hyp}\" "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-s {wildcards.seg} "
        "-id {wildcards.id_code} "

#TODO change parameters since it is running inside a container
rule baysor_assign:
    input: 
        '{results}/{dataset}/segments_{seg}.tif',
        mol = lambda w: parsed.get_data_file(w.dataset, 'molecules')
    params:
        hyp = lambda w: get_params('baysor', int(w.id_code), 'p')
    output:
        touch('{results}/{dataset}/blank_{seg}_baysor-{id_code}.txt')
    container:
        "docker://vpetukhov/baysor:master"
    shell:
        "baysor"

rule baysor_segment:
    container:
        "docker://vpetukhov/baysor:master"
    input:
        mol = lambda w: parsed.get_data_file(w.dataset, 'molecules')
    params:
        hyp = lambda w: get_params('baysor', int(w.id_code), 'p')
    output:
        touch('{results}/{dataset}/blank_baysor-{id_code}.txt')
    shell:
        "baysor"

rule format_baysor:
    input:
        '{results}/{dataset}/blank_{prev}baysor-{id_code}.txt'
    output:
        '{results}/{dataset}/segments_{prev}baysor-{id_code}.tif',
        '{results}/{dataset}/assignments_{prev}baysor-{id_code}.csv',
        '{results}/{dataset}/areas_{prev}baysor-{id_code}.csv'
    wildcard_constraints:
        prev=".*"


rule total_norm:
    conda:
        "envs/base-env.yaml"
    input:
        '{results}/{dataset}/assignments_{assign}.csv'
    params:
        hyp = lambda w: get_params('area', int(w.id_code), 'p')
    output:
        '{results}/{dataset}/counts_{assign}_total-{id_code}.h5ad'
    shell:
        "python3 scripts/gen_counts.py "
        "-as {wildcards.assign} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-n total "
        "-id {wildcards.id_code} "
        "-p \"{params.hyp}\" "

rule alpha_area:
    threads: 4
    conda:
        "envs/area-env.yaml"
    input:
        assign = '{results}/{dataset}/assignments_{method}.csv'
    params:
        hyp = lambda w: get_params('area', int(w.id_code), 'p')
    output:
        '{results}/{dataset}/counts_{method}_area-{id_code}.h5ad'
    shell:
        "python3 scripts/gen_counts.py "
        "-as {wildcards.method} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-n area "
        "-id {wildcards.id_code} "
        "-p \"{params.hyp}\" "

rule metric:
    conda:
        "envs/base-env.yaml"
    input:
        '{results}/{dataset}/counts_{methods}.h5ad',
        scd = lambda w: parsed.get_data_file(w.dataset, 'sc_data')
    output:
        '{results}/{dataset}/metrics_{methods}.txt'
    shell:
        "python3 scripts/calc_metrics.py "
        "-m {wildcards.methods} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-sc {input.scd} "
