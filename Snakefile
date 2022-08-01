from TxsimConfig import *

parsed = ParsedConfig('configs/config.yaml')
final_files = parsed.gen_file_names()

#TODO Type hints
def get_hyperparams(method, id_code):
    if parsed.get_method_params(method, id_code) is not None:
        return parsed.get_method_params(method, id_code).get('p')
    return None

#Main rule (will always be run)
rule all:
    input:
        final_files

#Rules corresponding to each method
rule watershed:
    output:
        '{results}/{dataset}/segments_watershed-{shp}.tif',
        '{results}/{dataset}/areas_watershed-{shp}.csv'
    wildcard_constraints:
        shp="\d+"
    params:
        img = lambda w: parsed.get_data_file(w.dataset, 'image'),
        exp = lambda w: parsed.get_method_params('watershed', int(w.shp)).get('expand'),
        bry = lambda w: "-b " if parsed.get_method_params('watershed', int(w.shp))['binary'] else ""
    shell:
        "python scripts/segment_image.py "
        "-i {params.img} "
        "-o {wildcards.results}/{wildcards.dataset} "
        "-e {params.exp} "
        "-s watershed "
        "-id {wildcards.shp} "
        "{params.bry}"

rule cellpose:
    output:
        '{results}/{dataset}/segments_cellpose-{shp}.tif',
        '{results}/{dataset}/areas_cellpose-{shp}.csv'
    wildcard_constraints:
        shp="\d+"
    params:
        img = lambda w: parsed.get_data_file(w.dataset, 'image'),
        exp = lambda w: parsed.get_method_params('cellpose', int(w.shp)).get('expand'),
        bry = lambda w: "-b " if parsed.get_method_params('cellpose', int(w.shp))['binary'] else ""
    shell:
        "python scripts/segment_image.py "
        "-i {params.img} "
        "-o {wildcards.results}/{wildcards.dataset} "
        "-e {params.exp} "
        "-s cellpose "
        "-id {wildcards.shp} "
        "{params.bry}"

rule pciSeq:
    input:
        '{results}/{dataset}/segments_{seg}.tif'
    params:
        mol = lambda w: parsed.get_data_file(w.dataset, 'molecules'),
        scd = lambda w: parsed.get_data_file(w.dataset, 'sc_data'),
        hyp = lambda w: get_hyperparams('pciSeq', int(w.ahp))
    output:
        '{results}/{dataset}/assignments_{seg}_pciSeq-{ahp}.csv',
        '{results}/{dataset}/celltypes_{seg}_pciSeq-{ahp}.csv'
    wildcard_constraints:
        ahp="\d+"
    shell:
        "python scripts/run_pciseq.py "
        "-m {params.mol} "
        "-p \"{params.hyp}\" "
        "-sc {params.scd} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-s {wildcards.seg} "
        "-id {wildcards.ahp} "

rule basic_assign:
    input:
        '{results}/{dataset}/segments_{seg}.tif'
    params:
        mol = lambda w: parsed.get_data_file(w.dataset, 'molecules'),
        hyp = lambda w: get_hyperparams('pciSeq', int(w.ahp))
    output:
        '{results}/{dataset}/assignments_{seg}_basic-{ahp}.csv'
    wildcard_constraints:
        ahp="\d+"
    shell:
        "python scripts/basic_assignment.py "
        "-m {params.mol} "
        "-p \"{params.hyp}\" "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-s {wildcards.seg} "
        "-id {wildcards.ahp} "

rule baysor_assign:
    input: 
        '{results}/{dataset}/segments_{seg}.tif'
    params:
        mol = lambda w: parsed.get_data_file(w.dataset, 'molecules'),
        hyp = lambda w: get_hyperparams('baysor', int(w.ahp))
    output:
        touch('{results}/{dataset}/blank_{seg}_baysor-{ahp}.txt')
    wildcard_constraints:
        ahp="\d+"
    container:
        "docker://vpetukhov/baysor:master"
    shell:
        "baysor"

rule baysor_segment:
    container:
        "docker://vpetukhov/baysor:master"
    params:
        mol = lambda w: parsed.get_data_file(w.dataset, 'molecules'),
        hyp = lambda w: get_hyperparams('baysor', int(w.ahp))
    output:
        touch('{results}/{dataset}/blank_baysor-{ahp}.txt')
    wildcard_constraints:
        ahp="\d+"
    shell:
        "baysor"

rule format_baysor:
    input:
        '{results}/{dataset}/blank_{prev}baysor-{ahp}.txt'
    output:
        '{results}/{dataset}/segments_{prev}baysor-{ahp}.tif',
        '{results}/{dataset}/assignments_{prev}baysor-{ahp}.csv',
        '{results}/{dataset}/areas_{prev}baysor-{ahp}.csv'
    wildcard_constraints:
        prev=".*",
        ahp="\d"


rule total_norm:
    input:
        '{results}/{dataset}/assignments_{assign}.csv'
    params:
        hyp = lambda w: get_hyperparams('area', int(w.nhp))
    output:
        '{results}/{dataset}/counts_{assign}_total-{nhp}.h5ad'
    wildcard_constraints:
        nhp="\d+"
    shell:
        "python scripts/gen_counts.py "
        "-as {wildcards.assign} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-n total "
        "-id {wildcards.nhp} "
        "-p \"{params.hyp}\" "

rule alpha_area:
    input:
        assign = '{results}/{dataset}/assignments_{method}.csv'
    params:
        hyp = lambda w: get_hyperparams('area', int(w.nhp))
    output:
        '{results}/{dataset}/counts_{method}_area-{nhp}.h5ad'
    wildcard_constraints:
        nhp="\d+"
    shell:
        "python scripts/gen_counts.py "
        "-as {wildcards.method} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-n area "
        "-id {wildcards.nhp} "
        "-p \"{params.hyp}\" "

rule metric:
    input:
        '{results}/{dataset}/counts_{methods}.h5ad'
    params:
        scd = lambda w: parsed.get_data_file(w.dataset, 'sc_data'),
    output:
        '{results}/{dataset}/metrics_{methods}.txt'
    shell:
        "python scripts/calc_metrics.py "
        "-m {wildcards.methods} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-sc {params.scd} "
