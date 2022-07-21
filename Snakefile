import yaml

cfg = yaml.load(open('configs/config.yaml', 'r'), Loader=yaml.FullLoader)
final_files = []
for run in cfg['PREPROCESSING']:
    final_files.extend(
            expand(
                cfg['ROOT']+'/{results}/{run}/{dataset}/metrics_{seg}_{assign}_{norm}.txt',
                results=cfg['RESULTS'],
                run=run,
                dataset=cfg['PREPROCESSING'][run]['dataset'],
                seg=cfg['PREPROCESSING'][run]['segmentation'],
                assign=cfg['PREPROCESSING'][run]['assignment'],
                norm=cfg['PREPROCESSING'][run]['normalize']
            )
        )

#TODO Add file name option for different method params

rule all:
    input:
        final_files

rule label:
    output:
        '{results}/{run}/{dataset}/label_{seg}.tif'
    params:
        img = lambda w: cfg['ROOT'] + '/' + cfg['DATA_SCENARIOS'][w.dataset]['image'],
        exp = lambda w: cfg['PREPROCESSING'][w.run]['seg_params']['expand'],
        bry = lambda w: "-b " if cfg['PREPROCESSING'][w.run]['seg_params']['binary'] else ""
    shell:
        "python load_images.py "
        "-i {params.img} "
        "-o {wildcards.results}/{wildcards.run}/{wildcards.dataset} "
        "-e {params.exp} "
        "-s {wildcards.seg} "
        "{params.bry}"

rule assign:
    input:
        '{results}/{run}/{dataset}/label_{seg}.tif'
    params:
        mol = lambda w: cfg['ROOT'] + '/' + cfg['DATA_SCENARIOS'][w.dataset]['molecules'],
        scd = lambda w: cfg['ROOT'] + '/' + cfg['DATA_SCENARIOS'][w.dataset]['sc_data'],
        opt = lambda w: cfg['PREPROCESSING'][w.run]['assign_params']['opts'],
    output:
        '{results}/{run}/{dataset}/assignments_{seg}_{assign}.csv'
    shell:
        "python run_{wildcards.assign}.py "
        "-m {params.mol} "
        "-o \"{params.opt}\" "
        "-sc {params.scd} "
        "-d {wildcards.results}/{wildcards.run}/{wildcards.dataset} "
        "-s {wildcards.seg} "

#TODO add params for area method
rule counts:
    input:
        '{results}/{run}/{dataset}/assignments_{seg}_{assign}.csv'
    output:
        '{results}/{run}/{dataset}/counts_{seg}_{assign}_{norm}.h5ad'
    shell:
        "python gen_counts.py "
        "-a {wildcards.assign} "
        "-d {wildcards.results}/{wildcards.run}/{wildcards.dataset} "
        "-s {wildcards.seg} "
        "-n {wildcards.norm} "

rule metric:
    input:
        '{results}/{run}/{dataset}/counts_{seg}_{assign}_{norm}.h5ad'
    params:
        scd = lambda w: cfg['ROOT'] + '/' + cfg['DATA_SCENARIOS'][w.dataset]['sc_data']
    output:
        '{results}/{run}/{dataset}/metrics_{seg}_{assign}_{norm}.txt'
    shell:
        "python calc_metrics.py "
        "-a {wildcards.assign} "
        "-d {wildcards.results}/{wildcards.run}/{wildcards.dataset} "
        "-s {wildcards.seg} "
        "-n {wildcards.norm} "
        "-sc {params.scd} "
