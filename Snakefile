import pandas as pd
from pathlib import PurePath
configfile:"config.yaml"

out_dir = config['out_dir']

rule all:
    input:
        expand(out_dir + '/CPM{CPMcutoff}_mincell{mincell}/pval{pval}_qval{qval}_logfc{logfc}/GSEAresult/{sample}_{cluster}/{sample}_{cluster}_{direction}_clusterprofiler_result.csv',
        sample = config['sampleLIST'],
        direction = config['directionLIST'], 
        CPMcutoff = config['cpmLIST'], 
        mincell = config['mincellLIST'], 
        pval = [f'{x:.2f}' for x in config['pvalLIST']], 
        qval = [f'{x:.2f}' for x in config['qvalLIST']], 
        logfc = [f'{x:.2f}' for x in config['logfcLIST']], 
        cluster = config['clusterLIST'],
        )

rule diffxpy:
    params:
        CPMcutoff = '{CPMcutoff}',
        mincell = '{mincell}',
        script_dir = config['script_dir'],
        outpath = out_dir + '/CPM{CPMcutoff}_mincell{mincell}'
    input:
        h5ad = config['h5ad'],
        cpmFile = config['CPMfile']
    output:
        out_dir + '/CPM{CPMcutoff}_mincell{mincell}/Diffxpy.csv'
    shell:
         '''
         export PATH=/public/home/lizw/anaconda3/envs/diffxpy_v3/bin/:$PATH
         export PATH=/public/home/lizw/anaconda3/envs/sc_py_221214/bin/:$PATH
         python {params.script_dir}/diffxpy_song.py --h5ad {input.h5ad} --cpmfile {input.cpmFile} --output {output} --cpmcutoff {params.CPMcutoff} --mincell {params.mincell}
         '''
        
rule DEGfilter:
    params:
        logfc = '{logfc}', 
        pval = '{pval}',
        qval = '{qval}',
        script_dir = config['script_dir'],
        outpath = out_dir +'/CPM{CPMcutoff}_mincell{mincell}',
        cpm_b417 = config['cpm_b417_cluster'],
        cpm_gmi = config['cpm_gmi_cluster'],
        cpm_bf = config['cpm_bf_cluster'],
    input:
        out_dir + '/CPM{CPMcutoff}_mincell{mincell}/Diffxpy.csv',
    output:
        out_dir + '/CPM{CPMcutoff}_mincell{mincell}/pval{pval}_qval{qval}_logfc{logfc}/{sample}_{direction}_Diffxpy.csv'
    shell:
        '''
        export PATH=/public/home/lizw/anaconda3/envs/diffxpy_v3/bin/:$PATH
        export PATH=/public/home/lizw/anaconda3/envs/sc_py_221214/bin/:$PATH
        python {params.script_dir}/diffxpy_filter.py --qval {params.qval} --pval {params.pval} --logfc {params.logfc} --cpm_b417file {params.cpm_b417} --cpm_gmifile {params.cpm_gmi} --cpm_bffile {params.cpm_bf} --diffraw {input} --outdir {params.outpath}        
        '''

rule GSEA:
    input:
         out_dir + '/CPM{CPMcutoff}_mincell{mincell}/pval{pval}_qval{qval}_logfc{logfc}/{sample}_{direction}_Diffxpy.csv',
    params:
         outpath = out_dir + '/CPM{CPMcutoff}_mincell{mincell}/pval{pval}_qval{qval}_logfc{logfc}',
         database = config['database'],
         script_dir = config['script_dir'],
         cluster = '{cluster}',
         direction = '{direction}'
    output:
         out_dir + '/CPM{CPMcutoff}_mincell{mincell}/pval{pval}_qval{qval}_logfc{logfc}/GSEAresult/{sample}_{cluster}/{sample}_{cluster}_{direction}_clusterprofiler_result.csv'
    shell:
         '''
         export PATH=/public/home/lizw/anaconda3/envs/diffxpy_v3/bin/:$PATH
         Rscript {params.script_dir}/clusterprofiler_G1-G3.R --DEGfile {input} --database {params.database} --outdir {params.outpath} --cluster "{params.cluster}" --direction {params.direction}
         ''' 
