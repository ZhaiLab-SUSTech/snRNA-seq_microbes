import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import click
@click.command()
@click.option('--qval',type = float, help ='min qval for DEG filter')
@click.option('--pval',type = float, help ='min pval for DEG filter')
@click.option('--logfc',type = float, help ='min abs(logfc) for DEG filter')
@click.option('--cpm_b417file',type = click.File('rb'), help ='Cluster_Cpm_for b417')
@click.option('--cpm_gmifile',type = click.File('rb'), help ='Cluster_Cpm_for gmi')
@click.option('--cpm_bffile',type = click.File('rb'), help ='Cluster_Cpm_for bf')
@click.option('--diffraw',type = str, help ='input diffxpy raw result')
@click.option('--outdir',type = str, help ='outdir')
def main(qval,pval,logfc,cpm_b417file,cpm_gmifile,cpm_bffile,diffraw,outdir):
    cpm_b417 = sc.read_h5ad(cpm_b417file).to_df().unstack().reset_index()
    cpm_gmi =  sc.read_h5ad(cpm_gmifile).to_df().unstack().reset_index()
    cpm_bf = sc.read_h5ad(cpm_bffile).to_df().unstack().reset_index()
    cpm_b417.columns = ['gene','cluster','CPM']
    cpm_gmi.columns = ['gene','cluster','CPM']
    cpm_bf.columns = ['gene','cluster','CPM_Bg']
    diff_result = pd.read_csv(diffraw,sep='\t',dtype={'cluster':'str'})
    diff_result = diff_result.query('sampleBg == "BF-6h"')
    diff_up_index = diff_result.query('(pval < @pval ) & (qval < @qval) & (log2fc > @logfc)').index
    diff_down_index = diff_result.query('(pval < @pval ) & (qval < @qval) & (log2fc < -@logfc)').index
    diff_result['direction'] = 'Fail'
    diff_result.loc[diff_up_index,'direction'] = 'Up'
    diff_result.loc[diff_down_index,'direction'] = 'Down'

    order_list = ['0',
    '2',
    '6',
    '8',
    '13',
    '16',
    '18',
    '20',
    'Not root',
    'Meristematic zone',
    'Pericycle',
    'Procambium',
    'Cortex',
    'Endodermis',
    'Hair cell',
    'No hair cell',
    'Phloem',
    'Xylem']
    order_dict = {order_list[x]:x for x in range(len(order_list))}

    import os
    folder = os.path.exists(f'{outdir}/pval{pval:.2f}_qval{qval:.2f}_logfc{logfc:.2f}')
    if not folder:
        os.makedirs(f'{outdir}/pval{pval:.2f}_qval{qval:.2f}_logfc{logfc:.2f}')

    for sample in ['B417-6h','GMI-6h']:
        diff_sample = diff_result.query('sample == @sample')
        cpm_sample = cpm_b417 if sample == 'B417-6h' else cpm_gmi
        diff_sample_withCPM = pd.merge(diff_sample,cpm_sample,on=['gene','cluster'],how='left')
        diff_sample_withCPM = pd.merge(diff_sample_withCPM,cpm_bf,on=['gene','cluster'],how='left')
        diff_sample_withCPM_UP = diff_sample_withCPM.query('(direction == "Up")')
        diff_sample_withCPM_DOWN = diff_sample_withCPM.query('(direction == "Down")')
        
        diff_sample_withCPM_UP.to_csv(f'{outdir}/pval{pval:.2f}_qval{qval:.2f}_logfc{logfc:.2f}/{sample}_Up_Diffxpy.csv')
        diff_sample_withCPM_DOWN.to_csv(f'{outdir}/pval{pval:.2f}_qval{qval:.2f}_logfc{logfc:.2f}/{sample}_Down_Diffxpy.csv')

        #draw
        font_color = '#525252'
        hfont = {'fontname':'Calibri'}
        color_red = '#ff7f0e'
        color_blue = '#1f77b4'
        column0 = diff_sample_withCPM_DOWN.groupby(
            'cluster').count(
            )['gene'].sort_index(key = lambda x:x.map(order_dict)
            )
        column1 = diff_sample_withCPM_UP.groupby(
            'cluster').count(
            )['gene'].sort_index(key = lambda x:x.map(order_dict)
            )
        title0 = 'Down'
        title1 = 'Up'
        fig, axes = plt.subplots(figsize=(10,5), ncols=2, sharey=True)
        fig.tight_layout()
        axes[0].barh(column0.index, column0, align='center', color=color_blue, zorder=10, )
        axes[0].set_title(title0, fontsize=18, pad=15, color=color_blue, **hfont)
        axes[1].barh(column1.index, column1, align='center', color=color_red, zorder=10)
        axes[1].set_title(title1, fontsize=18, pad=15, color=color_red, **hfont)

        axes[0].invert_xaxis()
        plt.gca().invert_yaxis()
        axes[0].set(yticks=column1.index, yticklabels=column1.index)
        axes[0].yaxis.tick_left()

        plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
        plt.savefig(f'{outdir}/pval{pval:.2f}_qval{qval:.2f}_logfc{logfc:.2f}/{sample}_DEG_barplot.png',format='png',bbox_inches='tight',dpi=100)

if __name__ == '__main__':
    main()
