import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.ipython.html
from jpy_tools.rTools import py2r, r2py, r_inline_plot, rHelp, trl, rSet, rGet, ad2so, so2ad
from jpy_tools import loadPkl, toPkl
rBase = importr('base')
rUtils = importr('utils')
dplyr = importr('dplyr')
R = ro.r
R("options(expressions = 5e5)")
R("options(browser='firefox', shiny.port=6511)")

import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.font_manager as font_manager
# plt.rcParams['figure.dpi'] = 150
font_dirs = ["/public/home/mowp/test/fonts/"]
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    font_manager.fontManager.addfont(font_file)
plt.rcParams["font.family"] = "Arial"
sns.despine(top=True, right=True)
from itertools import product
from functools import reduce
import patchworklib as pw
def show_figure(fig):
    dummy = plt.figure()
    new_manager = dummy.canvas.manager
    new_manager.canvas.figure = fig
    fig.set_canvas(new_manager.canvas)
    fig.show()
    plt.show()
pw.show = show_figure
from jpy_tools.otherTools import pwRecoverSeaborn
fc_recoverSns = pwRecoverSeaborn()

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as ss
from tqdm import tqdm
from cool import F
from loguru import logger

from jpy_tools import singleCellTools

import diffxpy.api as de
import click
@click.command()
@click.option('--h5ad',type = click.File('rb'),help='h5ad file for 10xscRNAseq')
@click.option('--cpmfile',help='sample cpm file')
@click.option('--cpmcutoff',type = float, help='cpm cutoff')
@click.option('--mincell',type = int, help='mincell express gene for the comparing samles.')
@click.option('--output',type = click.File('wb'),help='result')
def main(h5ad,cpmfile,cpmcutoff,mincell,output):
    ad = sc.read_h5ad(h5ad)
    list(ad.obs['Sample'].cat.categories)
    ls_sample = list(ad.obs['Sample'].cat.categories)
    lsDf = []

    cpm_h5ad = sc.read_h5ad(cpmfile)
    cpm_table = cpm_h5ad.to_df().T
    bf_index = cpm_table[cpm_table['BF-6h'] >=cpmcutoff].index
    gmi_index = cpm_table[cpm_table['GMI-6h'] >=cpmcutoff].index
    b417_index = cpm_table[cpm_table['B417-6h'] >=cpmcutoff].index

    for cellType, _ad in singleCellTools.basic.splitAdata(ad, 'leiden_r0.6_anno_v2', needName=True):
        for sample in ['B417-6h','GMI-6h']:
            _ls = [sample, 'BF-6h']
            _ad_forDiffxpy = _ad[_ad.obs.eval("Sample in @_ls ")].copy()
            _ad_forDiffxpy.X = _ad_forDiffxpy.layers['raw'].A
            sc.pp.filter_genes(_ad_forDiffxpy, min_cells=mincell)
            #gene filter for sample
            if sample == 'B417-6h':
                cpm_filter_index = bf_index|b417_index
            elif sample == 'GMI-6h':
                cpm_filter_index = bf_index|gmi_index
            _ad_forDiffxpy = _ad_forDiffxpy[:,_ad_forDiffxpy.var.index.isin(cpm_filter_index)]

            _ad_forDiffxpy.obs = _ad_forDiffxpy.obs.assign(
                diffxpy_temp=lambda df: np.where(df["Sample"] == sample, 1, 0)
            ).assign(
                diffxpy_temp=lambda df: df["diffxpy_temp"]
                .astype("category")
                .cat.set_categories([0, 1])
            )

            diffxpyTestResult = de.test.wald(
                _ad_forDiffxpy,
                formula_loc="~1 + diffxpy_temp",
                coef_to_test="diffxpy_temp[T.1]",
                noise_model="nb",
                quick_scale=True,
            )
            df_diffxpyResult = diffxpyTestResult.summary().assign(cluster = cellType, sample=sample, sampleBg='BF-6h')
            lsDf.append(df_diffxpyResult)
            del(_ad_forDiffxpy)
        del(_ad)
    df_diffxpyResultOneVsOne = pd.concat(lsDf)
    df_diffxpyResultOneVsOne.to_csv(output,sep='\t')

if __name__ == '__main__':
    main()
