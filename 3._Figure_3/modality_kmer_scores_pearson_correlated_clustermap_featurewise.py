#!/usr/bin/env python

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

sns.set(style='ticks', context='paper', rc={'font.sans-serif':'Arial', 'pdf.fonttype': 42})


import flotilla
flotilla_dir = '/projects/ps-yeolab/obotvinnik/flotilla_projects'

study = flotilla.embark('singlecell_pnm_figure2_modalities_bayesian_kmers', flotilla_dir=flotilla_dir)

corr = study.supplemental.kmer_zscores.fillna(0).T.corr()
print corr.shape
corr = corr.dropna(how='all', axis=1).dropna(how='all', axis=0)
print corr.shape

folder = '/home/obotvinnik/Dropbox/figures2/singlecell_pnm/figure2_modalities/bayesian'
figure_folder = '{}/kmer_counting'.format(folder)

g = sns.clustermap(corr)
g.savefig('{}/modality_kmer_scores_pearson_correlated_clustermap_featurewise.pdf'.format(figure_folder))