
import gffutils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import re
from outrigger.junctions_to_events import consolidate_junction_events, get_isoform_transcripts
from outrigger.psi import calculate_psi 

import seaborn as sns
sns.set(style='ticks', context='talk')
folder = '/projects/ps-yeolab/obotvinnik/singlecell_pnms'
csv_folder = '/projects/ps-yeolab/obotvinnik/singlecell_pnms/csvs_for_paper'


v19db_filename = '/projects/ps-yeolab/genomes/hg19/gencode/v19/gencode.v19.annotation.gtf.db'
v19db = gffutils.FeatureDB(v19db_filename)

mxe_consolidated = pd.read_csv('{}/mutually_exclusive_exon_consolidated_events.csv'.format(csv_folder), index_col=0)
se_consolidated = pd.read_csv('{}/skipped_exon_consolidated_events.csv'.format(csv_folder), index_col=0)

sj = pd.read_csv('{}/sj_raw.csv'.format(csv_folder))
    
sj = sj.set_index(['intron_location', 'sample_id'])
sj =  sj.sort_index()



inclusion_cols = ['junction_exon12', 'junction_exon23']
exclusion_cols = ['junction_exon13']
junction_cols = inclusion_cols + exclusion_cols

se_psi = calculate_psi(se_consolidated, sj, reads_col='unique_junction_reads',
                       isoform2_junctions=inclusion_cols, isoform1_junctions=exclusion_cols)

## MXE Psi
isoform1_cols = ['junction13', 'junction34']
isoform2_cols = ['junction12', 'junction24']
mxe_psi = calculate_psi(mxe_consolidated, sj, reads_col='unique_junction_reads',
                       isoform1_junctions=isoform1_cols, isoform2_junctions=isoform2_cols)
psi = pd.concat([se_psi_consolidated, mxe_psi], axis=1)
print psi.shape

fig, ax = plt.subplots()
bins = np.linspace(0, 1, 50)
sns.distplot(psi.values.flat, bins=bins, kde=False)
fig.savefig('{}/psi_values_raw.pdf'.format(figure0))


psi.to_csv("{}/psi_unfiltered.csv".format(csv_folder))

notnull = psi.notnull()

constitutively0 = (psi == 0)[notnull].all()

constitutively1 = (psi == 1)[notnull].all()
alternative = psi.columns[~constitutively0 & ~constitutively1]
print len(alternative)


constitutively0 = constitutively0[constitutively0].index
constitutively1 = constitutively1[constitutively1].index

print len(constitutively0)
print len(constitutively1)

constitutive = constitutively0 | constitutively1
print len(constitutive)

psi_constitutive1 = psi[constitutively1]
psi_constitutive1.to_csv('{}/psi_constitutively1.csv'.format(csv_folder))
psi_constitutive1.head()

psi_constitutive0 = psi[constitutively0]
psi_constitutive0.to_csv('{}/psi_constitutively0.csv'.format(csv_folder))
psi_constitutive0.head()

psi_alternative = psi[alternative]
psi_alternative.to_csv('{}/psi.csv'.format(csv_folder))
psi_alternative.to_csv('{}/splicing.csv'.format(csv_folder))
psi_alternative.head()