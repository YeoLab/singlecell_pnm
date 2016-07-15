import glob
import os
import numpy as np
import pandas as pd
import pybedtools
import HTSeq
from flotilla.util import timestamp

import cPickle as pickle



folder = '/projects/ps-yeolab/obotvinnik/singlecell_pnms/'

csv_folder = '{}/csvs_for_paper'.format(folder)
conservation_folder = '{}/conservation'.format(csv_folder)

splicing_feature_folder = '{}/splicing_feature_data'.format(csv_folder)
alternative_feature_folder = '{}/alternative'.format(splicing_feature_folder)
constitutive_feature_folder = '{}/constitutive'.format(splicing_feature_folder)

alt_exons_bedfile = '{}/exons.bed'.format(alternative_feature_folder)
constitutive_bedfile = '{}/exons.bed'.format(constitutive_feature_folder)


memmap_dir = '/home/obotvinnik/projects/singlecell_pnms/analysis/htseq_memmap'
filename = '{}/hg19_phastcons_placental_mammal_htseq.pickle'.format(memmap_dir)
with open(filename) as f:
    conservation = pickle.load(f)


bedfiles = {'alternative': alt_exons_bedfile, 'constitutive': constitutive_bedfile}

directions = 'upstream', 'downstream'

nt = 200


for exon_type, bedfile in bedfiles.items():
    bed = pybedtools.BedTool(bedfile)

    basename = os.path.basename(bedfile)
    prefix = basename.split('.bed')[0]
    
    for direction in directions:
        if direction == 'downstream':
            # Get downstream intron
            intron = bed.flank(l=0, r=nt, s=True, g=pybedtools.chromsizes('hg19'))
        else:
            intron = bed.flank(l=nt, r=0, s=True, g=pybedtools.chromsizes('hg19'))
        
        # Get just unique upstream,/downstream
        intron = pybedtools.BedTool(list(set(x for x in intron)))
        nrow = len(intron)
        ncol = nt
        array = np.zeros(shape=(nrow, ncol), dtype=float)
        junction_ids = pd.Series(index=np.arange(nrow))

        print 'Iterating over {} intervals in {} {}nt of {} ...'.format(nrow, direction, nt, basename)
        for i, interval in enumerate(intron):
            if (i+1) % 10000 == 0:
                print '\t{}\t{}/{}'.format(timestamp(), i+1, nrow)
            junction_ids[i] = interval.name
            region = conservation[HTSeq.GenomicInterval(str(interval.chrom), interval.start, interval.stop)]
            count = sum(1 for _ in region.values())
            subset = np.fromiter(region.values(),
                                 dtype=float, count=count)
            if interval.strand == '-':
                subset = subset[::-1]
                j = nt - count
                array[i][j:] = subset
            else:
                j = count
                array[i][:j] = subset
        intron_conservation = pd.DataFrame(array, index=junction_ids.values)
        filename = '{}/{}_{}_placental_mammal_conservation.csv'.format(splicing_feature_folder, prefix, direction, nt)
        print '\t', filename
        intron_conservation.to_csv(filename)