# Single-cell alternative splicing

This is the README.md file for code accompanying the paper, "" by Song and Botvinnik, et al, *Journal* (Year). The focus of this study was to investigate alternative splicing at the single-cell level, in the system of motor neuron differentiation and applied to a model of neurodegenerative disease. 

Some abbreviations:

* The celltypes are:
    * Induced pluripotent stem cells (iPSCs)
    * Neural progenitor cells (NPCs)
    * Motor neurons (MNs)
    * Motor neurons subjected to sodium arsenite, a cellular stress. Aka stressed motor neurons (sMNs)
* `hg19` - The genome build of the human genome used in this project
* [Splice types/modes](https://en.wikipedia.org/wiki/Alternative_splicing#Modes)
    * Skipped exons (SE)
        * Isoform 1: `exon1-exon3`
        * Isoform 2: `exon1-exon2-exon3`
    * Mutually exclusive exons (MXE)
        * Isoform 1: `exon1-exon3-exon4`
        * Isoform 2: `exon1-exon2-exon4`
* Percent spliced-in (Psi/$\Psi$)
    * A value between 0 and 1 (mathematically, [0, 1]) which represents what percentage of transcripts use this particular isoform, versus another isoform, e.g. isoform1 vs isoform2.
    * For example, for skipped exon events, with `exon1-exon2-exon3` as the three possible exons, this value represents what percentage of transcripts use `exon1-exon2-exon3` versus ones that use either `exon1-exon2-exon3` AND `exon1-exon3`

## 0. Data collection and aggregation

Notebooks whose numbers begin with "0" are collecting data *before* any figures are created. Notebooks whose numbers begin 1-5 correspond to figures or data manipulation to do with that particular figure.

### 0.0 Software versions

To ensure posterity of results, this document shows the versions of all software used in this paper: Python, Perl, and command line tools.

#### Note for Olga: Submit `poshsplice`, `sj2psi`, `flotilla`, `bonvoyage`, `modish` to pypi

### 0.1 Combine Data

This notebook reads the `sailfish` and STAR output files to create gene expression and mapping stats matrices, and splice junction files. This also subsets the data separately on iPSC, NPC, and MN samples (Figures 1-4) and saves the stressed motor neurons (sMNs for Figure 5) separately.

### 0.2 Splicing annotations

A key aspect of this study is the annotation of alternative splicing events. We calcluate alternative splicing of skipped exon and mutually exclusive exons, based on read which uniquely mapped to a junction between exons. 

As this is a large endeavor, the annotation steps are broken down into the following notebooks (single digits are prefixed by "0" for sorting purposes):

1. Get annotated exons corresponding to splice junctions
    1. Identify all annotated (Gencode v19 - all levels) exons which are exactly upstream and downstream of all junctions.
    2. Create (junction, direction, exon) triples for storage in a graph database using [`graphlite`](http://eugene-eeo.github.io/graphlite/)
2. Create SE, MXE, and constitutive exon annotations *de novo* from (1)
    1. Aggregate exons and junctions into *de novo* splicing annotations of skipped exon (SE) and mutually exclusive exon (MXE) events. 
    2. If multiple flanking exons are available, use the shortest CDS sequence. If no CDS, the shortest sequence.
3. Calculate percent spliced-in (Psi) scores for MXE and SE events from (2)
    1. Calculate percent spliced-in (Psi/$\Psi$) of SE and MXE events, based on junction reads.
4. Splicing gene id, name annotations with event definitions from (2)
    1. Start a splicing event annotation table.
    2. Get gene ids corresponding to SE and MXE event annotations
    3. Use these ids to intersect with previously created gene annotations.
5. Extract constitutive and alternative exons, write to bed file, submit job to calculate conservation
    1. Create bed files of alternative and constitutive exons.
    2. Submit a compute job to calculate mean placental mammal conservation of exon bodies (`bigWigAverageOverBed`).
6. Calculate exon properties (conservation, splice site strength)
    2. Use bed files from (5) to calculate splice site strength ([Yeo et al, *J Comput Biol* 2004](http://www.ncbi.nlm.nih.gov/pubmed/15285897)
    3. Read calcualted conservation from (5)
    4. Add this information back to the splicing event annotation table.
7. Ancient alternatively spliced exons (Merkin et al, 2012)
    1. Obtain ancient alternatively spliced exon table from [Merkin et al, *Science* 2012](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3568499/)
    2. Use `liftOver` to convert rhesus macaque (`rheMac2`) coordinates to human (`hg19`) coordinates
    2. Overlap with identified alternative exons from (2).
    3. Add this information back to the splicing event annotation table.
8. Save conservation wiggle as a memory-mapped Genomic Array
    1. Use `HTSeq` to save a 79 gigabyte memory-mapped file of the placental mammal conservation at every genomic base of the human genome.
    2. This is done because it takes a long time to populate the array (overnight on a compute node) and it's easier to read in the pickle file, which takes about 30 min.
9. Create CSVs of base-wise conservation of alternative and constitutive introns
    1. Read the pickled memory-mapped conservation file from (8)
    2. Get base-wise conservation of introns flanking alternative and constitutive exons. This will be used in Figure 2 when describing intron conservation of alternative splicing modalities.
10. Transcribe isoforms to get RNA sequence
11. Submit miRNA finding compute job on sequences from (10)
    1. Use first 17bp of microRNAs from [mirbase](http://www.mirbase.org/), since the ends of them are bound by Ago
    2. Use `RNAhybrid` to estimate which microRNAs may bind an isoform's transcript
12. Properties of the isoforms' RNA sequence from (10)
    1. Use BioPython to calculate...
        1. GC content
        2. maybe: [Estiamte RNA structure (single or double stranded)]
    3. Add this information back to the splicing event annotation table.
13. Translate isoforms to protein products
14. Submit compute job for `hmmscan` on Pfam-A, calculate disordered scores on translations from (13)
    1. Use `hmmscan` to find protein domains on the isoform translations from (13)
    2. Use `IUPred` to calculate average disordered region score of each isoform
15. Calculate properties of isoform protein products (isoelectric point, etc) on translations from (13)
    1. Properties calculated:
        1. Molecular weight
        2. Isoelectric point
        3. ...
    2. Read disordered scores from (14)
    3. Add this information back to the splicing event annotation table.
16. Read `hmmscan` output from (14) and get domain switching events
    1. Aggregate domains on pfam clans to see whether the entire concept of a protein's function has changed


### 0.3 Gene feature annotations

The gene annotations are less comprehensive and can be found in a single notebook. The steps are:

1. Get gene id, name, level, biotype from `gtf` file
2. Calculate maximum number of exons per gene
3. Intersect with annotations from [Gerstberger et al, *Nat Rev Genet* (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25365966)
    1. Transcription factors (TFs) 
        1. Really defined as any RNA biologist would define TFs: "stuff that binds DNA" - i.e. proteins which induce transcription versus chromatin remodelers vs topoisomerases are all included in this group.
    2. RNA binding proteins, which target ...
        1. mRNA
        2. snRNA
        3. tRNA
        4. ribosome
4. Intersect with annotations from [Animal TFDB](http://www.bioguo.org/AnimalTFDB/)
    1. Transcription factors
    2. Chromatin remodelers
    3. Transcription co-factors
5. Intersect with Phylostratum annotations from [Domazet-Loso and Tautz, *Mol Biol Evol* (2008)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2582983/)
    1. This gives the "age" of the gene, relative to the origin of life
6. Intersect with genes associated with the cell cycle from Gene Ontology
7. Cell surface markers [from ???????????????]

### 0.4 Filter out non-alternative splicing events and create a *flotilla* `Study`

## 1. Figure 1

### 1.0 Figure 1, part 1

### 1.1 Supplementary Figure 1

### 1.2 Figure 1, part 2

## 2. Figure 2

### 2.0 Figure 2, part 1

1. Get splicing events with at least 20 single cell observations per celltype
2. Estimate per-cellytpe modality for each splicing event using [`modish`](https://github.com/olgabot/modish)

### 2.1 Supplementary Figure 2

### 2.2 Figure 2, part 2

1. Investigate enrichment of properties from (0.2.\*) in each modality

### 2.3 Figure 2, part 3

1. Subset the bed files created in 0.2.5 to get each (phenotype, modality) subset
2. Submit a compute job to find enrichment of intron motifs in (phenotype, modality) vs (phenotype, other modalities) using [HOMER](http://homer.salk.edu/homer/index.html)

### 2.4 Figure 2, part 4

1. Change in modalities between celltypes
1. Venn diagrams of modalities
2. Heatmap/checkerboard of how modalities change from one celltype to the next

### 2.5 Figure 2, part 5

1. Gene Ontology (GO) enrichment of modaliites
    1. (phenotype, modality) vs (phenotype, other modalities)
    2. (phenotype, modality) vs (other phenotype, modality)
    3. (modality in any phenotype) vs (other modality in any phenotype)
    
### 2.6 Figure 2, part 6

1. Basewise conservation of flanking introns of modality exons.
    1. Use basewise conservation CSVs created in (0.2.09)
    
## 3. Figure 3

Changes in alternative splicing using change in modality and principal component analysis (PCA) on splicing events residing in constitutively expressed genes.

## 4. Figure 4

Introduces a new approach to visualizing and finding changing splicing events using `bonvoyage` to create a "Voyage Space"

## 5. Figure 5

Application of changing splicing events from modality finding using `modish`, voyaging events from `bon voyage`, and PCA on constitutively expressed, to the model of neurodegenerative disease of cellularly stressed motor neurons.