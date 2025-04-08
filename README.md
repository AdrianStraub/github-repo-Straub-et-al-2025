Data analysis:

IMMUNITY_01_Preprocessing_Demultiplexing displays data loading from cellranger output of raw sequencing data alongside demultiplexing via HashSolo and CITESeq of individual mouse donors.\Single cells and genes are filtered via calculation of QC metrics (mitochondrial genes and gene counts), as recommened by current single-cell best practrices (Heumos, L., Schaar, A.C., Lance, C. et al. Best practices for single-cell analysis across modalities. Nat Rev Genet (2023)) (see https://www.sc-best-practices.org/preamble.html)\
\
IMMUNITY_01_Preprocessing_SoupX_R displays further proprocessing of the single-cell data in R running the SoupX (Matthew D Young and Sam Behjati. SoupX removes ambient termRNA contamination from droplet-based single-cell termRNA sequencing data. GigaScience, December 2020)\ package to remove ambient mRNA from analysis as recommended by current single-cell best practrices (see https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html)\
\
IMMUNITY_02_Concatenation_and_feature_selection displays concatenation of processed single-cell data, normalizatzion of gene expression, feature selection and batch correction.\
\
IMMUNITY_03_figures displays data analysis as presented in the submitted figures and phenotype enrichment score calculations as described by Abdullah et al. (Abdullah L, Emiliani FE, Vaidya CM, Stuart H, Musial SC, Kolling FW, Obar JJ, Rosato PC, Ackerman ME, Song L, McKenna A, Huang YH. The endogenous antigen-specific CD8+ T cell repertoire is composed of unbiased and biased clonotypes with differential fate commitments. Immunity. 2025 Mar 11;58(3):601-615.e9.)
