OB# scerevisiae_cbf5

Some RNASeq analyses in yeast!

# List of files contained herein:

*  Metadata
    *  [Sample sheet](sample_sheets/all_samples.xlsx):  Metadata defining the samples used.
*  Experiment 1
    *  [E1 raw reads and plots](E1/E1_written_20180606.xlsx):  Raw reads, raw plots,
        normalized reads, normalized plots for the first experiment.
    *  [E1 differential expression](E1/E1_de_20180606.xlsx):  Differentially
       expressed genes between cbf5-D95A and WT, from first experiment.
    * [E1 significant by DESeq2](E1/E1_sig_20180606.xlsx):  Significantly
       differentially expressed genes according to DESeq2, first experiment.
    * [E1 gprofiler up genes](E1/E1_gprofiler_up_20180606.xlsx):  gprofiler results
       using differentially expressed genes according to DESeq2, first experiment.
    * [E1 gprofiler down genes](E1/E1_gprofiler_down_20180606.xlsx):  gprofiler results
       using differentially expressed genes according to DESeq2, first experiment.
* Experiment 1 and Experiment 2
    * [E1E2 raw reads and plots](E1E2/E1E2_written_20180212.xlsx):  Raw reads, raw plots,
       normalized reads, normalized plots for both experiments.
    * [E1E2 differential expression](E1E2/E1E2_de_20180212.xlsx):  Differentially
       expressed genes between cbf5-D95A and WT, from both experiments.
    * [E1E2 significant by DESeq2](E1E2/E1E2_sig_20180212.xlsx):  Significantly
       differentially expressed genes according to DESeq2, both experiments.
    * [E1E2 gprofiler up genes](E1E2/E1E2_gprofiler_up_20180606.xlsx):  gprofiler results
       using differentially expressed genes according to DESeq2, first experiment.
    * [E1E2 gprofiler down genes](E1E2/E1E2_gprofiler_down_20180606.xlsx):  gprofiler results
       using differentially expressed genes according to DESeq2, first experiment.
* Worksheets
    * [Annotation R markdown](01_annotation-v20180606.Rmd): R worksheet used to
       gather annotation data for Saccharomyces cerevisiae.
    * [Annotation html](01_annotation-v20180606.html): html transcript of Ibid.
    * [Sample estimation markdown](02_sample_estimation_merged-v20180606.Rmd):
       R worksheet used to examine the raw/normalized data.
    * [Sample estimation html](02_sample_estimation_merged-v20180606.html): html transcript of Ibid.
    * [Differential expression markdown](03_differential_expression_merged-v20180606.Rmd):
       R worksheet used to perform differential expression analyses.
    * [Differential expression html](03_differential_expression_merged-v20180606.html):
       html transcript of Ibid.
    * [gProfileR ontology searches](04_gene_ontology_merged-v20180606.Rmd):
       R worksheet used to perform ontology searches using gProfileR and goseq.
    * [Differential expression html](04_gene_ontology_merged-v20180606.html): html transcript of Ibid.
