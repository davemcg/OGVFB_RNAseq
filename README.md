# OGVFB_RNAseq
Generic RNA-seq quant for OGVFB group

# Usage (bulk)
1. Copy config.yaml to your directory and create a sample <-> fastq file table (tab separated, no columns)
2. Run: `sbatch --time=8:00:00 ~/git/OGVFB_RNAseq/Snakemake.wrapper.sh config.yaml`
3. [Generic analysis template](https://github.com/davemcg/OGVFB_RNAseq/blob/master/analysis_template.Rmd)
4. 
# Usage (single cell)
1. Copy (and edit) configSC.yaml to your directory and create [sample_info.csv](https://github.com/davemcg/OGVFB_RNAseq/blob/master/sample_info_SC_example.csv)
2. fastq files in the `fastq` folder
3. Run: `sbatch --time=16:00:00 ~/git/OGVFB_RNAseq/SnakemakeSC.wrapper.sh configSC.yaml ~/git/OGVFB_RNAseq/clusterSC.json`
4. [Generic analysis template](https://github.com/davemcg/OGVFB_RNAseq/blob/master/scAnalysis/analysis_template_SC.Rmd)
