# SMS-6026-22-kidney

This repo hosts the code and report for SMS Project #6026.

# How to rerun the analysis?

To just rerun the analysis, clone this repo to your local computer and move at least the files *analysis.nf*, *nextflow.config* and the *modules* directory to Bianca. Make sure the parameters at the top of the anaysis script are up-to-date and reachable (i.e the metadata, transcript fastas and the genome fasta). Also explicitly set the outdir parameter and workDir in the nextflow.config file.

```
ml load bioinfo-tools Nextflow
nextflow run analysis.nf
```

# How to rerun the report?

A Docker container running quarto is available on Dockerhub. This need to be pulled by Singularity to your local computer and uploaded to Bianca. This container can be used to render the .qmd file in this repo.

```
singularity pull docker://vvanhoef/quarto:latest

# Upload to Bianca
scp ...

# On Bianca, first load R_packages
ml load R_packages/4.1.1

# Then instruct the container to render the report
singularity exec quarto_latest.sif quarto render index.qmd
```



