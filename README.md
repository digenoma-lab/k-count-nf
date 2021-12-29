# k-count-nf
A nextflow pipeline to count k-mers and estimate genome size from WGS data

## usage

```
#using latest github code
nextflow run  k-count-nf/main.nf -profile singularity  --reads 'reads/*R{1,2}_001.fastq.gz' 
#using latest version
nextflow run  digenoma-lab/k-count-nf -r 1.0 -profile singularity  --reads 'reads/*R{1,2}_001.fastq.gz' 
```


