## RNA-Seq Exp Workflow

### Softwares

| Name         | Version                                                  |
| ------------ | -------------------------------------------------------- |
| hisat2       | 2.2.1                                                    |
| samtools     | 1.14                                                     |
| stringtie    | 2.2.1                                                    |
| ballgown     | 2.26.0                                                   |
| count        | Compatiable<br> Link to prepDE.py from stringtie package |
| qualimap     | 2.2.2d                                                   |
| fastq-screen | 0.15.2                                                   |
| fastqc       | 0.11.9                                                   |
| fastp        | 0.23.2                                                   |

### Installation and Testing

```
mamba create -n rnaseq-workflow hisat2=2.2.1 samtools=1.14 bioconductor-ballgown=2.26.0 bioconductor-genefilter=1.76.0 qualimap=2.2.2d fastq-screen=0.15.2 fastqc=0.11.9 fastp=0.23.2 stringtie=2.2.1
conda activate rnaseq-workflow
```