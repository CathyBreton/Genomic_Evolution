![Logo](images/Musa_NGS_Suite.png)

# *Musa_NGS_Suite_Replicate*

Purpose of Musa_NGS_Suite_Replicate
------------------------------------------

<div align="justify">
Musa NGS Suite, is a collection of scripts that are used in different projects to analyse banana genomes using NGS datasets. 
This suite comprises two main parts. The first part starts from the cleaning of reads to SNP calling. This procedure applied in Elan et al, crop science 2020.
The second part adds additional steps to study the structure of chromosomes that we applied for the detection of large insertion and deletions to check the genetic integrity of banana samples maintained in collection
</div>


Dependencies
------------
The tools are developed in Perl, bash, Python3, Java and work on the Linux system and require:

| Tools  | Website | Version |
| ------ | ------- | ------- |
| Bamtools      | https://github.com/pezmaster31/bamtools                         | bamtools/2.4.0 |
| BWA           | http://bio-bwa.sourceforge.net                                  | bwa/0.7.12 |
| Cutadapt      | https://cutadapt.readthedocs.io/en/stable/                      | cutadapt/2.10  |
| FastQC        | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/      | FastQC/0.11.7 |
| GATK V4       | https://software.broadinstitute.org/gatk/                       | GenomeAnalysisTK/4.1.9.0 |
| GATK V3       | https://software.broadinstitute.org/gatk/                       | GenomeAnalysisTK/3.7-0   |
| Picard Tools  | https://broadinstitute.github.io/picard/                        | picard-tools/2.7.0   |
| sambamba      | https://lomereiter.github.io/sambamba/                          | sambamba/0.6.6 |
| Samtools      | https://github.com/samtools/samtools                            | samtools/1.2  |
| STAR          | https://github.com/alexdobin/STAR                               | STAR/2.5.0b |
| VCFHunter     | https://github.com/SouthGreenPlatform/VcfHunter                 |  |
| Vcftools      | https://vcftools.github.io/index.html                           | vcftools/0.1.14  |

For the installation of the differents software, in is under construction.
Be carefull, those pipeline are part of a cluster architecture, so change the path of the software inside the scipt, this is a futur improvement of the pipeline. 


Sun Grid Engine (SGE) and SLURM job scheduler concepts are quite similar.
All the script are developped to use them in SGE systeme HPC. 
In order to adapt those script to SLURM , one script is added, but not finish in the SLURM folder.

