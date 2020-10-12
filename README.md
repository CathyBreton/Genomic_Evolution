# Genomic_Evolution
This repesitory is a collection script used to analyse NGS data in Evolution.



Skip to content
Pull requests
Issues
Marketplace
Explore
@CathyBreton
CathyBreton /
Genomic_Evolution

0
0

    0

Code
Issues
Pull requests
Actions
Projects
Wiki
Security
Insights

    Settings

Genomic_Evolution/Musa_NGS_Suite/

1

![Logo](images/Musa_NGS_Suite.png)

2

​

3

# *Musa_NGS_Suite*

4

​

5

Purpose of Musa_NGS_Suite

6

--------------------------

7

​

8

<div align="justify">

9

Musa NGS Suite, is a collection of scripts that are used in different projects to analyse banana genomes using NGS datasets. 

10

This suite comprises two main parts. The first part starts from the cleaning of reads to SNP calling. This procedure applied in Elan et al, crop science 2020.

11

The second part adds additional steps to study the structure of chromosomes that we applied for the detection of large insertion and deletions to check the genetic integrity of banana samples maintained in collection

12

</div>

13

​

14

​

15

Dependencies

16

------------

17

The tools are developed in Perl, bash, Python3, Java and work on the Linux system and require:

18

​

19

| Tools  | Website | Version |

20

| ------ | ------- | ------- |

21

| Bamtools      | https://github.com/pezmaster31/bamtools                         | bamtools/2.4.0 |

22

| BWA           | http://bio-bwa.sourceforge.net                                  | bwa/0.7.12 |

23

| Cutadapt      | https://cutadapt.readthedocs.io/en/stable/                      | cutadapt/2.10  |

24

| FastQC        | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/      | FastQC/0.11.7 |

25

| GATK V4       | https://software.broadinstitute.org/gatk/                       | GenomeAnalysisTK/4.0.5.2 |

26

| GATK V3       | https://software.broadinstitute.org/gatk/                       | GenomeAnalysisTK/3.7-0   |

27

| Picard Tools  | https://broadinstitute.github.io/picard/                        | picard-tools/2.7.0   |

28

| sambamba      | https://lomereiter.github.io/sambamba/                          | sambamba/0.6.6 |

29

| Samtools      | https://github.com/samtools/samtools                            | samtools/1.2  |

30

| STAR          | https://github.com/alexdobin/STAR                               | STAR/2.5.0b |

31

| VCFHunter     | https://github.com/SouthGreenPlatform/VcfHunter                 |  |

32

| Vcftools      | https://vcftools.github.io/index.html                           | vcftools/0.1.14  |

33

​

34

​

35

​

36

<details>

37

<summary>Table of content</summary>

38

​

39

## Table of contents

40

​

41

- [**How to cite**](#How-to-cite)

42

- [**Introduction**](#Introduction)

43

  - Genomic Complexity Reduction

44

  - DARTseq

@CathyBreton
Commit changes
Commit summary
Optional extended description
Commit directly to the master branch.
Create a new branch for this commit and start a pull request. Learn more about pull requests.
© 2020 GitHub, Inc.
Terms
Privacy

    Security
    Status
    Help
    Contact GitHub
    Pricing
    API
    Training
    Blog
    About


