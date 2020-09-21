#! /usr/bin/perl -w
#$ -q bigmem.q
#$ -cwd
#$ -V
#$ -t 1-6
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1
#$ -l mem_free=100G

#charge module
#use "bioinfo/bwa/0.7.2";
#use "system/java/jre8";
#use "bioinfo/FastQC/0.11.3";
#use "compiler/gcc/4.9.2";
#use "bioinfo/bwa/0.7.12";
#use "bioinfo/GATK/4.0.5.2" 



my $gatk3_path = "/usr/local/bioinfo/GenomeAnalysisTK/3.7-0";
my $gatk_path = "/usr/local/bioinfo/GenomeAnalysisTK/4.0.5.2";
my $java_path = "/usr/local/java/jre8/bin";
my $bwa_path = "/usr/local/bioinfo/bwa/0.7.12";
my $picard_path = "/usr/local/bioinfo/picard-tools/1.130";
my $fastqc_path = "/usr/local/bioinfo/FastQC/0.11.3";
my $cutadapt_path = "/usr/local/bioinfo/python/2.7.9_build2/bin";
my $samtools_path = "/usr/local/bioinfo/samtools/1.2/bin";
my $star_path = "/usr/local/bioinfo/STAR/2.5.0b/bin";

###################
#
#  Licencied under GNU General Public License : GPL-3.0-or-later 
#
#  Intellectual property belongs to Bioversity
#
#  Written  Catherine Breton and Yann Huebert
#
#  Contact : Catherine Breton c.breton@cgiar.org
#
#  Script Name : Rnaseq_single_end_fastq_to_Vcf_gvcf_jobarray_Total_GATK4.pl
#
#  Usage : perl Radseq_single_end_fastq_to_Vcf_gvcf_jobarray_Total_GATK4.pl -r <reference fasta file assembly> -x <extension gvcf> -cu <cultivars>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses/> or
#  write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
#################### 


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;


############### MAIN ##############


if (@ARGV < 6)
{
        print "\n\nUsage: perl $0 -g genome_dir_path -r reference_fasta -x extension_file_to_treat -cu cultivar\n\n\n";
        exit 1;
}


#Variables
my ($extension, $reference, $genome_indexes_dir, $cultivar);

#Getoptions
GetOptions("x|extension=s" => \$extension, "r|reference=s" => \$reference, "g|genomeindexesdir=s" => \$genome_indexes_dir, "c|cultivar=s" => \$cultivar) or die ("Error in command line arguments\n");



#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");  
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


#relative file path 
my $sge_id = `echo \$SGE_TASK_ID` or die ("Cannot catch SGE_TASK_ID\n");
chomp($sge_id);
$sge_id--; #perl is in basis 0
my $input_se_fastq = $files[$sge_id]; #se means single end


#Filename without extension
my $filenm_root = fileparse($input_se_fastq, qr/\.$extension$/);


#Programs 

#Discard filtered reads 
#my $discard_filtered_reads = "zcat $input_se_fastq | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v \"^--\$\" | gzip > $filenm_root.RPF.fastq.gz"; #RPF = reads primary filtered
#print ("Filter out bad reads on file $input_se_fastq\n");
#system("$discard_filtered_reads");

#Cutadapt
my $rm_adapter = "$cutadapt_path/cutadapt -b AGATCGGAAGAGC -O 7 -m 30 -q 20,20 $input_se_fastq | gzip > $filenm_root.cutadapt.fastq.gz"; #Remove adaptator, Trim on quality, discard read with length inf to 30 
print ("cutadapt on $filenm_root.fastq.gz\n");
system("$rm_adapter");

#fastqc
my $fastqc = "$fastqc_path/fastqc -j $java_path/java $filenm_root.cutadapt.fastq.gz";
print ("fastqc on $filenm_root.cutadapt.fastq.gz\n");
system("$fastqc");

#STAR mapping
my $star_mapping = "$star_path/STAR --genomeDir $genome_indexes_dir --readFilesIn $filenm_root.cutadapt.fastq.gz --readFilesCommand zcat --outSAMunmapped Within --outFileNamePrefix $filenm_root. --outSAMmapqUnique 255 --twopassMode Basic";
print ("Mapping with STAR (2-pass mode) on file $filenm_root.cutadapt.fastq.gz...\n");
system("$star_mapping");

#Add read group, and sort bam picardtools
my $condition = $filenm_root;
my $add_rg = "$java_path/java -Xmx1g -jar $picard_path/picard.jar AddOrReplaceReadGroups I=$filenm_root.Aligned.out.sam O=$filenm_root.RG.sorted.bam SO=coordinate RGLB=$cultivar RGPL=illumina RGPU=run RGSM=$condition RGID=$condition";
print ("Adding Read Group and sorting on file $filenm_root.Aligned.out.sam...\n");
system("$add_rg");

#Markduplicates and create index picardtools
my $markduplicates = "$java_path/java -Xmx1g -jar $picard_path/picard.jar MarkDuplicates I=$filenm_root.RG.sorted.bam O=$filenm_root.dedupped.RG.sorted.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$filenm_root.output.mectrics";
print ("Marking duplicates and creating index on file $filenm_root.RG.sorted.bam...\n");
system("$markduplicates");


#=================================================================================================================================
#SplitNCigarReads with GATK3 ######## add -U ALLOW_N_CIGAR_READS in case... For Old data
#my $split_and_trim = "$java_path/java -Xmx10g -XX:ParallelGCThreads=1 -jar $gatk3_path/GenomeAnalysisTK.jar -T SplitNCigarReads -R $reference -I $filenm_root.dedupped.RG.sorted.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -o $filenm_root.split.dedupped.RG.sorted.bam";
#print ("Split N CigarReads on file $filenm_root.dedupped.RG.sorted.bam...\n");
#system("$split_and_trim");

#SplitNCigarReads with GATK4 ########  (--skip-mapping-quality-transform / -skip-mq-transform  pour les données anciennes)
my $split_and_trim = "$gatk_path/gatk --java-options -Xmx1G SplitNCigarReads -R $reference -I $filenm_root.dedupped.RG.sorted.bam -O $filenm_root.split.dedupped.RG.sorted.bam";
print ("Split N CigarReads on file $filenm_root.dedupped.RG.sorted.bam...\n");
system("$split_and_trim");

#=================================================================================================================================

##### Indel Realignement, 
# This step is not present in GATK4, not necessary to do for SNP calling but for Indel , need to do so this step is done with GATK3

#RealignerTargetCreator with GATK3 ######Add --fix_misencoded_quality_scores option if needed
my $realigner_target = "$java_path/java -Xmx4g -XX:ParallelGCThreads=1 -jar $gatk3_path/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $filenm_root.split.dedupped.RG.sorted.bam --filter_reads_with_N_cigar -o $filenm_root.forIndelRealigner.intervals";
print ("Create RealignTargetCreator intervals with $filenm_root.split.dedupped.RG.sorted.bam...\n");
system("$realigner_target");

#IndelRealigner ######Add --fix_misencoded_quality_scores option if needed
my $indelrealigner = "$java_path/java -Xmx4g -XX:ParallelGCThreads=1 -jar $gatk3_path/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $filenm_root.split.dedupped.RG.sorted.bam --filter_reads_with_N_cigar -targetIntervals $filenm_root.forIndelRealigner.intervals -o $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam";
print ("IndelRealigning on file $filenm_root.split.dedupped.RG.sorted.bam...\n");
system("$indelrealigner");

#=================================================================================================================================


#=================================================================================================================================
#
##### First HaplotypeCaller to calibrate the bam.
#
#=================================================================================================================================


#HaplotypeCaller GATK4 (pour BaseRecalibrator)
my $haplotype_caller = "$gatk_path/gatk --java-options -Xmx4g HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam -ploidy 2 --dont-use-soft-clipped-bases -O $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.vcf";
print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$haplotype_caller");


#HaplotypeCaller (for BaseRecalibrator)  ###### Pour une ancienne version GATK3.7
#my $haplotype_caller = "$java_path/java -Xmx4g -XX:ParallelGCThreads=4 -jar $gatk_path/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference -I $filenm_root.dedupped.RG.sorted.bam -ploidy 2 -O $filenm_root.RG.sorted.HC.vcf";
#print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam...\n");
#system("$haplotype_caller");


#================================================ VariantAnnotator GATK4 is Beta version so use Gatk3
############VariantAnnotator (add annotation: MappingQualityZero)
#my $variant_annotation_First = "$gatk_path/gatk --java-options -Xmx3G VariantAnnotator -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.vcf -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam -A MappingQualityZero -O $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.vcf";
#print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.vcf...\n");
#system("$variant_annotation_First");

#VariantAnnotator GATK3 (add annotation: MappingQualityZero)
my $variant_annotation_First = "$java_path/java -Xmx4g -jar $gatk3_path/GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.vcf -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam -A MappingQualityZero -o $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.vcf";
print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.vcf...\n");
system("$variant_annotation_First");
#================================================ VariantAnnotator


#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
my $variant_filtration_First  = "$gatk_path/gatk --java-options -Xmx4G VariantFiltration -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' --mask-extension 0 --filter-name HARD_TO_VALIDATE --filter-name QD_FILTER --filter-name DP_FILTER --mask-name Mask -O $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.vcf";
print ("Variant filtration --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' on $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.vcf...\n");
system("$variant_filtration_First");

#SelectVariant to filter out filtered variants
my $select_variant_First = "$gatk_path/gatk --java-options -Xmx4g SelectVariants -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -select-type SNP -O $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.filtered.vcf";
print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.vcf...\n");
system("$select_variant_First");


#================================================  Base Recalibration 
#BaseRecalibrator
my $base_recalibrator = "$gatk_path/gatk --java-options -Xmx4g BaseRecalibrator -R $reference -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam --known-sites $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.filtered.vcf -O $filenm_root.recal_data.table";
print ("Create a recalibration table file with BaseRecalibrator with knownSites in file $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.filtered.vcf with BAM file $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam...\n");
system("$base_recalibrator");

#Apply recalibration with ApplyBQSR not PrintReads
my $apply_recalibration = "$gatk_path/gatk --java-options -Xmx4g ApplyBQSR -R $reference -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam -bqsr $filenm_root.recal_data.table -O $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.bam";
print ("Apply recalibration with ApplyBQSR on file $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam...\n");
system("$apply_recalibration");
#================================================  Base Recalibration 

#=================================================================================================================================
#
##### Second HaplotypeCaller on recalibrated Bam.
#
#=================================================================================================================================


##### Second HaplotypeCaller on recalibrated Bam.
#HaplotypeCaller GATK4 (pour BaseRecalibrator)
my $haplotype_caller_vcf = "$gatk_path/gatk --java-options -Xmx4g HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.bam -ploidy 2 -stand-call-conf 20.0 --dont-use-soft-clipped-bases -O $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.vcf";
print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam ...\n");
system("$haplotype_caller_vcf");


#================================================ VariantAnnotator
#VariantAnnotator GATK4 BETA (add annotation: MappingQualityZero)
#my $variant_annotation_Second = "$gatk_path/gatk --java-options -Xmx4G VariantAnnotator -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.vcf -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.bam -A MappingQualityZero -O $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.ann.vcf";
#print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.vcf...\n");
#system("$variant_annotation_Second");

#VariantAnnotator GATK3 (add annotation: MappingQualityZero)
my $variant_annotation_Second = "$java_path/java -Xmx4g -jar $gatk3_path/GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.vcf -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.bam -A MappingQualityZero -o $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.ann.vcf";
print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.vcf...\n");
system("$variant_annotation_Second");
#================================================ VariantAnnotator


#VariantFiltration (Filter SNPs after applying BaseRecalibrator)
my $variant_filtration_Second  = "$gatk_path/gatk --java-options -Xmx4g VariantFiltration -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.ann.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' --mask-extension 0 --filter-name HARD_TO_VALIDATE --filter-name QD_FILTER --filter-name DP_FILTER --mask-name Mask -O $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.ann.VF.vcf";
print ("Variant filtration --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' on $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.ann.vcf...\n");
system("$variant_filtration_Second");

#SelectVariant to filter out filtered variants
my $select_variant_Second = "$gatk_path/gatk --java-options -Xmx4g SelectVariants -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -select-type SNP -O $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.ann.VF.filtered.vcf";
print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -select-type SNP on file $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.HC.ann.VF.vcf...\n");
system("$select_variant_Second");

#=================================================================================================================================


print "Finished\n";


exit 1;
