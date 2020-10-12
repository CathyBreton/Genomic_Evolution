#! /usr/bin/perl -w
#$ -q normal.q
#$ -cwd
#$ -V
#$ -t 1-100
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1
#$ -l mem_free=120G

#charge modules
#use "bioinfo/bwa/0.7.12";
#use "system/java/jre8";
#use "bioinfo/FastQC/0.11.3";
#use "compiler/gcc/4.9.2";
#use "bioinfo/bwa/0.7.12";
#use "bioinfo/GATK/4.0.5.2";



my $gatk3_path = "/usr/local/bioinfo/GenomeAnalysisTK/3.7-0";
my $gatk_path = "/usr/local/bioinfo/GenomeAnalysisTK/4.0.5.2";
my $java_path = "/usr/local/java/jre8/bin";
my $bwa_path = "/usr/local/bioinfo/bwa/0.7.12";
my $picard_path = "/usr/local/bioinfo/picard-tools/2.7.0";
my $fastqc_path = "/usr/local/bioinfo/FastQC/0.11.7";
my $cutadapt_path = "/usr/local/bioinfo/python/2.7.9_build2/bin";
my $samtools_path = "/usr/local/bioinfo/samtools/1.9/bin";



###################
#
#  Licencied under GNU General Public License : GPL-3.0-or-later 
#
#  Intellectual property belongs to Bioversity
#
#  Written  Catherine Breton
#
#  Contact : Catherine Breton c.breton@cgiar.org
#
#  Script Name : Radseq_single_end_fastq_to_Vcf_gvcf_jobarray_Total_GATK4.pl
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
        print "\n\nUsage: perl $0 -r reference_file_path -x extension_file_to_treat -cu cultivar\n\n\n";
        exit 1;
}


#Variables
my ($extension, $reference, $cultivar);

#Getoptions
GetOptions("x|extension=s" => \$extension, "r|reference=s" => \$reference, "c|cultivar=s" => \$cultivar) or die ("Error in command line arguments\n");



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

#=================================================================================================================================





#=================================================================================================================================
#Programs 
#=================================================================================================================================

#Cutadapt
#my $cutadapt = "$cutadapt_path/cutadapt -b AGATCGGAAGAGC -O 7 -m 30 -q 20,20 $input_se_fastq | gzip > $filenm_root.cutadapt.fastq.gz"; #Remove adaptator, Trim on quality, discard read with length inf to 30 
#print ("cutadapt on $filenm_root.fastq.gz\n");
#system("$cutadapt");

#fastqc
my $fastqc = "$fastqc_path/fastqc -j $java_path/java $filenm_root.cutadapt.fastq.gz";
print ("fastqc on $filenm_root.cutadapt.fastq.gz\n");
system("$fastqc");


# bwa mapping on the reference choosen in argument -r
my $bwa_mapping = "$bwa_path/bwa mem $reference $filenm_root.cutadapt.fastq.gz > $filenm_root.sam";
print ("Mapping with BWA on file $filenm_root.cutadapt.fastq.gz...\n");
system("$bwa_mapping");

#Add read group, and sort bam , with the tools picardtools
my $condition = $filenm_root;
my $add_rg = "$java_path/java -Xmx1g -jar $picard_path/picard.jar AddOrReplaceReadGroups I=$filenm_root.sam O=$filenm_root.RG.sorted.bam SO=coordinate RGLB=$cultivar RGPL=illumina RGPU=run RGSM=$condition RGID=$condition";
print ("Adding Read Group and sorting on file $filenm_root.sam...\n");
system("$add_rg");

#Index bam, with the tools picardtools
my $index_bam = "$samtools_path/samtools index $filenm_root.RG.sorted.bam";
print ("Indexing on bam file $filenm_root.RG.sorted.bam...\n");
system("$index_bam");

#SplitNCigarReads with GATK4 ######## add -U ALLOW_N_CIGAR_READS in case... becarefull not in GATK4
####my $split_and_trim = "$gatk_path/gatk --java-options -Xmx4G SplitNCigarReads -R $reference -I $filenm_root.RG.sorted.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -O $filenm_root.Split.bam";
####print ("Split N CigarReads on file $filenm_root.RG.sorted.bam...\n");
####system("$split_and_trim");


#=================================================================================================================================

##### Indel Realignement, 
# This step is not present in GATK4, not necessary to do for SNP calling but for Indel , need to do so this step is done with GATK3


#RealignerTargetCreator ######Add --fix_misencoded_quality_scores option if needed --fix_misencoded_quality_scores 
my $realigner_target = "$java_path/java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $filenm_root.RG.sorted.bam -o $filenm_root.forIndelRealigner.intervals";
print ("Create RealignTargetCreator intervals with $filenm_root.Split.bam...\n");
system("$realigner_target");

#IndelRealigner ######Add --fix_misencoded_quality_scores option if needed
my $indelrealigner = "$java_path/java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $filenm_root.RG.sorted.bam -targetIntervals $filenm_root.forIndelRealigner.intervals -o $filenm_root.indelrealigned.RG.sorted.bam";
print ("IndelRealigning on file $filenm_root.Split.bamq...\n");
system("$indelrealigner");


#=================================================================================================================================


#=================================================================================================================================
#
##### First HaplotypeCaller to calibrate the bam.
#
#=================================================================================================================================


#HaplotypeCaller GATK4 (pour BaseRecalibrator)
my $haplotype_caller = "$gatk_path/gatk --java-options -Xmx1G HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -ploidy 2 -O $filenm_root.indelrealigned.RG.sorted.HC.vcf";
print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$haplotype_caller");


###Obtain a gVCF with HaplotypeCaller
####my $haplotype_caller_gvcf = "$gatk_path/gatk --java-options -Xmx4G HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -O $filenm_root.indelrealigned.RG.sorted.HC.gvcf";
#####print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.bam and output a gvcf (--emitRefConfidence GVCF)...\n");
#####system("$haplotype_caller_gvcf");

#VariantAnnotator (add annotation: MappingQualityZero)
my $variant_annotation_First = "$gatk_path/gatk --java-options -Xmx1G VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.bam -A MappingQualityZero -O $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf";
print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.HC.vcf...\n");
system("$variant_annotation_First");

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
my $variant_filtration_First  = "$gatk_path/gatk --java-options -Xmx1G VariantFiltration -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' --mask-extension 0 --filter-name HARD_TO_VALIDATE --filter-name QD_FILTER --filter-name DP_FILTER --mask-name Mask -O $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf";
print ("Variant filtration --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --filterExpression 'DP < 15' on $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf...\n");
system("$variant_filtration_First");

#SelectVariant to filter out filtered variants
my $select_variant_First = "$gatk_path/gatk --java-options -Xmx1G SelectVariants -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -select-type SNP -O $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf";
print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf...\n");
system("$select_variant_First");

#BaseRecalibrator
my $base_recalibrator = "$gatk_path/gatk --java-options -Xmx1G BaseRecalibrator -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam --known-sites $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf -O $filenm_root.recal_data.table";
print ("Create a recalibration table file with BaseRecalibrator with knownSites in file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf with BAM file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$base_recalibrator");

#Apply recalibration with ApplyBQSR (exPrintreads)
my $apply_recalibration = "$gatk_path/gatk --java-options -Xmx1G ApplyBQSR -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -bqsr $filenm_root.recal_data.table -O $filenm_root.indelrealigned.RG.sorted.recalibrated.bam";
print ("Apply recalibration with PrintReads on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$apply_recalibration");

#=================================================================================================================================
#
##### Second HaplotypeCaller for a gvcf on recalibrated Bam. gVCF Format
#
#=================================================================================================================================


##### Second HaplotypeCaller on recalibrated Bam.
#HaplotypeCaller GATK4 
#Obtain a gVCF with HaplotypeCaller
#my $haplotype_caller_gvcf = "$gatk_path/gatk --java-options -Xmx4G HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -ploidy 2 --emit-ref-confidence GVCF -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.gvcf";
#print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.recalibrated.bam and output a gvcf (--emitRefConfidence GVCF)...\n");
#system("$haplotype_caller_gvcf");


#========================
#VariantAnnotator GATK4 BETA (add annotation: MappingQualityZero)
#my $variant_annotation = "$gatk_path/gatk --java-options -Xmx4G VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.gvcf -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -A MappingQualityZero -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.gvcf";
#print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf...\n");
#system("$variant_annotation");

#VariantAnnotator GATK3 (add annotation: MappingQualityZero)
#my $variant_annotation = "$java_path/java -Xmx4g -jar $gatk3_path/GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.gvcf -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -A MappingQualityZero -o $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.gvcf";
#print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.gvcf...\n");
#system("$variant_annotation");
#========================

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
#my $variant_filtration  = "$gatk_path/gatk --java-options -Xmx4G VariantFiltration -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.gvcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' --mask-extension 0 --filter-name HARD_TO_VALIDATE --filter-name QD_FILTER --filter-name DP_FILTER --mask-name Mask -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.gvcf";
#print ("Variant filtration --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' on $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.gvcf...\n");
#system("$variant_filtration");

#SelectVariant to filter out filtered variants
#my $select_variant = "$gatk_path/gatk --java-options -Xmx4G SelectVariants -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.gvcf -select 'vc.isNotFiltered()' -select-type SNP -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.filtered.gvcf";
#print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.gvcf...\n");
#system("$select_variant");

#=================================================================================================================================


#=================================================================================================================================
#
##### Second HaplotypeCaller for a vcf on recalibrated Bam.  VCF Format
#
#=================================================================================================================================

#HaplotypeCaller GATK4 
#Obtain a gVCF with HaplotypeCaller
my $haplotype_caller_gvcf = "$gatk_path/gatk --java-options -Xmx1G HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -ploidy 2 -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf";
print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.recalibrated.bam ...\n");
system("$haplotype_caller_gvcf");


#========================
#VariantAnnotator GATK4 BETA (add annotation: MappingQualityZero)
my $variant_annotation_Second = "$gatk_path/gatk --java-options -Xmx1G VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -A MappingQualityZero -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf";
print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf...\n");
system("$variant_annotation_Second");

#VariantAnnotator GATK3 (add annotation: MappingQualityZero)
#my $variant_annotation = "$java_path/java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -A MappingQualityZero -o $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf";
#print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf...\n");
#system("$variant_annotation");
#========================

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
my $variant_filtration_Second  = "$gatk_path/gatk --java-options -Xmx1G VariantFiltration -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' --mask-extension 0 --filter-name HARD_TO_VALIDATE --filter-name QD_FILTER --filter-name DP_FILTER --mask-name Mask -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf";
print ("Variant filtration --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' on $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf...\n");
system("$variant_filtration_Second");

#SelectVariant to filter out filtered variants
my $select_variant_Second = "$gatk_path/gatk --java-options -Xmx1G SelectVariants -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -select-type SNP -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.filtered.vcf";
print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -select-type SNP on file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf...\n");
system("$select_variant_Second");

#=================================================================================================================================





print "Finished\n";

exit 1;
