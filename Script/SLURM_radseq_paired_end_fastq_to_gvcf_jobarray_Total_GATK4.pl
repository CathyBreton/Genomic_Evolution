#! /usr/bin/perl -w






###################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
#
# Intellectual property belongs to Bioversity
#
# Written by  Catherine Breton
#
# Script Name : SLURM_radseq_paired_end_fastq_to_gvcf_jobarray_Total_GATK4.pl
#
#################### 


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;
use Data::Dumper;
use Carp qw (cluck confess croak);
use Pod::Usage;




my $cutadapt_path = "/nfs/work/agap_id-bin/img/cutadapt/3.1/cutadapt.3.1.img";
############### MAIN ##############

print "Job starting\n";


if (@ARGV < 6)
{
        print "\n\nUsage: perl $0 -r reference_file_path -x extension_file_to_treat -d directory -c cultivar\n\n\n";
        exit 1;
}




#Variables
my ($help, $man, $resume, $reference, $extension, $directory, $cultivar);

#Getoptions
GetOptions 
(   "x|extension=s"         => \$extension, 
    "r|reference=s"         => \$reference, 
    "d|directorty=s"        => \$directory,
    "c|cultivar=s" 			=> \$cultivar
) or die ("Error in command line arguments\n");

#Variables Environnement
foreach (sort keys %ENV) { 
  print "$_  =  $ENV{$_}\n"; 
}




#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");  
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


my %hash_file;
my $i=0;
my $j=1;
foreach my $file (sort( {$a cmp $b} @files)) {
	if ($j%2 != 0) {
		$hash_file{$i}=$file;
		} else {
		$hash_file{$i}.= "-" . $file;	
		$i++;
		}
	$j++;

}
print (Dumper([%hash_file]) . "\n");	

#relative file path 

#relative file path 
#my $slurm_id = `echo $SLURM_ARRAY_TASK_ID` or die ("Cannot catch SLURM_ARRAY_TASK_ID\n");
my $slurm_id = $ENV{'SLURM_ARRAY_TASK_ID'};

chomp($slurm_id);
print(" Debugg : " . $slurm_id ."\n");
#$slurm_id--; #perl is in basis 0

my $files_to_treat = $hash_file{$slurm_id};


print($files_to_treat);
my @files_to_treat_split = split(m/-/,$files_to_treat);
my $forward = $files_to_treat_split[0];
print "$forward\n";
my $forward_without_ext = fileparse($forward, qr/.$extension$/);
print "$forward_without_ext\n";
my $reverse = $files_to_treat_split[1];
print "$reverse\n";
my $reverse_without_ext = fileparse($reverse, qr/.$extension$/);
print "$reverse_without_ext\n";

#Filename without extension
my $filenm_root = fileparse($forward, qr/_\d\.$extension$/);
print ("$filenm_root\n");



#=================================================================================================================================

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#Programs 

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#=================================================================================================================================



#Check homogenieity of fastq files forward and reverse : Cutadapt



my $cutadapt_cmd = "singularity exec $cutadapt_path cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -O 7 -m 30 -q 20,20 -o $forward_without_ext.cutadapt.fastq -p $reverse_without_ext.cutadapt.fastq $directory/$forward $directory/$reverse";
#print ("cutadapt -a AGATCGGAAGAGC -O 7 -m 30 -q 20,20 -o $forward_without_ext.cutadapt.fastq -p $reverse_without_ext.cutadapt.fastq $forward $reverse\n");
print "$cutadapt_cmd\n"; #+debug

if (system("$cutadapt_cmd"))
         {
                # system() returned an error code, execution failed.
                warn "Failed to execute: $!\n";
         } 


#foreach my $file ( ("$forward_without_ext.cutadapt.fastq", "$reverse_without_ext.cutadapt.fastq") ) {

	#fastqc
        #my $fastqc = "$fastqc_path/fastqc -j $java_path/java $file";
        #print ("fastqc on $file (two of two)\n");
        #system("$fastqc");
#}

# bwa mapping on the reference choosen in argument -r
#my $bwa_mapping = "$bwa_path/bwa mem $reference $forward_without_ext.cutadapt.fastq $reverse_without_ext.cutadapt.fastq > $filenm_root.sam";
#print ("Mapping with BWA on file $forward_without_ext.cutadapt.fastq and $reverse_without_ext.cutadapt.fastq...\n");
#system("$bwa_mapping");

#Add read group, and sort bam , with the tools picardtools
#my $condition = $filenm_root;
#my $add_rg = "$java_path/java -Xmx1g -jar $picard_path/picard.jar AddOrReplaceReadGroups I=$filenm_root.sam O=$filenm_root.RG.sorted.bam SO=coordinate RGLB=$cultivar RGPL=illumina RGPU=run RGSM=$condition RGID=$condition";
#print ("Adding Read Group and sorting on file $filenm_root.sam...\n");
#system("$add_rg");

#Index bam, with the tools picardtools
#my $index_bam = "$samtools_path/samtools index $filenm_root.RG.sorted.bam";
#print ("Indexing on bam file $filenm_root.RG.sorted.bam...\n");
#system("$index_bam");

#SplitNCigarReads with GATK4 ######## add -U ALLOW_N_CIGAR_READS in case... becarefull not in GATK4
####my $split_and_trim = "$gatk_path/gatk --java-options -Xmx4G SplitNCigarReads -R $reference -I $filenm_root.RG.sorted.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -O $filenm_root.Split.bam";
####print ("Split N CigarReads on file $filenm_root.RG.sorted.bam...\n");
####system("$split_and_trim");


#=================================================================================================================================

##### Indel Realignement, 
# This step is not present in GATK4, not necessary to do for SNP calling but for Indel , need to do so this step is done with GATK3


#RealignerTargetCreator ######Add --fix_misencoded_quality_scores option if needed --fix_misencoded_quality_scores 
#my $realigner_target = "$java_path/java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $filenm_root.RG.sorted.bam -o $filenm_root.forIndelRealigner.intervals";
#print ("Create RealignTargetCreator intervals with $filenm_root.Split.bam...\n");
#system("$realigner_target");

#IndelRealigner ######Add --fix_misencoded_quality_scores option if needed
#my $indelrealigner = "$java_path/java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $filenm_root.RG.sorted.bam -targetIntervals $filenm_root.forIndelRealigner.intervals -o $filenm_root.indelrealigned.RG.sorted.bam";
#print ("IndelRealigning on file $filenm_root.Split.bamq...\n");
#system("$indelrealigner");


#=================================================================================================================================


#=================================================================================================================================
#
##### First HaplotypeCaller to calibrate the bam.
#
#=================================================================================================================================


#HaplotypeCaller GATK4 (pour BaseRecalibrator)
#my $haplotype_caller = "$gatk_path/gatk --java-options -Xmx1G HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -ploidy 2 -O $filenm_root.indelrealigned.RG.sorted.HC.vcf";
#print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
#system("$haplotype_caller");


#VariantAnnotator (add annotation: MappingQualityZero)
#my $variant_annotation_First = "$gatk_path/gatk --java-options -Xmx1G VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.bam -A MappingQualityZero -O $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf";
#print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.HC.vcf...\n");
#system("$variant_annotation_First");

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
#my $variant_filtration_First  = "$gatk_path/gatk --java-options -Xmx1G VariantFiltration -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' --mask-extension 0 --filter-name HARD_TO_VALIDATE --filter-name QD_FILTER --filter-name DP_FILTER --mask-name Mask -O $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf";
#print ("Variant filtration --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --filterExpression 'DP < 15' on $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf...\n");
#system("$variant_filtration_First");

#SelectVariant to filter out filtered variants
#my $select_variant_First = "$gatk_path/gatk --java-options -Xmx1G SelectVariants -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -select-type SNP -O $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf";
#print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf...\n");
#system("$select_variant_First");

#BaseRecalibrator
#my $base_recalibrator = "$gatk_path/gatk --java-options -Xmx1G BaseRecalibrator -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam --known-sites $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf -O $filenm_root.recal_data.table";
#print ("Create a recalibration table file with BaseRecalibrator with knownSites in file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf with BAM file $filenm_root.indelrealigned.RG.sorted.bam...\n");
#system("$base_recalibrator");

#Apply recalibration with ApplyBQSR (exPrintreads)
#my $apply_recalibration = "$gatk_path/gatk --java-options -Xmx1G ApplyBQSR -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -bqsr $filenm_root.recal_data.table -O $filenm_root.indelrealigned.RG.sorted.recalibrated.bam";
#print ("Apply recalibration with PrintReads on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
#system("$apply_recalibration");

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


#=================================================================================================================================


#=================================================================================================================================
#
##### Second HaplotypeCaller for a vcf on recalibrated Bam.  VCF Format
#
#=================================================================================================================================

#HaplotypeCaller GATK4 
#Obtain a VCF with HaplotypeCaller
#my $haplotype_caller_vcf = "$gatk_path/gatk --java-options -Xmx1G HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -ploidy 2 -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf";
#print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.recalibrated.bam ...\n");
#system("$haplotype_caller_vcf");


#========================
#VariantAnnotator GATK4 BETA (add annotation: MappingQualityZero)
#my $variant_annotation_Second = "$gatk_path/gatk --java-options -Xmx1G VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -A MappingQualityZero -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf";
#print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf...\n");
#system("$variant_annotation_Second");

#VariantAnnotator GATK3 (add annotation: MappingQualityZero)
#my $variant_annotation = "$java_path/java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -A MappingQualityZero -o $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf";
#print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf...\n");
#system("$variant_annotation");
#========================

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
#my $variant_filtration_Second  = "$gatk_path/gatk --java-options -Xmx1G VariantFiltration -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' --mask-extension 0 --filter-name HARD_TO_VALIDATE --filter-name QD_FILTER --filter-name DP_FILTER --mask-name Mask -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf";
#print ("Variant filtration --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' on $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf...\n");
#system("$variant_filtration_Second");

#SelectVariant to filter out filtered variants
#my $select_variant_Second = "$gatk_path/gatk --java-options -Xmx1G SelectVariants -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -select-type SNP -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.filtered.vcf";
#print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -select-type SNP on file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf...\n");
#system("$select_variant_Second");

#=================================================================================================================================





print "Finished\n";

exit 1;
