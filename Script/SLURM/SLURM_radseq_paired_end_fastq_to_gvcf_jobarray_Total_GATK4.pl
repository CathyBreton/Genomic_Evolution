#!/usr/bin/perl -w






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
my $fastp_path    = "/nfs/work/agap_id-bin/img/fastp/0.20.1/fastp.0.20.1.img";
my $fastqc_path   = "/nfs/work/agap_id-bin/img/FastQC/0.11.7";
my $java_path     = "/nfs/work/agap_id-bin/img/java/jre1.8.0_31/bin";
my $bwa_path      = "/nfs/work/agap_id-bin/img/bwa/0.7.17/bwa.0.7.17.img";
##my $bwa_meme2_path="/nfs/work/agap_id-bin/img/bwa-mem2/2.0" 
my $gatk3_path = "/nfs/work/agap_id-bin/img/GATK/3.6";
my $gatk_path = "/nfs/work/agap_id-bin/img/GATK/4.2.0.0";
my $picard_path_2_24 = "/home/bretonc/work_agap_id-bin/img/picard-tools/2.24.0/picard-tools.2.24.0.img";
my $picard_path_2_72 = "/home/bretonc/work_agap_id-bin/img/picard/2.7.2";

my $samtools_path = "/nfs/work/agap_id-bin/img/samtools/1.2/samtools.1.2.img";
my $sambamba_path = "/nfs/work/agap_id-bin/img/sambamba/0.8.0/sambamba.0.8.0.img";



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



#========================================================================================================================================================

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#Programs 

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#=======================================================================================================================================================

print "#---------------------------------------------------------------------------------#\n";
print "Check homogenieity of fastq files forward and reverse : Cutadapt ok avec cutadapt  \n";
print "#---------------------------------------------------------------------------------#\n";

#Check homogenieity of fastq files forward and reverse : Cutadapt


#Check homogenieity (Cutadapt)
if (!-e "$directory/$forward_without_ext.cutadapt.fastq" && "$directory/$reverse_without_ext.cutadapt.fastq")
{
        #my $cutadapt_cmd = "module load cutadapt/3.1;cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -O 7 -m 30 -q 20,20 -o $forward_without_ext.cutadapt.fastq -p $reverse_without_ext.cutadapt.fastq $directory/$forward $directory/$reverse";
        #my $cutadapt_cmd = "cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -O 7 -m 30 -q 20,20 -o $forward_without_ext.cutadapt.fastq -p $reverse_without_ext.cutadapt.fastq $directory/$forward $directory/$reverse";
        my $cutadapt_cmd = "singularity exec $cutadapt_path cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -O 7 -m 30 -q 20,20 -o $forward_without_ext.cutadapt.fastq -p $reverse_without_ext.cutadapt.fastq $directory/$forward $directory/$reverse";
        #print ("cutadapt -a AGATCGGAAGAGC -O 7 -m 30 -q 20,20 -o $forward_without_ext.cutadapt.fastq -p $reverse_without_ext.cutadapt.fastq $forward $reverse\n");
        print "$cutadapt_cmd\n"; #+debug

        if (system("$cutadapt_cmd"))#
        {
                # system() returned an error code, execution failed.
                warn "Failed to execute: $!\n";
        } 
}
else
{
    print("$directory/$forward_without_ext.cutadapt.fastq and $directory/$reverse_without_ext.cutadapt.fastq exists! the file will not be overwritten : skip combine step \n\n");
}




print "#----------------------------------------------------------------#\n";
print "Check homogenieity of fastq files forward and reverse : fastp OK  avec singularity    \n";
print "#----------------------------------------------------------------#\n";


if (!-e "$directory/$forward_without_ext.fastp.fastq.gz" && "$directory/$reverse_without_ext.fastp.fastq.gz")
{
        my $fastp_cmd = "singularity exec $fastp_path fastp --detect_adapter_for_pe"
        #my $fastp_cmd = "fastp --detect_adapter_for_pe"
                . ' ' . "--overrepresentation_analysis"
                . ' ' . "--correction --cut_right --thread 2"
                . ' ' . "--html $directory/anc.fastp.html --json $directory/anc.fastp.json"
                . ' ' . "-i $directory/$forward -I $directory/$reverse"
                . ' ' . "-o $directory/$forward_without_ext.fastp.fastq.gz -O $directory/$reverse_without_ext.fastp.fastq.gz"
                . ' ' . "\n"; 
        print "$fastp_cmd\n"; #+debug
        if (system("$fastp_cmd"))#
                {
                        # system() returned an error code, execution failed.
                       warn "Failed to execute: $!\n";
                } 
}
else
{
    print("$directory/$forward_without_ext.fastp.fastq.gz and $directory/$reverse_without_ext.fastp.fastq.gz exists! the file will not be overwritten : skip combine step \n\n");
}


print "#----------------------------------------------------------------#\n";
print "Check homogenieity of fastq files forward and reverse : fastqc  ok module   \n";
print "#----------------------------------------------------------------#\n";

if (!-e "$directory/$forward_without_ext.cutadapt.fastq" && "$directory/$reverse_without_ext.cutadapt.fastq")
{
        foreach my $file ( ("$directory/$forward_without_ext.cutadapt.fastq", "$directory/$reverse_without_ext.cutadapt.fastq") ) {
	        #fastqc
                my $fastqc_cmd = "fastqc -j java $file";
                #my $fastqc_cmd = "$fastqc_path/fastqc -j $java_path/java $file";
                print ("fastqc on $file (two of two)\n");
                if (system("$fastqc_cmd"))#
                {
                        # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }     
        }
}       
else
{
    print("$directory/$forward_without_ext.cutadapt.fastq and $directory/$reverse_without_ext.cutadapt.fastq exists! the file will not be overwritten : skip combine step \n\n");
}



print "#---------------------------------------------------#\n";
print "bwa mapping on the reference choosen in argument -r Ok with singularity \n";
print "#---------------------------------------------------#\n";

# bwa mapping on the reference choosen in argument -r

if (!-e "$directory/$filenm_root.sam")
{
        my $bwa_mapping_cmd = "singularity exec $bwa_path bwa mem $reference $forward_without_ext.cutadapt.fastq $reverse_without_ext.cutadapt.fastq > $filenm_root.sam";
        #my $bwa_mapping_cmd = "bwa mem $reference $forward_without_ext.cutadapt.fastq $reverse_without_ext.cutadapt.fastq > $filenm_root.sam";
        print ("Mapping with BWA on file $forward_without_ext.cutadapt.fastq and $reverse_without_ext.cutadapt.fastq...\n");
        if (system("$bwa_mapping_cmd"))#
                {
                        # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }
}
else
{
    print("$directory/$filenm_root.sam exists! the file will not be overwritten : skip bwa mapping step \n\n");
}



# bwa mapping on the reference choosen in argument -r
#my $bwa_mapping = "$bwa_path/bwa-mem2 mem $reference $forward_without_ext.cutadapt.fastq $reverse_without_ext.cutadapt.fastq > $filenm_root.sam";
#print ("Mapping with BWA on file $forward_without_ext.cutadapt.fastq and $reverse_without_ext.cutadapt.fastq...\n");
#system("$bwa_mapping");


print "#---------------------------------------------------------#\n";
print "Add read group, and sort bam , with the tools picardtools ok java jar  \n";
print "#---------------------------------------------------------#\n";


#Add read group, and sort bam , with the tools picardtools
my $condition = $filenm_root;

if (!-e "$directory/$filenm_root.RG.sorted.bam")
{
        #my $add_rg_cmd = "singularity exec $picard_path_2_24 picard AddOrReplaceReadGroups I=$filenm_root.sam O=$filenm_root.RG.sorted.bam SO=coordinate RGLB=$cultivar RGPL=illumina RGPU=run RGSM=$condition RGID=$condition";
        my $add_rg_cmd = "java -Xmx1g -jar $picard_path_2_72/picard.jar AddOrReplaceReadGroups I=$filenm_root.sam O=$filenm_root.RG.sorted.bam SO=coordinate RGLB=$cultivar RGPL=illumina RGPU=run RGSM=$condition RGID=$condition";
        print ("Adding Read Group and sorting on file $filenm_root.sam...\n");
        if (system("$add_rg_cmd"))#
                {
                        # system() returned an error code, execution failed.
                       warn "Failed to execute: $!\n";
                }
}
else
{
    print("$directory/$filenm_root.RG.sorted.bam exists! the file will not be overwritten : skip Add read group step \n\n");
}




print "#-------------------------------------------------------------------------#\n";
print "Index bam, with the tools samtools       ok avec singularity               \n";
print "#-------------------------------------------------------------------------#\n";


#Index bam, with the tools picardtools
if (!-e "$directory/$filenm_root.RG.sorted.bai")
{
        my $index_bam_cmd = "singularity exec $samtools_path samtools index $filenm_root.RG.sorted.bam";
        #my $index_bam_cmd = "samtools index $filenm_root.RG.sorted.bam";
        print ("Indexing on bam file $filenm_root.RG.sorted.bam...\n");
        if (system("$index_bam_cmd"))#
                {
                        # system() returned an error code, execution failed.
                       warn "Failed to execute: $!\n";
                }
}
else
{
    print("$directory/$filenm_root.RG.sorted.bai exists! the file will not be overwritten : skip bwa mapping step \n\n");
}


#SplitNCigarReads with GATK4 ######## add -U ALLOW_N_CIGAR_READS in case... becarefull not in GATK4
####my $split_and_trim = "$gatk_path/gatk --java-options -Xmx4G SplitNCigarReads -R $reference -I $filenm_root.RG.sorted.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -O $filenm_root.Split.bam";
####print ("Split N CigarReads on file $filenm_root.RG.sorted.bam...\n");
####system("$split_and_trim");


#=================================================================================================================================

print "#---------------------------------------------------------#\n";
print "Indel Realignement with gatk 3  \n";
print "#---------------------------------------------------------#\n";




##### Indel Realignement, 
# This step is not present in GATK4, not necessary to do for SNP calling but for Indel , need to do so this step is done with GATK3

#RealignerTargetCreator ######Add --fix_misencoded_quality_scores option if needed --fix_misencoded_quality_scores 
if (!-e "$directory/$filenm_root.forIndelRealigner.intervals")
{    
        #my $realigner_target_cmd = "$java_path/java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $filenm_root.RG.sorted.bam -o $filenm_root.forIndelRealigner.intervals";
        my $realigner_target_cmd = "java -Xmx1g -jar $gatk3_path/GenomeAnalysisTK.jar -T RealignerTargetCreator --fix_misencoded_quality_scores -R $reference -I $filenm_root.RG.sorted.bam -o $filenm_root.forIndelRealigner.intervals";
        print ("Create RealignTargetCreator intervals with $filenm_root.Split.bam...\n");
        if (system("$realigner_target_cmd"))#
                {
                        # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }
}
else
{
    print("$directory/$filenm_root.forIndelRealigner.intervals exists! the file will not be overwritten : skip bwa mapping step \n\n");
}



#IndelRealigner ######Add --fix_misencoded_quality_scores option if needed
if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.bam")
{    
        my $indelrealigner_cmd = "java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T IndelRealigner --fix_misencoded_quality_scores -R $reference -I $filenm_root.RG.sorted.bam -targetIntervals $filenm_root.forIndelRealigner.intervals -o $filenm_root.indelrealigned.RG.sorted.bam";
        #my $indelrealigner_cmd = "$java_path/java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $filenm_root.RG.sorted.bam -targetIntervals $filenm_root.forIndelRealigner.intervals -o $filenm_root.indelrealigned.RG.sorted.bam";
        print ("IndelRealigning on file $filenm_root.Split.bamq...\n");
        if (system("$indelrealigner_cmd"))#
                {
                        # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }        
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.bam exists! the file will not be overwritten : skip IndelRealigner step \n\n");
}

#=================================================================================================================================


#=================================================================================================================================
#
##### First HaplotypeCaller to calibrate the bam.
#
#=================================================================================================================================


print "#---------------------------------------------------------#\n";
print "HaplotypeCaller GATK4                                      \n";
print "#---------------------------------------------------------#\n";


#HaplotypeCaller GATK4 (pour BaseRecalibrator)
if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.HC.vcf")
{    
        my $haplotype_caller_cmd = "$gatk_path/gatk --java-options -Xmx1G HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -ploidy 2 -O $filenm_root.indelrealigned.RG.sorted.HC.vcf";
        print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
        if (system("$haplotype_caller_cmd"))#
                {
                       # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.HC.vcf exists! the file will not be overwritten : skip HaplotypeCaller step \n\n");
}


print "#---------------------------------------------------------#\n";
print "VariantAnnotator GATK4                                     \n";
print "#---------------------------------------------------------#\n";

#VariantAnnotator (add annotation: MappingQualityZero)
if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.HC.ann.vcf")
{    
        my $variant_annotation_First_cmd = "$gatk_path/gatk --java-options -Xmx1G VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.bam -A MappingQualityZero -O $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf";
        print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.HC.vcf...\n");
        if (system("$variant_annotation_First_cmd"))#
                {
                       # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }        
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.HC.ann.vcf exists! the file will not be overwritten : skip VariantAnnotator step \n\n");
}



print "#---------------------------------------------------------#\n";
print "VariantFiltration GATK4                                    \n";
print "#---------------------------------------------------------#\n";

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf")
{    
        my $variant_filtration_First_cmd  = "$gatk_path/gatk --java-options -Xmx1G VariantFiltration -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' --mask-extension 0 --filter-name HARD_TO_VALIDATE --filter-name QD_FILTER --filter-name DP_FILTER --mask-name Mask -O $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf";
        print ("Variant filtration --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --filterExpression 'DP < 15' on $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf...\n");
        if (system("$variant_filtration_First_cmd"))#
                {
                        # system() returned an error code, execution failed.
                       warn "Failed to execute: $!\n";
                }        
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf exists! the file will not be overwritten : skip VariantFiltration step \n\n");
}



print "#---------------------------------------------------------#\n";
print "SelectVariant GATK4                                        \n";
print "#---------------------------------------------------------#\n";


#SelectVariant to filter out filtered variants
if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf")
{    
        my $select_variant_First_cmd = "$gatk_path/gatk --java-options -Xmx1G SelectVariants -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -select-type SNP -O $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf";
        print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf...\n");
        if (system("$select_variant_First_cmd"))#
                {
                        # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }        
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf exists! the file will not be overwritten : skip SelectVariant step \n\n");
}





print "#---------------------------------------------------------#\n";
print "BaseRecalibrator GATK4                                        \n";
print "#---------------------------------------------------------#\n";

#BaseRecalibrator

if (!-e "$directory/$filenm_root.recal_data.table")
{    
        my $base_recalibrator_cmd = "$gatk_path/gatk --java-options -Xmx1G BaseRecalibrator -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam --known-sites $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf -O $filenm_root.recal_data.table";
        print ("Create a recalibration table file with BaseRecalibrator with knownSites in file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf with BAM file $filenm_root.indelrealigned.RG.sorted.bam...\n");
        if (system("$base_recalibrator_cmd"))#
                {
                        # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }        
}
else
{
    print("$directory/$filenm_root.recal_data.table exists! the file will not be overwritten : skip BaseRecalibrator step \n\n");
}



print "#---------------------------------------------------------#\n";
print "Apply recalibration with ApplyBQSR GATK4                   \n";
print "#---------------------------------------------------------#\n";


#Apply recalibration with ApplyBQSR (exPrintreads)

if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.bam")
{    
        my $apply_recalibration_cmd = "$gatk_path/gatk --java-options -Xmx1G ApplyBQSR -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -bqsr $filenm_root.recal_data.table -O $filenm_root.indelrealigned.RG.sorted.recalibrated.bam";
        print ("Apply recalibration with PrintReads on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
        if (system("$apply_recalibration_cmd"))#
                {
                        # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }        
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.bam exists! the file will not be overwritten : skip recalibration step \n\n");
}



#=================================================================================================================================
#
##### Second HaplotypeCaller for a gvcf on recalibrated Bam. gVCF Format
#
#=================================================================================================================================

print "#---------------------------------------------------------#\n";
print "Second HaplotypeCaller on recalibrated Bam GATK4   Gvcf    \n";
print "#---------------------------------------------------------#\n"; 


##### Second HaplotypeCaller on recalibrated Bam.
#HaplotypeCaller GATK4 
#Obtain a gVCF with HaplotypeCaller

if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.gvcf")
{    
        my $haplotype_caller_gvcf_cmd = "$gatk_path/gatk --java-options -Xmx4G HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -ploidy 2 --emit-ref-confidence GVCF -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.gvcf";
        print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.recalibrated.bam and output a gvcf (--emitRefConfidence GVCF)...\n");
        if (system("$haplotype_caller_gvcf_cmd"))#
                {
                        # system() returned an error code, execution failed.
                        warn "Failed to execute: $!\n";
                }        
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.gvcf exists! the file will not be overwritten : skip Second HaplotypeCaller step \n\n");
}

#=================================================================================================================================


#=================================================================================================================================
#
##### Second HaplotypeCaller for a vcf on recalibrated Bam.  VCF Format
#
#=================================================================================================================================

#HaplotypeCaller GATK4 
#Obtain a VCF with HaplotypeCaller

print "#-------------------------------------------------------------#\n";
print "Second HaplotypeCaller on recalibrated Bam GATK4  VCF          \n";
print "#-------------------------------------------------------------#\n";  


if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf")
{    
        #my $haplotype_caller_vcf_cmd = "$gatk_path/gatk --java-options -Xmx1G HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -ploidy 2 -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf";
        #print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.recalibrated.bam ...\n");
        #if (system("$haplotype_caller_vcf_cmd"))#
          #       {
                        # system() returned an error code, execution failed.
            #            warn "Failed to execute: $!\n";
          #       }       
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf exists! the file will not be overwritten : skip HaplotypeCaller step \n\n");
}

#========================

#VariantAnnotator GATK4 BETA (add annotation: MappingQualityZero)

print "#-------------------------------------------------------------#\n";
print "Second VariantAnnotator on recalibrated Bam GATK4  VCF         \n";
print "#-------------------------------------------------------------#\n";  


if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf")
{    
        #my $variant_annotation_Second_cmd = "$gatk_path/gatk --java-options -Xmx1G VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -A MappingQualityZero -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf";
        #print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf...\n");
        #if (system("$variant_annotation_Second_cmd"))#
          #       {
                        # system() returned an error code, execution failed.
          #              warn "Failed to execute: $!\n";
          ##       }       
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf exists! the file will not be overwritten : skip VariantAnnotator step \n\n");
}

#VariantAnnotator GATK3 (add annotation: MappingQualityZero)

print "#-------------------------------------------------------------#\n";
print "Second VariantAnnotator on recalibrated Bam GATK4  VCF         \n";
print "#-------------------------------------------------------------#\n";  



if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf")
{    
        #my $variant_annotation_cmd = "$java_path/java -Xmx1G -jar $gatk3_path/GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam -A MappingQualityZero -o $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf";
        #print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.vcf...\n");
        #if (system("$variant_annotation_cmd"))#
            #     {
                        # system() returned an error code, execution failed.
             #           warn "Failed to execute: $!\n";
            #     }        
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf exists! the file will not be overwritten : skip bwa mapping step \n\n");
}
#========================

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)

if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf")
{    
#my $variant_filtration_Second_cmd  = "$gatk_path/gatk --java-options -Xmx1G VariantFiltration -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' --mask-extension 0 --filter-name HARD_TO_VALIDATE --filter-name QD_FILTER --filter-name DP_FILTER --mask-name Mask -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf";
#print ("Variant filtration --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-expression 'QD < 1.5' --filter-expression 'DP < 15' on $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.vcf...\n");
#if (system("$variant_filtration_Second_cmd"))#
  #       {
                # system() returned an error code, execution failed.
    #            warn "Failed to execute: $!\n";
   #      }        
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf exists! the file will not be overwritten : skip VariantFiltration step \n\n");
}

#SelectVariant to filter out filtered variants

if (!-e "$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.filtered.vcf")
{    
#my $select_variant_Second_cmd = "$gatk_path/gatk --java-options -Xmx1G SelectVariants -R $reference -V $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -select-type SNP -O $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.filtered.vcf";
#print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -select-type SNP on file $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.vcf...\n");
#if (system("$select_variant_Second_cmd"))#
  #       {
                # system() returned an error code, execution failed.
   #             warn "Failed to execute: $!\n";
  #       }        
}
else
{
    print("$directory/$filenm_root.indelrealigned.RG.sorted.recalibrated.HC.ann.VF.filtered.vcf exists! the file will not be overwritten : skip bwa mapping step \n\n");
}

#=================================================================================================================================





print "Finished\n";

exit 1;
