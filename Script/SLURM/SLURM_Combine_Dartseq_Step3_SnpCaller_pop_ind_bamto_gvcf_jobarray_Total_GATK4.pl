#! /usr/bin/perl -w





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
#  Script Name : Dartseq_single_end_fastq_to_bam_Total_GATK4.pl
#
#
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
use Data::Dumper;
use Carp qw (cluck confess croak);
use Pod::Usage;




my $cutadapt_path    = "/nfs/work/agap_id-bin/img/cutadapt/3.1/cutadapt.3.1.img";
my $fastp_path       = "/nfs/work/agap_id-bin/img/fastp/0.20.1/fastp.0.20.1.img";
my $fastqc_path      = "/nfs/work/agap_id-bin/img/FastQC/0.11.7";
my $java_path        = "/nfs/work/agap_id-bin/img/java/jre1.8.0_31/bin";
my $bwa_path         = "/nfs/work/agap_id-bin/img/bwa/0.7.17/bwa.0.7.17.img";
##my $bwa_meme2_path="/nfs/work/agap_id-bin/img/bwa-mem2/2.0" 
my $gatk3_path       = "/nfs/work/agap_id-bin/img/GATK/3.6";
##my $gatk_path        = "/home/bretonc/work_agap_id-bin/img/GATK/4.2.0.0";
my $gatk_path        = "/home/bretonc/Software/GATK/4.2.0.0";
my $picard_path_2_24 = "/nfs/work/agap_id-bin/img/picard-tools/2.24.0/picard-tools.2.24.0.img";
my $picard_path_2_72 = "/nfs/work/agap_id-bin/img/picard/2.7.2";
my $samtools_path    = "/nfs/work/agap_id-bin/img/samtools/1.2/samtools.1.2.img";
my $sambamba_path    = "/nfs/work/agap_id-bin/img/sambamba/0.8.0/sambamba.0.8.0.img";






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
    "c|cultivar=s" 	        => \$cultivar
) or die ("Error in command line arguments\n");




#List of all file in directory with specific extension
opendir(CWD,"$directory") or die("Cannot open current directory\n");  
my @files = grep { m/\.reheader.Merged.sorted\.$extension$/ } readdir(CWD);
closedir(CWD);

#List of all file in directory with specific extension
opendir(CWD,"$directory") or die("Cannot open current directory\n");  
my @names = grep { m/\_.+\.$extension$/ } readdir(CWD);
closedir(CWD);


#relative file path 
#my $slurm_id = `echo $SLURM_ARRAY_TASK_ID` or die ("Cannot catch SLURM_ARRAY_TASK_ID\n");
my $slurm_id = $ENV{'SLURM_ARRAY_TASK_ID'};

chomp($slurm_id);
print(" Debugg : " . $slurm_id ."\n");
$slurm_id--; #perl is in basis 0
my $input_bam = $files[$slurm_id]; 


#Filename without extension
my $filenm_root = fileparse($input_bam, qr/\.$extension$/);
#print ("$filenm_root\n");



#=================================================================================================================================

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#Programs 

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#=================================================================================================================================






#=================================================================================================================================
#
##### Step  Combining sorting and Indexing the replicated bam files
#
#=================================================================================================================================


   print "#============================================================#\n";
   


   
#=================================================================================================================================
#
##### Second HaplotypeCaller for a gvcf on recalibrated Bam. gVCF Format
#
#=================================================================================================================================


##### Second HaplotypeCaller on recalibrated Bam.
#HaplotypeCaller GATK4 
#Obtain a gVCF with HaplotypeCaller

   

#my $haplotype_caller_gvcf = "$gatk_path/gatk --java-options -Xmx1G HaplotypeCaller -R $reference -I $prefix.R1.indelrealigned.RG.sorted.recalibrated.reheader.Merged.bam -ploidy 2 --emit-ref-confidence GVCF -O  $prefix.R1.indelrealigned.RG.sorted.recalibrated.reheader.Merged.HC.gvcf";



   print "#--------------------------------#\n";
   print "Step 5 Haplotype Caller with GATK \n";  
   print "#--------------------------------#\n";


if (!-e "$directory/$filenm_root.HC.gvcf")
{    
    my $haplotype_caller_gvcf_cmd =
      "$gatk_path/gatk --java-options -Xmx10G HaplotypeCaller"
      . ' ' . "-R $reference"
      . ' ' . "-I $directory/$filenm_root.bam"
      . ' ' . "-ploidy 2"
      . ' ' . "--emit-ref-confidence GVCF"
      . ' ' . "-O $directory/$filenm_root.HC.gvcf\n";
 
   print ("SnpCalling with HaplotypeCaller of GATK on files $filenm_root.bam...\n");
   print "$haplotype_caller_gvcf_cmd\n"; #+debug
   if (system("$haplotype_caller_gvcf_cmd"))
     {
         # system() returned an error code, execution failed.
         warn "Failed to execute: $!\n";
     } 
}
else
{
    print("$directory/$filenm_root.HC.gvcf exists! the file will not be overwritten : skip bwa mapping step \n\n");
}
   

   print "#---------------------------------------#\n";
   print "Step 6 Haplotype Caller with OCTOPUS \n";  
   print "#---------------------------------------#\n";



#if (!-e "$directory/$filenm_root.sorted.OCTOPUS")
#{    
    #my $haplotype_caller_Octopus_cmd =
#      "$octopus_path/octopus"
#      . ' ' . "-R $reference"
#      . ' ' . "-I $directory/$filenm_root.bam"
#      . ' ' . "--threads 4"
#      . ' ' . "--emit-ref-confidence GVCF"
#      . ' ' . "-O $directory/$filenm_root.sorted.OCTOPUS\n";
 
#   print ("SnpCalling with HaplotypeCaller of GAOCTOPUSTK on files $filenm_root.bam...\n");
   #print "  $haplotype_caller_Octopus\n"; #+debug
   #if (system("$haplotype_caller_Octopus"))
    # {
         # system() returned an error code, execution failed.
    #     warn "Failed to execute: $!\n";
    # } 
#}
#else
#{
#    print("$directory/$filenm_root.sorted.OCTOPUS exists! the file will not be overwritten : skip bwa mapping step \n\n");
#}


    print "#============================================================#\n";




print "Finished\n";

exit 1;
