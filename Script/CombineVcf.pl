#! /usr/bin/perl -w
#$ -q bigmem.q
#$ -cwd
#$ -V
#$ -S /usr/bin/perl
#$ -pe parallel_fill 50
#$ -l mem_free=100G

#charge modules
#use "bioinfo/bwa/0.7.12";
#use "system/java/jre8";
#use "bioinfo/FastQC/0.11.3";
#use "compiler/gcc/4.9.2";
#use "bioinfo/bwa/0.7.12";

#my $gatk_path = "/usr/local/bioinfo/GenomeAnalysisTK/3.4-46";
#my $gatk_path = "/homedir/cbreton/Software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef";
my $gatk_path = "/usr/local/bioinfo/GenomeAnalysisTK/3.7-0";
my $java_path = "/usr/local/java/jre8/bin";
my $picard_path = "/usr/local/bioinfo/picard-tools/1.130";
my $samtools_path = "/usr/local/bioinfo/samtools/1.2/bin";
my $bcftools_path = "/usr/local/bioinfo/bcftools/1.9/bin";




###################
#
#  Licencied under GNU General Public License : GPL-3.0-or-later 
#
#  Intellectual property belongs to Bioversity
#
#  Written by Yann Hueber modified by Catherine Breton
#
#  Contact : Catherine Breton c.breton@cgiar.org
#
#  Script Name : CombineVcf.pl
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


############### MAIN ##############


if (@ARGV < 4)
{
        print "\n\nUsage: perl $0 -r reference_fasta -p output_prefix -x extension_file_to_treat\n\n\n";
        exit 1;
}


#Variables
my ($extension, $reference, $output_prefix);

#Getoptions
GetOptions
(       "x|extension=s" => \$extension, 
        "r|reference=s" => \$reference, 
        "p|outputprefix=s" => \$output_prefix) or die ("Error in command line arguments\n");



#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


#Concat vcf filenames
my $ivcfs = "";
foreach (@files) {
	$ivcfs .= " -V $_ ";
}

my $vcfgzs = "";
foreach (@files) {
	$vcfgzs .= " $_.gz ";
}



#-V:control,vcf my_control_samples.vcf


#Programs 



#MergeVCFs
foreach my $ivcfs (@files) {
	my $filenm_root = fileparse($ivcfs, qr/\.$extension$/);	
	print ("$ivcfs\n");
    print ("$filenm_root\n"); 		
	
	#ConvertVCFgz
	my $Convert_vcf = "$bcftools_path/bcftools view $ivcfs -Oz -o $filenm_root.vcf.gz";
	print ("Convert_vcf vcf  ...\n");
	system("$Convert_vcf");
	
	my $Index_vcf = "$bcftools_path/bcftools index $filenm_root.vcf.gz";
	print ("Index vcf  ...\n");
	system("$Index_vcf");
	
}


my $filenm_root = fileparse($vcfgzs, qr/\.$extension$/);
print ("$vcfgzs\n"); 
my $merge_vcf = "$bcftools_path/bcftools merge $vcfgzs --threads 5 -o $output_prefix.vcf.gz ";
print ("Merge vcf  ...\n");
system("$merge_vcf");


#RemoveVariant
#my $RemoveVariant = "sed 's/.variant[0-9]//' $output_prefix.vcf > Test.vcf";
#print ("Remove the variant caracter on output $output_prefix.vcf...\n");
#system("$RemoveVariant");


#SelectVariants
#my $SelectVariants = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V $output_prefix.vcf --restrictAllelesTo BIALLELIC -o $output_prefix._biallelic.vcf";
#print ("Select biallelic SNP with SelectVariants on files $output_prefix.vcf...\n");
#system("$SelectVariants");


print "Finished\n";

exit 1;
