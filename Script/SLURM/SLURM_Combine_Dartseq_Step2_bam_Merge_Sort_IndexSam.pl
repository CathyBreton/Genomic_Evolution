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

#Variables Environnement
#foreach (sort keys %ENV) { 
#  print "$_  =  $ENV{$_}\n"; 
#}




#List of all file in directory with specific extension
opendir(CWD,"$directory") or die("Cannot open current directory\n");  
my @files = grep { m/\.$extension$/ } readdir(CWD);
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




#List of the bam files 
foreach my $input_bam (@files) {
   #Filename without extension
   #my $filenm_root = fileparse($input_bam, qr/\.$extension$/);
   my $filenm_prefix = fileparse($input_bam, qr/\_.+\.$extension$/);
      #print ("$filenm_root\n");
      #print ("$filenm_prefix\n");
   #print ("$file\n");
}



#List of all prefix in directory with specific extension
opendir(CWD,"$directory") or die("Cannot open current directory\n");  
my @prefixe = grep { m/\_.+\.$extension$/ } readdir(CWD);
closedir(CWD);





#========================================="
#List of the prefix files 


# Function to get the uniq prefix.
sub get_unique {
   my %seen;
   grep !$seen{$_}++, @_;
} 
#========================================="



### Table with individual
#my @unique = ();
my @uniq = ();
foreach my $prefix (@prefixe) 
{
   my $file_prefix = fileparse($prefix, qr/\_.+\.$extension$/);
   #print ("$file_prefix\n");
    #if we get here, we have not seen it before
         push(@uniq, $file_prefix);
}
#my @unique_inds = get_unique(@uniq);
##print "@unique_inds\n"; 
#print "@unique\n"; 



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
##### Step 1 check the replicate on the prefix . 
#
#=================================================================================================================================


my @bamfile_rep = (); # Correspond to the bam replicates
my $sample_bam = (); 
my %individuals; ## Faire un aray vide
foreach my $file (@files)
{
   #if ($file =~ m/^(\d+)_(.*)\.fastq$/) 

   if ($file =~ m/^(.\d.)_(.*)\.bam/) ## regex for all the sample H20 or 654
   #if ($file =~ m/^(\d+)_(.*)\.bam/) ## regex for the sample 654
   #if ($file =~ m/^(.*)_(.*)\.bam$/) 
      {
         $individuals{$1} ||= [];
         push(@{$individuals{$1}}, $2);
         #print "matched: '$file'\n";
      }
   else
      {
         print "Unmatched: '$file'\n";
      }

}


# Find the replicate 
my $individual =();
foreach my $individual (keys(%individuals)) 
      {
         foreach my $position (@{$individuals{$individual}})
               {
                  #print $individual . '_' . $position . ".bam\n";
                  my $number_of_files = @{$individuals{$individual}};
                  #print "$individual\n" ;
                  #print "$number_of_files\n";
                  #print scalar "@{$individuals{$individual}}\n";
                  #print $individual . '_' . $position . ".fastq\n";
                  my $sample_bam = $individual . '_' . $position . ".bam";
                  

                     if (($number_of_files > 1) )
                        {
                           #print $individual . '_' . $position . ".bam\n";
                           my $sample_bam = $individual . '_' . $position . ".bam";
                           my $out_bam = $individual . '_' . $position ;
                           #print "Matches for $individual:$sample_bam \n";
                           #print "$individual \n";
                           #print "$sample_bam \n";

                           push(@bamfile_rep, $sample_bam);  

                           #Add read group, and sort bam , with the tools picardtools

                           #my $add_rg_rep = "$java_path/java -Xmx1g -jar $picard_path/picard.jar AddOrReplaceReadGroups"
                           #. ' ' ."I=$sample_bam"
                           #. ' ' . "O=$out_bam.reheader.bam"
                           #. ' ' . "RGLB=$cultivar"
                           #. ' ' . "RGPL=illumina"
                           #. ' ' . "RGPU=BforBB"
                           #. ' ' . "RGSM=$individual"
                           #. ' ' . "\n";

                           #print ("Adding Read Group and sorting on file $sample_bam...\n");
                           #print ("$add_rg_rep\n");
                           #system("$add_rg_rep");
                           #if (system("$add_rg_rep"))
                           #{
                              #system() returned an error code, execution failed. 
                              #warn "Failed to execute: $!\n";
                           #} 



                        }
               }
      }


### check the replicates.
#print("@bamfile_rep\n"); 



#print "@bamfile_rep \n";    
my %bam_by_pre_reheader=();## Faire un aray vide 
foreach my $bam (@bamfile_rep)
{


#=================================================================================================================================
#
##### Step  Keep the prefix to be able to use it as Sample Name
#
#=================================================================================================================================



    # Récupérer le préfixe numérique (pour un autre type de préfixe, il faudra
    # ajusetr la regex).
    # Cette regex match tout ce qui commence par au moins 1 chiffre suivi d'un
    # "_" ou d'un "." et le capture dans $1 (pas de capture sur le second groupe
    # de parenthèses car il y a le modificateur "?:").
    #if ($bam =~ m/^(\d+)(?:_|\.)/)
    if ($bam =~ m/^(.\d.)_(.*)\.reheader.bam/) ## regEx to find everything with a letter and a number example H20 or 923
    #if ($bam =~ m/^(\d+)_(.*)\.reheader.bam/) ## regEx to find a number example 923
    #if ($bam =~ m/^(.+)_(.*)\.reheader.bam/)
    {
        my $prefix = $1;
        my $description = $2;
        #print "$bam --> $prefix\n";
        #print "$bam --> $description\n";
        # Initialise à un array vide si pas encore initialisé, sinon laisse comme c'était.
        $bam_by_pre_reheader{$prefix} ||= [];
        # Ajoute le nouveau nom de fichier à la liste.
        push(@{$bam_by_pre_reheader{$prefix}}, $bam);
        #print "@{$bam_by_pre{$prefix}}\n";

   }
}    


#=================================================================================================================================
#
##### Step  Combining sorting and Indexing the replicated bam files
#
#=================================================================================================================================


# Maintenant que les fichiers sont triés, nous allons procéder par préfix.
foreach my $prefix (keys(%bam_by_pre_reheader))
{

   print "#============================================================#\n";
   print "Processing $prefix files\n";

   print "#--------------------------------#\n";
   print "Step 2 Combine Bam\n";
   print "#--------------------------------#\n";

   my $combine_bams_cmd =
      "singularity exec $sambamba_path sambamba merge"
      . " $directory/$prefix.indelrealigned.RG.sorted.recalibrated.reheader.Merged.bam"
      . ' ' . join(' ', @{$bam_by_pre_reheader{$prefix}})
      . "\n";
   print ("Merge Bam with samtools on files $prefix.*.bam...\n");
   print "  $combine_bams_cmd\n"; #+debug
   #system($combine_bams_cmd) and die "cannot run sambamba : $!";
   if (system("$combine_bams_cmd"))
   {
         # system() returned an error code, execution failed.
         warn "Failed to execute: $!\n";
   } 
   


   print "#--------------------------------#\n";
   print "Step 3 Sort Bam with Sambamba\n";
   print "#--------------------------------#\n";

   my $SortBamSambamba_cmd =
      "singularity exec $sambamba_path sambamba sort"
      . ' ' . "--tmpdir tmp"
      . ' ' . "-o $directory/$prefix.indelrealigned.RG.sorted.recalibrated.reheader.Merged.sorted.bam"
      . ' ' . "$directory/$prefix.indelrealigned.RG.sorted.recalibrated.reheader.Merged.bam"
      . "\n";
   print ("Sort Bam with sambamba on files $prefix.indelrealigned.RG.sorted.recalibrated.reheader.Merged.bam...\n");
   print "  $SortBamSambamba_cmd\n"; #+debug
   if (system("$SortBamSambamba_cmd"))
   {
       # system() returned an error code, execution failed.
         warn "Failed to execute: $!\n";
   } 




   print "#--------------------------------#\n";
   print "Step 4 IndexBam with Samtools \n"; 
   print "#--------------------------------#\n"; 

   my $IndexBam_cmd =
      "singularity exec $samtools_path samtools index"
      . ' ' . "$directory/$prefix.indelrealigned.RG.sorted.recalibrated.reheader.Merged.sorted.bam"
      . "\n";
   print ("Index Bam with samtools on files $prefix.indelrealigned.RG.sorted.recalibrated.reheader.Merged.sorted.bam...\n");
   print "  $IndexBam_cmd\n"; #+debug
   if (system("$IndexBam_cmd"))
     {
         # system() returned an error code, execution failed.
         warn "Failed to execute: $!\n";
     } 



    
    print "#============================================================#\n";

}


print "Finished\n";

exit 1;
