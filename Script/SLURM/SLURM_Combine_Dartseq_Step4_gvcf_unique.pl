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
#  Script Name : 
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


if (@ARGV < 5)
{
        print "\n\nUsage: perl $0 -x extension_file_to_treat -d directory_in -do directory_out\n\n\n";
        exit 1;
}




#Variables
my ($help, $man, $resume, $extension, $directory, $directory_out, $cultivar);

#Getoptions
GetOptions 
(   "x|extension=s"           => \$extension, 
    "d|directory_in=s"        => \$directory,
    "do|directory_out=s"      => \$directory_out 
) or die ("Error in command line arguments\n");

#Variables Environnement
#foreach (sort keys %ENV) { 
#  print "$_  =  $ENV{$_}\n"; 
#}

if (!-e($directory))
{
        print "ERROR: no  directory to process with option -d !\n\n";
        pod2usage(2);
}

if (!-e($directory_out))
{
        print "ERROR: no  directory to process with option -do !\n\n";
        
        my $create_directory_cmd = "mkdir $directory_out";
        print ("Create the directory...\n");
        #if (system("$create_directory_cmd"))#
        #        {
        #                warn "Failed to execute: $!\n"; # system() returned an error code, execution failed.
        #        }   
        pod2usage(2);
}


if (!($extension))
{
        print "Warning : no extension was provided with -x, so we will use gvcf!\n\n";
        $extension = 'gvcf';
}


#List of all file in directory with specific extension
opendir(CWD,"$directory") or die("Cannot open current directory\n");  
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


#List of all file in directory with specific extension
opendir(CWD,"$directory") or die("Cannot open current directory\n");  
my @names = grep { m/\_.+\.$extension$/ } readdir(CWD);
closedir(CWD);



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
my @unique_inds = get_unique(@uniq);
#print "@unique_inds\n"; 
#print "@unique\n"; 




#=================================================================================================================================

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#Programs 

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#=================================================================================================================================


# In sge the job array is work , not with slurm : has to be improve.


#=================================================================================================================================
#
##### Step 1 check the replicate on the prefix . 
#
#=================================================================================================================================


my @gvcffile_rep = (); # Correspond to the gvcf replicates
my $sample_gvcf = (); 
my %individuals; ## Faire un aray vide
foreach my $file (@files)
{
  
   #if ($file =~ m/^(.\d.)_(.*)\.gvcf/) ## regex for all the sample H20 or 654
   #if ($file =~ m/^(\d+)_(.*)\.gvcf/) ## regex for the sample 654  for BforBB data
   if ($file =~ m/^(.*)_(\d+.*)\.gvcf/) ## regex for the RTB Data 
   
   #if ($file =~ m/^(.*)_(.*)\.gvcf$/) 
      {
         $individuals{$1} ||= [];
         push(@{$individuals{$1}}, $2);
         print "matched: '$file'\n";
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
                  print $individual . '_' . $position . ".gvcf\n";
                  my $number_of_files = @{$individuals{$individual}};
                  print "$individual\n" ;
                  print "$number_of_files\n";
                  print scalar "@{$individuals{$individual}}\n";
                  print $individual . '_' . $position . ".gvcf\n";
                  my $sample_gvcf = $individual . '_' . $position . ".gvcf";
                  

                     if (($number_of_files == 1) )
                        {
                           print $individual . '_' . $position . ".gvcf\n";
                           #my $sample_gvcf = $individual . '_' . $position . ".gvcf";
                           #my $out_gvcf = $individual . '_' . $position ;
                           #print "Matches for $individual:$sample_gvcf \n";
                           #print "$individual \n";
                           #print "$sample_gvcf \n";
                           push(@gvcffile_rep, $sample_gvcf);                         
                        }
               }
      }





my %gvcf_by_pre=();## Faire un aray vide 
foreach my $gvcf (@gvcffile_rep)
{
#print "$gvcf \n";

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
    #if ($gvcf =~ m/^(\d+)(?:_|\.)/)

   #if ($gvcf =~ m/^(\d+)_(.*)\.gvcf/) ## regEx to find a number example 923 for BforBB 

   if ($gvcf =~ m/^(.*)_(\d+.*)\.gvcf/) ## regEx to find a number example 923 for RTB 
   
   #if ($gvcf =~ m/^(.\d.)_(.*)\.gvcf/) ## regEx to find everything with a letter and a number example H20 or 923
   
   #if ($gvcf =~ m/^(.*)_(.*)\.gvcf/) ## regEx to find H20_5   or 923_5
      {
        my $prefix = $1;
        my $description = $2;
        print "$gvcf --> $prefix\n";
        print "$gvcf --> $description\n";
        #print $prefix . '_' . $description;

        my $name_file = $prefix . '_' . $description;
        print "$name_file\n";
        # Initialise à un array vide si pas encore initialisé, sinon laisse comme c'était.
        $gvcf_by_pre{$prefix} ||= [];
        # Ajoute le nouveau nom de fichier à la liste.
        push(@{$gvcf_by_pre{$prefix}}, $gvcf);
        #print "@{$gvcf_by_pre{$prefix}}\n";


#=================================================================================================================================
#
##### Step  Add read groupe for each replicate with the same Sampl Name. 
#
#=================================================================================================================================


#VariantAnnotator GATK4 BETA (add annotation: MappingQualityZero)




         print "#============================================================#\n";
         print "Processing $prefix files\n";


         print "#--------------------------------#\n";
         print "Step 4 Transfer gvcf              \n";
         print "#--------------------------------#\n";

         
         my $transfert_gvcf_cmd = " cp $gvcf $directory_out/ ";
         
         print ("Transfer files $gvcf uniq ...\n");
         print "  $transfert_gvcf_cmd\n"; #+debug
         if (system("$transfert_gvcf_cmd"))
         {   
               warn "Failed to execute: $!\n";  # system() returned an error code, execution failed.
         } 
         print "#============================================================#\n";

   
      }    



}

print "Finished\n";

exit 1;




