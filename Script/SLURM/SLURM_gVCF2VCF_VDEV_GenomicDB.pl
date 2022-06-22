#! /usr/bin/perl



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
#  Script Name : SLURM_gVCF2VCF_VDEV_GenomicDB.pl
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


=pod

=head1 NAME

SLURM_gVCF2VCF_VDEV.pl

=head1 SYNOPSIS

perl SLURM_gVCF2VCF_VDEV_GenomicDB.pl -r <reference fasta file assembly> -p <file prefix> -v <version number eg. 4.2> -x <extension gvcf> -d <path directory>


=head1 REQUIRES

Perl5

=head1 DESCRIPTION

convert format with GATK

=cut


#use strict;
use warnings;
use Carp qw (cluck confess croak);
use Getopt::Long;
use File::Basename;
use File::Path;
use Pod::Usage;
use File::Path qw(mkpath);

our $stdout = \*STDOUT;
our $stderr = \*STDERR;
our $debug = 1;

my $java_path           = "/nfs/work/agap_id-bin/img/java/jre1.8.0_31/bin";
my $picard_path_2_24    = "/nfs/work/agap_id-bin/img/picard-tools/2.24.0/picard-tools.2.24.0.img";
my $picard_path_2_72    = "/nfs/work/agap_id-bin/img/picard/2.7.2";
my $samtools_path       = "/nfs/work/agap_id-bin/img/samtools/1.2/samtools.1.2.img";
my $tabix_path          = "/home/bretonc/work_agap_id-bin/img/tabix/0.2.6";

my %gatk_version = (
    '3.6' => "/nfs/work/agap_id-bin/img/GATK/3.6",
    '4.2' => "/nfs/work/agap_id-bin/img/GATK/4.2.0.0" 
);



#Variables
my ($help, $man, $resume, $extension, $reference, $output_prefix, $directory, $version);

#Getoptions
GetOptions 
(   'help|?'                => \$help,
    'debug'                 => \$debug,
    'resume=s'              => \$resume,
    "x|extension=s"         => \$extension, 
    "r|reference=s"         => \$reference, 
    "d|directorty=s"        => \$directory,
    "v|version=s"           => \$version,
    "p|outputprefix=s"      => \$output_prefix
) 
or pod2usage(0);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}


if (!-e($reference))
{
        print "ERROR: please indicate a reference genome to process with option -r!\n\n";
        pod2usage(2);
}

if (!-e($directory))
{
        print "ERROR: no  directory to process with option -d !\n\n";
        pod2usage(2);
}

if (!($extension))
{
        print "Warning : no extension was provided with -x, so we will use gvcf!\n\n";
        $extension = 'gvcf';
}

my $gatk_path;
if ($version && exists $gatk_version{'$version'})
{
    $gatk_path = $gatk_version{'$version'}; 
}
else
{
    print "Warning: unknown version or no version specified! we will use version 4.2 \n\n";
    $gatk_path = $gatk_version{'4.2'};  
}


#Create the GenomicBD folder
my $dir_GenomicDB = "$directory/$output_prefix\_GenomicDB";
# if dir not exists create it
unless (-d "$dir_GenomicDB") {mkpath("$dir_GenomicDB");}


#Create the VCF folder
my $dir_VCF = "$directory/$output_prefix\_VCF";
# if dir not exists create it
unless (-d "$dir_VCF") {mkpath("$dir_VCF");}



#List of all file in directory with specific extension
my $igvcfs;
opendir(CWD,"$directory") or die("Cannot open current directory: $directory\n");

my @files = grep { m/\.$extension$/ } readdir(CWD);

#Concat igvcfs filenames
my $i;
foreach (@files) {
	$igvcfs .= "-V $directory/$_ ";
	$i++;
}
closedir(CWD);


# The List of chromosome in an array
my @ChromArray = ('chr01' , 'chr02' , 'chr03' , 'chr04' , 'chr05' , 'chr06', 'chr07' , 'chr08' , 'chr09' , 'chr10' , 'chr11' , 'chloro', 'chrUn_random', 'mito1' , 'mito2' , 'mito3' , 'mito4' , 'mito5' , 'mito6' , 'mito7' , 'mito8' , 'mito9' , 'mito10' , 'mito11' , 'mito12' );
        ####=================================================
        ##manually determine which chromosomes to work with
        @chromosomes=();
        push (@chromosomes, @ChromArray);


#========================================================================================================================================================

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#Programs 

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#=======================================================================================================================================================






foreach $chr (@chromosomes) 
{
    print ( "$chr\n" );
       
    #GenomicDB_gvcf_cmd (GATK4)
           

    if (!-e "$dir_GenomicDB/$chr\_$output_prefix\_gdb")
        {
        my $GenomicDB_gvcf_cmd = "$gatk_path/gatk --java-options -Xmx2G GenomicsDBImport"
                . ' ' . "-R $reference"
        . ' ' . "--genomicsdb-workspace-path $dir_GenomicDB/$chr\_$output_prefix\_gdb"
        . ' ' . "$igvcfs"
        . ' ' . "--tmp-dir $dir_GenomicDB/"
        . ' ' . "--max-num-intervals-to-import-in-parallel 3"
        . ' ' . "--intervals $chr"
        . ' ' . "\n"; 
        print ("Generate a Genomic Database with GenomicsDBImport on file $igvcfs\n with $GenomicDB_gvcf_cmd ...\n");

        if (system($GenomicDB_gvcf_cmd))#
            {
            # system() returned an error code, execution failed.
            warn "Failed to execute: $!\n";
            }        
    }
    else
    {
        print("$dir_GenomicDB/$chr\_$output_prefix\_gdb exists! the file will not be overwritten : skip combine step \n\n");
    }



  

    #GenotypeGVCFs (GATK4) use --allow-old-rms-mapping-quality-annotation-data for older version than 4.1.6  , 

    
    if (!-e "$dir_VCF/$chr\_$output_prefix.vcf.gz")
        {
        my $GenotypeGVCFs_cmd = "$gatk_path/gatk --java-options -Xmx2G GenotypeGVCFs"
        . ' ' . "-R $reference"
        . ' ' . "-V gendb://$output_prefix\_GenomicDB/$chr\_$output_prefix\_gdb"
        . ' ' . "-O $dir_VCF/$chr\_$output_prefix.vcf.gz"
        . ' ' . "--genomicsdb-use-bcf-codec"
        . ' ' . "--tmp-dir $dir_VCF/"
        . ' ' . "\n"; 
        print ("Genotype combination with GenotypeGVCFs file $chr\_$output_prefix\_gdb\n with  $GenotypeGVCFs_cmd...\n");

            if (system($GenotypeGVCFs_cmd))#
                {
                 #system() returned an error code, execution failed.
                warn "Failed to execute: $!\n";
                }      
        }      
    else
        {
            print("$dir_VCF/$chr\_$output_prefix.vcf.gz ! the file will not be overwritten : skip Genotype step \n\n");
        }
 

}

#####=========================================










print "Finished\n";




=pod

=head2 printDebug

B<Description>: display a debug message on STDERR if debug mode is ON.

B<ArgsCount>: 1

=over 4

=item $message: (string) (R)

debug message to display.

=back

B<Example>:

    printDebug("Starting the process using paramaters:!\na=25\nb=50\nc=3");

=cut

sub printDebug
{
    my ($message) = @_;
    # check arguments
    if (1 != @_)
    {
        confess "usage: printDebug(message);";
    }

    # check if debugging is disabled
    if (!$debug)
    {
        #debugging disabled, nothing to do!
        return;
    }

    # remove trailing invisible characters
    $message =~ s/\s+$//s;
    print {$stderr} "\n";
    # check for multi-lines warning
    if ($message =~ m/[\r\n]/)
    {
        # multi-lines
        my @message_lines = split(/(?:\n|\r\n|\r)/, $message);
        print {$stderr} "DEBUG: " . shift(@message_lines) . "\n       ";
        print {$stderr} join("\n       ", @message_lines);
        print {$stderr} "\n";
    }
    else
    {
        # single line
        print {$stderr} "DEBUG: " . $message . "\n";
    }
    return;
}

exit(0);

# CODE END
###########

=pod

=head1 AUTHORS

Catherine Breton c.breton@cgiar.org


###################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
#
#
# Written by Catherine Breton
#
#################### 
