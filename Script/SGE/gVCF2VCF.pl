#! /usr/bin/perl

#$ -q normal.q
#$ -cwd
#$ -V
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1
#$ -l mem_free=10G


=pod

=head1 NAME

gVCF2VCF.pl

=head1 SYNOPSIS

qsub perl gVCF2VCF.pl -r <reference fasta file assembly> -p <file prefix> -v <version number eg. 4.0, 4.1> -x <extension gvcf> -d <path directory>


=head1 REQUIRES

Perl5
module load bioinfo/GATK/4.1.9.0
module load bioinfo/picard-tools/1.130
module load bioinfo/samtools/1.9
module load bioinfo/tabix/0.2.6
module load bioinfo/vcftools

=head1 DESCRIPTION

convert format with samtools


=cut


use strict;
use warnings;
use Carp qw (cluck confess croak);
use Getopt::Long;
use File::Basename;
use File::Path;
use Pod::Usage;


our $stdout = \*STDOUT;
our $stderr = \*STDERR;
our $debug = 1;

my $java_path       = "/usr/local/java/jre8/bin";
my $picard_path     = "/usr/local/bioinfo/picard-tools/1.130";
my $samtools_path   = "/usr/local/bioinfo/samtools/1.9/bin/";
my $tabix_path      = "/usr/local/bioinfo/tabix/0.2.6";

my %gatk_version = (
    '3.4' => "/usr/local/bioinfo/GenomeAnalysisTK/3.4-46",
    '3.7' => "/usr/local/bioinfo/GenomeAnalysisTK/3.7-0",
    '4.0' => "/usr/local/bioinfo/GenomeAnalysisTK/4.0.5.2",
    '4.1' => "/usr/local/bioinfo/GenomeAnalysisTK/4.1.9.0"
);



#Variables
my ($help, $man, $resume, $extension, $reference, $output_prefix, $directory, $version);

#Getoptions
GetOptions 
(   'help|?'             => \$help,
    'debug'              => \$debug,
    'resume=s'            => \$resume,
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
    print "Warning: unknown version or no version specified! we will use version 4.0 \n\n";
    $gatk_path = $gatk_version{'4.0'};  
}


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



#Programs 

#Indexgvcfs
#my $Indexgvcfs = "$tabix_path/tabix -p vcf $gvcfs";
#print ("Index gVcf on files $gvcfs...\n");
#system("$Indexgvcfs");


#CombineGVCFs (GATK4)
if (!-e "$directory/$output_prefix.gvcf")
{
    my $combinegvcfs_cmd = "$gatk_path/gatk --java-options -Xmx5g CombineGVCFs -R $reference $igvcfs -O $directory/$output_prefix.gvcf";
    print("Call SNPs with CombineGVCFs on $i files...\n");
    printDebug($combinegvcfs_cmd);
    system($combinegvcfs_cmd) and die "cannot run gatk : $!";
}
else
{
    print("$directory/$output_prefix.gvcf exists! the file will not be overwritten : skip combine step \n\n");
}


#GenotypeGVCFs (GATK4)
if (-e "$directory/$output_prefix.gvcf")
{
    if (!-e "$directory/$output_prefix.Genotype.vcf")
    {
        my $GenotypeGVCFs_cmd = "$gatk_path/gatk --java-options -Xmx5g GenotypeGVCFs -R $reference -V $directory/$output_prefix.gvcf -O $directory/$output_prefix.Genotype.vcf";
        print ("Call SNPs with GenotypeGVCFs on files $directory/$output_prefix.gvcf...\n");
        printDebug($GenotypeGVCFs_cmd);
        system($GenotypeGVCFs_cmd) and die "cannot run gatk : $!";
    }
    else
    {
        print("$directory/$output_prefix.Genotype.vcf! the file will not be overwritten : skip Genotype step \n\n");
    }
}
else {
    print("combined gVCF file missing. please check prevoius step!");
}


#VariantFiltration (GATK4)
if (-e "$directory/$output_prefix.Genotype.vcf")
{
    my $variant_filtration_cmd  = "$gatk_path/gatk --java-options -Xmx4G VariantFiltration -R $reference -V $directory/$output_prefix.Genotype.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression 'QD < 1.5' --filter-name QD_FILTER --mask-extension 0 --mask-name Mask -O $directory/$output_prefix.VF.vcf";
    print ("Variant filtration --filterExpression 'QD < 1.5' on $output_prefix.vcf...\n...\n");
    printDebug($variant_filtration_cmd);
    system($variant_filtration_cmd) and die "cannot run gatk : $!";
}
else 
{
    print ("genotype VCF file missing. please check prevoius step!");
}

#SelectVariant to filter out filtered variants (GATK4)
if (-e "$directory/$output_prefix.VF.vcf")
{
    my $select_variant_cmd = "$gatk_path/gatk --java-options -Xmx4G SelectVariants -R $reference -V $directory/$output_prefix.VF.vcf -select 'vc.isNotFiltered()' -select-type SNP -O $directory/$output_prefix.VF.filtered.vcf";
    print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $output_prefix.VF.vcf...\n");
    printDebug($select_variant_cmd);
    system($select_variant_cmd)and die "cannot run gatk : $!";
}
else 
{
    print ("filtered genotype VCF file missing. please check prevoius step!");
}



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
Mathieu ROUARD  m.rouard@cgiar.org

###################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
#
#
# Written by Catherine Breton
#
#################### 
