#! /usr/bin/perl


=pod

=head1 NAME

gVCF2VCF.pl

=head1 SYNOPSIS

sbatch perl SLURM_gVCF2VCF_gz.pl -r <reference fasta file assembly> -p <file prefix> -v <version number eg. 4.0, 4.2> -x <extension gvcf> -d <path directory> -m <memory>


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




if (@ARGV < 7)
{
        print "\n\nUsage: perl $0 -r reference_file_path -p output_prefix -v version -x extension_file_to_treat -d directory -m memory\n\n\n";
        exit 1;
}


#Variables
my ($help, $man, $resume, $extension, $reference, $output_prefix, $directory, $version, $memory);

#Getoptions
GetOptions 
(   'help|?'                => \$help,
    'debug'                 => \$debug,
    'resume=s'              => \$resume,
    "x|extension=s"         => \$extension, 
    "r|reference=s"         => \$reference, 
    "d|directorty=s"        => \$directory,
    "v|version=s"           => \$version,
    "p|outputprefix=s"      => \$output_prefix,
    "m|memory=s"            => \$memory  
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



#========================================================================================================================================================

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#Programs 

###---------------------------------------------------------------------------------------------------------------------------------------------------###

#=======================================================================================================================================================



## If the index is needed use this code.
#Indexgvcfs
#my $Indexgvcfs = "$tabix_path/tabix -p vcf $gvcfs";
#print ("Index gVcf on files $gvcfs...\n");
#system("$Indexgvcfs");


#CombineGVCFs (GATK4)
if (!-e "$directory/$output_prefix.gvcf.gz")
{
    my $combinegvcfs_cmd = "$gatk_path/gatk --java-options -$memory CombineGVCFs -R $reference $igvcfs -O $directory/$output_prefix.gvcf.gz";
    print("Call SNPs with CombineGVCFs on $i files...\n");
    printDebug($combinegvcfs_cmd);
    system($combinegvcfs_cmd) and die "cannot run gatk : $!";
}
else
{
    print("$directory/$output_prefix.gvcf exists! the file will not be overwritten : skip combine step \n\n");
}


#GenotypeGVCFs (GATK4) use --allow-old-rms-mapping-quality-annotation-data for older version than 4.1.6  , 
if (-e "$directory/$output_prefix.gvcf.gz")
{
    if (!-e "$directory/$output_prefix.Genotype.vcf.gz")
    {
        my $GenotypeGVCFs_cmd = "$gatk_path/gatk --java-options -$memory GenotypeGVCFs --allow-old-rms-mapping-quality-annotation-data -R $reference -V $directory/$output_prefix.gvcf.gz -O $directory/$output_prefix.Genotype.vcf.gz";
        print ("Call SNPs with GenotypeGVCFs on files $directory/$output_prefix.gvcf.gz...\n");
        printDebug($GenotypeGVCFs_cmd);
        system($GenotypeGVCFs_cmd) and die "cannot run gatk : $!";
    }
    else
    {
        print("$directory/$output_prefix.Genotype.vcf.gz! the file will not be overwritten : skip Genotype step \n\n");
    }
}
else {
    print("combined gVCF file missing. please check prevoius step!");
}


#VariantFiltration (GATK4)
if (-e "$directory/$output_prefix.Genotype.vcf.gz")
{
    my $variant_filtration_cmd  = "$gatk_path/gatk --java-options -$memory VariantFiltration -R $reference -V $directory/$output_prefix.Genotype.vcf.gz --cluster-size 3 --cluster-window-size 10 --filter-expression 'QD < 1.5' --filter-name QD_FILTER --mask-extension 0 --mask-name Mask -O $directory/$output_prefix.VF.vcf.gz";
    print ("Variant filtration --filterExpression 'QD < 1.5' on $output_prefix.vcf.gz...\n...\n");
    printDebug($variant_filtration_cmd);
    system($variant_filtration_cmd) and die "cannot run gatk : $!";
}
else 
{
    print ("genotype VCF file missing. please check prevoius step!");
}

#SelectVariant to filter out filtered variants (GATK4)
if (-e "$directory/$output_prefix.VF.vcf.gz")
{
    my $select_variant_cmd = "$gatk_path/gatk --java-options -$memory SelectVariants -R $reference -V $directory/$output_prefix.VF.vcf.gz -select 'vc.isNotFiltered()' -select-type SNP -O $directory/$output_prefix.VF.filtered.vcf.gz";
    print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $output_prefix.VF.vcf.gz...\n");
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




