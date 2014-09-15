#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Copy;
use File::Basename;
use Cwd;

our $VERSION = 			'$Revision: 518 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2013-02-08 11:10:50 -0800 (Fri, 08 Feb 2013) $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $buildver, $genefile, $omimfile);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'buildver=s'=>\$buildver, 'genefile=s'=>\$genefile,
	'omimfile=s'=>\$omimfile) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

my $path=cwd();
#$path and $ENV{PATH} = "$path:$ENV{PATH}";		#set up the system executable path to include the path where this program is located in

($queryfile, $dbloc) = @ARGV;

$outfile ||= "$queryfile.table";



if (not defined $buildver) {
	$buildver = 'hg18';
	print STDERR "NOTICE: the --buildver argument is set as 'hg18' by default\n";
}
$buildver eq 'hg18' or $buildver eq 'hg19' or pod2usage ("Error in argument: the --buildver argument must be 'hg18' or 'hg19'");

my $sc;
$sc = "perl $path/../../bin/convert2annovar.pl -format bed $queryfile > $outfile.avinput";
system ($sc) and die "Error: cannot execute system command $sc\n";

$sc = "perl $path/../../bin/annotate_variation.pl -geneanno -buildver $buildver -outfile $outfile $outfile.avinput $dbloc";
system ($sc) and die "Error: cannot execute system command $sc\n";

my ($countregion, $countexonic) = qw/0 0/;
my (%hitgene, $totallen);
open (VF, "$outfile.variant_function") or die "Error: cannot read $outfile.variant_function file\n";
while (<VF>) {
	my @field = split (/\t/, $_);
	$countregion++;
	$field[0] =~ m/exonic/ or next;
	$countexonic++;
	my @gene = split (/,/, $field[1]);
	
	for my $gene (@gene) {
		$hitgene{$gene} = "$field[2]:$field[3]-$field[4]";
	}
	$totallen += ($field[4]-$field[3]+1);
}

print STDERR "NOTICE: Among $countregion BED regions ($totallen base pairs), $countexonic disrupt exons, and ", scalar keys %hitgene, " genes are affected\n";


my (%omim);
if (defined $omimfile) {
	print STDERR "NOTICE: reading OMIM to fill in more information\n";
	open (OMIM, $omimfile) or die "Error: cannot read from omimfile $omimfile\n";
	while (<OMIM>) {
		s/[\r\n]+$//;
		my @field = split (/\|/, $_);
		my ($cyto, $symbol, $title1, $title2, $disorder1, $disorder2, $disorder3) = @field[4,5,7,8,13,14,15];
		$symbol =~ m/(\w+),/ and $symbol = $1;		#take the first symbol when synonyms are present
		$omim{$symbol} = join ("\t", $cyto, $title1 . $title2, $disorder1 . $disorder2 . $disorder3);
	}
}
	
my (%candgene);
if (defined $genefile) {
	my @genefile = split (/,/, $genefile);
	for my $nextfile (@genefile) {
		print STDERR "NOTICE: the list of genes that overlap with user-supplied genefile are below:\n";
		open (GENE, $nextfile) or die "Error: cannot read from genefile $nextfile\n";
		while (<GENE>) {
			m/^(\w+)/ and $candgene{$1}++;
		}
		close (GENE);
	}
	for my $key (keys %hitgene) {
		$candgene{$key} or delete $hitgene{$key};
	}
}

open (OUT, ">$outfile") or die "Error: cannot write to output file $outfile.table: $!\n";
for my $key (keys %hitgene) {
	print OUT "$hitgene{$key}\t$key", defined $omim{$key}?"\t$omim{$key}":"\t\t\t", "\n";
}
close (OUT);
print STDERR "NOTICE: ", scalar (keys %hitgene), " candidate genes are identified\n";

exit (0);

=head1 SYNOPSIS

 bed2gene.pl [arguments] <query-file> <database-location>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --outfile <string>		output file name prefix
            --buildver <string>		genome build version (default: hg18)
            

 Function: automatically run a pipeline on a list of variants (potentially 
 whole-genome SNPs from a patient with Mendelian disease) and identify a small 
 subset that are most likely causal for Mendelian diseases
 
 Example: bed2gene.pl infile humandb/ -buildver hg19
                  
 Version: $LastChangedDate: 2013-02-08 11:10:50 -0800 (Fri, 08 Feb 2013) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--outfile>

the prefix of output file names

=item B<--buildver>

specify the genome build version

=item B<--remove>

remove all temporary files. By default, all temporary files will be kept for 
user inspection, but this will easily clutter the directory.

=item B<--genetype>

specify the gene definition, such as refgene (default), ucsc known gene, ensembl 
gene and gencode gene.

=item B<--maf_threshold>

specify the MAF threshold for allele frequency databases. This argument works 
for 1000 Genomes Project, ESP database and CG (complete genomics) database.

=item B<--checkfile>

the program will check if all required database files exist before execution of annotation

=item B<--genericdbfile>

specify the genericdb file used in -dbtype generic

=item B<--ljb_sift_threshold>

specify the LJB_SIFT threshold for filter operation (default: 0.95)

=item B<--ljb_pp2_threshold>

specify the LJB_PP2 threshold for filter operation (default: 0.85)

=back

=head1 DESCRIPTION

ANNOVAR is a software tool that can be used to functionally annotate a list of 
genetic variants, possibly generated from next-generation sequencing 
experiments. For example, given a whole-genome resequencing data set for a human 
with specific diseases, typically around 3 million SNPs and around half million 
insertions/deletions will be identified. Given this massive amounts of data (and 
candidate disease- causing variants), it is necessary to have a fast algorithm 
that scans the data and identify a prioritized subset of variants that are most 
likely functional for follow-up Sanger sequencing studies and functional assays.

by default, for hg18, the arguments are

variants_reduction.pl x1.avinput humandb -protocol nonsyn_splicing,1000g2010jul_ceu,1000g2010jul_jptchb,1000g2010jul_yri,snp132,esp5400_ea,esp5400_aa,recessive -operation g,f,f,f,f,f,f,m

for hg19, the arguments are

variants_reduction.pl x1.avinput humandb -protocol nonsyn_splicing,1000g2012feb_all,snp135,esp5400_ea,esp5400_aa,recessive -operation g,f,f,f,m


ANNOVAR is freely available to the community for non-commercial use. For 
questions or comments, please contact kai@openbioinformatics.org.


=cut