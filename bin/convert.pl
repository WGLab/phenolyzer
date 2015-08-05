#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


our ($verbose, $help, $man);
our ($variantfile);
our ($outfile, $format, $includeinfo, $snpqual, $snppvalue, $coverage, $maxcoverage, $chr, $chrmt, $altcov, $allelicfrac, $fraction, $species, 
	$filterword, $confraction, $allallele, $withzyg, $comment);

our %iupac = (R=>'AG', Y=>'CT', S=>'CG', W=>'AT', K=>'GT', M=>'AC', A=>'AA', C=>'CC', G=>'GG', T=>'TT', B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT', '.'=>'-', '-'=>'-'); ### <<< FOR 5500SOLiD LifeScope ( S=>'GC' is replaced by S=>'CG')
our %iupacrev = reverse %iupac; ### <<< FOR 5500SOLiD LifeScope

GetOptions('verbose'=>\$verbose, 'help|h'=>\$help, 'man'=>\$man, 'outfile=s'=>\$outfile, 'format=s'=>\$format, 'includeinfo'=>\$includeinfo,
	'snpqual=f'=>\$snpqual, 'snppvalue=f'=>\$snppvalue, 'coverage=i'=>\$coverage, 'maxcoverage=i'=>\$maxcoverage, 'chr=s'=>\$chr, 'chrmt=s'=>\$chrmt, 
	'fraction=f'=>\$fraction, 'altcov=i'=>\$altcov, 'allelicfrac'=>\$allelicfrac,
	'species'=>\$species, 'filter=s'=>\$filterword, 'confraction=f'=>\$confraction, 'allallele!'=>\$allallele, 'withzyg'=>\$withzyg,
	'comment'=>\$comment) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 1 or pod2usage ("Syntax error");

($variantfile) = @ARGV;

$chrmt ||= 'M';

if (not $format) {
	$format = 'pileup';
	print STDERR "NOTICE: the default --format argument is set as 'pileup'\n";
}

if (defined $outfile) {
	open (STDOUT, ">$outfile") or die "Error: cannot write to output file $outfile: $!\n";
}

defined $snpqual and $format eq 'pileup' || $format eq 'vcf4' || pod2usage ("Error in argument: the --snpqual is supported only for the 'pileup' or 'vcf4' format");
defined $snppvalue and $format eq 'gff3-solid' || pod2usage ("Error in argument: the --snppvalue is supported only for the 'gff3-solid' format");
if (not defined $snpqual and $format eq 'pileup') {
	$snpqual = 20;
	print STDERR "NOTICE: the default --snpqual argument for pileup format is set as 20\n";
}

if (not defined $snppvalue) {
	$snppvalue = 1;		#by default, do not use any of the P-value cutoff in filtering out GFF3-SOLID files (this is differnt from handling pileup files)
}

if (not defined $coverage) {
	$coverage = 0;
}

if (defined $fraction) {
	$format eq 'pileup' or $format eq 'vcf4' or pod2usage ("Error in argument: the '--fraction' argument is supported for the pileup or vcf4 format only");
	$format eq 'vcf4' and print STDERR "NOTICE: the --fraction argument works ONLY on indels for vcf4 format\n";
	$fraction >= 0 and $fraction <=1 or pod2suage ("Error in argument: the --fraction argument must be between 0 and 1 inclusive");
} else {
	$fraction = 0;
}

if (defined $confraction) {
	$format eq 'vcf4' and print STDERR "NOTICE: the --confraction argument works ONLY on indels for vcf4 format\n";
	$confraction >= 0 and $fraction <=1 or pod2suage ("Error in argument: the --confraction argument must be between 0 and 1 inclusive");
} else {
	$confraction = 0;
}

if (defined $altcov) {
	$format eq 'pileup' or pod2usage ("Error in argument: the '--altcov' argument is supported for the '--format pileup' only");
	$altcov < $coverage or pod2suage ("Error in argument: the --altcov argument must be less than --coverage");
	$altcov > 0 or pod2suage ("Error in argument: the --altcov argument must be a positive integer");
}

if (defined $species) {
	$format eq 'gff3-solid' or pod2usage ("Error in argument: the '--species' argument is only necessary for the '--format gff3-solid'");
}

if ($allallele) {
	$format eq 'vcf4' or pod2usage ("Error in argument: the '--allallele' argument is only supported for the '--format vcf4'");
}

if ($format eq 'pileup') {
	convertPileup ($variantfile);
} elsif ($format eq 'cg') {
	convertCG ($variantfile);
} elsif ($format eq 'cgmastervar') {
	convertCGMasterVar ($variantfile);
} elsif ($format eq 'gff3-solid') {
	convertGFF3SolidSNP ($variantfile);
} elsif ($format eq 'soap') {
	print STDERR "WARNING: the support for '--format soap' is not well developed yet and may contain bugs for indel analysis.\n";
	convertSOAP ($variantfile);
} elsif ($format eq 'maq') {
	print STDERR "WARNING: the support for '--format maq' is not well developed yet and may contain bugs.\n";
	convertMAQSNP ($variantfile);
} elsif ($format eq 'casava') {
	if (not defined $chr) {
		pod2usage ("Error in argument: please supply --chr argument for the '--format casava'");
	}
	convertCASAVA ($variantfile, $chr);
} elsif ($format eq 'vcf4') {
	convertVCF4 ($variantfile);
} elsif ($format eq 'annovar') {
	convertANNOVAR ($variantfile);
} elsif ($format eq 'bed') {
	convertBED ($variantfile);
} else {
	pod2usage ("Error in argument: the --format $format is not currently supported. Please contact ANNOVAR developer for adding the support");
}


sub convertBED {
	my ($variantfile) = @_;
	
	my ($foundheader, $countline, $countvar) = qw/1 0 0/;
	my ($prechr, $prestart, $preend, $prevartype, $preref, $preobs, $prescore, $prexref) = qw/0 0 0 0 0 0 0 0/;

	if ($variantfile eq 'stdin') {
		*VAR = *STDIN;
	} elsif ($variantfile =~ m/^\.gz$/) {
		open (VAR, "gunzip -c $variantfile |") or die "Error: cannot read from STDIN uncompressing variant file $variantfile: $!\n";
	} else {
		open (VAR, $variantfile) or die "Error: cannot read from variant file $variantfile: $!\n";
	}


	print STDERR "NOTICE: Converting variants from $variantfile\n";
	while (<VAR>) {
		s/[\r\n]+$//;
		$countline++;
		m/^#/ and next;		#comment lines
		if (m/^track/) {
			$foundheader++;
			next;
		}
		if (not $foundheader) {
			$countline > 10 and die "Error: invalid BED file format for $variantfile (track record is not found within the first 10 lines)\n";
		}
		my ($chrom, $start, $end, @otherinfo) = split (/\t/, $_);
		$chrom =~ s/^chr//;
		
		print join ("\t", $chrom, $start+1, $end, 0, 0, @otherinfo), "\n";
		$countvar++;
		
		
	}
	print STDERR "NOTICE: Done with $countline lines and $countvar variants\n";
}


