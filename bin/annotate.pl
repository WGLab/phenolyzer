#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
use Cwd;

our $REVISION = '$Revision$';
our $DATE =	'$Date$';  
our $AUTHOR =	'$Author$';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $separate, $batchsize, $dbtype, $neargene, $genomebinsize, $geneanno, $regionanno, $filter, $downdb, $buildver, $score_threshold, $normscore_threshold, $minqueryfrac, $expandbin, $splicing_threshold,
	$maf_threshold, $chromosome, $zerostart, $rawscore, $memfree, $memtotal, $sift_threshold, $gff3dbfile, $genericdbfile, $vcfdbfile, $time, $wget, $precedence,
	$webfrom, $colsWanted, $comment, $scorecolumn, $transfun, $exonsort, $avcolumn, $bedfile, $hgvs, $reverse, $indexfilter_threshold, $otherinfo, $infoasscore,
	$firstcodondel, $aamatrixfile, $gff3attr, $exonicsplicing, $infosep);
our $seq_padding;			# If set, create a new file with cDNA, aa sequence padded by this much on either side.
our $indel_splicing_threshold;		# If set, use this value for allowed indel size for splicing variants, else, use splicing_threshold;
our (%valichr, $dbtype1);
our (@precedence, @colsWanted, @avcolumn);
sub printerr;				#declare a subroutine

our %codon1 = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"*", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codon3 = (TTT=>"Phe", TTC=>"Phe", TCT=>"Ser", TCC=>"Ser", TAT=>"Tyr", TAC=>"Tyr", TGT=>"Cys", TGC=>"Cys", TTA=>"Leu", TCA=>"Ser", TAA=>"*", TGA=>"*", TTG=>"Leu", TCG=>"Ser", TAG=>"*", TGG=>"Trp", CTT=>"Leu", CTC=>"Leu", CCT=>"Pro", CCC=>"Pro", CAT=>"His", CAC=>"His", CGT=>"Arg", CGC=>"Arg", CTA=>"Leu", CTG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", ATT=>"Ile", ATC=>"Ile", ACT=>"Thr", ACC=>"Thr", AAT=>"Asn", AAC=>"Asn", AGT=>"Ser", AGC=>"Ser", ATA=>"Ile", ACA=>"Thr", AAA=>"Lys", AGA=>"Arg", ATG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"Arg", GTT=>"Val", GTC=>"Val", GCT=>"Ala", GCC=>"Ala", GAT=>"Asp", GAC=>"Asp", GGT=>"Gly", GGC=>"Gly", GTA=>"Val", GTG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonfull = (TTT=>"Phenylalanine", TTC=>"Phenylalanine", TCT=>"Serine", TCC=>"Serine", TAT=>"Tyrosine", TAC=>"Tyrosine", TGT=>"Cysteine", TGC=>"Cysteine", TTA=>"Leucine", TCA=>"Serine", TAA=>"Stop", TGA=>"Stop", TTG=>"Leucine", TCG=>"Serine", TAG=>"Stop", TGG=>"Tryptophan", CTT=>"Leucine", CTC=>"Leucine", CCT=>"Proline", CCC=>"Proline", CAT=>"Histidine", CAC=>"Histidine", CGT=>"Arginine", CGC=>"Arginine", CTA=>"Leucine", CTG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", ATT=>"Isoleucine", ATC=>"Isoleucine", ACT=>"Threonine", ACC=>"Threonine", AAT=>"Asparagine", AAC=>"Asparagine", AGT=>"Serine", AGC=>"Serine", ATA=>"Isoleucine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Arginine", ATG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Arginine", GTT=>"Valine", GTC=>"Valine", GCT=>"Alanine", GCC=>"Alanine", GAT=>"Aspartic acid", GAC=>"Aspartic acid", GGT=>"Glycine", GGC=>"Glycine", GTA=>"Valine", GTG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine");
our %codonr1 = (UUU=>"F", UUC=>"F", UCU=>"S", UCC=>"S", UAU=>"Y", UAC=>"Y", UGU=>"C", UGC=>"C", UUA=>"L", UCA=>"S", UAA=>"*", UGA=>"*", UUG=>"L", UCG=>"S", UAG=>"*", UGG=>"W", CUU=>"L", CUC=>"L", CCU=>"P", CCC=>"P", CAU=>"H", CAC=>"H", CGU=>"R", CGC=>"R", CUA=>"L", CUG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", AUU=>"I", AUC=>"I", ACU=>"T", ACC=>"T", AAU=>"N", AAC=>"N", AGU=>"S", AGC=>"S", AUA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", AUG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GUU=>"V", GUC=>"V", GCU=>"A", GCC=>"A", GAU=>"D", GAC=>"D", GGU=>"G", GGC=>"G", GUA=>"V", GUG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codonr3 = (UUU=>"Phe", UUC=>"Phe", UCU=>"Ser", UCC=>"Ser", UAU=>"Tyr", UAC=>"Tyr", UGU=>"Cys", UGC=>"Cys", UUA=>"Leu", UCA=>"Ser", UAA=>"*", UGA=>"*", UUG=>"Leu", UCG=>"Ser", UAG=>"*", UGG=>"Trp", CUU=>"Leu", CUC=>"Leu", CCU=>"Pro", CCC=>"Pro", CAU=>"His", CAC=>"His", CGU=>"Arg", CGC=>"Arg", CUA=>"Leu", CUG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", AUU=>"Ile", AUC=>"Ile", ACU=>"Thr", ACC=>"Thr", AAU=>"Asn", AAC=>"Asn", AGU=>"Ser", AGC=>"Ser", AUA=>"Ile", ACA=>"Thr", AAA=>"Lys", AGA=>"Arg", AUG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"Arg", GUU=>"Val", GUC=>"Val", GCU=>"Ala", GCC=>"Ala", GAU=>"Asp", GAC=>"Asp", GGU=>"Gly", GGC=>"Gly", GUA=>"Val", GUG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonrfull = (UUU=>"Phenylalanine", UUC=>"Phenylalanine", UCU=>"Serine", UCC=>"Serine", UAU=>"Tyrosine", UAC=>"Tyrosine", UGU=>"Cysteine", UGC=>"Cysteine", UUA=>"Leucine", UCA=>"Serine", UAA=>"Stop", UGA=>"Stop", UUG=>"Leucine", UCG=>"Serine", UAG=>"Stop", UGG=>"Tryptophan", CUU=>"Leucine", CUC=>"Leucine", CCU=>"Proline", CCC=>"Proline", CAU=>"Histidine", CAC=>"Histidine", CGU=>"Arginine", CGC=>"Arginine", CUA=>"Leucine", CUG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", AUU=>"Isoleucine", AUC=>"Isoleucine", ACU=>"Threonine", ACC=>"Threonine", AAU=>"Asparagine", AAC=>"Asparagine", AGU=>"Serine", AGC=>"Serine", AUA=>"Isoleucine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Arginine", AUG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Arginine", GUU=>"Valine", GUC=>"Valine", GCU=>"Alanine", GCC=>"Alanine", GAU=>"Aspartic acid", GAC=>"Aspartic acid", GGU=>"Glycine", GGC=>"Glycine", GUA=>"Valine", GUG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine");

our %codon1m = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"W", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"M", ACA=>"T", AAA=>"K", AGA=>"*", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"*", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codon3m = (TTT=>"Phe", TTC=>"Phe", TCT=>"Ser", TCC=>"Ser", TAT=>"Tyr", TAC=>"Tyr", TGT=>"Cys", TGC=>"Cys", TTA=>"Leu", TCA=>"Ser", TAA=>"*", TGA=>"Trp", TTG=>"Leu", TCG=>"Ser", TAG=>"*", TGG=>"Trp", CTT=>"Leu", CTC=>"Leu", CCT=>"Pro", CCC=>"Pro", CAT=>"His", CAC=>"His", CGT=>"Arg", CGC=>"Arg", CTA=>"Leu", CTG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", ATT=>"Ile", ATC=>"Ile", ACT=>"Thr", ACC=>"Thr", AAT=>"Asn", AAC=>"Asn", AGT=>"Ser", AGC=>"Ser", ATA=>"Met", ACA=>"Thr", AAA=>"Lys", AGA=>"*", ATG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"*", GTT=>"Val", GTC=>"Val", GCT=>"Ala", GCC=>"Ala", GAT=>"Asp", GAC=>"Asp", GGT=>"Gly", GGC=>"Gly", GTA=>"Val", GTG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonfullm = (TTT=>"Phenylalanine", TTC=>"Phenylalanine", TCT=>"Serine", TCC=>"Serine", TAT=>"Tyrosine", TAC=>"Tyrosine", TGT=>"Cysteine", TGC=>"Cysteine", TTA=>"Leucine", TCA=>"Serine", TAA=>"Stop", TGA=>"Tryptophan", TTG=>"Leucine", TCG=>"Serine", TAG=>"Stop", TGG=>"Tryptophan", CTT=>"Leucine", CTC=>"Leucine", CCT=>"Proline", CCC=>"Proline", CAT=>"Histidine", CAC=>"Histidine", CGT=>"Arginine", CGC=>"Arginine", CTA=>"Leucine", CTG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", ATT=>"Isoleucine", ATC=>"Isoleucine", ACT=>"Threonine", ACC=>"Threonine", AAT=>"Asparagine", AAC=>"Asparagine", AGT=>"Serine", AGC=>"Serine", ATA=>"Methionine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Stop", ATG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Stop", GTT=>"Valine", GTC=>"Valine", GCT=>"Alanine", GCC=>"Alanine", GAT=>"Aspartic acid", GAC=>"Aspartic acid", GGT=>"Glycine", GGC=>"Glycine", GTA=>"Valine", GTG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine");
our %codonr1m = (UUU=>"F", UUC=>"F", UCU=>"S", UCC=>"S", UAU=>"Y", UAC=>"Y", UGU=>"C", UGC=>"C", UUA=>"L", UCA=>"S", UAA=>"*", UGA=>"W", UUG=>"L", UCG=>"S", UAG=>"*", UGG=>"W", CUU=>"L", CUC=>"L", CCU=>"P", CCC=>"P", CAU=>"H", CAC=>"H", CGU=>"R", CGC=>"R", CUA=>"L", CUG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", AUU=>"I", AUC=>"I", ACU=>"T", ACC=>"T", AAU=>"N", AAC=>"N", AGU=>"S", AGC=>"S", AUA=>"M", ACA=>"T", AAA=>"K", AGA=>"*", AUG=>"M", ACG=>"T", AAG=>"K", AGG=>"*", GUU=>"V", GUC=>"V", GCU=>"A", GCC=>"A", GAU=>"D", GAC=>"D", GGU=>"G", GGC=>"G", GUA=>"V", GUG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codonr3m = (UUU=>"Phe", UUC=>"Phe", UCU=>"Ser", UCC=>"Ser", UAU=>"Tyr", UAC=>"Tyr", UGU=>"Cys", UGC=>"Cys", UUA=>"Leu", UCA=>"Ser", UAA=>"*", UGA=>"Trp", UUG=>"Leu", UCG=>"Ser", UAG=>"*", UGG=>"Trp", CUU=>"Leu", CUC=>"Leu", CCU=>"Pro", CCC=>"Pro", CAU=>"His", CAC=>"His", CGU=>"Arg", CGC=>"Arg", CUA=>"Leu", CUG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", AUU=>"Ile", AUC=>"Ile", ACU=>"Thr", ACC=>"Thr", AAU=>"Asn", AAC=>"Asn", AGU=>"Ser", AGC=>"Ser", AUA=>"Met", ACA=>"Thr", AAA=>"Lys", AGA=>"*", AUG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"*", GUU=>"Val", GUC=>"Val", GCU=>"Ala", GCC=>"Ala", GAU=>"Asp", GAC=>"Asp", GGU=>"Gly", GGC=>"Gly", GUA=>"Val", GUG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonrfullm = (UUU=>"Phenylalanine", UUC=>"Phenylalanine", UCU=>"Serine", UCC=>"Serine", UAU=>"Tyrosine", UAC=>"Tyrosine", UGU=>"Cysteine", UGC=>"Cysteine", UUA=>"Leucine", UCA=>"Serine", UAA=>"Stop", UGA=>"Tryptophan", UUG=>"Leucine", UCG=>"Serine", UAG=>"Stop", UGG=>"Tryptophan", CUU=>"Leucine", CUC=>"Leucine", CCU=>"Proline", CCC=>"Proline", CAU=>"Histidine", CAC=>"Histidine", CGU=>"Arginine", CGC=>"Arginine", CUA=>"Leucine", CUG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", AUU=>"Isoleucine", AUC=>"Isoleucine", ACU=>"Threonine", ACC=>"Threonine", AAU=>"Asparagine", AAC=>"Asparagine", AGU=>"Serine", AGC=>"Serine", AUA=>"Methionine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Stop", AUG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Stop", GUU=>"Valine", GUC=>"Valine", GCU=>"Alanine", GCC=>"Alanine", GAU=>"Aspartic acid", GAC=>"Aspartic acid", GGU=>"Glycine", GGC=>"Glycine", GUA=>"Valine", GUG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine");

our %iupac = (R=>'AG', Y=>'CT', S=>'GC', W=>'AT', K=>'GT', M=>'AC', A=>'AA', C=>'CC', G=>'GG', T=>'TT', B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT', '.'=>'-', '-'=>'-');
our $aamatrix;				#amino acid substitution matrix


processArguments ();			#process program arguments, set up default values, check for errors, check for existence of db files
    
# Padded output
my $pad_fh; # Write seq pad here
my $cDNA_pad;
if ($seq_padding) {
    open $pad_fh, ">$outfile.seqpad" or die "Error: cannot write to output file $outfile.seqpad: $!\n";
    $cDNA_pad = $seq_padding * 3;
}

if ($geneanno) {
	annotateQueryByGene ();		#generate gene-based annoations (classify variants into intergenic, introgenic, non-synonymous, synonymous, UTR, frameshift, etc)
}  elsif ($downdb) {
	downloadDB ();			#download annotation databases from Internet
}


sub processArguments {
	my @command_line = @ARGV;		#command line argument
	GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'separate'=>\$separate,
	'batchsize=s'=>\$batchsize, 'dbtype=s'=>\$dbtype, 'neargene=i'=>\$neargene, 'genomebinsize=s'=>\$genomebinsize,
	'geneanno'=>\$geneanno, 'regionanno'=>\$regionanno, , 'filter'=>\$filter, 'downdb'=>\$downdb, 'buildver=s'=>\$buildver, 'score_threshold=f'=>\$score_threshold, 
	'normscore_threshold=i'=>\$normscore_threshold,	'minqueryfrac=f'=>\$minqueryfrac, 'expandbin=i'=>\$expandbin, 'splicing_threshold=i'=>\$splicing_threshold,
	'maf_threshold=f'=>\$maf_threshold, 'chromosome=s'=>\$chromosome, 'zerostart'=>\$zerostart, 'rawscore'=>\$rawscore, 'memfree=i'=>\$memfree, 
	'memtotal=i'=>\$memtotal, 'sift_threshold=f'=>\$sift_threshold, 'gff3dbfile=s'=>\$gff3dbfile, 'genericdbfile=s'=>\$genericdbfile, 'vcfdbfile=s'=>\$vcfdbfile,
	'time'=>\$time, 'wget!'=>\$wget, 'precedence=s'=>\$precedence, 'webfrom=s'=>\$webfrom, 'colsWanted=s'=>\$colsWanted, 'comment'=>\$comment,
	'scorecolumn=i'=>\$scorecolumn, 'transcript_function'=>\$transfun, 'exonsort'=>\$exonsort, 'avcolumn=s'=>\$avcolumn, 'bedfile=s'=>\$bedfile,
	'hgvs'=>\$hgvs, 'reverse'=>\$reverse, 'indexfilter_threshold=f'=>\$indexfilter_threshold, 'otherinfo'=>\$otherinfo, 
	'seq_padding=i'=>\$seq_padding, 'indel_splicing_threshold=i'=>\$indel_splicing_threshold, 'infoasscore'=>\$infoasscore, 'firstcodondel!'=>\$firstcodondel,
	'aamatrixfile=s'=>\$aamatrixfile, 'gff3attribute'=>\$gff3attr, 'exonicsplicing'=>\$exonicsplicing, 'infosep'=>\$infosep) or pod2usage ();
	
	$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
	$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
	@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
	@ARGV == 2 or pod2usage ("Syntax error");

	($queryfile, $dbloc) = @ARGV;

	$dbloc =~ s/[\\\/]$//;			#delete the trailing / or \ sign as part of the directory name
	if (defined $batchsize) {
		$batchsize =~ s/k$/000/;
		$batchsize =~ s/m$/000000/;
		$batchsize =~ m/^\d+$/ or pod2usage ("Error: the --batchsize argument must be a positive integer (suffix of k or m is okay)");
	} else {
		$batchsize = 5_000_000;
	}
	if (defined $genomebinsize) {
		$genomebinsize =~ s/k$/000/;
		$genomebinsize =~ s/m$/000000/;
		$genomebinsize =~ m/^\d+$/ or pod2usage ("Error: the --genomebinsize argument must be a positive integer (suffix of k or m is okay)");
		$genomebinsize > 1000 or pod2suage ("Error: the --genomebinsize argument must be larger than 1000");
	} else {
		if ($geneanno) {
			$genomebinsize = 100_000;		#gene usually span large genomic regions
		} else {
			$genomebinsize = 10_000;		#MCE, TFBS, miRNA, etc are small genomic regions
		}
	}

	$verbose ||= 0;			#when it is not specified, it is zero
	$neargene ||= 1_000;		#for upstream/downstream annotation of variants, specify the distance threshold between variants and genes
	$expandbin ||= int(2_000_000/$genomebinsize);		#for gene-based annotations, when intergenic variants are found, expand to specified number of nearby bins to find closest genes
	$outfile ||= $queryfile;	#specify the prefix of output file names

	#set up log file
	if ($downdb) {
		if (not -d $dbloc) {
			mkdir ($dbloc) or die "Error: the directory $dbloc does not exist and cannot be created\n";
		}
		my $errfile = File::Spec->catfile ($dbloc, "annovar_downdb.log");
		open (LOG, ">$errfile") or die "Error: cannot write LOG information to log file $errfile: $!\n";
	} else {
		open (LOG, ">$outfile.log") or die "Error: cannot write LOG information to log file $outfile.log: $!\n";
	}
	print LOG "ANNOVAR Version:\n\t", q/$Date$/, "\n";
	print LOG "ANNOVAR Information:\n\tFor questions, comments, documentation, bug reports and program update, please visit http://www.openbioinformatics.org/annovar/\n";
	print LOG "ANNOVAR Command:\n\t$0 @command_line\n";
	print LOG "ANNOVAR Started:\n\t", scalar (localtime), "\n";
	
	my $num = 0;
	$geneanno and $num++;
	$downdb and $num++;
	$filter and $num++;
	$regionanno and $num++;
	$num <= 1 or pod2usage ("Error in argument: please specify only one of --geneanno, -regionanno, --downdb, --filter");
	if (not $num) {
		$geneanno++;
		printerr "NOTICE: The --geneanno operation is set to ON by default\n";
	}
	
	my %dbtype1 = ('gene'=>'refGene', 'refgene'=>'refGene', 'knowngene'=>'knownGene', 'ensgene'=>'ensGene', 'band'=>'cytoBand', 'cytoband'=>'cytoBand', 'tfbs'=>'tfbsConsSites', 'mirna'=>'wgRna',
			'mirnatarget'=>'targetScanS', 'segdup'=>'genomicSuperDups', 'omimgene'=>'omimGene', 'gwascatalog'=>'gwasCatalog', 
			'1000g_ceu'=>'CEU.sites.2009_04', '1000g_yri'=>'YRI.sites.2009_04', '1000g_jptchb'=>'JPTCHB.sites.2009_04', 
			'1000g2010_ceu'=>'CEU.sites.2010_03', '1000g2010_yri'=>'YRI.sites.2010_03', '1000g2010_jptchb'=>'JPTCHB.sites.2010_03',
			'1000g2010jul_ceu'=>'CEU.sites.2010_07', '1000g2010jul_yri'=>'YRI.sites.2010_07', '1000g2010jul_jptchb'=>'JPTCHB.sites.2010_07',
			'1000g2010nov_all'=>'ALL.sites.2010_11', '1000g2011may_all'=>'ALL.sites.2011_05'
			);

		
	if ($geneanno) {
		$dbtype ||= 'refGene';
		$dbtype1 = $dbtype1{$dbtype} || $dbtype;
	} elsif ($regionanno) {
		defined $dbtype or pod2usage ("Error in argument: please specify --dbtype (required for the --regionanno operation)");
		$dbtype1 = $dbtype1{$dbtype} || $dbtype;
		if ($dbtype =~ m/^mce(\d+)way/) {			#added 2010Feb16
			$dbtype1 = "phastConsElements$1way";
		}
		if ($dbtype1 eq 'gff3') {
			defined $gff3dbfile or pod2usage ("Error in argument: please specify --gff3dbfile for the --dbtype of 'gff3'");
		}
	} elsif ($filter) {
		defined $dbtype or pod2usage ("Error in argument: please specify --dbtype (required for the --filter operation)");
		#as of Feb 2012, I no longer check the validity of the database name for -filter operation, to give users the maximum amount of flexibility in designing and using their own favorite databases
		$dbtype =~ m/^avsift|generic|1000g_(ceu|yri|jptchb)|1000g2010_(ceu|yri|jptchb)|1000g20\d\d[a-z]{3}_[a-z]+|snp\d+\w+?|vcf|(ljb_[\w\+]+)|esp\d+_[\w]+$/ or print STDERR "NOTICE: the --dbtype $dbtype is assumed to be in generic ANNOVAR database format\n";
		
		$dbtype1 = $dbtype1{$dbtype} || $dbtype;
		
		if ($dbtype1 =~ m/^1000g(20\d\d)([a-z]{3})_([a-z]+)$/) {
			my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
			$dbtype1 = uc ($3) . '.sites.' . $1 . '_' . $monthhash{$2};
		}
		
		if ($dbtype1 eq 'generic') {
			defined $genericdbfile or pod2usage ("Error in argument: please specify --genericdbfile for the --dbtype of 'generic'");
		}
		if ($dbtype eq 'vcf') {
			defined $vcfdbfile or pod2usage ("Error in argument: please specify --vcfdbfile for the --dbtype of 'vcf'");
		}
	} elsif ($downdb) {
		defined $dbtype and pod2usage ("Error in argument: please do not specify --dbtype for the --downdb operation");
		$dbtype1 = $dbtype1{$queryfile} || $queryfile;
		if ($queryfile =~ m/^mce(\d+)way/) {			#added 2013may08
			$dbtype1 = "phastConsElements$1way";
		}
	}
	
	if (not $buildver) {
		$buildver = 'hg18';
		printerr "NOTICE: The --buildver is set as 'hg18' by default\n";
	}
	
	if (defined $score_threshold) {
		#$score_threshold >= 0 or pod2usage ("Error in argument: the --score_threshold must be a positive number or zero (you specified $score_threshold)");	#20130208: score_threshold can be anything now as long as it is a number
		$geneanno || $downdb and pod2usage ("Error in argument: the --score_threshold is not useful for --geneanno or --downdb operations");
	}
	if ($normscore_threshold) {
		$normscore_threshold >= 0 and $normscore_threshold <= 1000 or pod2usage ("Error in argument: the --normscore_threshold must be between 0 and 1000 (you specified $normscore_threshold)");
		$regionanno or pod2usage ("Error in argument: the --normscore_threshold is supported only for the --regionanno operation");
	}
	
	if ($zerostart) {
		pod2usage ("Error: the -zerostart argument is now obselete and will no longer be supported in ANNOVAR");
	}
	
	if (defined $sift_threshold) {
		$filter or pod2usage ("Error in argument: the --sift_threshold is supported only for the --filter operation");
		$dbtype1 eq 'avsift' or pod2usage ("Error in argument: the --sift_threshold argument can be used only if '--dbtype avsift' is used");
		$sift_threshold >= 0 and $sift_threshold <= 1 or pod2usage ("Error in argument: the --sift_threshold must be between 0 and 1 inclusive");
	} else {
		$dbtype1 eq 'avsift' and printerr "NOTICE: The --sift_threshold is set as 0.05 by default\n";
		$sift_threshold = 0.05;
	}
	
	if (defined $indexfilter_threshold) {
		$filter or pod2usage ("Error in argument: the --indexfilter_threshold is supported only for the --filter operation");
		$indexfilter_threshold >= 0 and $indexfilter_threshold <= 1 or pod2usage ("Error in argument: the --indexfilter_threshold must be between 0 and 1 inclusive");
	} else {
		$indexfilter_threshold = 0.9;
	}
	
	#operation-specific argument
	if (defined $splicing_threshold) {
		$geneanno or pod2usage ("Error in argument: the --splicing_threshold is supported only for the --geneanno operation");
	} else {
		$splicing_threshold = 2;	#for splicing annotation, specify the distance threshold between variants and exon/intron boundaries
	}
	if (defined $indel_splicing_threshold) {
		$geneanno or pod2usage ("Error: the --indel_splicing_threshold is supported only for the --geneanno operation");
	}
	else {
		$indel_splicing_threshold = $splicing_threshold;    #if not set, preserve original behavior;
	}
	if (defined $maf_threshold) {
		$filter or pod2usage ("Error in argument: the --maf_threshold is supported only for the --filter operation");
		$dbtype =~ m/^1000g/ or pod2usage ("Error in argument: the --maf_threshold is supported only for 1000 Genomes Project data set (try -score_threshold instead)");
	} else {
		$maf_threshold = 0;		#for filter-based annotations on 1000 Genomes Project data, specify the MAF threshold to be used in filtering
	}
	if (defined $minqueryfrac) {
		$regionanno or pod2usage ("Error in argument: the --minqueryfrac is supported only for the --regionanno operation");
	} else {
		$minqueryfrac = 0;		#minimum query overlap to declare a "match" with database records
	}
	if (defined $gff3dbfile) {
		$dbtype eq 'gff3' or pod2usage ("Error in argument: the --gff3dbfile argument can be used only if '--dbtype gff3' is used");
		$geneanno or $regionanno or pod2usage ("Error in argument: the --gff3dbfile argument is supported only for the --geneanno or --regionanno operation");
	}
	if (defined $bedfile) {
		$dbtype eq 'bed' or pod2usage ("Error in argument: the --bedfile argument can be used only if '--dbtype bed' is used");
		$regionanno or pod2usage ("Error in argument: the --bedfile argument is supported only for the --regionanno operation");
	}
	if (defined $genericdbfile) {
		$filter or pod2usage ("Error in argument: the --genericdbfile argument is supported only for the --filter operation");
	}
	if (defined $wget) {
		$downdb or pod2usage ("Error in argument: the --wget argument is supported only for the --downdb operation");
	} else {
		$wget = 1;			#by default, use wget for downloading files from Internet
	}
	if (defined $precedence) {
		$geneanno or pod2usage ("Error in argument: the --precedence argument is supported only for the --geneanno operation");
		@precedence = split (/,/, $precedence);
		@precedence >= 2 or pod2usage ("Error in argument: the --precedence argument should be comma delimited");
		for my $i (0 .. @precedence-1) {
			$precedence[$i] =~ m/^(exonic|intronic|splicing|utr5|utr3|upstream|downstream|splicing|ncrna)$/ or pod2usage ("Error in argument: the --precedence argument contains invalid keywords (valid ones are exonic|intronic|splicing|utr5|utr3|upstream|downstream|splicing)");
		}
	}
	
	if (defined $colsWanted) {
		$regionanno or $filter or pod2usage ("Error in argument: the --colWanted argument is supported only for the --regionanno and --filter operation");
		if (lc $colsWanted eq 'all') {
			@colsWanted = ('all');
		} elsif (lc $colsWanted eq 'none') {
			@colsWanted = ('none');
		} else {
			@colsWanted = split (/,/, $colsWanted);
			for my $i (0 .. @colsWanted-1) {
				$colsWanted[$i]=~m/^\d+$/ or pod2usage ("Error in argument: the --colsWanted argument ($colsWanted) must be a list of comma delimited numbers or be 'all' or be 'none'");
			}
		}
	}
	
	if (defined $scorecolumn) {
		$regionanno or pod2usage ("Error in argument: the --scorecolumn argument is supported only for the --regionanno operation");
	}
	
	if ($exonsort) {
		$geneanno or pod2usage ("Error in argument: the --exonsort argument is supported only for the --geneanno operation");
	}
	
	if (defined $avcolumn) {
		$avcolumn =~ m/^\d+,\d+,\d+,\d+,\d+$/ or pod2usage ("Error in argument: the --avcolumn argument must be five integer numbers separated by comma");
		@avcolumn = split (/,/, $avcolumn);
		@avcolumn = map {$_-1} @avcolumn;
	} else {
		@avcolumn = (0..4);		#by default, the first five columns are the required AVINPUT information
	}
	
	if (defined $webfrom) {
		if ($webfrom ne 'ucsc' and $webfrom ne 'annovar') {
			$webfrom =~ m#^(http://|ftp://)# or pod2usage ("Error: the --webfrom argument needs to be 'ucsc', 'annovar', or a URL");
		}
	}
	
	$maf_threshold >= 0 and $maf_threshold <= 0.5 or pod2usage ("Error in argument: the --maf_threshold must be between 0 and 0.5 (you specified $maf_threshold)");
	$minqueryfrac >= 0 and $minqueryfrac <= 1 or pod2usage ("Error in argument: the --minqueryfrac must be between 0 and 1 (you specified $minqueryfrac)");
	$memfree and $memfree >= 100_000 || pod2usage ("Error in argument: the --memfree argument must be at least 100000 (in the order of kilobytes)");
	$memtotal and $memtotal >= 100_000 || pod2usage ("Error in argument: the --memtotal argument must be at least 100000 (in the order of kilobytes)");
	
	if ($chromosome) {
		my @chr = split (/,/, $chromosome);
		for my $i (0 .. @chr-1) {
			if ($chr[$i] =~ m/^(\d+)-(\d+)$/) {
				for my $j ($1 .. $2) {
					$valichr{$j}++;
				}
			} else {
				$valichr{$chr[$i]}++;
			}
		}
		printerr "NOTICE: These chromosomes in database will be examined: ", join (",", sort keys %valichr), "\n";
	}
	
	if (not defined $firstcodondel) {
		$firstcodondel = 1;
	}
	
	if (defined $aamatrixfile) {
		$geneanno or pod2usage ("Error in argument: the --aamatrix argument can be used only for gene-based annotation");
		$aamatrix = readAAMatrixFile ($aamatrixfile);
	}
	
	if ($gff3attr) {
		$dbtype eq 'gff3' or pod2usage ("Error in argument: the --gff3attr argument can be used only if '--dbtype gff3' is used");
	}
}

sub annotateQueryByGene {
	my ($queryfh);							#query file handle
	my ($totalquerycount, $totalinvalidcount, $batchcount) = qw/0 0 1/;
	open ($queryfh, $queryfile) or die "Error: cannot read from --queryfile ($queryfile): $!\n";
	
	open (OUT, ">$outfile.variant_function") or die "Error: cannot write to output file $outfile.variant_function: $!\n";
	open (EXONIC, ">$outfile.exonic_variant_function") or die "Error: cannot write to output file $outfile.exonic_variant_function: $!\n";
	open (INVALID, ">$outfile.invalid_input") or die "Error: cannot write to output file $outfile.invalid_input: $!\n";

	my ($genedb, $geneidmap, $cdslen, $mrnalen) = readUCSCGeneAnnotation ($dbloc);
	
	$time and printerr "NOTICE: Current time (before examining variants) is ", scalar (localtime), "\n";
	while (1) {
		my ($linecount, $invalidcount) = newprocessNextQueryBatchByGene ($queryfh, $batchcount, $batchsize, $genedb, $geneidmap, $cdslen, $mrnalen);
		$totalquerycount += $linecount;
		$totalinvalidcount += $invalidcount;
		$linecount == $batchsize or last;
		$batchcount++;
		printerr "NOTICE: Begin processing batch $batchcount (each batch contains $batchsize variants)\n";
	}
	close (INVALID);
	close (EXONIC);
	close (OUT);
	close ($queryfh);
	$time and printerr "NOTICE: Current time (after examining variants) is ", scalar (localtime), "\n";

	$totalinvalidcount or unlink ("$outfile.invalid_input");	#delete the file as it is empty
	printerr "NOTICE: Finished gene-based annotation on $totalquerycount genetic variants in $queryfile";
	$totalinvalidcount and printerr " (including $totalinvalidcount with invalid format written to $outfile.invalid_input)";
	printerr "\n";
	printerr "NOTICE: Output files were written to $outfile.variant_function, $outfile.exonic_variant_function\n";
}

sub newprocessNextQueryBatchByGene {
	my ($queryfh, $batchcount, $batchsize, $genedb, $geneidmap, $cdslen, $mrnalen) = @_;
	my (%refseqvar);
	
	my ($chr, $start, $end, $ref, $obs);
	my ($name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend, $name2);
	my ($invalid);
	my ($linecount, $invalidcount) = qw/0 0/;
	
	for my $i (1 .. $batchsize) {					#process up to batchsize variants
		my $nextline = <$queryfh>;				#read the next line in variant file
		defined $nextline or last;
		$nextline =~ s/[\r\n]+$//;
		
		if ($nextline =~ m/^#/ and $comment) {			#comment line start with #, do not include this is $linecount
			print OUT "#comment\t$comment\t$nextline\n";
			next;
		}
		
		$linecount++;						#linecount does not include the comment line
		$invalid = 0;
		
		my @nextline = split (/\s+/, $nextline);
		($chr, $start, $end, $ref, $obs) = @nextline[@avcolumn];
		if ( not (defined $chr and defined $start and defined $end and defined $ref and defined $obs)) {
			$invalid++;
		} else {
			($ref, $obs) = (uc $ref, uc $obs);
			$zerostart and $start++;
			$chr =~ s/^chr//;
			if ($chr =~ m/[^\w\.]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {		#chr name could contain . (example: GL000212.1, or Zv9_NA###
				$invalid++;
			} elsif ($ref eq '-' and $obs eq '-' 		#both are empty allele
				or $ref =~ m/[^ACTG0\-]/ 		#non-standard nucleotide code
				or $obs =~ m/[^ACGT0\-]/ 		#non-standard nucleotide code
				or $start =~ m/[^\d]/ 			#start is not a number
				or $end =~ m/[^\d]/ 			#end is not a number
				or $start > $end			#start is more than end
				or $ref ne '0' and $end-$start+1 != length ($ref) 	#length mismatch with ref
				or $ref eq '-' and $start != $end	#length mismatch for insertion
				) {
				$invalid++;
			}
		}



		if ($invalid) {
			print INVALID $nextline, "\n";			#invalid record found
			$invalidcount++;
			next;
		}
		
		my (%intronic, %utr5, %utr3, %exonic, %upstream, %downstream, %ncrna, %intergenic, %splicing, %splicing_anno);
		my $foundgenic;						#variant found in genic region (between start and end position of a gene in genome)
		my ($distl, $distr, $genel, $gener);			#for intergenic variant, the distance and gene name to the left and right side of gene
		my $bin1 = int ($start/$genomebinsize)-1;		#start bin
		$bin1 < 0 and $bin1=0;
		my $bin2 = int ($end/$genomebinsize)+1;			#end bin (usually same as start bin, unless the query is really big that spans multiple megabases)
		
		while (not exists $genedb->{$chr, $bin1} and $bin1 > int ($start/$genomebinsize)-$expandbin) {		#examine at least 5 bins (by default 5Mb) to the left to make sure that a gene is found in the bin
			$bin1 > 0 or last;
			$bin1--;
		}
		
		while (not exists $genedb->{$chr, $bin2} and $bin2 < int ($end/$genomebinsize)+$expandbin) {		#examine at least 5 bins (by default 5Mb) to the right to make sure that a gene is found in the bin
			$bin2++;
		}

		my (%seen);
		for my $nextbin ($bin1 .. $bin2) {
			exists $genedb->{$chr, $nextbin} or next;		#this genome bin has no annotated gene (a complete intergenic region)
			for my $nextgene (@{$genedb->{$chr, $nextbin}}) {	#when $genedb->{$chr, $nextbin} is undefined, this automatically create an array!!!
				($name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend, $name2) = @$nextgene;
				defined $name2 or printerr "WARNING: name2 field is not provided for transcript $name (start=$txstart end=$txend)\n" and $name2='';
				$seen{$name, $txstart} and next;		#name and txstart uniquely identify a transcript and chromosome position (sometimes same transcript may map to two nearby positions, such as nearby segmental duplications)
				$seen{$name, $txstart}++;			#a transcript may be in two adjacent bins, so if one is already scanned, there is no need to work on it again
				my $current_ncRNA;
				
				if ($transfun) {	#variant_function output contains transcript name, rather than gene name
					$name2 = $name;
				}
				
				if (not $foundgenic) {				#this variant has not hit a genic region yet
					if ($start > $txend) {
						defined $distl or $distl = $start-$txend and $genel=$name2;
						$distl > $start-$txend and $distl = $start-$txend and $genel=$name2;	#identify left closest gene
					}
		
					if ($end < $txstart) {
						defined $distr or $distr = $txstart-$end and $gener=$name2;
						$distr > $txstart-$end and $distr = $txstart-$end and $gener=$name2;	#identify right closest gene
					}
				}
				
				if ($end < $txstart) {
					#query ---
					#gene		<-*----*->
					$foundgenic and last;					#if found a genic annotation already, end the search of the bins
					if ($end > $txstart - $neargene) {
						if ($dbstrand eq '+') {
							$upstream{$name2}++;
						} else {
							$downstream{$name2}++;
						}
					} else {
						last;						#if transcript is too far away from end, end the search of the bins
					}
				} elsif ($start > $txend) {
					#query            ---
					#gene  <-*----*->
					if (not $foundgenic and $start < $txend + $neargene) {
						if ($dbstrand eq '+') {
							$downstream{$name2}++;
						} else {
							$upstream{$name2}++;
						}
					}
				#} elsif ($cdsstart == $cdsend+1) {				#non-coding RNA (could be microRNA, or could be due to lack of CDS annotation for mRNA such as NR_026730 or BC039000). Previously we already did cdsstart++ so here the cdsstart is more than cdsend
				#	if ($start >= $txstart and $start <= $txend or $end >= $txstart and $end <= $txend or $start <= $txstart and $end >= $txend) {
				#		$ncrna{$name2}++;
				#		$foundgenic++;
				#	}
				} else {							#query overlaps with coding region of gene
					
					#change in 2011jul24: handle ncRNA and protein coding gene together but with the ncRNA flag when printing out results in the future
					if ($cdsstart == $cdsend+1) {				#non-coding RNA (could be microRNA, or could be due to lack of CDS annotation for mRNA such as NR_026730 or BC039000). Previously we already did cdsstart++ so here the cdsstart is more than cdsend
						if ($start >= $txstart and $start <= $txend or $end >= $txstart and $end <= $txend or $start <= $txstart and $end >= $txend) {
							$ncrna{$name2}++;
							$foundgenic++;
						}
						
						#now treat this ncRNA as if it is a protein-coding gene
						($cdsstart, $cdsend) = ($txstart, $txend);
						$current_ncRNA++;		#current transcript is a noncoding transcript
					}
					
					
					
					my ($lenintron, $lenexon) = (0, 0);			#cumulative intron/exon length at a given exon
					my ($rcdsstart, $rvarstart, $rvarend);			#start of coding and variant in reference mRNA sequence
					my @exonstart = @$exonstart;
					my @exonend = @$exonend;
					my $foundexonic;
					if ($dbstrand eq '+') {					#forward strand, search from left to right (first exon to last exon)
						for my $k (0 .. @exonstart-1) {
							$k and $lenintron += ($exonstart[$k]-$exonend[$k-1]-1);		#calculate cumulative intron length
							$lenexon += ($exonend[$k]-$exonstart[$k]+1);
							if ($cdsstart >= $exonstart[$k]) {				#calculate CDS start accurately by considering intron length
								$rcdsstart = $cdsstart-$txstart-$lenintron+1;
								
								if ($cdsstart <= $exonend[$k]) {	#CDS start is within this exon
									$lenexon = ($exonend[$k]-$cdsstart+1);
								} else {				#CDS start is in previous exon
									#$lenexon += ($exonend[$k]-$exonstart[$k]+1);
								}
							}
							
							if ($exonicsplicing) {			#when -exonicsplicing argument is set, exonic variants near exon/intron boundary is reported as "exonic;splicing" variant
								#splicing calculation (changed 2012may24)
								if (@exonstart != 1) {
									if ( $k == 0 and $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1 or $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query start site is close to exon start or exon end
										}
									}
									
									if ($k == 0 and $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1 or $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query end site is close to exon start or exon end
										}
									}
									
									if ($k == 0 and $start <= $exonend[$k] and $end >= $exonend[$k]) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $start <= $exonstart[$k] and $end>=$exonstart[$k]) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($start <= $exonstart[$k] and $end>=$exonstart[$k] or $start <= $exonend[$k] and $end >= $exonend[$k]) {
											$splicing{$name2}++;		#when query encompass the exon/intron boundary
										}
									}
								}
							} else {
								#splicing calculation (changed 2013feb10)
								#splicing calculation (changed again 2013feb21 per Mitsuhiro Komura suggestion)
								if (@exonstart != 1) {
									if ( $k == 0 and $start >= $exonend[$k]+1 and $start <= $exonend[$k]+$splicing_threshold) {			#first exon
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]-1) {	#last exon
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {										#middle exon
										if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]-1 or $start >= $exonend[$k]+1 and $start <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query start site is close to exon start or exon end
										}
									}
								}
							}
							
							
							
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if this splicing leads to amino acid change (rather than UTR change)
							if ($splicing{$name2} and $start==$end and $start>=$cdsstart) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$lenexon -= ($exonend[$k]-$exonstart[$k]);		#formerly "$exonend[$k]-$exonstart[$k]+1"; changed this line 2011oct01 given German's comments. The coding portion for acceptor site should be added by one, but for donor site should be fine
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:c.$lenexon-" . ($exonstart[$k]-$start) . "$ref>$obs,";
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:c.$lenexon+" . ($start-$exonend[$k]) . "$ref>$obs,";
								}
							}
							
							if ($start < $exonstart[$k]) {
								if ($end >= $exonstart[$k]) {	#exonic 
									$rvarstart = $exonstart[$k]-$txstart-$lenintron+1;
									
									for my $m ($k .. @exonstart-1) {
										$m > $k and $lenintron += ($exonstart[$m]-$exonend[$m-1]-1);
										if ($end < $exonstart[$m]) {
											#query           --------
											#gene     <--**---******---****---->
											$rvarend = $exonend[$m-1]-$txstart-$lenintron+1 + ($exonstart[$m]-$exonend[$m-1]-1);
											last;
										} elsif ($end <= $exonend[$m]) {
											#query           -----------
											#gene     <--**---******---****---->
											$rvarend = $end-$txstart-$lenintron+1;
											last;
										}
									}
									if (not defined $rvarend) {
										$rvarend = $txend-$txstart-$lenintron+1;		#if this value is longer than transcript length, it suggest whole gene deletion
									}
									
									#here the trick begins to differentiate UTR versus coding exonic
									if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
										#query  ----
										#gene     <--*---*->
										$utr5{$name2}++;		#positive strand for UTR5
									} elsif ($start > $cdsend) {
										#query             ----
										#gene     <--*---*->
										$utr3{$name2}++;		#positive strand for UTR3
									} else {									
										$exonic{$name2}++;
										not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '+', $i, $k+1, $nextline];	#refseq CDS start, refseq variant start. obs is non-zero (obs is specified by user)
									}
									$foundgenic++;
									last;
								} elsif ($k and $start > $exonend[$k-1]) {	#intronic
									$intronic{$name2}++;
									$foundgenic++;
									last;
								}
							} elsif ($start <= $exonend[$k]) {	#exonic
								$rvarstart = $start-$txstart-$lenintron+1;
								
								for my $m ($k .. @exonstart-1) {
									$m > $k and $lenintron += ($exonstart[$m]-$exonend[$m-1]-1);
									if ($end < $exonstart[$m]) {
										#query              ------
										#gene     <--**---******---****---->
										$rvarend = $exonend[$m-1]-$txstart-$lenintron+1 + ($exonstart[$m]-$exonend[$m-1]-1);
										last;
									} elsif ($end <= $exonend[$m]) {
										#query           -----------
										#gene     <--**---******---****---->
										$rvarend = $end-$txstart-$lenintron+1;
										last;
									}
								}
								if (not defined $rvarend) {
									$rvarend = $txend-$txstart-$lenintron+1;		#if this value is longer than transcript length, it suggest whole gene deletion
								}
								
								#here is the trick begins to differentiate UTR versus coding exonic
								if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
									#query  ----
									#gene     <--*---*->
									$utr5{$name2}++;		#positive strand for UTR5
								} elsif ($start > $cdsend) {
									#query             ----
									#gene     <--*---*->
									$utr3{$name2}++;		#positive strand for UTR3
								} else {
									$exonic{$name2}++;
									not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '+', $i, $k+1, $nextline];		#queryindex, refseq CDS start, refseq variant start
								}
								$foundgenic++;
								last;
							}
						}
					} elsif ($dbstrand eq '-') {		#process negative strand (in the future, this should be fused to the paragraph above for positive strands; for now, I keep them separate for easier debugging)
						for (my $k = @exonstart-1; $k>=0; $k--) {
							$k < @exonstart-1 and $lenintron += ($exonstart[$k+1]-$exonend[$k]-1);
							$lenexon += ($exonend[$k]-$exonstart[$k]+1);
							if ($cdsend <= $exonend[$k]) {		#calculate CDS start accurately by considering intron length
								$rcdsstart = $txend-$cdsend-$lenintron+1;
								
								if ($cdsend >= $exonstart[$k]) {	#CDS start within this exon
									$lenexon = ($cdsend-$exonstart[$k]+1);
								} else {				#CDS start in prevous exon
									#$lenexon += ($exonend[$k]-$exonstart[$k]+1);
								}
							}

							if ($exonicsplicing) {			#when -exonicsplicing argument is set, exonic variants near exon/intron boundary is reported as "exonic;splicing" variant
								#splicing calculation (changed 2012may24)
								if (@exonstart != 1) {
									if ( $k == 0 and $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1 or $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query start site is close to exon start or exon end
										}
									}
									
									if ($k == 0 and $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1 or $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query end site is close to exon start or exon end
										}
									}
									
									if ($k == 0 and $start <= $exonend[$k] and $end >= $exonend[$k]) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $start <= $exonstart[$k] and $end>=$exonstart[$k]) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($start <= $exonstart[$k] and $end>=$exonstart[$k] or $start <= $exonend[$k] and $end >= $exonend[$k]) {
											$splicing{$name2}++;		#when query encompass the exon/intron boundary
										}
									}
								}
							} else {
								#splicing calculation (changed 2013feb10)
								#splicing calculation (changed again 2013feb21 per Mitsuhiro Komura suggestion)
								if (@exonstart != 1) {
									if ( $k == 0 and $start >= $exonend[$k]+1 and $start <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]-1 or $start >= $exonend[$k]+1 and $start <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query start site is close to exon start or exon end
										}
									}
								}
							}
							
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if this splicing leads to amino acid change (rather than UTR change)
							if ($splicing{$name2} and $start==$end and $start<=$cdsend) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-$k+1)}:c.$lenexon+" . ($exonstart[$k]-$start) . revcom($ref) . '>' . revcom ($obs) . ',';
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$lenexon -= ($exonend[$k]-$exonstart[$k]);	#formerly "$exonend[$k]-$exonstart[$k]+1"; changed this line 2011oct01 given German's comments. The coding portion for acceptor site should be added by one, but for donor site should be fine
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-$k+1)}:c.$lenexon-" . ($start-$exonend[$k]) . revcom($ref) . '>' . revcom($obs) . ',';
								}
							}
							
							if ($end > $exonend[$k]) {
								if ($start <= $exonend[$k]) {
									$rvarstart = $txend-$exonend[$k]-$lenintron+1;
									
									for (my $m = $k; $m >= 0; $m--) {
										$m < $k and $lenintron += ($exonstart[$m+1]-$exonend[$m]-1);
										if ($start > $exonend[$m]) {
											#query           --------
											#gene     <--**---******---****---->
											#$rvarend = $txend-$exonstart[$m]-$lenintron+1 - ($exonstart[$m+1]-$exonend[$m]-1);	#commented out 2011feb18
											$rvarend = $txend-$exonstart[$m+1]+1-$lenintron + ($exonstart[$m+1]-$exonend[$m]-1);	#fixed this 2011feb18
											last;		#finsih the cycle!!!!!!!!!!!!!!!!!!!
										} elsif ($start >= $exonstart[$m]) {		#start within exons
											#query               ----
											#gene     <--**---******---****---->
											$rvarend = $txend-$start-$lenintron+1;
											last;
										}
									}
									if (not defined $rvarend) {				#if rvarend is not found, then the whole tail of gene is covered
										$rvarend = $txend-$txstart-$lenintron+1;
									}
									
									#here is the trick begins to differentiate UTR versus coding exonic
									if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
										#query  ----
										#gene     <--*---*->
										$utr3{$name2}++;		#negative strand for UTR5
									} elsif ($start > $cdsend) {
										#query             ----
										#gene     <--*---*->
										$utr5{$name2}++;		#negative strand for UTR3
									} else {
										$exonic{$name2}++;
										not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '-', $i, @exonstart-$k, $nextline];
									}
									$foundgenic++;
									last;
								} elsif ($k < @exonstart-1 and $end < $exonstart[$k+1]) {
									$intronic{$name2}++;
									$foundgenic++;
									last;
								}
							} elsif ($end >= $exonstart[$k]) {
								$rvarstart = $txend-$end-$lenintron+1;		#all the rvarstart, rvarend are with respect to the cDNA sequence (so rvarstart corresponds to end of variants)
								
								for (my $m = $k; $m >= 0; $m--) {
									$m < $k and $lenintron += ($exonstart[$m+1]-$exonend[$m]-1);
									if ($start > $exonend[$m]) {
										#query           ----
										#gene     <--**---******---****---->
										#$rvarend = $txend-$exonstart[$m]-$lenintron+1 - ($exonstart[$m+1]-$exonend[$m]-1);		#commented out 2011feb18 due to bug (10 42244567 42244600 CACCTTTGCTTGATATGATAATATAGTGCCAAGG - hetero)
										$rvarend = $txend-$exonstart[$m+1]+1 - $lenintron + ($exonstart[$m+1]-$exonend[$m]-1);		#fixed this 2011feb18
										last;			#finish the circle of counting exons!!!!!
									} elsif ($start >= $exonstart[$m]) {			#the start is right located within exon
										#query        -------
										#gene     <--**---******---****---->
										$rvarend = $txend-$start-$lenintron+1;
										last;						#finish the cycle
									}
								}
								if (not defined $rvarend) {					#if rvarend is not found, then the whole tail of gene is covered
									$rvarend = $txend-$txstart-$lenintron+1;
								}
								
								#here the trick begins to differentiate UTR versus coding exonic
								if ($end < $cdsstart) {			#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
									#query  ----
									#gene     <--*---*->
									$utr3{$name2}++;		#negative strand for UTR5
								} elsif ($start > $cdsend or $start==$cdsend&&$ref eq '-') {	#insertions at the first base in negative strand is not a real insertion since it is not translated
									#query             ----
									#gene     <--*---*->
									$utr5{$name2}++;		#negative strand for UTR3
								} else {
									$exonic{$name2}++;
									not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '-', $i, @exonstart-$k, $nextline];
								}
								$foundgenic++;
								last;
							}
						}
					}
				}
			}
		}
		$foundgenic or $intergenic{''}++;		#changed $name2 to '' on 20110924
		$i =~ m/000000$/ and printerr "NOTICE: Finished analyzing $i query variants\n";

	
		my (@txname, %genename);
		my (%newsplicing);
		
		#process splicing annotation (change %splicing hash to %newsplicing, where gene name is replaced by gene name plus splicing annotation)
		if (%splicing) {
			if ($end-$start+1<=$indel_splicing_threshold) {		#make sure that long indel are not considered here
				for my $tempname (keys %splicing) {
					if ($splicing_anno{$tempname}) {
						$splicing_anno{$tempname} =~ s/,$//;	#remove the trailing comma
						$tempname .= "($splicing_anno{$tempname})";
					}
					$newsplicing{$tempname}++;
				}
			} else {
				%newsplicing = %splicing;
			}
		}
		
		
		if ($separate) {		#separately print out each effect on one line
			if (%exonic or %splicing or %intronic or %utr5 or %utr3 or %ncrna or %upstream or %downstream) {
				#if (%ncrna) {
				#	%exonic and print OUT "ncRNA_exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#	%splicing and $end-$start+1<=$splicing_threshold and print OUT "ncRNA_splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#	%intronic and print OUT "ncRNA_intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
				#	
				#	for my $key (keys %ncrna) {
				#		delete $exonic{$key};
				#		delete $splicing{$key};
				#		delete $intronic{$key};
				#	}
				#}
				if (%exonic) {
					my (@coding, @noncoding);
					for my $key (keys %exonic) {
						if ($ncrna{$key}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
					@coding and print OUT "exonic\t", join(",", sort @coding), "\t", $nextline, "\n";
					@noncoding and print OUT "ncRNA_exonic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
				#changed per Oscar 20131109
				#if (%splicing and $end-$start+1<=$indel_splicing_threshold) {
				#	my (@coding, @noncoding);
				#	for my $key (keys %splicing) {
				#		if ($ncrna{$key}) {
				if (%splicing) {
					my (@coding, @noncoding);
					for my $key (keys %newsplicing) {
						$key =~ m/^([^\(]+)/;
						if ($ncrna{$1}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
					@coding and print OUT "splicing\t", join(",", sort @coding), "\t", $nextline, "\n";
					@noncoding and print OUT "ncRNA_splicing\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
				if (%intronic) {
					my (@coding, @noncoding);
					for my $key (keys %intronic) {
						if ($ncrna{$key}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
					@coding and print OUT "intronic\t", join(",", sort @coding), "\t", $nextline, "\n";
					@noncoding and print OUT "ncRNA_intronic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
				
				#the following paragraph is commented out on 2011oct02
				#%exonic and print OUT "exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#%splicing and $end-$start+1<=$splicing_threshold and print OUT "splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#%intronic and print OUT "intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
				%utr5 and print OUT "UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
				%utr3 and print OUT "UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
				#if (%ncrna) {
				#	if (%exonic) {
				#		print OUT "ncRNA_exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#	}
				#	if (%splicing and $end-$start+1<=$splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
				#		print OUT "ncRNA_splicing\t",join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				#	}
				#	if (%utr5) {	#ncRNA should not have UTR5 or UTR3. If such an output exists, then there is a bug that should be reported and debugged!!!!!!!!!!!!!!!!
				#		print OUT "ncRNA_UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
				#	}
				#	if (%utr3) {	#ncRNA should not have UTR5 or UTR3. If such an output exists, then there is a bug that should be reported and debugged!!!!!!!!!!!!!!!!
				#		print OUT "ncRNA_UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
				#	}
				#	if (%intronic) {
				#		print OUT "ncRNA_intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
				#	}
				#}
				%upstream and print OUT "upstream\t", join(",", sort keys %upstream), "\t", $nextline, "\n";
				%downstream and print OUT "downstream\t", join(",", sort keys %downstream), "\t", $nextline, "\n";

			} elsif (%intergenic) {
				$genel ||= "NONE";
				$gener ||= "NONE";
				$distl ||= "NONE";
				$distr ||= "NONE";
				print OUT "intergenic\t", "$genel(dist=$distl),$gener(dist=$distr)", "\t", $nextline, "\n";
			} else {
				die "FATAL ERROR: please report bug to ANNOVAR author with your input file\n";
			}
		} else {			
			if (@precedence) {
				my $foundmatch;
				for my $i (0 .. @precedence-2) {
					$precedence[$i] eq 'exonic' and %exonic and $foundmatch++;
					$precedence[$i] eq 'splicing' and %splicing and $foundmatch++;
					$precedence[$i] eq 'intronic' and %intronic and $foundmatch++;
					$precedence[$i] eq 'utr5' and %utr5 and $foundmatch++;
					$precedence[$i] eq 'utr3' and %utr3 and $foundmatch++;
					$precedence[$i] eq 'ncrna' and %ncrna and $foundmatch++;
					$precedence[$i] eq 'upstream' and %upstream and $foundmatch++;
					$precedence[$i] eq 'downstream' and %downstream and $foundmatch++;
					$precedence[$i] eq 'intergenic' and %intergenic and $foundmatch++;
					if ($foundmatch) {
						for my $j ($i+1 .. @precedence-1) {
							$precedence[$j] eq 'exonic' and %exonic = ();
							$precedence[$j] eq 'splicing' and %splicing = ();
							$precedence[$j] eq 'intronic' and %intronic = ();
							$precedence[$j] eq 'utr5' and %utr5 = ();
							$precedence[$j] eq 'utr3' and %utr3 = ();
							$precedence[$j] eq 'ncrna' and %ncrna = ();
							$precedence[$j] eq 'upstream' and %upstream = ();
							$precedence[$j] eq 'downstream' and %downstream = ();
							$precedence[$j] eq 'intergenic' and %intergenic = ();
						}
						last;
					}
				}
			}
			
		
				
			if (%exonic) {
				my (@coding, @noncoding);
				for my $key (keys %exonic) {
					if ($ncrna{$key}) {
						push @noncoding, $key;
					} else {
						push @coding, $key;
					}
				}
				if (@coding and %splicing and $end-$start+1<=$indel_splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
					print OUT "exonic;splicing\t", join(",", sort @coding), ";", join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				} elsif (@coding) {
					print OUT "exonic\t", join(",", sort @coding), "\t", $nextline, "\n";
				} elsif (@noncoding) {
					print OUT "ncRNA_exonic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
			} elsif (%splicing) {
				my (@coding, @noncoding);
				for my $key (keys %newsplicing) {
					$key =~ m/^([^\(]+)/;
					if ($ncrna{$1}) {
						push @noncoding, $key;
					} else {
						push @coding, $key;
					}
				}
				if (@coding) {
					print OUT "splicing\t", join (',', sort @coding), "\t", $nextline, "\n";
				} elsif (@noncoding) {
					print OUT "ncRNA_splicing\t", join (',', sort @noncoding), "\t", $nextline, "\n";
				}
			} elsif (%ncrna) {
				my (@coding, @noncoding);
				for my $key (keys %intronic) {
					if ($ncrna{$key}) {
						push @noncoding, $key;
					} else {
						push @coding, $key;
					}
				}
				#print OUT "ncRNA\t", join(",", sort keys %ncrna), "\t", $nextline, "\n";
				
				#if (%exonic) {
				#	if (%splicing and $end-$start+1<=$splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
				#		print OUT "ncRNA_exonic;ncRNA_splicing\t", join(",", sort keys %exonic), ";", join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				#	} else {
				#		print OUT "ncRNA_exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#	}
				#} elsif (%splicing) {
				#	print OUT "ncRNA_splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#} elsif (%utr5 or %utr3) {		#ncRNA should not have UTR5 or UTR3. If such an output exists, then there is a bug that should be reported and debugged!!!!!!!!!!!!!!!!
				if (%utr5 or %utr3) {
					if (%utr5 and %utr3) {
						print OUT "ncRNA_UTR5;ncRNA_UTR3\t", join(",", sort keys %utr5), ";", join(",", sort keys %utr3), "\t", $nextline, "\n";		#use ";" to separate UTR5 and UTR3 genes
					} elsif (%utr5) {
						print OUT "ncRNA_UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
					} else {
						print OUT "ncRNA_UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
					}
				} elsif (@noncoding) {
					print OUT "ncRNA_intronic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				} else {
					die "FATAL ERROR: please report bug to ANNOVAR author with your input line <$nextline>\n";
				}
			} elsif (%utr5 or %utr3) {
				if (%utr5 and %utr3) {
					print OUT "UTR5;UTR3\t", join(",", sort keys %utr5), ";", join(",", sort keys %utr3), "\t", $nextline, "\n";		#use ";" to separate UTR5 and UTR3 genes
				} elsif (%utr5) {
					print OUT "UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
				} else {
					print OUT "UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
				}
			} elsif (%intronic) {
				print OUT "intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
			} elsif (%upstream or %downstream) {
				if (%upstream and %downstream) {
					print OUT "upstream;downstream\t", join(",", sort keys %upstream), ";", join(",", sort keys %downstream), "\t", $nextline, "\n";
				} elsif (%upstream) {
					print OUT "upstream\t", join(",", sort keys %upstream), "\t", $nextline, "\n";
				} else {
					print OUT "downstream\t", join(",", sort keys %downstream), "\t", $nextline, "\n";
				}
			} elsif (%intergenic) {
				$genel ||= "NONE";
				$gener ||= "NONE";
				$distl ||= "NONE";
				$distr ||= "NONE";
				print OUT "intergenic\t", "$genel(dist=$distl),$gener(dist=$distr)", "\t", $nextline, "\n";
			} else {
				die "FATAL ERROR: please report bug to ANNOVAR author with your input line <$nextline>\n";
			}
		}
	} 
	%refseqvar and annotateExonicVariants (\%refseqvar, $geneidmap, $cdslen, $mrnalen, $batchcount, $batchsize);

	return ($linecount, $invalidcount);
}

sub filterQuery {
	open (FIL, ">$outfile.${buildver}_${dbtype1}_filtered") or die "Error: cannot write to output file $outfile.${buildver}_${dbtype1}_filtered: $!\n"; 
	open (DROPPED, ">$outfile.${buildver}_${dbtype1}_dropped") or die "Error: cannot write to output file $outfile.${buildver}_${dbtype1}_dropped: $!\n";
	open (INVALID, ">$outfile.invalid_input") or die "Error: cannot write to output file $outfile.invalid_input: $!\n";
	
	printerr "NOTICE: Variants matching filtering criteria are written to $outfile.${buildver}_${dbtype1}_dropped, other variants are written to $outfile.${buildver}_${dbtype1}_filtered\n";
	
	open (QUERY, $queryfile) or die "Error: cannot read from query file $queryfile: $!\n";
	
	my (%variant, $filedone, $batchdone);
	my ($linecount, $batchlinecount, $invalid, $invalidcount) = (0, 0);
	my ($chr, $start, $end, $ref, $obs, $info);
	while (1) {
		$_ = <QUERY>;
		if (not defined $_) {
			$filedone++;
		} else {
			s/[\r\n]+$//;
			
			if (m/^#/ and $comment) {				#comment line start with #, do not include this is $linecount
				print FIL "$_\n";
				print DROPPED "#comment\t#comment\t$_\n";
				next;
			}
			
			$linecount++;
			$batchlinecount++;
			if ($batchlinecount == $batchsize) {
				$batchdone++;
			}
			
			if ($memfree or $memtotal) {		#if these arguments are specified
				if ($linecount =~ m/00000$/) {						#about 40Mb memory per 10k lines for a typical input dataset
					my ($availmem, $allmem) = currentAvailMemory();
					$verbose and printerr "NOTICE: Current available system memory is $availmem kb (this program uses $allmem bytes memory), after reading $linecount query\n";
					if ($availmem and $availmem <= $memfree+50_000) {		#some subsequent steps may take ~50Mb memory, so here we try to allocate some more memory
						$batchdone++;
					}
					if ($memtotal and $allmem >= $memtotal-50_000) {	#when --memtotal is specified, ensure that program use less memory
						$batchdone++;
					}
				}
			}
	
			$invalid = 0;						#reset invalid status
			
			my @nextline = split (/\s+/, $_);
			($chr, $start, $end, $ref, $obs) = @nextline[@avcolumn];
			if ( not (defined $chr and defined $start and defined $end and defined $ref and defined $obs)) {
				$invalid++;
			} else {
				($ref, $obs) = (uc $ref, uc $obs);
				$zerostart and $start++;
				$chr =~ s/^chr//;
				if ($chr =~ m/[^\w\.]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {
					$invalid++;
				} elsif ($ref eq '-' and $obs eq '-' 		#both are empty allele
					or $ref =~ m/[^ACTG0\-]/ 		#non-standard nucleotide code
					or $obs =~ m/[^ACGT0\-]/ 		#non-standard nucleotide code
					or $start =~ m/[^\d]/ 			#start is not a number
					or $end =~ m/[^\d]/ 			#end is not a number
					or $start > $end			#start is more than end
					or $ref ne '0' and $end-$start+1 != length ($ref) 	#length mismatch with ref
					or $ref eq '-' and $start != $end	#length mismatch for insertion
					) {		
					$invalid++;
				}
			}				
			
			if ($invalid) {
				print INVALID $_, "\n";	#invalid record found
				$invalidcount++;
				next;
			}
			
			if ($start == $end and $ref eq '-') {	#insertion
				$obs = "0$obs";
			} elsif ($obs eq '-') {			#deletion
				$obs = $end-$start+1;
			} elsif ($end>$start or $start==$end and length($obs)>1) {	#block substitution	#fixed the bug here 2011feb19
				$obs = ($end-$start+1) . $obs;
			}
			
			if (exists $variant{$chr, $start, $obs}) {
				$variant{$chr, $start, $obs} .= "\n$_";
			} else {
				$variant{$chr, $start, $obs} = "$ref\n$_";	#store the information for the variant in newline-separated strings, with the first field being ref allele
			}
		}
		
		if ($filedone or $batchdone) {
			printerr "NOTICE: Processing next batch with ${\(scalar keys %variant)} unique variants in $batchlinecount input lines\n";
			filterNextBatch (\%variant);			#filter the next batch of variants (print output to FIL, DROPPED), then clear %variant hash for next batch
			%variant = ();
			$batchlinecount = 0;				#reset the line count for this batch
			$batchdone = 0;
		}
		if ($filedone) {
			last;
		}
	}
	close (INVALID); close (DROPPED); close (FIL);
	if ($invalidcount) {
		printerr "NOTICE: Variants with invalid input format were written to $outfile.invalid_input\n";
	} else {
		unlink ("$outfile.invalid_input");
	}
}

sub revcom {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/acgtACGT/tgcaTGCA/;
	return ($seq);
}  

sub readUCSCGeneAnnotation {			#read RefGene annotation database from the UCSC Genome Browser, convert 0-based coordinates to 1-based coordinates
	my ($dbloc) = @_;
	my ($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes);
	my (%genedb, %geneidmap, %name2count, %cdslen, %mrnalen);
	my ($genecount, $ncgenecount) = (0, 0);
	
	my $dbfile;
	my $kgxref;
	my %iscoding;		#this gene name is a coding gene (if it has coding and noncoding transcripts, ignore all noncoding transcripts)
	
	if ($dbtype1 eq 'refGene') {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");
	} elsif ($dbtype1 eq 'knownGene') {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");
		my $kgxreffile = File::Spec->catfile($dbloc, $buildver . "_kgXref.txt");
		-f $kgxreffile or die "Error: the knownGene cross-reference file $kgxreffile does not exist. Please use 'annotate_variation.pl --downdb knownGene $dbloc' to download the database.\n";
		$kgxref = readKgXref ($kgxreffile);
	} elsif ($dbtype1 eq 'ensGene') {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");
	} else {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");		#added 2011feb18
		#die "FATAL ERROR: the dbype $dbtype1 is not supported in the readUCSCGeneAnnotation() subroutine.\n";		#commented 2011feb18
	}
	-f $dbfile or die "Error: The gene annotation database $dbfile does not exist. Please use 'annotate_variation.pl --downdb $dbtype $dbloc -build $buildver' to download the database.\n";

	open (GENEDB, $dbfile) or die "Error: cannot read from gene annotaion database $dbfile: $!\n";
	printerr "NOTICE: Reading gene annotation from $dbfile ... ";
	while (<GENEDB>) {
		s/[\r\n]+$//;							#deleting the newline characters
		my @record = split (/\t/, $_);

		if ($dbtype1 eq 'refGene') {
			@record == 16 or die "Error: invalid record in $dbfile (expecting 16 tab-delimited fields in refGene file): <$_>\n";
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes) = @record[1..15];		#human hg18, mouse
		} elsif ($dbtype1 eq 'knownGene') {
			@record >= 11 or die "Error: invalid record in $dbfile (>=11 fields expected in knownGene file): <$_>\n";	#mm8=11, hg18=hg19=12
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend) = @record[0..9];
			$name2 = $kgxref->{$name} || $name;
		} elsif ($dbtype1 eq 'ensGene') {
			@record == 16 or die "Error: invalid record in $dbfile (expecting 16 fields in ensGene file): <$_>\n";
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes) = @record[1..15];
		} else {
			@record >= 11 or die "Error: invalid record in $dbfile (>=11 fields expected in $dbtype1 gene definition file): <$_>\n";
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes) = @record[1..15];
			$name2 or $name2=$name;		#changed on 20130821 (deleted "defined $name2" so that ccdsGene can be handled, by relacing gene name by name2 field)
			#die "FATAL ERROR: the --dbtype $dbtype is not supported in readUCSCGeneAnnotation() subroutine.\n";		#commented 2011feb18
		}
	
		#handle situations where the same transcript is mapped to several chromosomes or regions (for example, NM_019105 is mapped to chr6, chr6_cox_hap1, chr6_qbl_hap2; NM_002538 is mapped to chr5 positive and negative strand and also in chr5_h2_hap1)
		if ($chr =~ m/hap\d+$/) {
			next;			#this is a temporary solution on 2011feb19, to ignore alternative haplotype chromosomes
		}
	
		#$chr =~ s/^chr// or die "Error: invalid record found in $dbfile (chrom field not found): <$_>\n";						#UCSC always prefix "chr" to the chromosome identifier, so this is a good check to make sure that the file is the correct file
		$chr =~ s/^chr//;			#some genomes like zebrafish does not start with chr in their chromosome names.
		
		$dbstrand eq '+' or $dbstrand eq '-' or die "Error: invalid dbstrand information found in $dbfile (dbstrand has to be + or -): <$_>\n";		#dbstrand is important to know and cannot be optional
		my @exonstart = split (/,/, $exonstart); 			#remove trailing comma
		my @exonend = split (/,/, $exonend);				#remove trailing comma
		$exoncount == @exonstart or die "Error: invalid record found in $dbfile (exoncount discordance): <$exoncount vs ${\(scalar @exonstart)}>\n";
		@exonstart == @exonend or die "Error: invalid record found in $dbfile (exonstart and exonend count discordance): <${\(scalar @exonstart)} vs ${\(scalar @exonend)}>\n";
		$txstart++; $cdsstart++; map {$_++} @exonstart;			#convert 0-based coordinate to 1-based coordinate

		#LOGIC here:
		#first calcluate mRNA length, and if the transcript maps to multiple locations with discordant mRNA length, only consider the leftmost chromosome and leftmost coordinate (because the FASTA file is sorted in this manner)

		my $cdslength = 0;
		my $mrnalength = 0;
		for my $i (0 .. @exonstart-1) {
			$mrnalength += $exonend[$i]-$exonstart[$i]+1;
		}
		for my $i (0 .. @exonstart-1) {					#this calculation is valid regardless of strand
			#$mrnalength += $exonend[$i]-$exonstart[$i]+1;
			if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
				if ($cdsend <= $exonend[$i]) {
					$cdslength = $cdsend-$cdsstart+1;
					last;
				} else {
					$cdslength += $exonend[$i]-$cdsstart+1;
					next;
				}
			}
			if ($cdslength and $cdsend < $exonstart[$i]) {
				die "FATAL ERROR: impossible scenario for $name in $dbfile (cdsend is less than exon start)";
			} elsif ($cdslength and $cdsend <= $exonend[$i]) {
				$cdslength += $cdsend-$exonstart[$i]+1;
				last;
			} elsif ($cdslength and $cdsend > $exonend[$i]) {
				$cdslength += $exonend[$i]-$exonstart[$i]+1;
			}
		}
		
		if ($cdsstart != $cdsend+1) {		#coding gene
			if (defined $mrnalen{$name} and $mrnalen{$name} != $mrnalength) {
				$verbose and printerr "WARNING: $name occurs more than once in $dbfile with different mRNA length. The first occurences with identical mRNA length will be uesd in analysis.\n";
				next;
			}
			
			if (defined $cdslen{$name} and $cdslen{$name} != $cdslength) {
				$verbose and printerr "WARNING: $name occurs more than once in $dbfile with different CDS length. The first occurences with identical CDS length will be uesd in analysis.\n";
				next;
			}
			
			$iscoding{$name2}++;		#name2 is a coding gene, and if there is a noncoding transcript, ignore such transcripts in future analysis
		} else {		#noncoding gene
			1;
		}
		
		$cdslen{$name} = $cdslength;
		$mrnalen{$name} = $mrnalength;
				
		my ($bin1, $bin2) = (int(($txstart - $neargene)/$genomebinsize), int(($txend + $neargene)/$genomebinsize));
		for my $nextbin ($bin1 .. $bin2) {
			push @{$genedb{$chr, $nextbin}}, [$name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, [@exonstart], [@exonend], $name2];
		}
		$geneidmap{$name} = $name2;
		$genecount++;
		$name2count{$name2}++;
		$cdsstart == $cdsend+1 and $ncgenecount++;			#non-coding gene has the same start and end site
	} 
	close (GENEDB);
	
	my %badgene;
	for my $key (keys %genedb) {
		my @newgenedb;
		for my $geneinfo (@{$genedb{$key}}) {
			if (not $cdslen{$geneinfo->[0]} and $iscoding{$geneinfo->[8]}) {
				$badgene{$geneinfo->[0]}++;
				$verbose and printerr "WARNING: $geneinfo->[0] will be ignored in analysis, because it is a non-coding transcript but the associated gene has another coding transcript\n";
			} else {
				push @newgenedb, $geneinfo;
			}
		}
		@{$genedb{$key}} = @newgenedb;
	}
	
	for my $key (keys %genedb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$genedb{$key}} = sort {$a->[2] <=> $b->[2]} @{$genedb{$key}};
	}
	printerr "Done with $genecount transcripts (including $ncgenecount without coding sequence annotation) for ", scalar (keys %name2count), " unique genes\n";
	$verbose and %badgene and printerr "NOTICE: ", scalar (keys %badgene), " noncoding transcripts will be ignored, because their associated genes have annotated coding transcript\n";
	return (\%genedb, \%geneidmap, \%cdslen, \%mrnalen);
}

sub downloadDB {
	my ($cwd, $msg, $sc);
	
	$cwd = Cwd::cwd();		#save current directory information to go back later
	
	-w $dbloc or die "Error: the directory $dbloc is not writable by the current user\n";
	chdir ($dbloc) or die "Error: the directory $dbloc cannot be accessed\n";
	
	my (@urlin, @filein, @fileout, %fail);		#the fail hash contains index of files that fail to be downloaded
	my $count_success;
	my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
	if ($dbtype1 eq 'refGene') {
		if ($webfrom eq 'annovar') {
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_refGene.txt.gz";
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_refLink.txt.gz";
		} else {
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/refGene.txt.gz";
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/refLink.txt.gz";
		}
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_refGeneMrna.fa.gz";
	} elsif ($dbtype1 eq 'knownGene') {
		if ($webfrom eq 'annovar') {
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_knownGene.txt.gz";
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_kgXref.txt.gz";
		} else {
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/knownGene.txt.gz";
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/kgXref.txt.gz";
		}
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_knownGeneMrna.fa.gz";
	} elsif ($dbtype1  eq 'ensGene') {
		if ($webfrom eq 'annovar') {
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_ensGene.txt.gz";
		} else {
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/ensGene.txt.gz";
		}
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_ensGeneMrna.fa.gz";
	} elsif ($dbtype1 eq 'seq') {
		push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/chromFa.zip";		#example: hg18, hg19
		push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/chromFa.tar.gz";	#example: panTro2
		push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/$buildver.fa.gz";	#example: bosTau4
	} elsif ($dbtype1 eq '1000g' or $dbtype1 eq '1000g2010' or $dbtype1 eq '1000g2010jul') {
		$buildver eq 'hg18' or die "Error: currently the --dbtype of '$dbtype1' only support --buildver of 'hg18'\n";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_$dbtype1.zip";
	} elsif ($dbtype1 =~ m/^1000g(\d{4})(\w{3})$/) {	#1000g2010nov, 1000g2011may, 1000g2012feb, 1000g2012apr, etc
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_$dbtype1.zip";
	} elsif ($dbtype1 eq 'null') {
		1;
	} else {
		$webfrom ||= 'ucsc';
		if ($webfrom eq 'annovar') {
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_$dbtype1.txt.gz";
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_$dbtype1.txt.idx.gz";
		} elsif ($webfrom eq 'ucsc') {
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/$dbtype1.txt.gz";
		} else {
			push @urlin, "$webfrom/$dbtype1.txt.gz";
		}
	}
	
	
	@filein = @urlin;
	map {s/.+\///} @filein;
	@fileout = @filein;
	map {s/\.gz$//; s/\.zip$//} @fileout;
	
	if ($wget) {
		$msg = qx/wget --help 2>&1/ || '';		#collect the output of the system command
	} else {
		$msg = '';					#when --nowget is specified, do not use wget to retrieve files from Internet
	}
	if ($msg =~ m/Usage/) {
		#checkProgramUpdate ("wget");
		for my $i (0 .. @urlin-1) {
			printerr "NOTICE: Downloading annotation database $urlin[$i] ... ";
			if ($verbose) {
				$sc = "wget -t 1 -T 30 -O $filein[$i] $urlin[$i]";
			} else {
				$sc = "wget -t 1 -T 30 -q -O $filein[$i] $urlin[$i]";
			}
			if (system ($sc)) {	#time-out is 30 seconds, with 1 retry attempt
				printerr "Failed\n";
				$verbose and print "WARNING: unable to execute system command: <$sc>\n";
				unlink ($filein[$i]);		#delete the temporary files generated by wget
				$fail{$i}++;
			} else {
				printerr "OK\n";
				$count_success++;
			}
		}
	} else {
		eval {
			require Net::FTP;
			require LWP::UserAgent;
		};
		if ($@) {
			printerr "WARNING: cannot retrieve remote files automatically (by 'wget' command or by standard Net::FTP/LWP::UserAgent Perl module).\n";
			printerr "Please manually download the following file, uncompress the files to $dbloc directory, then add a ${buildver}_ prefix to the file names.\n";
			printerr join ("\n", @urlin), "\n";
			exit (100);
		}
		
		#checkProgramUpdate ("lwp");
		my ($http, $ftp);
		for my $i (0 .. @urlin-1) {
			printerr "NOTICE: Downloading annotation database $urlin[$i] ... ";
			if ($urlin[$i] =~ m/^http/) {
				$http = LWP::UserAgent->new (timeout=>10, show_progress=>$verbose);
				$http->env_proxy;
				
				my $response = $http->get ($urlin[$i], ':content_file'=>$filein[$i]);
				if ($response->is_success) {
					printerr "Done\n";
					$count_success++;
				} else {
					printerr "Failed\n";
					$verbose and printerr "WARNING: cannot retrieve remote files ($urlin[$i]) via LWP::UserAgent Perl module: ", $response->status_line, "\n";
					$fail{$i}++;
				}
			} elsif ($urlin[$i] =~ m#^ftp://([^\\\/]+)#) {		#for hgdownload.cse.ucsc.edu, ftp-trace.ncbi.nih.gov, ftp.ensembl.org, etc
				my $urlroot = $1;
				if ($ftp = Net::FTP->new($urlroot, Timeout=>10, Debug=>$verbose)) {
					$ftp->login("anonymous", 'anonymous@');
					$ftp->binary();
					my $url = $urlin[$i];
					$url =~ s#ftp://[\w\.\-]+/##;		#remove the URL root
					if (not $ftp->get($url)) {
						printerr "Failed\n";
						$verbose and printerr "WARNING: cannot retrieve remote file ($url) in FTP server $urlroot\n";
						$fail{$i}++;
					} else {
						printerr "Done\n";
						$count_success++;
					}
				} else {
					printerr "Failed\n";
					$verbose and printerr "WARNING: cannot retrieve remote file ($urlin[$i]) via Net::FTP Perl module\n";
					$fail{$i}++;
				}
				
			} else {
				die "Error: The URL $urlin[$i] uses an unsupported protocol. Download cannot continue\n";
			}
		}
	}
	
	$count_success and printerr "NOTICE: Uncompressing downloaded files\n";
	for my $i (0 .. @filein-1) {
		$fail{$i} and next;
		if ($filein[$i] =~ m/\.zip$/) {
			$msg = qx/unzip --help 2>&1/ || '';		#collect the output of the system command
			if ($msg =~ m/Usage/i) {
				if ($verbose) {
					system ("unzip -o $filein[$i]");
				} else {
					system ("unzip -o -q $filein[$i]");
				}
			} else {
				printerr "ERROR: unzip is not installed in your system.\nPlease manually uncompress the files (@filein) at the $dbloc directory", $dbtype1 eq 'seq'?".\n":", and rename them by adding ${buildver}_ prefix to the file names.\n";
				exit (101);
			}
			unlink ($filein[$i]);				#delete the ZIP file
		} elsif ($filein[$i] =~ m/\.tar\.gz$/) {		#panTro2 FASTA sequence is stored as tar.gz rather than zip
			$msg = qx/tar --help 2>&1/ || '';		#collect the output of the system command
			if ($msg =~ m/Usage/i or $msg =~ m/Options/i) {	#BSD-derived version of tar on Mac OS does not list "Usage" information
				system ("tar -x -z -f $filein[$i]");
			} else {
				printerr "ERROR: tar/gunzip is not installed in your system.\nPlease manually uncompress the files (@filein) at the $dbloc directory", $dbtype1 eq 'seq'?".\n":", and rename them by adding ${buildver}_ prefix to the file names.\n";
				exit (102);
			}
		} elsif ($filein[$i] =~ m/\.gz$/) {
			$msg = qx/gunzip --help 2>&1/ || '';		#collect the output of the system command
			if ($msg =~ m/Usage/i) {
				system ("gunzip -f $filein[$i]");
			} else {
				printerr "ERROR: gunzip is not installed in your system.\nPlease manually uncompress the files (@filein) at the $dbloc directory", $dbtype1 eq 'seq'?".\n":", and rename them by adding ${buildver}_ prefix to the file names.\n";
				exit (103);
			}
		}
	}

	for my $i (0 .. @fileout-1) {
		$fail{$i} and next;				#skip the file that failed to be downloaded
		my $fileout = $fileout[$i];
		$dbtype1 eq 'seq' and next;			#the zip file contains dozens of FASTA files so cannot rename them automatically
		if (not $fileout =~ m/^${buildver}_/) {		#if the buildver is not the prefix of the files
			rename ($fileout, "${buildver}_$fileout") or die "Error: cannot rename $fileout to ${buildver}_$fileout\n";
			$fileout = "${buildver}_$fileout";
		}
		if (not $fileout =~ m/\.txt$/ and not $fileout =~ m/\.fa$/ and not $fileout =~ m/\.idx$/) {
			rename ($fileout, "$fileout.txt");
		}
	}
	
	$count_success and printerr "NOTICE: Finished downloading annotation files for $buildver build version, with files saved at the '$dbloc' directory\n";
	$cwd and chdir ($cwd);
	if (%fail) {
		my @failindex = keys %fail;
		if ($dbtype1 eq 'seq' and @failindex == 1) {	#not really a fail, because for seq, ANNOVAR attempts on tar.gz and zip file
			1;
		} else {
			printerr "WARNING: Some files cannot be downloaded, including ", join (', ', @urlin[@failindex]), "\n";
		}
		
		for my $index (@failindex) {
			if ($urlin[$index] =~ m#^http://www\.openbioinformatics\.org.+Mrna.fa.gz$#) {
				printerr "---------------------------ADDITIONAL PROCEDURE---------------------------\n";
				printerr "--------------------------------------------------------------------------\n";
				printerr "NOTICE: the FASTA file $urlin[$index] is not available to download but can be generated by the ANNOVAR software. ";
				printerr "PLEASE RUN THE FOLLOWING TWO COMMANDS CONSECUTIVELY TO GENERATE THE FASTA FILES:\n\n";
				printerr "\tannotate_variation.pl --buildver $buildver --downdb seq $dbloc/${buildver}_seq\n";
				printerr "\tretrieve_seq_from_fasta.pl $dbloc/${buildver}_$dbtype1.txt -seqdir $dbloc/${buildver}_seq -format $dbtype1 -outfile $dbloc/${buildver}_${dbtype1}Mrna.fa\n";
				printerr "--------------------------------------------------------------------------\n";
				printerr "--------------------------------------------------------------------------\n";
			}
		}
	}	
}

sub printerr {
	print STDERR @_;
	print LOG @_;
}

