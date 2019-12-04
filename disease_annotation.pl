#!/usr/bin/env perl

######################################################################
## LIBRARIES

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use Cwd;
use File::Basename;
use Bio::OntologyIO;
use Graph::Directed;
use File::Copy;
use JSON;

## LIBRARIES CONDITIONALLY LOADED
## Shown here for documentation purposes. 
# Parallel::ForkManager;

##
######################################################################

######################################################################
## INITIALISATION

our (@out_gene_scores_file);

# Declare help options
our (
		$verbose,
		$help,
		$man,
		$buildver,
		$bedfile
	);

# Declare query options
our (
		$query_diseases,
		$if_file,
		$if_exact_match,
		$prediction,
		$is_phenotype,
		$if_wordcloud,
		$if_use_precalc
	);

# Declare test options
our (
		$if_hi,
		$if_rvis
	);

# Declare
our (
		$out,
		$database_directory,
		$precalc_db_directory,
		$if_logistic_regression,
		$user_nproc
	);

# Declare weight options
our (
		$HPRD_WEIGHT,
		$BIOSYSTEM_WEIGHT,
		$GENE_FAMILY_WEIGHT,
		%GENE_WEIGHT,
		$HTRI_WEIGHT,
		$GENE_DISEASE_WEIGHT,
		$HPO_GENE_WEIGHT,
		$INTERCEPT
	);

# Declare files options
our  (
		$ctd_disease_file,
		$hprd_file,
		$biosystem_file,
		$disease_count_file,
		$gene_disease_score_file,
		$hpo_annotation_file,
		$gene_annotation_file,
		$omim_disease_id_file,
		$biosystem_to_info_file,
		$gene_family_file,
		$htri_file,
		$omim_description_file,
		$hi_gene_score_file,
		$rvis_gene_score_file,
		$addon_gene_disease_score_file,
		$addon_gene_gene_score_file,
		$addon_gene_disease_weight,
		$addon_gene_gene_weight,
		$genelist,
		@genes,
		%gene_hash,
		%gene_id,
		%hi_score,
		%rvis_score,
		$path,
		$work_path,
		$phenotype_to_gene_file
		);

# %gene_hash ( $gene => "Not Annotated" or "Chr:Pos1 - Pos2")
# %omim_disease (lc $first_disease => join(";",$each_disease)  )

# Declare options
		GetOptions(
				'verbose|v'             =>\$verbose,
				'help|h'                =>\$help,
				'man|m'                 =>\$man,
				'file|f'                =>\$if_file,
				'directory|d=s'         =>\$database_directory,
				'use_precalc'           =>\$if_use_precalc,
				'work_directory|w=s'    =>\$work_path,
				'out=s'                 =>\$out,
				'prediction|p'          =>\$prediction,
				'buildver=s'            =>\$buildver,
				'bedfile=s'             =>\$bedfile,
				'gene=s'                =>\$genelist,
				'phenotype|ph'          =>\$is_phenotype,
				'exact'                 =>\$if_exact_match,
				'logistic'              =>\$if_logistic_regression,
				'haploinsufficiency|hi' =>\$if_hi,
				'intolerance|it'        =>\$if_rvis,
				'addon=s'               =>\$addon_gene_disease_score_file,
				'addon_gg=s'            =>\$addon_gene_gene_score_file,
				'addon_weight=s'        =>\$addon_gene_disease_weight,
				'addon_gg_weight=s'     =>\$addon_gene_gene_weight,
				'hprd_weight=s'         =>\$HPRD_WEIGHT,
				'biosystem_weight=s'    =>\$BIOSYSTEM_WEIGHT,
				'gene_family_weight=s'  =>\$GENE_FAMILY_WEIGHT,
				'htri_weight=s'         =>\$HTRI_WEIGHT,
				'hpo_gene_weight=s'		=>\$HPO_GENE_WEIGHT,
				'gwas_weight=s'         =>\$GENE_WEIGHT{"GWAS"},
				'gene_reviews_weight=s' =>\$GENE_WEIGHT{"GENE_REVIEWS"},
				'clinvar_weight=s'      =>\$GENE_WEIGHT{"CLINVAR"},
				'omim_weight=s'         =>\$GENE_WEIGHT{"OMIM"},
				'orphanet_weight=s'     =>\$GENE_WEIGHT{"ORPHANET"},
				'wordcloud'             =>\$if_wordcloud,
				'nproc=i'               =>\$user_nproc
				) or pod2usage ();

# Launch help or man
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage   ("      ERROR: Please enter disease names!!");
@ARGV == 1 or pod2usage ("       ERROR: too many input arguments");

# Initialise logistics regression
if($if_logistic_regression) {
	print STDOUT "NOTICE: The logistic regression model was used!!!\n";
	defined $GENE_DISEASE_WEIGHT or $GENE_DISEASE_WEIGHT = 9.5331966;			#20191029: users can change these parameters through command line
	defined $HPRD_WEIGHT         or $HPRD_WEIGHT         = 0.8335866;
	defined $BIOSYSTEM_WEIGHT    or $BIOSYSTEM_WEIGHT    = 0.1755904 ;
	defined $GENE_FAMILY_WEIGHT  or $GENE_FAMILY_WEIGHT  = 0.3561601 ;
	defined $HTRI_WEIGHT         or $HTRI_WEIGHT         = 4.1003533 ;

	$GENE_DISEASE_WEIGHT = 24.62841;
	$HPRD_WEIGHT         = 0.61063;
	$BIOSYSTEM_WEIGHT    = 0.17723;
	$GENE_FAMILY_WEIGHT  = 0.52574;
	$HTRI_WEIGHT         = 6.41266;
}

# Initialise variables
$GENE_DISEASE_WEIGHT         = 1.0  unless (defined $GENE_DISEASE_WEIGHT);
$HPRD_WEIGHT                 = 0.1  unless (defined $HPRD_WEIGHT);
$BIOSYSTEM_WEIGHT            = 0.05 unless (defined $BIOSYSTEM_WEIGHT);
$GENE_FAMILY_WEIGHT          = 0.05 unless (defined $GENE_FAMILY_WEIGHT);
$HTRI_WEIGHT                 = 0.05 unless (defined $HTRI_WEIGHT );
$GENE_WEIGHT{"GWAS"}         = 1.0  unless (defined $GENE_WEIGHT{"GWAS"});
$GENE_WEIGHT{"GENE_REVIEWS"} = 1.0  unless (defined $GENE_WEIGHT{"GENE_REVIEWS"});
$GENE_WEIGHT{"CLINVAR"}      = 1.0  unless (defined $GENE_WEIGHT{"CLINVAR"} );
$GENE_WEIGHT{"OMIM"}         = 1.0  unless (defined $GENE_WEIGHT{"OMIM"} );
$GENE_WEIGHT{"ORPHANET"}     = 1.0  unless (defined $GENE_WEIGHT{"ORPHANET"} );
$GENE_WEIGHT{"HPO_PHENOTYPE_GENE"}     = 1.0  unless (defined $GENE_WEIGHT{"HPO_PHENOTYPE_GENE"} );
$addon_gene_disease_weight   = 1.0  unless (defined $addon_gene_disease_weight);
$addon_gene_gene_weight      = 1.0  unless (defined $addon_gene_gene_weight);
$user_nproc                  = 1    unless (defined $user_nproc);

$HPO_GENE_WEIGHT             = 0.1	unless (defined $HPO_GENE_WEIGHT);		#In recent version of Phenolyzer, we add HPO annotation (from JAX) to supplement the known phenotype-gene relationships

# Test input values
$user_nproc >= 1 or pod2usage ("       ERROR: number of requested processes is zero or negative");

##
######################################################################

######################################################################
## Main program

setup_variables();

annotate_bedfile() if (defined $bedfile);

# output_gene_prioritization() works by: (1) process_terms() to generate output file for input phenotype terms. It uses pre-computed HPO score optionally, and can support multiple threads, and it uses process_individual_terms() subroutine to score each phenotype term. The key subroutine include score_genes() which counts for each term,  (2) it generates merge_gene_scores, annotated_gene_scores and seed_gene_list file, to combine scores from multiple phenotype terms (3) find max/min score for all genes, and then normalize so that max is 1.
# MERGE file has: MERGE "{\"Name\":\"$gene\",\"Id\":\"ID:$gene_id{$gene}\",\"Score\":$normalized_score,\"Diseases\":[$content_json]}";
# ANNOTATED file has: $gene."\tID:$gene_id{$gene} ".$gene_hash{$gene}."\t$normalized_score\n".$content."\n"
# PREDICTED file has: $gene."\t"."ID:$gene_id{$gene} - $status\t".$score."\tNormalized score: $normalized_score\n".$content."\n";

output_gene_prioritization();

##
######################################################################

######################################################################
## Subroutines

sub print_array
{
	foreach (@_){
		print "$_\n";
	}
}

sub process_individual_term
{
	my $individual_term = $_[0];
	my $raw_individual_term = $individual_term;

	if ($individual_term =~ /^\W*$/){ return; } # next if there is no word 
		$individual_term=TextStandardize($individual_term);
# For a normal term, disease name extension is needed
	if ($individual_term !~/^all([\s_\-]diseases?)?$/i) {
# Expand the disease term and save them into files
# @diseases are lower case disease array, %disease_hash keeps the original disease name
		my $diseases_reference = disease_extension($individual_term);
		my %disease_hash = %$diseases_reference;
		my @diseases;
		for my $disease_key (keys %disease_hash) {
			my  @records =split("\n", $disease_hash{$disease_key});
			for my $record (@records) {
				return if(not $record);
				my ($disease_line, $source) = split("\t", $record);
				my @disease_terms = split(";", $disease_line);
				for  (@disease_terms){
					my $each = $_;
					$each=TextStandardize($each);
					my $change2= $each =~ s/\btype //ig;
					if($change2){
						push @diseases,$each;
						$disease_hash{$disease_key} = $each.';'.$disease_hash{$disease_key};
					}
				}
				push @diseases,@disease_terms;
			}
		}
		my ($hash, %disease_score_hash, @hpo_ids);
		if ($is_phenotype) {		#this is controlled by the --phenotype argument in command line; we will perform phenotype extension on the input phenotype, so @diseases array becomes larger
			($hash, @hpo_ids) = phenotype_extension($individual_term);
			%disease_score_hash = %$hash;
			for (@diseases) {
				my $disease_key = lc $_;
				delete $disease_score_hash{$disease_key} if($disease_score_hash{$disease_key});
				$disease_key = TextStandardize($disease_key);
				delete $disease_score_hash{$disease_key} if($disease_score_hash{$disease_key});
			}
			for (keys %disease_score_hash) {
				my $disease_score = join ("\t", ($disease_score_hash{$_}[1], $disease_score_hash{$_}[0]) );
				push (@diseases, lc $disease_score);
			}
		}

#The non-word characters are changed into '_'
		$individual_term=~s/\W+/_/g;
		if(@hpo_ids) {
			open(OUT_PHENOTYPE, ">$out"."_$individual_term"."_hpo") or die;
			print OUT_PHENOTYPE $_."\n" for @hpo_ids;
			close(OUT_PHENOTYPE);
		}

		open (OUT_DISEASE,">$out"."_$individual_term"."_diseases") or die;
		for (keys %disease_hash) {
			my @lines=split("\n", $disease_hash{$_});
			for my $line (@lines) {
				return if(not $line);
				my @words = split("\t", $line);
				my @diseases = split(";", $words[0]);
				@diseases = Unique(@diseases);
				my $disease_line = join(";",@diseases);
				my $out_line = join("\t",($disease_line,$words[1]));
				print OUT_DISEASE $out_line."\n";
			}
		}

		if ($is_phenotype) {
			print OUT_DISEASE join ("\t", ($disease_score_hash{$_}[1], $disease_score_hash{$_}[0]) )."\n"
				for (keys %disease_score_hash);
		}
		close (OUT_DISEASE);

		generate_wordcloud($individual_term, \%disease_hash, \%disease_score_hash) if($if_wordcloud);
		if(@diseases==0){
			print STDERR "ERROR: The input term -----$individual_term------ has no corresponding names in the disease database, please check your spelling!!!\n";
			return;
		}

		my $i=0;
# Output the gene_score files
# $item = {  $gene => [$score, $information_string] }
# $information_string = "ID (SOURCE)	DISEASE_NAME RAW_SCORE
		@diseases = map {
			my @words = split("\t");
			$words[0] = TextStandardize($words[0]);
			$words[0] = lc $words[0];
			join("\t", @words);
		} @diseases;

		@diseases = Unique(@diseases);
		my ($item,$count) = score_genes (\@diseases, \@hpo_ids, $raw_individual_term);		#item is hash (key=gene value=[score, information])
		my %output=();
		@{$output{$_}} = @{$item->{$_}} for keys %$item;

		if($count==0){
			print STDERR "ERROR: The input term -----$individual_term------ has no results!!!\n";
			return;
		}

		open( OUT_GENE_SCORE,">$out"."_$individual_term"."_gene_scores") or die "can't write to "."$out"."_$individual_term"."_gene_scores";
		print OUT_GENE_SCORE "Tuple number in the gene_disease database for the term $individual_term: $count\n";

		for my $gene (sort{ $output{$b}[0] <=> $output{$a}[0] }keys %output){
#The probability of the gene when the disease is given
			my $p=$output{$gene}[0]/$count;
			print OUT_GENE_SCORE $gene."\t"."Normalized score: $p\tRaw Score: $output{$gene}[0]\n".$output{$gene}[1]."\n";
		}
		print STDOUT "$raw_individual_term: -------------------------------------------------- \n";
		close (OUT_GENE_SCORE);
	}
# For the term 'all disease(s)'(case insensitive)
	else {
		$individual_term=~s/\W+/_/g;
		open(OUT_GENE_SCORE,">$out"."_$individual_term"."_gene_scores") or die "can't write to "."$out"."_$individual_term"."_gene_scores";
		my ($item,$count)=score_all_genes();
		my %output=();
		@{$output{$_}} = @{$item->{$_}} for keys %$item;

		print OUT_GENE_SCORE "Tuple number in the gene_disease database for the term $individual_term: $count\n";
		for my $gene(sort{ $output{$b}[0] <=> $output{$a}[0] }keys %output){
#The probability of the gene when the disease is given
			my $p=$output{$gene}[0]/$count;
			print OUT_GENE_SCORE $gene."\t"."Normalized score: $p\tRaw Score: $output{$gene}[0]\n".$output{$gene}[1]."\n";
		}
		print STDOUT "------------------------------------------------------------------------ \n";
		close (OUT_GENE_SCORE);
	}
}

sub get_terms_precalc_db
{
	open(DB_MASTER, "$precalc_db_directory/hpos.txt") or die "Can't open $precalc_db_directory/hpos.txt !";
	chomp(my @lines = <DB_MASTER>);
	close(DB_MASTER);
	return @lines
}

sub process_terms
{
	my @disease_input = sort @_;
## Process each individual term
# Find pre-processed terms from the input
	my @hpo_list = grep /hp:[0-9]*/, @disease_input;
	# Search HPO DB only if there are HPO inputs
	if ($if_use_precalc && @hpo_list > 0) {
		my @found_list;
		# Read master table
		my @db_terms = get_terms_precalc_db();
		my $i = 0;
		for my $hpo(@hpo_list) {
			(my $num = $hpo) =~ s/hp:*//;
			# This assumes the HPO inputs are ordered and so are the DB terms
			for(; $i < @db_terms; $i++) {
				if( $num == $db_terms[$i] ){
					push (@found_list, $hpo);
					$i++;
					last;
				} elsif($num < $db_terms[$i]) {
					last;
				}
			}
		}
		my $exe_ext = "";
		unless ($^O =~ /Unix|Linux|unix|linux/) {
			$exe_ext = ".exe";
		}

		# Copy preprocessed terms to directory
		for my $hpo(@found_list) {
			print("$hpo: Using HPO expansion found in pre-expanded database ...\n");
			(my $mod_hpo = $hpo) =~ s/:/_/;
			my $basename = "precalc_db_$mod_hpo";
			my $zipfile = "$precalc_db_directory/db/$basename";
			my $tempdir = "$precalc_db_directory/temp";
			my $dst = "$out"."_$mod_hpo";
			foreach (qw(_hpo _diseases _gene_scores)) {
				# Extract from DB
				if (system("unzip$exe_ext -q $zipfile $basename$_ -d $tempdir")) {
					if (system("7z$exe_ext -e -y $zipfile $basename$_ -o $tempdir")) {
						die "Failed unzipping file from DB. Install unzip or 7z command line program. File: $zipfile to $basename$_";
					}
				}
				if (system("mv $tempdir/$basename$_ $dst$_")) {
					die "Failed moving file from temp DB: $tempdir/$basename$_ -> $dst$_";
				}
			}
			# Push to list of gene_scores files
			push @out_gene_scores_file, "$dst"."_gene_scores";
		}
		# Continue with remaining terms
		# @disease_input = @remaining_list;

		my %temp;
		@temp{ @disease_input } = ();
		delete @temp{ @found_list };
		@disease_input = keys %temp;
		if (@disease_input == 0) {
			return;
		}
	}
# Determine number of parallel subprocesses
	my $n_procs = $user_nproc <= @disease_input ? $user_nproc : @disease_input;
	print STDOUT "NOTICE: Processing $n_procs phenotypes at a time!\n";
	if ($n_procs > 1) {
# Run in parallel after loading ForkManager. Test if the module is available.
		eval {require Parallel::ForkManager};
		if ($@) {
			pod2usage ("Error in argument: you need to install Parallel::ForkManager module before parallelizing with the -nproc argument");
		}
		my $parallel_mngr = Parallel::ForkManager->new($n_procs);
		for my $individual_term(@disease_input)
		{
			$parallel_mngr->start and next;
			process_individual_term($individual_term);
			$parallel_mngr->finish;
		}
		$parallel_mngr->wait_all_children();
	} else {
# Process one input at a time.
		for my $individual_term(@disease_input)
		{
			process_individual_term($individual_term);
		}
	}
## Fill @out_gene_scores_file array
	for my $individual_term(@disease_input) 
	{
		#(my $individual_term = $individual_term) =~ s/:/_/;	#not sure what exactly this is doing (probably change hp:1234 to hp_1234)
		$individual_term=TextStandardize($individual_term);		#added 20191029, to be identical as process_individual_term
		$individual_term=~s/\W+/_/g;	#The non-word characters are changed into '_'	(added 20191029)
		my $dst = "$out"."_$individual_term";
		if (-f "$dst"."_gene_scores") {
			push @out_gene_scores_file, "$dst"."_gene_scores";
		}
	}
}

# The main sub to output prioritized genelist
sub output_gene_prioritization 
{
	my @disease_input = split (qr/[^ _,\w\.\-'\(\)\[\]\{\}:]+/,lc $query_diseases);		#' the $query_diseases include all phenotype terms separated by tab or other characters, convert it to lower case and split into separate phrases
	s/^\s+|\s+$// for @disease_input;		#remove heading and trailing spaces



# Process individual phenotype terms, possibly with parallel processing and with pre-computed scores where available. Generate $out_${individual_term}_gene_scores file
	process_terms(@disease_input);

# Finish processing individual terms
# Merge the gene_score files
# %output  ( gene =>[score,  content])
	my ($item, $count)=merge_result();
	my %output=();
	my ($max_score, $min_score) = (0, 1);

	@{$output{$_}} = @{$item->{$_}} for keys %$item;

	open (MERGE,    ">$out.merge_gene_scores") or die;
	# open (MERGE_LEGACY,    ">$out.merge_gene_scores_legacy") or die;
	open (ANNOTATED,">$out.annotated_gene_scores") if (%gene_hash and not $prediction);
	open (my $seed_fh, ">$out.seed_gene_list");

	my $annotated_seed_fh;
#### header for seed gene list ####
	printHeader($seed_fh, 0);
	if (%gene_hash and not $prediction) {
		open ($annotated_seed_fh, ">$out.annotated_gene_list");
		printHeader($annotated_seed_fh,0);
	}

	# print MERGE_LEGACY     "{\n\"header\": \"Tuple number in the gene_disease database for all the terms: $count\" \n";
	print MERGE     "[";
	print ANNOTATED "Tuple number in the gene_disease databse for all the terms: $count \n" if (%gene_hash and not $prediction);

#Find the max and min score
	for my $gene (keys %output) {
		$max_score = $output{$gene}[0] if ($max_score < $output{$gene}[0]);
		$min_score = $output{$gene}[0] if ($min_score > $output{$gene}[0]);
	}

	my $rank = 0;
	my $annotated_rank = 0;
	my (%hi_out, %rvis_out);

	for my $gene(sort{ $output{$b}[0] <=> $output{$a}[0] }keys %output) {
		$rank++;
		my ($score, $content) = ($output{$gene}[0], $output{$gene}[1]);
#$score = ($score - $min_score)/$diff;
		my $normalized_score = $score/$max_score;
                #my $normalized_score = $score;

#normalize scores of each detail
		chomp($content);
		my @content_lines = split("\n", $content);
		my @content_lines_output;
		my @content_lines_next;
		my %content_json_dict;
		for(@content_lines) {
			my @words = split("\t");
			my $detail_score = $words[3];
			$detail_score = $detail_score/$max_score;
			my $new_score = $detail_score * $GENE_DISEASE_WEIGHT;
			my $line = join("\t", (@words[0,1,2],$detail_score));
			my $new_line = join("\t", (@words[0,1,2],$new_score));
			push (@content_lines_output,$line);
			push (@content_lines_next,  $new_line);
			# Build disease hash, organized by origin (OMIM, ORPHANET...: words[0])
			my ($db, $temp) = split(':', $words[0]);
			my ($db_ids, $db_note) = $words[0] =~ /(.*) \((.*)\)/; # first group captures all
			$db_note =~ s/(\r)//g; # remove "\r" if present
			$words[2] =~ s/(\D+)//g; # remove "hp "
			if (!exists($content_json_dict{$db})) {
				$content_json_dict{"$db"}=[];
			}
			my %temp = (Id => $db_ids, IdNote => $db_note, Condition => $words[1], Hp => $words[2], Score => $detail_score);
			push(@{ $content_json_dict{$db} }, \%temp);
		}
		my $content_json = encode_json \%content_json_dict;
		$content = join("\n", @content_lines_output)."\n";
		my $new_content = join("\n", @content_lines_next)."\n";

#Normalized the input for prediction
		$item->{$gene}[0] = $normalized_score * $GENE_DISEASE_WEIGHT;
		$item->{$gene}[1] = $new_content;
		$item->{$gene}[2] = $normalized_score;

#print out results
		if($rank != 1) {
			print MERGE ","
		}
		print MERGE "{\"Name\":\"$gene\",\"Id\":\"ID:$gene_id{$gene}\",\"Score\":$normalized_score,\"Diseases\":[$content_json]}";
		# print MERGE_LEGACY $gene."\t"."ID:$gene_id{$gene} -\t$normalized_score\n".$content."\n";
		print ANNOTATED $gene."\tID:$gene_id{$gene} ".$gene_hash{$gene}."\t$normalized_score\n".$content."\n"
			if ($gene_hash{$gene} and not $prediction);

#Normalized score for the genelist
		$normalized_score = sprintf('%.4g', $normalized_score);
		print $seed_fh $rank."\t".$gene."\t".$gene_id{$gene}."\t".$normalized_score;
		if($if_hi) {
			$hi_score{$gene} = 0 if not defined $hi_score{$gene};
			print $seed_fh "\t$hi_score{$gene}";
		}
		if($if_rvis) {
			$rvis_score{$gene} = 0 if not defined $rvis_score{$gene};
			print $seed_fh "\t$rvis_score{$gene}";
		}
		print $seed_fh "\n";

#Normalized score for the annotation list
		if ($gene_hash{$gene} and not $prediction) {
			print $annotated_seed_fh join("\t",(++$annotated_rank, $gene, $gene_id{$gene}, $normalized_score));
			if($if_hi) {
				$hi_score{$gene} = 0 if not defined $hi_score{$gene};
				print $annotated_seed_fh "\t$hi_score{$gene}";
			}
			if($if_rvis){
				$rvis_score{$gene} = 0 if not defined $rvis_score{$gene};
				print $annotated_seed_fh "\t$rvis_score{$gene}";
			}
			print $annotated_seed_fh "\n";
		}
	}
	print MERGE "]";
	close (MERGE);
	# close (MERGE_LEGACY);
	close ($seed_fh);
	close (ANNOTATED) if (%gene_hash and not $prediction);
	close (ANNOTATED_GENE_LIST) if (%gene_hash and not $prediction);

# Integrate scores from relation databases
	if($prediction) {
		($max_score, $min_score) = (0, 1);
		my $predicted_item = predict_genes($item);		#this is a hash with key=gene and value=[predicted_score, evidence]. The predicted_score was calculated in predict_genes() with weights already applied there
		my %predicted_output = ();
		my $annotated_fh;
		@{$predicted_output{$_}} = @{$predicted_item->{$_}} for keys %$predicted_item;

		if(%gene_hash) {
			open (ANNOTATED,">$out.annotated_gene_scores") ;
			open ($annotated_fh, ">$out.annotated_gene_list");
			printHeader($annotated_fh,1);
		}

		open (PREDICTED, ">$out.predicted_gene_scores");
		open (my $final_fh,">$out.final_gene_list");
		printHeader($final_fh, 1);
		print  PREDICTED   "Tuple number in the gene_disease database for all the terms: $count \n";
		print  ANNOTATED   "Tuple number in the gene_disease database for all the terms: $count \n"
			if %gene_hash;

#Find the max and min score
		for my $gene (keys %predicted_output) {
			$max_score = $predicted_output{$gene}[0] if ($max_score < $predicted_output{$gene}[0]);
			$min_score = $predicted_output{$gene}[0] if ($min_score > $predicted_output{$gene}[0]);
		}

# my $diff = $max_score - $min_score;
		my $rank = 0;
		my $annotated_rank = 0;
		for my $gene (sort{ $predicted_output{$b}[0] <=> $predicted_output{$a}[0] } keys %predicted_output) {
			$rank++;
			$predicted_output{$gene}[1] =~ /^.*?\(  (.+?)  \).*?\t/x;
			my $source = $1;
			my $status;
			if(
					($source eq "HPRD") or
					($source eq "BIOSYSTEM") or
					($source eq "GENE_FAMILY") or
					($source eq "HTRI") or
					($source eq "ADDON_GENE_GENE")
			  ) {
				$status = "Predicted";
			} else {
				$status = "SeedGene";
			}

			my ($score, $content) = ($predicted_output{$gene}[0], $predicted_output{$gene}[1]);

			my $normalized_score = $predicted_output{$gene}[0]/$max_score;
                        #my $normalized_score = $predicted_output{$gene}[0];

			print PREDICTED $gene."\t"."ID:$gene_id{$gene} - $status\t".$score."\tNormalized score: $normalized_score\n".$content."\n";
			print ANNOTATED $gene."\tID:$gene_id{$gene} $gene_hash{$gene} $status"."\t".$score."\tNormalized score: $normalized_score\n".$content."\n"
				if ($gene_hash{$gene});

#Normalize score for the genelist
			$normalized_score = sprintf('%.4g', $normalized_score);
			print $final_fh $rank."\t".$gene."\t".$gene_id{$gene}."\t".$normalized_score."\t".$status;
			if($if_hi) {
				$hi_score{$gene} = 0 if not defined $hi_score{$gene};
				print $final_fh "\t$hi_score{$gene}";
			}

			if($if_rvis) {
				$rvis_score{$gene} = 0 if not defined $rvis_score{$gene};
				print $final_fh "\t$rvis_score{$gene}";
			}
			print $final_fh "\n";

			if ($gene_hash{$gene}) {
				print $annotated_fh  join("\t",(++$annotated_rank,$gene,$gene_id{$gene},$normalized_score,$status)) ;
				if($if_hi) {
					$hi_score{$gene} = 0 if not defined $hi_score{$gene};
					print $annotated_fh "\t$hi_score{$gene}";
				}
				if($if_rvis) {
					$rvis_score{$gene} = 0 if not defined $rvis_score{$gene};
					print $annotated_fh "\t$rvis_score{$gene}";
				}
				print $annotated_fh "\n";
			}
		}

		close (PREDICTED);
		close ($final_fh);
		close (ANNOTATED) if %gene_hash;
		close (ANNOTATED_GENE_LIST) if (%gene_hash);
	}
}

# Input some disease terms and return all its extended diseases
sub disease_extension{
	@_==1 or die "The input should be only one string!";

	my $input_term=$_[0];
	print STDOUT "$input_term: The journey to find all related disease names of your query starts!\n";
	$input_term =~ s/[\W_]+/ /g;

	-f "${path}/$disease_count_file" or die "Could not open ${path}/$disease_count_file";
	-f "${path}/${ctd_disease_file}" or die "Could not open ${path}/${ctd_disease_file}";

	open(DISEASE,"${path}/$disease_count_file") or die "can't open ${path}/$disease_count_file";
	open(CTD_DISEASE,"${path}/${ctd_disease_file}") or die "can't open ${path}/${ctd_disease_file}";
	open(OMIM_DISEASE_ID, "${path}/$omim_disease_id_file") or die "can't open ${path}/$omim_disease_id_file";

	my %disease_extend=();
	my @disease_occur=<DISEASE>;
	my @disease_ctd=<CTD_DISEASE>;

	print STDOUT "$input_term: The item was queried in the databases!! \n";
	print STDOUT "$input_term: The exact match (case non-sensitive) was used for disease/phenotype name match!! \n"
		if($if_exact_match);

	for (<OMIM_DISEASE_ID>) {
		chomp();
		my ($id, $disease_line) = split("\t");
		next if ($id eq "OMIM_ID");

# When compare, get rid of '-' if it is not after a number
		my $disease_line_key = $disease_line;
		$disease_line_key =~ s/\bs\b//g;
		$disease_line_key =~ s/\W+/ /g;

		my $query_term = $input_term;
		$query_term =~ s/\bs\b//g;
		$query_term =~ s/\W+/ /g;

		if ($disease_line_key =~/\b$query_term\b/i or $query_term eq $id) {
# If exact match
			next if($if_exact_match and $disease_line_key !~ /(^|;)$query_term($|;)/i);
			my @diseases = split(";",$disease_line);
			my $disease_key =lc $diseases[0];
			$disease_extend{$disease_key} = $disease_line;
		}
	}

	print STDOUT "$input_term: The word matching in OMIM disease synonyms file has been done!! \n";

	$input_term = lc $input_term;
#Query disease in the compiled list from gene_disease relations
	for my $term (@disease_occur) {
		chomp($term);

		my @words=split('\t',$term);
		my $disease=$words[0];
		my $id = $words[1];
		my ($id_source,$id_num) = split (":", $id);

		next if(not $id_source);

		my $disease_key = lc $disease;
		$disease_key =~ s/\bs\b//g;
		$disease_key =~ s/\W+/ /g;

		my $query_term = $input_term;
		$query_term =~ s/\bs\b//g;
		$query_term =~ s/\W+/ /g;

# If the term matches
		if(
				$disease_key =~ /\b$query_term\b/i or
				(
				 $id_num and
				 $query_term eq $id_num and
				 $id_source eq "OMIM"
				)
		  ) {
# If exact match
			next if($if_exact_match and $disease_key !~ /(^|;)$query_term($|;)/i);
			if ($disease_extend{$disease_key}) {
				$disease_extend{$disease_key}.=";".$disease;
			} else {
#Save the disease name into hash, key is the disease name in lower case form withougt "-"
				$disease_extend{$disease_key}=$disease;
			}
		}
	}

	$disease_extend{$_} .= "\tGENE_DISEASE\n" for keys %disease_extend;
	my @tree_number=();

	print STDOUT "$input_term: The word matching search in the compiled disease databases for gene_disease relations has been done!\n";

	for my $term(@disease_ctd) {
		chomp($term);

		my @words = split('\t',$term);
		my $disease = $words[0];

		my $disease_key = lc $disease;
		$disease_key =~ s/\bs\b//g;
		$disease_key =~ s/\W+/ /g;

		my $synonym_key = $words[1];
		$synonym_key =~ s/\bs\b//g;
		$synonym_key =~ s/\W+/ /g;

		my $query_term = $input_term;
		$query_term =~ s/\bs\b//g;
		$query_term =~ s/\W+/ /g;

# First push all the matched disease names or synonyms in
		if($disease_key=~/\b$query_term\b/i or $synonym_key=~/\b$query_term\b/i) {
#If exact match
			next if($if_exact_match and $disease_key.$words[1] !~ /(^|\|)$input_term($|\|)/i);

# Retrieve all the synonyms
			my @synonyms = split('\|',$words[1]);

#Record the tree_number of each term and trace all their children later
			push @tree_number,split('\|',$words[2]);
			$disease_extend{$disease_key} .= join(';', ($disease,@synonyms)) and
				$disease_extend{$disease_key} .= "\tCTD_DISEASE\n";
		}
	}

	print STDOUT "$input_term: The word matching search in the CTD (Medic) databases has been done! \n";

#Second find all children of the terms found in the first round
	for my $term(@disease_ctd) {
		chomp($term);

		my @words=split('\t',$term);
		my @synonyms=split('\|',$words[1]);
		my $disease = $words[0];
		my $disease_key = lc $disease;

		for my $each_tree_num(@tree_number) {
			if($words[2] =~ qr/(^|\|)   $each_tree_num [^\|]+    /x) {
				$disease_extend{$disease_key} .=join (';', ($disease, @synonyms) )
					and $disease_extend{$disease_key} .= "\tCTD_DISEASE\n"
					if(not defined $disease_extend{$disease_key} or $disease_extend{$disease_key}!~/\bCTD_DISEASE\b/);
			}
		}
	}

	print STDOUT "$input_term: The descendants search in the CTD (Medic) databases has been done! \n";

	if (-f "$work_path/ontology_search.pl") {
		print STDERR "ERROR: The doio.obo file couldn't be found!!! The disease_ontology search wouldn't be conducted properly!! \n"
			if ( not -f "$path/doid.obo");
		$input_term =~ s/[^ \w-]s?//g;

		my $system_command = "perl $work_path/ontology_search.pl -o $path/doid.obo '$input_term' -f name,id";

		$system_command .= " -exact" if($if_exact_match);
		$system_command .= " 2>/dev/null";

		my $line = `$system_command`;
		my @ontology_diseases = split('DOID:\d*\n',$line);
		for (@ontology_diseases) {
			next if(not $_);
			my ($disease, @synonyms) = split('\n');
			my $disease_key = lc $disease;

			$disease_extend{$disease_key} .= join(';', ($disease, @synonyms));
			$disease_extend{$disease_key} .= "\tDISEASE_ONTOLOGY\n";
		}
		print STDOUT "$input_term: The descendants search in disease_ontology (DO) database has been done! \n";
	} else {
		print STDERR "ERROR: $input_term: The $work_path/ontology_search.pl file couldn't be found, so the disease_ontology (DO) database won't be used! \n";
	}
	return \%disease_extend;
}

sub phenotype_extension{
	@_==1 or die "ERROR: Only one phenotype term is accepted!!! ";

	my %hpo_score_system = (
			"very rare" => 0.01,
			"rare"      => 0.05,
			"occasional"=> 0.075,
			"frequent"  => 0.33,
			"typical"   => 0.5,
			"variable"  => 0.5,
			"common"    => 0.75,
			"hallmark"  => 0.90,
			"obligate"  => 1.00
			);

	my $input_term = $_[0];
	my $raw_input_term = $input_term;
	print STDOUT "$raw_input_term: The phenotype search and annotation process starts!! \n";
	$input_term =~s/[\W_]+/ /g;
	my %disease_hash;
	my @hpo_ids;

# %disease_hash( "disease_name_key" => [score, original_disease_name] )
	if( -f "$work_path/ontology_search.pl" ) {
		print STDERR "ERROR: The hpo.obo file couldn't be found!!! The phenotype_ontology search wouldn't be conducted properly!! \n"
			if ( not -f "$path/hpo.obo");
		print STDERR "ERROR: The $hpo_annotation_file couldn't be found!!! The phenotype_annotation couldn't be conducted properly!! \n"
			if ( not -f "$path/$hpo_annotation_file" );

		open (HPO_ANNOTATION, "$path/$hpo_annotation_file") or die "ERROR: Can't open $hpo_annotation_file!!! \n";
		open (OMIM_DESCRIPTION, "$path/$omim_description_file") or die "ERROR: Can't open $omim_description_file!!! \n";

		my $line = `perl $work_path/ontology_search.pl -o $path/hpo.obo -format id -p '$input_term' `;
		print STDOUT "$raw_input_term: executing ontology_search.pl on '$input_term' to expand phenotype terms\n";
		@hpo_ids = split("\n", $line);
		print STDOUT "$raw_input_term: Found ", scalar (@hpo_ids), " additional phenotype terms\n";
		
		my @hpo_annotation = <HPO_ANNOTATION>;
		shift @hpo_annotation;
		@hpo_ids = sort @hpo_ids;
		my ($i,$j) = (0,0);

		while($i<@hpo_ids and $j<@hpo_annotation) {
			$hpo_annotation[$j] =~s/[\r\n]+//g;

#[0]HPO_ID	[1]SOURCE	[2]HPO_DISEASE_NAME	[3]FREQUENCY
			my @words=split("\t",$hpo_annotation[$j]);

			if($hpo_ids[$i] eq "HP:".$words[0]) {
				my @diseases = split(";",$words[2]);
				for my $individual_disease(@diseases) {
					next if(not $individual_disease);

					$individual_disease = TextStandardize($individual_disease);
					my $individual_disease_key = lc $individual_disease;

					$disease_hash{$individual_disease_key}[0] = 0;
					$disease_hash{$individual_disease_key}[1] =lc $individual_disease if not $disease_hash{$individual_disease_key}[1];

					my $score;

					if (defined $words[3] and $hpo_score_system{$words[3]}) {
						$score = $hpo_score_system{$words[3]};
					} else {
						if($words[3] and $words[3] =~ /^(\d*\.?\d+)\%$/) {
							$score = $1 * 0.01;
						} else {
							$score = $hpo_score_system{"frequent"};
						}
					}
					$disease_hash{$individual_disease_key}[0] = $score if($disease_hash{$individual_disease_key}[0] < $score);
				}
				$j++;
			}
			if($hpo_ids[$i] lt "HP:".$words[0]) { $i++; next; }
			if($hpo_ids[$i] gt "HP:".$words[0]) { $j++; next; }
		}
	} else {
		print STDERR "NOTICE: The $work_path/ontology_search.pl file couldn't be found, so the HPO database won't be used! \n";
	}
	my %omim_description;
	for my $line (<OMIM_DESCRIPTION>) {
		next if($line =~ /^OMIM_ID/);
		chomp($line);
		my ($id, $disease, $description) = split("\t", $line);
		$disease = TextStandardize($disease);
		if($description =~ /\b$input_term\b/i) {
			$omim_description{$disease} = 1 if(not $omim_description{$disease});
			$omim_description{$disease}++ if($omim_description{$disease});
		}
	}
	my $total=0;

	for my $disease (keys %omim_description) {
		$total += $omim_description{$disease};
	}

	for my $disease (keys %omim_description) {
		my $score = ($omim_description{$disease}+0.0)/$total;
		my @diseases = split(";",$disease);
		for my $individual_disease(@diseases) {
			next if(not $individual_disease);
			my $individual_disease_key = lc $individual_disease;

			if(not $disease_hash{$individual_disease_key}) {
				$disease_hash{$individual_disease_key}[0] = $score;
				$disease_hash{$individual_disease_key}[1] = $individual_disease;
			} else {
				$disease_hash{$individual_disease_key}[0] = $score if($disease_hash{$individual_disease_key}[0] < $score);
			}
		}
	}
	return (\%disease_hash,@hpo_ids);
}

# Input the disease list and return all the genes and item count

sub get_phenotype_to_gene_hash
{
	open(PHENOTYPE_TO_GENE, "${path}/$phenotype_to_gene_file") or die "can't open file: ${path}/$phenotype_to_gene_file";
	my @lines = <PHENOTYPE_TO_GENE>;
	shift @lines;
	my $line;
	my %hpo_id_to_gene_hash;
	my $hpo_id;
	my $gene_names;
	foreach $line (@lines)
	{
		chomp $line;
		my @split_line = split('\t', $line);
		$hpo_id = $split_line[0];
		$gene_names = $split_line[2];
		$hpo_id_to_gene_hash{$hpo_id} = $gene_names;
	}
	close (PHENOTYPE_TO_GENE);
	return \%hpo_id_to_gene_hash;
}

sub get_hpo_id_to_name_hash
{

	open(PHENOTYPE_TO_GENE, "${path}/$phenotype_to_gene_file") or die "can't open file: ${path}/$phenotype_to_gene_file";
	my @lines = <PHENOTYPE_TO_GENE>;
	shift @lines;
	my $line;
	my %hpo_id_to_name_hash;
	my $hpo_id;
	my $hpo_name;
	foreach $line (@lines)
	{
		chomp $line;
		my @split_line = split('\t', $line);
		$hpo_id = $split_line[0];
		$hpo_name = $split_line[1];
		if (not exists($hpo_id_to_name_hash{$hpo_id} ) ){
			$hpo_id_to_name_hash{$hpo_id} = "$hpo_name";
		}
	}
	close (PHENOTYPE_TO_GENE);
	return \%hpo_id_to_name_hash;
}


sub score_genes{

	@_ == 3 or die "input should contain references to two arrays and a reference input term!";

# item is a hash, keys are gene names, values are an array reference and a the total score for the gene
	my ($ref1, $ref2, $individual_term) = @_;

	my @diseases = @{$ref1};
	my @hpo_ids  = @{$ref2};


	my $hpo_id_to_gene_hash_ref = get_phenotype_to_gene_hash();

	my %hpo_id_to_gene_hash = %$hpo_id_to_gene_hash_ref;

	my $hpo_id_to_name_hash_ref = get_hpo_id_to_name_hash();

	my %hpo_id_to_name_hash = %$hpo_id_to_name_hash_ref;

	my %phenotype_related_genes_hash; # key = gene name, value = hpo_id 

	foreach (@hpo_ids)
	{
		my $hpo_id = $_;
		if (not exists($hpo_id_to_gene_hash{$hpo_id})){
			next;
		}
		my @gene_array = split(",", $hpo_id_to_gene_hash{$hpo_id});	
		foreach (@gene_array) {
			$phenotype_related_genes_hash{$_} = $hpo_id;
		}
	}

	my %item=();

#my @diseases=@{$_[0]};

# $count will record how many tuples are retrived from the GENE_DISEASE_SCORE database

	my $count = 0.0;

	open(SCORE,"${path}/$gene_disease_score_file") or die "could not open ${path}/$gene_disease_score_file";

	my @disease_gene_score = <SCORE>;

	shift @disease_gene_score;

	my @addon_disease_gene_score;

	if($addon_gene_disease_score_file) 
	{
		my @addon_files = split(',', $addon_gene_disease_score_file);

		for my $each_file (@addon_files)
		{
			open(ADDON,"${path}/$each_file") or die "could not open ${path}/$each_file";
			push(@addon_disease_gene_score, <ADDON>);
			@addon_disease_gene_score = map {s/[\n\r]+//g;$_; } @addon_disease_gene_score;
			close(ADDON);
			print STDOUT "$individual_term: The ${path}/$each_file is used as addons!!!\n";
		}
		push (@disease_gene_score,@addon_disease_gene_score);
	}

	@disease_gene_score = sort {
		my @words1=split("\t",$a);
		my @words2=split("\t",$b);
		$words1[1] cmp  $words2[1];
	}

	map {
		my @words=split("\t");
		$words[1]=lc $words[1];
		$words[1]=TextStandardize($words[1]);
		join("\t",@words);
	} @disease_gene_score;

	@diseases = sort {
		my @words1=split("\t",$a);
		my @words2=split("\t",$b);
		$words1[0] cmp $words2[0];
	} @diseases;

	my ($i,$j)=(0,0);


	while($i<@diseases and $j<@disease_gene_score) {

		chomp($disease_gene_score[$j]);

		last if(not $disease_gene_score[$j]);

# @words:[0]GENE	[1]DISEASE	[2]DISEASE_ID	[3]SCORE	[4]SOURCE

		my @words = split ("\t",$disease_gene_score[$j]);

		my ($query_disease, $inference_score) = split("\t",$diseases[$i]);

		if ($query_disease eq $words[1]) {
			$count++;
			$inference_score = 1.0 if (not $inference_score);

			my @genes = split(",",$words[0]);
			my $gene = $genes[0];

			if( not $gene) {
				$verbose and print STDOUT "WARNING: gene name not found in disease_gene_score array <$disease_gene_score[$j]>\n";		#some entries miss gene name such as "FEBRILE SEIZURES FAMILIAL 1     umls:C1852577   0.3     DISGENET" in DB_DISGENET_GENE_DISEASE_SCORE
				$j++;
				next;
			}

			$GENE_WEIGHT{$words[4]} = $addon_gene_disease_weight if (not defined $GENE_WEIGHT{$words[4]});

			my $score = $words[3]*$inference_score*$GENE_WEIGHT{$words[4]};

			if ($score != 0)
			{
				if ($item{$gene}[0]) {
					$item{$gene}[0] += $score ;
				} else {
					$item{$gene}[0] = $score ;
				}
				$item{$gene}[1].=$words[2]." ($words[4])"."\t".$words[1]."\t".$score."\n";
			}
			$j++;
		}
		if($query_disease lt $words[1]) {$i++;}
		if($query_disease gt $words[1]) {$j++;}
	}

	for (keys %phenotype_related_genes_hash)
	{
		my $gene = $_;
		$count++;
		my $hpo_id = $phenotype_related_genes_hash{$gene};
		my $hpo_name = $hpo_id_to_name_hash{$hpo_id};

		if ( exists($item{$gene}) )
		{
			$item{$gene}[0] += $HPO_GENE_WEIGHT; 
			$item{$gene}[1] .= "$hpo_id (HPO_PHENOTYPE_GENE)\t$hpo_name\t$HPO_GENE_WEIGHT\n"; # !!!!!!!!!!! Li changed on 11/27/2019
		}else{
			$item{$gene}[0] = $HPO_GENE_WEIGHT; 
			$item{$gene}[1] = "$hpo_id (HPO_PHENOTYPE_GENE)\t$hpo_name\t$HPO_GENE_WEIGHT\n";  # !!!!!!!!!!! Li changed on 11/27/2019
		}
	}

	for (keys %item){
		delete  $item{$_}  if (not $gene_id{$_});
	}

	my @value_array;

	print STDOUT "$individual_term: The gene score process has been done!!\n";

	return (\%item,$count);
}

# Merge gene_scores for each term, return %item reference and $count for the sum of the tuple count
sub merge_result {
# %item (  $gene[0] => score, $gene[1] => "SOURCE_ID	 DISEASE_NAME	DISEASE_TERM"     )
	print STDOUT "NOTICE: Start to merge gene scores! \n";

	#my $dirname=dirname($out);
	#my $basename=basename($out);
	#my @filelist = split("\n",`ls $dirname`);
	print STDOUT "NOTICE: Reading a list of gene_scores files: <@out_gene_scores_file>\n";

#item will save the results for output
	my %item=();
	my $count=0.0;
	my $individual_count;

	for my $filename(@out_gene_scores_file) {
		if($filename=~/${out}_(\w+)_gene_scores$/){
			my $term=$1;
			$term=~s/_/ /g;

			open(GENE_SCORE, $filename) or die "ERROR: Can't open $filename: $!!!!\n";
			my $gene="";
			my $i=0;

			while (my $line = <GENE_SCORE>) {
				chomp($line);

				if($i==0) {
					$line=~/: (\d+)/;
					$individual_count = $1;
					$count+=$1;
					$i++;
					next;
				}

#Proress the first line of each term_score
				if($line=~/^\s*$/) {
					$gene="";
					next;
				}

#[0]SOURCE	[1]EVIDENCE	[2]RAW_SCORE
				my @words=split("\t",$line);

				if(not $gene){
					$gene=$words[0];
					$words[1]=~/Normalized score: (.*?)$/;
					my $score=$1;
$score eq 'alzheimer' and die "Die gene=$gene score=<$score> current line died at <$line> ";
					$item{$gene}[0]+=$score if defined $item{$gene};
					$item{$gene}[0] =$score if not defined $item{$gene};
				} else {
					my ($source, $disease, $individual_score)= @words;
					$item{$gene}[1].= $source."\t".$disease."\t".$term."\t".$individual_score/$individual_count."\n";
				}
			}
		} else {
			die "Unknown error in file name processing: $filename is not recognized as a valid gene_scores file\n";
		}
	}
	return (\%item,$count);
}

#GENE	DISEASE	DISEASE_ID	SCORE	SOURCE
sub score_all_genes{
	print STDOUT "NOTICE: The process to score all genes in the database starts!!\n";

#item is a hash, keys are gene names, values are the total score for the gene an array reference and
	my %item=();
# $count will record how many tuples are retrived from the GENE_DISEASE_SCORE database
	my $count=0.0;

	open(SCORE,"${path}/$gene_disease_score_file")  or die "could not open ${path}/$gene_disease_score_file";
	my @addon_disease_gene_score;
	my @disease_gene_score=<SCORE>;

	if($addon_gene_disease_score_file) {
		my @addon_files = split(',', $addon_gene_disease_score_file);
		for my $each_file (@addon_files) {
			open(ADDON,"${path}/$each_file") or die "could not open ${path}/$each_file";

			push(@addon_disease_gene_score, <ADDON>);
			@addon_disease_gene_score = map {
				s/[\n\r]+//g;
				$_;
			} @addon_disease_gene_score;

			close(ADDON);

			print STDOUT "NOTICE: The ${path}/$each_file is used as addons!!!\n";
		}
		push (@disease_gene_score,@addon_disease_gene_score);
	}
	my $i=0;

	for (@disease_gene_score) {
		chomp;

		my ($genes, $disease, $disease_id, $score, $source)=split("\t");
		my @genes = split(",", $genes);
		my $gene = $genes[0];

		next if (not $gene);

		if($i==0){
			$i++;
			next;
		}

		$count++;
		$GENE_WEIGHT{$source}= $addon_gene_disease_weight if (not defined $GENE_WEIGHT{$source});

		$score *= $GENE_WEIGHT{$source};

		if($score!=0) {
			if($item{$gene}[0]) {
				$item{$gene}[0] += $score ;
			} else  {
				$item{$gene}[0] = $score ;
			}
			$item{$gene}[1].=$disease_id." ($source)"."\t".$disease."\t".$score."\n";
		}

	}

	print STDOUT "NOTICE: The process to score all genes has been done!!\n";

	for (keys %item){
		delete $item{$_} if (not $gene_id{$_});
	}
	return (\%item,$count);
}

# This function is to convert CNV into genes
## Code borrowed from bed2gene.pl
sub annotate_bedfile {
	$buildver = 'hg19' if(not $buildver);

	print STDOUT "NOTICE: the --buildver argument is set as 'hg19' by default\n"

		unless defined $buildver;

	$buildver eq 'hg18' or $buildver eq 'hg19' or pod2usage ("Error in argument: the --buildver argument must be 'hg18' or 'hg19'");

	my $sc;

	$sc = "perl $path/../../bin/convert.pl -format bed $bedfile > $out.avinput";
	system ($sc) and die "Error: cannot execute system command $sc\n";

	$sc = "perl $path/../../bin/annotate.pl -geneanno -buildver $buildver -outfile $out $out.avinput $path/../humandb";
	system ($sc) and die "Error: cannot execute system command $sc\n";

	my ($countregion, $countexonic) = qw/0 0/;
	my ($totallen);

	open (VF, "$out.variant_function") or die "Error: cannot read $out.variant_function file\n";
	open (REGION_GENE, ">$out.genelist_from_region") or die "Error: can't write to genelist_from_region\n";

	while (<VF>) {
		my @field = split (/\t/, $_);
		$countregion++;
		$field[0] =~ m/exonic/ or next;
		$countexonic++;
		my @gene = split (/,/, $field[1]);

		for my $gene (@gene) {
			$gene_hash{$gene} = "$field[2]:$field[3]-$field[4]";
		}

		$totallen += ($field[4]-$field[3]+1);
	}

	print REGION_GENE $_."\n" for keys %gene_hash;
	close (VF);
	close (REGION_GENE);

	print STDOUT "NOTICE: Among $countregion BED regions ($totallen base pairs), $countexonic disrupt exons, and ", scalar keys %gene_hash, " genes are affected\n";
### Code borrowed from bed2gene.pl
}

sub setup_variables{
#Input argument setup
	$query_diseases = $ARGV[0];

	# Initialisation of variables
	#$out_directory = cwd();			#in the older version of Phenolyzer, the output directory is fixed to be current directory, which is not ideal
	my $dirname       = dirname(__FILE__);		#script file's location
	#chdir $dirname;

	$path = "$dirname/lib/compiled_database";
	$database_directory and $path=$database_directory;
	$path =~s/\/$//;
	$work_path ||= $dirname;

	if (defined $precalc_db_directory) {
		# pass
	} else {
		$precalc_db_directory = "$dirname/lib/precalculated_expansions";
	}
	
	if (defined $out) {
		$out =~ s/[\\\/]*$//;		#remove trailing "/" or "\" if there are any
		(-d dirname($out)) or system("mkdir ".dirname($out)) and die "ERROR: $out is not a value output, make sure that the directory exists and is writtable!!";
	} else {
		$out = "out";
	}

	$disease_count_file      = "DB_COMPILED_DISEASE_COUNT";
	$gene_disease_score_file = "DB_COMPILED_GENE_DISEASE_SCORE";
	$ctd_disease_file        = "DB_COMPILED_CTD_DISEASES";
	$hprd_file               = "DB_COMPILED_BINARY_PROTEIN_INTERACTION_SCORE";
	$biosystem_file          = "DB_COMPILED_BIOSYSTEM_SCORE";
	$hpo_annotation_file     = "DB_HPO_ANNOTATION";
	$gene_annotation_file    = "DB_HUMAN_GENE_ID";
	$omim_disease_id_file    = "DB_COMPILED_OMIM_ID_DISEASE";
	$gene_family_file        = "DB_HGNC_GENE_FAMILY";
	$htri_file               = "DB_HTRI_TRANSCRIPTION_INTERACTION";
	$omim_description_file   = "DB_OMIM_DESCRIPTION";
	$hi_gene_score_file      = "DB_HI_SCORE";
	$rvis_gene_score_file    = "DB_RVIS_SCORE";
	$phenotype_to_gene_file  = "DB_COMPILED_HPO_PHENOTYPE_GENE";

	my %gene_transform;

# The phenotype/disease input will come from file
	if (defined $if_file) {
		print STDOUT  "NOTICE: The file name option was used!! \n";

		open(INPUT_DISEASE,"$query_diseases") or die "can't open $query_diseases: $!";
		my @input=<INPUT_DISEASE>;

		$query_diseases = join("\t",@input);

		if(defined $genelist){
			open(INPUT_GENELIST,"$genelist") or die "can't open $genelist";
			@input = <INPUT_GENELIST>;
			$genelist = join("\t",@input);
		}
	}

#build the haploinsufficiency hash, by default from lib/compiled_database/DB_HI_SCORE
	if($if_hi){
		print STDOUT "NOTICE: The haploinsufficency score is used!\n";
		open (HI, "$path/$hi_gene_score_file") or die;
		for (<HI>){
			s/[\r\n+]//g;
			my ($gene, $score) = split("\t");
			$hi_score{$gene} = $score;
		}
		close(HI);
	}

#build the gene intolerance shash, by default from lib/compiled_database/DB_RVIS_SCORE
	if($if_rvis){
		print STDOUT "NOTICE: The gene intolerance score is used!\n";
		open(RVIS, "$path/$rvis_gene_score_file") or die;
		for (<RVIS>){
			s/[\r\n]+//g;
			my ($gene, $score) = split("\t");
			$rvis_score{$gene} = $score;
		}
	}

# Read information from gene annotation file (3 tab-delimited columns, see example below), by default lib/compiled_database/DB_HUMAN_GENE_ID
# %gene_transform maps Entrez ID or gene name synonym to upper case gene name
#1       A1BG    A1B|ABG|GAB|HYST2477
#2       A2M     A2MD|CPAMD5|FWP007|S863-7
#3       A2MP1   A2MP
#9       NAT1    AAC1|MNAT|NAT-1|NATI
#10      NAT2    AAC2|NAT-2|PNAT
	open (GENE_ID, "$path/$gene_annotation_file") or die;
	my $i=0;

	for my $line (<GENE_ID>) {
		if($i==0) {
			$i++;
			next;
		}

		chomp($line);

		my ($id, $gene, $synonyms) = split("\t", $line);
		$gene_transform{$gene} = uc $gene;
		$gene_transform{uc $gene} = uc $gene;
		$gene_transform{$id} = uc $gene;
	}

	seek GENE_ID, 0,0;
	$i=0;

	for my $line (<GENE_ID>) {
		if($i==0) {
			$i++;
			next;
		}

		chomp($line);

		my ($id, $gene, $synonyms) = split("\t", $line);

		if($synonyms eq "-") {
			next;
		} else {
			my @synonyms = split('\|', $synonyms);
			for my $each (@synonyms) {
				$gene_transform{uc $each} = uc $gene if(not defined $gene_transform{$each});
			}
		}
	}

# THE genelist will be saved into %gene_hash
	if(defined $genelist) {
		my @genes=split(qr/[^_\w\.\-]+/m,$genelist);
		for my $individual_term(@genes){
			if($individual_term =~ /^\W*$/){
				next;
			}

# Get rid of the whitespaces in the beginnning and end
			$individual_term =~ s/^\W*(.*?)\W*$/$1/;
			my $upper_gene = uc $individual_term;

			if($gene_transform{$upper_gene}){
				$gene_hash{$gene_transform{$upper_gene}} ="-" unless defined $gene_hash{$gene_transform{$upper_gene}};
			}
		}
	}

	$i=0;
	seek GENE_ID,0,0;
	for my $line (<GENE_ID>) {
		if($i==0) {
			$i++;
			next;
		}
		chomp($line);
		my ($id, $gene, $synonyms) = split("\t", $line);
		$gene_id{$gene} = $id;
		$gene_id{uc $gene} = $id;
	}
	close (GENE_ID);
}

sub predict_genes{
	print STDOUT "----------------------------------------------------------------------\n";
	print STDOUT "NOTICE: The prediction process starts!!!\n";

	@_ == 1 or die "You can only have one input argument!!!";

	my %item = ();
	my %output = ();
	%item = %{$_[0]};
	@{$output{$_}} = @{$item{$_}} for keys %item;
	my $gene_count = keys %item;
	my $i = 0;
	my %biosystem_id_type = ();
#%biosystem_id_type = (
#          $biosystem_id => $type
#              )
	my %biosystem = ();
# %biosystem = (
#     $system_num => {
#                 $gene => $score
#                     }
#               )
	my %interaction = ();
# %interaction = (
#         $gene1 => { $gene2 => [$score, $biosystem_id] }       #The score here is the maximal score of all the possible biosystems
#                  )
	my %gene_family = ();
# %gene_family = (
#     $gene_family_tag => {
#                 $gene => true of false
#                     }
#               )

	# open (HPRD, "$path/$hprd_file") or die "Can't open $path/$hprd_file !";
	# open (BIOSYSTEM, "$path/$biosystem_file") or die "Can't open $path/$biosystem_file !";
	# open (GENE_FAMILY, "$path/$gene_family_file") or die "Can't open $path/$gene_family_file!";
	# open (HTRI, "$path/$htri_file") or die "Can't open $path/$htri_file!";

	my @ggfiles;
	if($addon_gene_gene_score_file) {
		@ggfiles = split(",", $addon_gene_gene_score_file);
	}

# Predict genes based on Addon Gene-Gene relations
	for my $each_file (@ggfiles) {
		open(ADDON_GG,"$path/$each_file") or die "Can't open $path/$each_file!";
		for my $line (<ADDON_GG>) {
			if ($i==0) {
				$i++;
				next;
			}
			chomp ($line);

			my ($gene1, $gene2, $evidence, $score, $pubmed_id) = split ("\t", $line);
			my $individual_score;

			$gene1 = uc $gene1;
			$gene2 = uc $gene2;

			next if($gene1 =~ /^\W*$/ or $gene2 =~ /^\W*$/);

			$pubmed_id =~s/,/ /g;
			if($item{$gene1}[0] and ($gene1 ne $gene2)) {
				$individual_score = $score  * $addon_gene_gene_weight * $item{$gene1}[2];   #$item{$gene1}[2] saves the normalized score
					$output{$gene2}[0] = 0 if (not $output{$gene2}[0]);
				if( $individual_score !=0) {
					$output{$gene2}[0] +=  $individual_score;
					$output{$gene2}[1] .= "PUBMED:$pubmed_id (ADDON_GENE_GENE)\t".$evidence."\t"
						."With $gene1"."\t".$individual_score."\n";
				}
			}

			if($item{$gene2}[0] and ($gene1 ne $gene2)) {
				$individual_score = $score * $addon_gene_gene_weight * $item{$gene2}[2];
				$output{$gene1}[0] = 0 if (not $output{$gene1}[0]);
				if($individual_score !=0) {
					$output{$gene1}[0] += $individual_score;
					$output{$gene1}[1] .= "PUBMED:$pubmed_id (ADDON_GENE_GENE)\t".$evidence."\t"
						."With $gene2"."\t".$individual_score."\n";
				}
			}

		}
		close(ADDON_GG);
		print STDOUT "NOTICE: The Addon Database loaded !\n";
	}

	print STDOUT "NOTICE: Loading HPRD Database ... ";
	open (HPRD, "$path/$hprd_file") or die "Can't open $path/$hprd_file !";
#Predict genes based on Human Protein Interactions
	for my $line (<HPRD>) {
		if ($i==0) {
			$i++;
			next;
		}
		chomp ($line);

		my ($gene1, $gene2, $evidence, $score, $pubmed_id) = split ("\t", $line);

		$gene1 = uc $gene1;
		$gene2 = uc $gene2;

		my $individual_score;

		next if($gene1 eq '-' or $gene2 eq '-');

		$pubmed_id =~s/,/ /g;

		if($item{$gene1}[0] and ($gene1 ne $gene2) ) {
			$individual_score = $score * $HPRD_WEIGHT * $item{$gene1}[2];   #$item{$gene1}[2] saves the normalized score ($normalized_score * $GENE_DISEASE_WEIGHT;$new_content; $normalized_score;)
				$output{$gene2}[0] = 0 if (not $output{$gene2}[0]);
			if($individual_score != 0) {
				$output{$gene2}[0] +=  $individual_score;
				$output{$gene2}[1] .= "PUBMED:$pubmed_id (HPRD)\t".$evidence."\t"
					."With $gene1"."\t".$individual_score."\n";
			}
		}

		if($item{$gene2}[0] and ($gene1 ne $gene2)) {
			$individual_score = $score * $HPRD_WEIGHT * $item{$gene2}[2];
			$output{$gene1}[0] = 0 if (not $output{$gene1}[0]);
			if($individual_score != 0) {
				$output{$gene1}[0] += $individual_score;
				$output{$gene1}[1] .= "PUBMED:$pubmed_id (HPRD)\t".$evidence."\t"
					."With $gene2"."\t".$individual_score."\n";
			}
		}
	}
	close(HPRD);
	print STDOUT "Done!\n";

	print STDOUT "NOTICE: Loading HTRI Database ... ";
#Predict genes based on transcription interaction
	$i = 0;
	my $TF_PENALTY=4;
	open (HTRI, "$path/$htri_file") or die "Can't open $path/$htri_file!";
	#TF	TG	EVIDENCE	PUBMED	SCORE
	for my $line (<HTRI>) {
		if ($i==0) {
			$i++;
			next;
		}
		chomp ($line);

		my ($gene1, $gene2, $evidence, $pubmed_id, $score) = split ("\t", $line);

		$gene1 = uc $gene1;
		$gene2 = uc $gene2;

		my $individual_score;

		next if($gene1 eq '-' or $gene2 eq '-');

		$pubmed_id =~s/;/ /g;

		if($item{$gene1}[0] and ($gene1 ne $gene2) ) {
			$individual_score = $score * $HTRI_WEIGHT * $item{$gene1}[2];
			$output{$gene2}[0] = 0 if (not $output{$gene2}[0]);
			if( $individual_score !=0) {
				$output{$gene2}[0] +=  $individual_score;
				$output{$gene2}[1] .= "PUBMED:$pubmed_id (HTRI)\t".$evidence."\t"
					."Regulated by $gene1"."\t".$individual_score."\n";
			}
		}

		if($item{$gene2}[0] and ($gene1 ne $gene2) ) {
			$output{$gene1}[0] = 0 if (not $output{$gene1}[0]);
			$individual_score = $score * $HTRI_WEIGHT * $item{$gene2}[2];
			if($individual_score !=0) {
				$output{$gene1}[0] += $individual_score/$TF_PENALTY;
				$output{$gene1}[1] .= "PUBMED:$pubmed_id (HTRI)\t".$evidence."\t"
					."Regulates $gene2"."\t".$individual_score."\n";
			}
		}
	}
	close(HTRI);
	print STDOUT "Done!\n";

#Predict genes based on gene family
	open (GENE_FAMILY, "$path/$gene_family_file") or die "Can't open $path/$gene_family_file!";
	print STDOUT "NOTICE: Loading HGNC Gene Family Database ... ";
	$i = 0;
	my %in_gene_family = ();

	for my $line (<GENE_FAMILY>) {
		if ($i==0) {
			$i++;
			next;
		}
		chomp ($line);
		my ($gene, $tag, $description) = split ("\t", $line);
		$gene = uc $gene;
		$gene_family{$tag}{$gene} = $description;
	}

#Each gene family
	for my $tag (keys %gene_family) {
#Each gene
		for my $gene ( keys %{$gene_family{$tag}} ) {
			if ($item{$gene}[0]) {
#Each related gene in the gene_family
				for my $related_gene (keys %{$gene_family{$tag}}) {
					if ($related_gene ne $gene) {
						my $description = $gene_family{$tag}{$related_gene};
						$in_gene_family{$related_gene}{$gene}[0] .= $tag." ";
						$in_gene_family{$related_gene}{$gene}[1] .= $description." ";
					}
				}
			}
		}
	}

	for my $related_gene (keys %in_gene_family) {
		for my $gene (keys %{$in_gene_family{$related_gene}}) {
			my $individual_score = $GENE_FAMILY_WEIGHT * $item{$gene}[2];
			$output{$related_gene}[0] = 0 if (not$output{$related_gene}[0] );
			if ($individual_score !=0 ) {
				$output{$related_gene}[0] += $individual_score;

				my $tag_string = $in_gene_family{$related_gene}{$gene}[0];
				$tag_string = substr ( $tag_string, 0, length($tag_string)-1 );

				my $description_string = $in_gene_family{$related_gene}{$gene}[1];
				$description_string = substr ($description_string, 0, length($description_string)-1 );

				$output{$related_gene}[1] .= "NAME:".$tag_string." (GENE_FAMILY)\t"."In the family ($description_string)\t"."With $gene\t".$individual_score."\n";
			}
		}
	}
	close (GENE_FAMILY);
	print STDOUT "Done!\n";

	print STDOUT "NOTICE: Loading Biosystem Database ... ";
	open (BIOSYSTEM, "$path/$biosystem_file") or die "Can't open $path/$biosystem_file !";
#Predict genes based on Biosystem
	my %biosystem_id_name = ();
	$i = 0;

	for my $line (<BIOSYSTEM>) {
		if ($i==0) {$i++; next; }
		chomp ($line);
		my ($biosystem_id, $gene, $score, $biosystem_name) = split ("\t", $line);
		$gene = uc $gene;
		$biosystem{$biosystem_id}{$gene} = $score;
		$biosystem_id_name{$biosystem_id} = $biosystem_name if(not defined $biosystem_id_name{$biosystem_id}) ;
	}

#Each biosystem
	for my $biosystem_id (keys %biosystem) {
#Each gene
		for my $gene (keys %{$biosystem{$biosystem_id}}) {
			if ($item{$gene}[0]) {
#Each related gene in the same biosystem
				for my $related_gene (keys %{$biosystem{$biosystem_id}}) {
					if ($related_gene ne $gene) {
						my $score = $biosystem{$biosystem_id}{$related_gene};
						$interaction{$related_gene}{$gene}[0] = 0 if(not $interaction{$related_gene}{$gene}[0]);
						$interaction{$related_gene}{$gene}[0] = $score
							if ( $score > $interaction{$related_gene}{$gene}[0] );
						$interaction{$related_gene}{$gene}[1] .= $biosystem_id." ";
						$interaction{$related_gene}{$gene}[2] .= $biosystem_id_name{$biosystem_id}."; ";
					}
				}
			}
		}
	}

	for my $related_gene (keys %interaction) {
		for my $gene (keys %{$interaction{$related_gene}}) {
			my $individual_score = $BIOSYSTEM_WEIGHT * $interaction{$related_gene}{$gene}[0] * $item{$gene}[2];
			$output{$related_gene}[0] = 0 if (not $output{$related_gene}[0] );

			if ($individual_score !=0) {
				$output{$related_gene}[0] += $individual_score;

				my $id_string = $interaction{$related_gene}{$gene}[1];
				$id_string = substr ( $id_string, 0, length($id_string)-1 );

				my $type_string = $interaction{$related_gene}{$gene}[2];
				$type_string = substr ($type_string, 0, length($type_string)-2 );

				$output{$related_gene}[1] .= "BIOSYSTEM:".$id_string." (BIOSYSTEM)\t"."In the same ($type_string)\t"."With $gene\t".$individual_score."\n";
			}
		}
	}
	close (BIOSYSTEM);
	print STDOUT "Done!\n";

	for (keys %output){
		delete $output{$_} if (not $gene_id{$_});
	}

	return \%output;
}

sub generate_wordcloud{
	@_==3 or die "Error: generate_wordcloud only accept 3 variables!!";

	my $term = $_[0];
	my %disease_hash = %{$_[1]};
	my %disease_score_hash = %{$_[2]};

	open (WORD_CLOUD,">${out}_${term}_wordcloud") or die "Error: can't write into ${term}_wordcloud.txt!!";

	my %output =();

	for (keys %disease_hash) {
		chomp($disease_hash{$_});
		my @words=split("\t",$disease_hash{$_});
		my @diseases=split(/\W+/,$words[0]);
		@diseases = map {lc $_;} @diseases;
		length($_)>2 and $output{$_}++ for (@diseases);
	}

	for (keys %disease_score_hash) {
		my @words=@{$disease_score_hash{$_}};
		my @diseases=split(/\W+/,$words[1]);
		@diseases = map {lc $_;} @diseases;
		length($_)>2 and $output{$_}+= $words[0] for (@diseases);
	}

	for (sort { $output{$b} <=> $output{$a}  } keys %output) {
		print WORD_CLOUD $_."\t".$output{$_}."\n";
	}

	system("Rscript $work_path/wordcloud.R ${out}_${term}_wordcloud > ${out}_Rwordcloud.log");
}

sub Unique {
	my @words = @_;
	my %repeat_check = ();

	for my $each (@words) {
		if($repeat_check{$each}) {
			$repeat_check{$each}++;
		} else {
			$repeat_check{$each}=1;
		}
	}

	return keys %repeat_check;
}

sub TextStandardize {
	my $word=$_[0];
	
	$word=~s/^\W*(.*?)\W*$/$1/;
	$word=~s/'s\b//g;		#'
	$word=~s/\W+/ /g;
	$word=~s/\berthematosus\b/erythematosus/gi;
	$word=~s/\bshow all\b//ig;

	return $word;
}

sub printHeader{
	@_==2 or die "ERROR: printHeader() must take two oarguments";
	my ($fh, $if_predict) = @_;
	if (not $if_predict) {
		print $fh join("\t", qw(Rank Gene ID Score));
	} else {
		print $fh join("\t", qw(Rank Gene ID Score Status));
	}

	if($if_hi) {
		print $fh "\tHaploinsufficiencyScore";
	}

	if($if_rvis) {
		print $fh "\tGeneIntoleranceScore";
	}

	print $fh "\n";
}

=head1 SYNOPSIS

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
        -out <string>                   output file name prefix (default:out)
        -d, --directory                 compiled database directory (default is ./lib/compiled_database)
        -f, --file                      the input will be treated as file names(both diseases and genes)
        -p, --prediction                Use the Protein interaction and Biosystem database to predict unreported gene 
                                        disease relations (like HPRD human protein interaction, Biosystem database and so on)
        -ph, --phenotype                the input term is also treated as a phenotype, the HPO annotation and OMIM description would be used      
        -hi, --haploinsufficiency       use haploinsufficiency score as weight to prioritize dominant disease genes
        -it, --intolerance              use gene intolerance score as weight to prioritize severe disease genes
        --bedfile                       the bed file as a genomic region used for selection and annotation of the genes
        --buildver                      the build version (hg18 or hg19) to annotate the bedfile
        --wordcloud                     generates a wordcloud of the interpretated diseases if used (not working if you input 'all diseases')
        --logistic                      uses the weight based on the logistic modeling with four different complex diseases
        --gene                          the genes used to select the results (file name if -f command is used)    
        --exact                         choose if you want only exact match but not just a word match
        --addon                         the name of a user-defined add-on gene-disease mapping file (has to be in the ./lib/compiled_database)
        --addon_gg                      the name of user-defined add-on gene-gene mapping file (has to be in the ./lib/compiled_database)
        --addon_weight                  the weight of add-on gene-disease mapping
        --addon_gg_weight               the weight of add-on gene-gene mapping
        --hprd_weight                   the weight for genes found in HPRD
        --biosystem_weight              the weight for genes found in Ncbi Biosystem 
        --gene_family_weight            the weight for genes found in HGNC Gene Family
        --htri_weight                   the weight for genes found in HTRI Transcription Interaction Database
        --gwas_weight                   the weight for gene disease pairs in Gwas Catalog
        --gene_reviews_weight           the weight for gene disease pairs in Gene Reviews  
        --clinvar_weight                the weight for gene disease pairs in Clinvar
        --omim_weight                   the weight for gene disease pairs in OMIM
        --orphanet_weight               the weight for gene disease pairs in Orphanet    
        --hpo_gene_weight				the weight for gene phenotype pairs documented in HPO
        --nproc                         number of parallel processes (forks) requested by the user. The code uses as much parallelism as 
                                        allowed by the data. Setting this to 1 means no child processes are created.
        --use_precalc                   use HPO phenotypes expansion found in the precalculated database
  
 Function:       
          automatically expand the input disease term to a list of professional disease names, 
          get a prioritized genelist based on these disease names or phenotypes, score the genes.

 Notice: 
          If you input 'all diseases' for disease name, then every item in the gene_disease database
          will be used and no disease expansion will be conducted. 
          Addon Gene Gene file should be in the format "GENE A	GENE B	EVIDENCE	SCORE	PMID"
          Addon Gene Disease file should be in the format "GENE	DISEASE	DISEASE_ID SCORE	SOURCE"    
          
 Example:  
          disease_annotation.pl alzheimer -prediction -phenotype -logistic -out test1 -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
          disease_annotation.pl 'autism spectrum disorder' -prediction -logistic -out test2 -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25
          disease_annotation.pl HP:0001250 -prediction -phenotype -logistic -out test3
          disease_annotation.pl -file example_hpo.txt -phenotype -out test4
          disease_annotation.pl -file example_disease.txt -prediction -logistic -out test5

=cut
