use strict;
use Carp;
open(OMIM,"DB_OMIM_GENE_DISEASE");
open(CLINVAR,"DB_CLINVAR_GENE_DISEASE");
open(ORPHANET,"DB_ORPHANET_GENE_DISEASE");
open(GENEREVIEWS,"DB_COMPILED_GENEREVIEWS");
open(GWAS, "DB_GWAS_GENE_DISEASE");

open(OUT,">DB_COMPILED_GENE_DISEASE_SCORE");
open(GENE_ID, "DB_HUMAN_GENE_ID") or die;
open(TEMP,">TEMP");
my %D_SCORE_OMIM=(
	  "1"=>0.25,
	  "2"=>0.5,
	  "3"=>1.00,
	  "4"=>0.75
);
my %S_SCORE_OMIM=(
       "C"=>1.00,
       "P"=>0.75,
       "L"=>0.5,
       "I"=>0.25
);
my %D_SCORE_ORPHANET=(
   'Disease-causing germline mutation(s) in' => 1.00,
   'Disease-causing somatic mutation(s) in' =>  1.00,
    'Major susceptibility factor in'  =>  0.75,
    'Part of a fusion gene in'    =>  0.5,
    'Modifying germline mutation in'  => 0.5,
    'Modifying somatic mutation in'  =>  0.5,
    'Role in the phenotype of'  =>  0.25,         
     'Candidate gene tested in' =>  0.25
);
my $ORPHANET_SOURCE_MAX = 7.0;
my $GENE_REVIEW_SCORE = 1.0;
my $GWAS_SCORE_WEIGHT = 0.25;
my %repeat_check=();
my %repeat_no_score_check=();
print OUT join("\t",qw/GENE	DISEASE	DISEASE_ID SCORE SOURCE/)."\n";

my $i=0;

my %gene_hash;
for my $line (<GENE_ID>)
{
	if($i==0) { $i++;next;  }
	chomp($line);
	my ($id, $gene, $synonyms) = split("\t", $line);
	$gene = uc $gene;
	if($synonyms eq "-") {$synonyms = "";
		$gene_hash{$gene} = $gene;
	}
	else {
		$synonyms =~s/\|/,/g;
		$gene_hash{$gene} = uc "$gene,$synonyms";
	}
}
seek GENE_ID, 0,0;
$i=0;
for my $line (<GENE_ID>)
{
	if($i==0) { $i++;next;  }
	chomp($line);
	my ($id, $gene, $synonyms) = split("\t", $line);
	$gene = uc $gene;
	if($synonyms eq "-") {
		next;
	}
	else {
		my @synonyms = split('\|', $synonyms);
		$synonyms = join(",", @synonyms); 
		for my $each (@synonyms)
		{
				$gene_hash{$each} = uc "$gene,$synonyms" if(not defined $gene_hash{$each});
		};
	}
}


$i=0;
for my $line (<OMIM>){   #GENE	DISEASE	MIM_NUMBER	SOURCE_CODE	LINGKAGE_INFO
	if($i==0) {$i++;next;}
	chomp($line);
	my @words=split("\t",$line);
	$words[0]=~/^(\w+)/;
	my $gene=uc $1;
	my $disease_id="OMIM:".$words[2];
	my $disease = $words[1];
	   $disease = getRidOfSusceptibility($disease);
	   $disease = eliminateNonWords($disease);
	   $disease =~ s/^\s*//;
	my $gene_score=$S_SCORE_OMIM{$words[3]}*$D_SCORE_OMIM{$words[4]};                 
	# $all_info{$line2}=1 if(not defined $all_info{$line2});                #This hash is to eliminate repeated
	my $repeat_line = join("\t",($gene_hash{$gene},$disease_id,$gene_score,"OMIM"));
	my $output = join("\t",($gene_hash{$gene},$disease,$disease_id,$gene_score,"OMIM"));
	print TEMP $output."\n" if ($gene_hash{$gene} and $gene_score and not $repeat_check{$repeat_line});
	$repeat_check{$repeat_line} = 1; 
}
$i=0;
for my $line (<CLINVAR>){                   #GENE	DISEASE	MIM_NUMBER	SOURCE_COUNT
	 if($i==0){$i++;next;}
	 chomp($line);
	 my @words=split("\t",$line);
	    $words[0] = uc $words[0];
	 my $disease = $words[1];
	    $disease = getRidOfSusceptibility($disease);
	    $disease = eliminateNonWords($disease);
	 my $gene_score=clinvarToScore($words[3]);
	 my $disease_id;
	 if ($words[2])
	 {
	 $disease_id="OMIM:".$words[2];
	 }
	 else{
	 $disease_id="Not available";
	 }
	 my $repeat_line = join("\t",($gene_hash{$words[0]},$disease_id,$gene_score,"CLINVAR"));
	 my $output = join("\t",($gene_hash{$words[0]},$disease,$disease_id,$gene_score,"CLINVAR"));
	 print TEMP $output."\n" if ($gene_hash{$words[0]} and $gene_score and not $repeat_check{$repeat_line});
	 $repeat_check{$repeat_line} = 1;
}
$i=0;
for my $line(<ORPHANET>){   #GENE	DISEASE	ORPHANET_NUMBER	SOURCE_COUNT	LINKAGE_INFO
	if($i==0){$i++;next;}
	chomp($line);
	my @words=split("\t",$line);
	   $words[0] = uc $words[0];
	my $disease = $words[1];
       $disease = getRidOfSusceptibility($disease);
       $disease = eliminateNonWords($disease);
       $words[4]=~s/^\s*(.*?)\s*$/$1/;
	my $gene_score=orphanetToScore($D_SCORE_ORPHANET{$words[4]},$words[3]);
	my $disease_id="ORPHANET:".$words[2];
	my $repeat_line = join("\t",($gene_hash{$words[0]}, $disease_id, $gene_score,"ORPHANET"));
	my $output = join("\t",($gene_hash{$words[0]}, $disease, $disease_id, $gene_score,"ORPHANET"));
	print TEMP $output."\n" if ($gene_hash{$words[0]} and $gene_score and not $repeat_check{$repeat_line});
	$repeat_check{$repeat_line} = 1;
}
$i=0;
for my $line(<GENEREVIEWS>){     #GENE	DISEASE	OMIM_NUMBER
	if($i==0){$i++;next;}
	chomp($line);
	my @words=split("\t",$line);
	next if ($words[0] =~ /not applicable/i);
	$GENE_REVIEW_SCORE = $words[3];
	$words[0] = uc $words[0];
	my $disease_id="OMIM:".$words[2];
	my $disease = $words[1];
       $disease = getRidOfSusceptibility($disease);
       $disease = eliminateNonWords($disease);
    my $repeat_line = join("\t",($gene_hash{$words[0]}, $disease_id, "GENE_REVIEWS"));  
	print TEMP join("\t",($gene_hash{$words[0]},$disease,$disease_id,$GENE_REVIEW_SCORE,"GENE_REVIEWS"))."\n" 
	 if (not $repeat_check{$repeat_line} and $gene_hash{$words[0]});
	 $repeat_check{$repeat_line} = 1;
	
}

$i=0;

for my $line (<GWAS>){          #GENE DISEASE PUBMED_NUMBER RAW_SCORE          
	if ($i==0) {$i++; next;}
	chomp ($line);
	my @words = split("\t",$line);
	my $disease = $words[1];
       $disease = getRidOfSusceptibility($disease);
       $disease = eliminateNonWords($disease);
	next if($words[0] eq "Intergenic"   or 
	        $words[0] eq "other genes"  or 
	        $words[0] eq "NR"           or
	        not $words[3]);
    $words[0] =~ s/^(.*?)\s.+$/\1/;  
    $words[0] = uc $words[0];      
	my $disease_id = "PUBMED:".$words[2];
	my $gene_score = $GWAS_SCORE_WEIGHT * $words[3];
	#my $gene_score = $GWAS_SCORE_WEIGHT * calculateGwasScore($words[3]);
	my $repeat_line = join ("\t", ($gene_hash{$words[0]},$disease, $disease_id, "GWAS"));
    if($gene_hash{$words[0]} and $disease and $disease_id)
    {
    if(not $repeat_no_score_check{$repeat_line})
    {
    $repeat_no_score_check{$repeat_line} = $gene_score;
    }
    else
    {
    $repeat_no_score_check{$repeat_line} = $gene_score if($gene_score > $repeat_no_score_check{$repeat_line});
    }
    }
 }
for my $line (keys %repeat_no_score_check)
{
	my @words = split("\t",$line);
	print TEMP join("\t", (@words[0..2], $repeat_no_score_check{$line}, "GWAS"))."\n";
	
}


#Finally sort items and choose unique ones
system("cat TEMP >> DB_COMPILED_GENE_DISEASE_SCORE");
system ("rm TEMP");

=cut
sub calculateGwasScore{                        # The smallest score is about 5.0, to make the range smaller, I add a sqrt here                 
	    @_ == 1 or die "input illegal!!";
	    my $raw_score = $_[0];
        return  ($raw_score/$GWAS_MAX_RAW_SCORE) ** 0.5;	
	
}
=cut

sub orphanetToScore{       #  D_SCORE and Source_number as input
	@_==2 or die "input illegal!!";
	my ($d_score,$source_num)=@_;
	return $source_num/$ORPHANET_SOURCE_MAX;
	
}

sub clinvarToScore{
   @_==1 or die "input illegal!!";
    return 1.00*0.25 if($_[0]==4);
    return 0.75*0.25 if($_[0]==3);
    return 0.50*0.25 if($_[0]==2);
    return 0.25*0.25 if($_[0]==1);
}

sub getRidOfSusceptibility{
	@_==1 or die "input illegal!!";
	$_[0] =~ s/^(.*?)\W*susc?eptibi?lity( to)?,?.*/$1/i;
	$_[0] =~ s/autism \d+/autism/gi;
	$_[0] =~ s/\berthematosus\b/erythematosus/gi;
	return $_[0];
}
sub eliminateNonWords{
	@_==1 or die "input illegal!!";
	$_[0] =~ s/'s\b//g;
	$_[0] =~ s/\W+/ /g;
	$_[0] =~ s/^\W+//;
	$_[0] =~ s/\W+$//;
	#get rid of with (without)
	$_[0] =~ s/\bwith\b.*$//gi;
	$_[0] =~ s/\bwithout\b.*$//gi;
	return $_[0];
}














