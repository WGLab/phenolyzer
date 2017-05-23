use strict;
open (DIS_GEN_NET, "../DisGenNet_Curated.txt") or die;
open (GAD, "../GAD_gene_disease.txt") or die;
open (GENE_ID, "./DB_HUMAN_GENE_ID") or die;
open (DISGENET_OUT, ">DB_DISGENET_GENE_DISEASE_SCORE") or die;
open (GAD_OUT,  ">DB_GAD_GENE_DISEASE_SCORE") or die;
my $i=0;
my %gene_transform;
for my $line (<GENE_ID>)
{
	if($i==0) { $i++;next;  }
	chomp($line);
	my ($id, $gene, $synonyms) = split("\t", $line);
		$gene_transform{$gene} = uc $gene;
		$gene_transform{uc $gene} = uc $gene;
		$gene_transform{$id} = uc $gene;
}
seek GENE_ID, 0,0;
$i=0;
for my $line (<GENE_ID>)
{
	if($i==0) { $i++;next;  }
	chomp($line);
	my ($id, $gene, $synonyms) = split("\t", $line);
	if($synonyms eq "-") {
		next;
	}
	else {
		my @synonyms = split('\|', $synonyms);
		for my $each (@synonyms)
		{
				$gene_transform{uc $each} = uc $gene if(not defined $gene_transform{$each});
		};
	}
}


open (TESTING_INDEX, ">TESING_INDEX");
my %repeat_check = ();
my $i=0;
for my $line (<DIS_GEN_NET>){
	if($i==0) {$i++; next;}
	chomp($line);
	my @words = split("\t", $line);	
	my ($gene, $disease, $score, $disease_id) = @words[6,9,1,8];
	$disease =~s/"//g;
	 $disease = TextStandardize($disease);
	 $disease = GetRidOfSusceptibility($disease);
	$gene=$gene_transform{uc $gene};
	my $repeat_line = join ("\t", ($gene, $disease_id, $score, "DISGENET"));
    my $output = join ("\t", ($gene, $disease, $disease_id, $score, "DISGENET"));
    print DISGENET_OUT $output."\n" if ($score and not $repeat_check{$repeat_line});
    $repeat_check{$repeat_line} = 1;
}
$i=0;
my %conflict_check = ();
for my $line (<GAD>){
    if($i<=2) {$i++; next;}	
	chomp($line);
	my @words = split("\t", $line);
	my ($gene, $disease, $if_association,$pubmed_id,$omim_id) = @words[8,2,1,13,29];
	$gene=$gene_transform{uc $gene};
	$pubmed_id;
	$omim_id = "OMIM:".$omim_id;
	my @diseases = split('\|',$disease);
	$disease = $diseases[0];
	$disease =~s/"|('s?)//g;
	$disease =~s/\bII\b/2/g;
	$disease =~s/\bIII\b/3/g;
	$disease =~s/\bI\b/1/g;
	 $disease = TextStandardize($disease);
	 $disease = GetRidOfSusceptibility($disease);
	 $disease = lc $disease;
	if($if_association and $gene and $disease!~/^\W*$/ and $pubmed_id >0 )
	{	
		next if($if_association !~ /^Y$/);
		my $gene_disease = join("\t",($gene,$disease));
        $conflict_check{$gene_disease} = $pubmed_id  if (not $conflict_check{$gene_disease});
        $conflict_check{$gene_disease}.= " ".$pubmed_id  if ( $conflict_check{$gene_disease}); 
       		
	};
 } 
 my @output;
 my %disease_check = ();
 my $max_count = 0;
 for my $each (keys %conflict_check){ 
 	 my @pubmed_ids = split(" ", $conflict_check{$each});
 	 $max_count = scalar(@pubmed_ids) if(scalar(@pubmed_ids)>$max_count);
 }
 for my $each (keys %conflict_check){ 
 my @pubmed_ids = split(" ", $conflict_check{$each});	
 push @output, join("\t",($each, "PUBMED:".$conflict_check{$each}, 0.25*(@pubmed_ids)/$max_count, "GAD"));
 }
@output = sort {my @words1 = split("\t",$a);
 	            my @words2 = split("\t",$b);
 	             $words1[1].$words1[2] cmp $words2[1].$words2[2]; }  @output;
 	 print GAD_OUT $_."\n" for (@output);

 	 
sub TextStandardize {
	my $word=$_[0];
	$word=~s/^\W*(.*?)\W*$/$1/;
	$word=~s/'s\b//g;
	$word=~s/\W+/ /g;
	$word=~s/\berthematosus\b/erythematosus/gi;
	$word=~s/\bshow all\b//ig;
	#get rid of with (without)
	$word =~ s/\bwith\b.*$//gi;
	$word =~ s/\bwithout\b.*$//gi;
	return $word;
} 
sub GetRidOfSusceptibility{
	@_==1 or die "input illegal!!";
	$_[0] =~ s/^(.*?)\W*susc?eptibi?lity( to)?,?.*/$1/i;
	$_[0] =~ s/autism \d+/autism/gi;
	return $_[0];
}
 	 
 	 
 	 
 	 