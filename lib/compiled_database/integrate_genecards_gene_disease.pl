use strict;
open (GENECARDS, "../raw_gene_cards") or die;
open (GENE_ID, "./DB_HUMAN_GENE_ID") or die;
open (OUT, ">./DB_GENECARDS_GENE_DISEASE_SCORE") or die;
my %gene_transform = ();
my %output_hash =();
Generate_gene_transform();
for my $line (<GENECARDS>)
{
	chomp($line);
	my ($gene, $name,$locus,$disease) = split(/\t+/,$line);
    next if(not $gene_transform{uc $gene});
	$gene=$gene_transform{uc $gene};
	my @diseases= split(/, /,$disease);
	@diseases = map { $_=TextStandardize($_);$_=GetRidOfSusceptibility($_);lc $_;  } @diseases;
	$output_hash{join("\t",($gene, $_, "unknown", 0.25, "GENE_CARDS"))}=1 for @diseases;
}
print OUT $_."\n" for (keys %output_hash);

sub Generate_gene_transform{  
	my $i=0;
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

}
 	 
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
 	 
 	 
 	 
 	 