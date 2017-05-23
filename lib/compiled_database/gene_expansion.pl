use strict;
open(GENE_ID, "DB_HUMAN_GENE_ID") or die;
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
while(<>){
	s/[\r\n]+//g;
	my @words = split("\t");
	if (not defined $gene_hash{$words[0]}) { next; }
	print join("\t", ($gene_hash{$words[0]}, @words[1..$#words])  )."\n";
}

