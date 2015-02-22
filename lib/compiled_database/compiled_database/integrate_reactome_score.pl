use strict;
open (GENE_ID,"Ensemble_homo_sapians_geneID.txt");
open (REACTOME,"Reactome_homo_sapian.interactions.txt");
open (OUTPUT,">DB_COMPILDE_REACTOME_SCORE");
my %gene_hash=();
my %line_hash=();
my $i=0;
for my $line(<GENE_ID>){
	if($i==0){$i++;next;}
	chomp($line);
	my @words = split(",",$line);
	$gene_hash{$words[0]}=$words[2] if not defined $gene_hash{$words[0]};
}
print OUTPUT join("\t",qw/PROTEIN1 PROTEIN2 INTERACTION SCORE/)."\n";
$i=0;
for my $line(<REACTOME>){
	if($i==0){$i++;next;}
	chomp($line);
	my @words = split("\t",$line);
	my $gene1 = $words[1];
	my $gene2 = $words[4];
	my $reaction = $words[4];
	next if($gene1=~/^\s*$/);
	next if($gene2=~/^\s*$/);
	$gene1=~s/^ENSEMBL:(.+)$/$1/;
	$gene2=~s/^ENSEMBL:(.+)$/$1/;
	$gene1=$gene_hash{$gene1};
	$gene2=$gene_hash{$gene2};
	$line=join("\t",($gene1,$gene2,$reaction))."\n";
	next if(defined $line_hash{$line});
	$line_hash{$line}=1;
	$score
	print OUTPUT $line;
}