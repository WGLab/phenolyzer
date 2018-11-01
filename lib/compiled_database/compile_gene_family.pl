use strict;
open (GENE_FAMILY, "../entrez_id_gene_family.txt") or die $!;
open (ENTREZ, "DB_HUMAN_GENE_ID") or die $!;
open (OUTPUT, ">DB_HGNC_GENE_FAMILY") or die $!;
my %id_gene;
my $i = 0;
for my $line ( <ENTREZ>)
{
	if ($i==0){ $i++;next;}
    my ($id, $gene, $synonym) = split("\t", $line);
    $id_gene{$id} = $gene; 
}
$i=0;
print OUTPUT join("\t", qw/GENE GENE_FAMILY_TAG DESCRIPTION/)."\n";
for my $line (<GENE_FAMILY>)
{
	if ($i==0){ $i++;next;}
	chomp($line);
	my ($id, $tag, $description) = split("\t", $line);
	next if (not ($id and $tag and $description));
	$description =~ s/^"?(.*?)"?$/$1/;
	$description =~ s/[",]//g ;
	my @tags = split(",",$tag);
	$tag = $tags[0];
	print OUTPUT join("\t", ($id_gene{$id}, $tag, $description))."\n";
	
}
