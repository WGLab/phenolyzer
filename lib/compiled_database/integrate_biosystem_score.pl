use strict;
use Carp;
open(GENE_ID,"../Homo_sapiens.gene_info") or die "can't open Homo_sapiens.gene_info!!";
open(BIOSYSTEM,"../biosystems_gene") or die "can't open biosystems_gene!!";
open(BIOSYSTEM_INFO, "../biosystem_bsid_to_info.txt") or die "Can't open biosystem_info file!!";
open(OUTPUT,">DB_COMPILED_BIOSYSTEM_SCORE");
my %gene_hash=();
my %biosystem_id_type;
for (<BIOSYSTEM_INFO>)
	    {
	    	chomp;
	    	my @words = split("\t");
	    	my ($biosystem_id, $name) = ($words[0], $words[3]); 	
	    	   $biosystem_id_type{$biosystem_id} = $name 
	    	   if (not defined $biosystem_id_type{$biosystem_id}); 
	    	
	    }
my $i=0;
for my $line(<GENE_ID>){
	if($i==0){$i++;next;}
	chomp($line);
	my @words=split("\t",$line);
	$gene_hash{$words[1]}=$words[2];
}
print OUTPUT join("\t",qw/BIOSYSTEM_ID GENE SCORE NAME/)."\n";
my @biosystem = <BIOSYSTEM>;
my $MAX = 0;
for my $line (@biosystem){
	my @words=split("\t",$line);
	$MAX = $words[2] if ($MAX < $words[2]);
	}

for my $line(@biosystem){
	my @words=split("\t",$line);
	my $score=$words[2]/$MAX;
	print OUTPUT join("\t",($words[0],$gene_hash{$words[1]},$score, $biosystem_id_type{$words[0]}))."\n" if defined $gene_hash{$words[1]};
    
}
