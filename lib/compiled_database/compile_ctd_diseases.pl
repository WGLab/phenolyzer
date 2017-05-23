use strict;
use warnings;
open(CTD, "../CTD_diseases.txt") or die;
open(CTD_OUT, ">DB_COMPILED_CTD_DISEASES") or die;
print CTD_OUT join("\t", ("DISEASE","SYNONYMS","TREE_NUM","PARENT_TREE_NUM"))."\n";
while(<CTD>){
	next if(/^#/);
	s/[\r\n]+//g;
	my @words = split("\t");
	my ($disease, $synonym, $tree, $parent_tree) = @words[0, 7, 5, 6];
	$parent_tree = '' if not $parent_tree;
	$synonym = '' if not $synonym;
	my @synonyms = split(/\|/, $synonym);
	$disease =getRidOfSusceptibility($disease);
	$disease = eliminateNonWords($disease);
	my %synonyms_check;
	for (my $i=0;$i<@synonyms;$i++){
		$synonyms[$i] = getRidOfSusceptibility($synonyms[$i]);
		$synonyms[$i] = eliminateNonWords($synonyms[$i]);
	}
	$synonym = join('|', @synonyms);
	print CTD_OUT join("\t", ($disease, $synonym, $tree, $parent_tree))."\n";
	
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