use strict;
open(CTD,"CTD_diseases.txt");

open(CTD_OUT,">>DB_COMPILED_CTD_DISEASES");
print CTD_OUT join("\t",qw/DISEASE SYNONYMS TREE_NUM PARENT_TREE_NUM/)."\n";
while(<CTD>){
	if(/^#/){next;}
	chomp();
	my @words=split("\t");
	$words[0]=~s/^"*(.*?)"$/\1/;
	$words[7]=~s/^"(.*?)"$/\1/;
	print CTD_OUT join("\t",($words[0],$words[7],$words[5],$words[6]))."\n";
	}

#system("sort -k 3 -d temp >>DB_COMPILED_CTD_DISEASES");
