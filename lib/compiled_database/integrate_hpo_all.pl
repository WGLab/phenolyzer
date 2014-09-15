use strict;
open(FILE,"HPO_ontology_all.txt");
open(FILE2,">> DB_HPO_NUM_TERM");
my @words=();
my $line=join("\t",("HPO_ID","HPO_TERM"));
$line=$line."\n";
print FILE2 $line;
while (<FILE>){
	chomp;
	if (@words==2)
	{
	    print FILE2 join("\t",@words)."\n";
		@words=();
	}
	
	if(/^"?id:\s*HP:(\d*\w)/){
		push @words,$1;
		
	}
	if(/^"?name:\s*(.*\w)/){
		push @words,$1;
		
	}
}
