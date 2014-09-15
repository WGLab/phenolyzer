use strict;
open(INPUT,"DB_COMPILED_DISEASE_COUNT");
open(PHENOTYPE,"phenotype_terms.txt");
open(OUTPUT,">hot_disease_term.txt");
my %disease_count;
my $i = 0;
for(<INPUT>){
	if ($i==0) {$i++; next; }
	chomp;
	my($disease, $id, $count) = split("\t");
	$disease =~ s/^(.*?)((,)|([, ]+\()|(\s*-\s*\d+)).*$/\1/i;
	$disease =~ s/^(.*?)\s-\s.*$/\1/i;
	$disease =~ s/^([^\d]+)\s*-\s*/\1 /g;
	$disease_count{$disease} += $count;
	$disease =~ s/^([^\d]+?)-.*$/\1/i;
	$disease =~ s/^(.*?)(([, ]+syndrome)|([, ]+type)|([, ]+disease)).*$/\1/i;
	$disease =~ s/^(.*?)[, ]+\d+$/\1/i;
	$disease_count{$disease} += $count;
          }
$i=0;
for my $phenotype (<PHENOTYPE>)
{
	if ($i==0) {$i++; next;}
	chomp($phenotype);
	if($disease_count{$phenotype}) { $disease_count{$phenotype}++;next; }
	else {  $disease_count{$phenotype} = 1;  }
}          
     
     
	 for my $disease_term (sort {$disease_count{$a} <=> $disease_count{$b} } keys %disease_count){
	 	print OUTPUT $disease_term."\n"  if($disease_term!~/^[\W\dq]*$/);
	 }


	 
	 
	#system("cat TEMP DB_COMPILED_DISEASE_COUNT > TEMP2");
	