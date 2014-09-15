use strict;
open(FILE,"DB_HPO_NUM_TERM");
open(FILE2,"HPO_NUM_USED");
open(OUTPUT,">DB_HPO_NUM_TERM_USED");
my @list1=<FILE>;
my $i=0;
for my $line(<FILE2>){
	 chomp($line);
	 
	 for my $line2(@list1[$i..$#list1]){
	 	  chomp($line2);
	 	  my @words2=split("\t",$line2);
	 	  $i++;
	 	  if($line==$words2[0]){
	 	  	  print OUTPUT $line2."\n";
                       last;        }
	 	  		 	 }
	
	
}