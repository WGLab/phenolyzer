use strict;
system("cut -f 2,3 DB_COMPILED_GENE_DISEASE_SCORE >temp");   #GENE	DISEASE	DISEASE_ID	SCORE	SOURCE
system("sort temp|uniq -c|sort -n -k 1 >temp2");
open(DISEASE,"temp2");
open(DISEASE_OUT,">DB_COMPILED_DISEASE_COUNT");
print DISEASE_OUT join("\t",qw/DISEASE DISEASE_COUNT/)."\n";
while(<DISEASE>){
	 chomp();
	 $_=~/^\s*(\d+)\s(.*)$/;
	 if($2 ne "DISEASE"){
	 print DISEASE_OUT join("\t",($2,$1))."\n";
	 }
}
system("rm -f temp");
system("rm -f temp2");
