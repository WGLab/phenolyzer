use strict;

open(INPUT,"../binary_protein_interaction.txt");
open(OUTPUT,">DB_COMPILED_BINARY_PROTEIN_INTERACTION_SCORE");
print OUTPUT (join "\t",qw/PROTEIN_1 PROTEIN_2 EVIDENCE SCORE PUBMED_ID/)."\n"; 
for my $line(<INPUT>)
{
	   chomp($line);
	   my @words=split("\t",$line);
	   $words[0]=~s/\$//g;
	   $words[1]=~s/\$//g;
	   $words[2]=~s/(^|;)\s*(.*?)\s*(;|$)/$2$3/;
	   $words[2]=~s/VT/in vitro/g;
	   $words[2]=~s/VV/in vivo/g;
	   $words[2]=~s/Y2H/yeast 2-hybrid/g;
	   print $words[2]."\n";
	   my @evidence=split(";",$words[2]);
	   my %evidence_hash=();
	   my $score=0;
	   for my $term (@evidence){
	   	 $score+=1 and $evidence_hash{$term}=1 if(!$evidence_hash{$term} and $term=~/^in vivo$/);
	   	 $score+=0.5 and $evidence_hash{$term}=1 if(!$evidence_hash{$term} and $term=~/^in vitro$/);
	   	 $score+=0.25 and $evidence_hash{$term}=1 if(!$evidence_hash{$term} and $term=~/^yeast 2-hybrid$/);
	   }
	   $score/=1.75;
	   $words[2]=join ";",(keys %evidence_hash);
	   print OUTPUT (join "\t",($words[0],$words[1],$words[2],$score,$words[3]))."\n";
}