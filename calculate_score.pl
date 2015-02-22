use strict;
use warnings;
@ARGV == 2 or die("The arguments should be <phenolyzer_gene_list> <wANNOVAR query.exome/genome_summary.txt>");
open(PHENOLYZER,$ARGV[0]) or die("ERROR:Can't open $ARGV[0]!");
open(ANNOVAR,  $ARGV[1]) or die("ERROR:Can't open $ARGV[1]!");
*OUTPUT = *STDOUT;
my %phenolyzer = ();
my %index_hash = ();
my %output = ();
my $i=0;
while(<PHENOLYZER>)
{
     if($i==0) {$i++; next; }    
     chomp();	
	 my @words = split("\t");
	 my ($gene, $score) = @words[1,3];
	 $phenolyzer{$gene}=$score;
}
my $head = <ANNOVAR>;
my @words = split("\t", $head);
for(my $i=0;$i<@words;$i++)
{
      $index_hash{"gene"} = $i if($words[$i] eq "Gene.refgene");  
      $index_hash{"RadialSVM"} = $i if($words[$i]=~/(^RadialSVM_score$)|(^ljb\d+_metasvm$)/i);
}
my $min_svm = "+inf";
my $max_svm = "-inf";

while(<ANNOVAR>)
{
     chomp();
     my @words = split("\t");
     my $gene = $words[$index_hash{"gene"}];
     my $RadialSVM;
        if($words[$index_hash{"RadialSVM"}]){
        	 $RadialSVM = $words[$index_hash{"RadialSVM"}]
        }
        else{
        	$RadialSVM = 0;
        }
        
     my $pheno = $phenolyzer{$gene};
     $pheno = 0 if( not $pheno);
     my $variant = join("\t", (@words[0..4],$gene));
     my @RadialSVMs = split(",", $RadialSVM);
     $RadialSVM = $RadialSVMs[0];
     $output{$variant}{"RadialSVM"} = $RadialSVM if($RadialSVM);
     $output{$variant}{"Phenolyzer"} = $pheno;
     
     if($RadialSVM and $RadialSVM ne "."){
     $min_svm = $RadialSVM if($RadialSVM < $min_svm); 
     $max_svm = $RadialSVM if($RadialSVM > $max_svm); 
     } 
}
$min_svm-=0.001;
for  my $variant (keys %output){
     my $RadialSVM = $output{$variant}{"RadialSVM"};
     my $pheno = $output{$variant}{"Phenolyzer"};
     if($RadialSVM and $RadialSVM ne ".") {
         $RadialSVM = ($RadialSVM-$min_svm)/($max_svm-$min_svm);
     }
     else {  $RadialSVM=0; }
     $output{$variant}{"RadialSVM"} = $RadialSVM;
     $output{$variant}{"Score"} = $pheno*$RadialSVM;
}
print OUTPUT join("\t", qw/Chr Start End Ref Alt Gene RadialSVM_Score Phenolyzer_Score Integrated_Score/)."\n"; 
for (sort { $output{$b}{"Score"} <=>  $output{$a}{"Score"}  }keys %output)
{
	   next if(not $output{$_}{"Score"});
       print OUTPUT join("\t", ($_, $output{$_}{"RadialSVM"}, $output{$_}{"Phenolyzer"}, $output{$_}{"Score"}))."\n";
}

