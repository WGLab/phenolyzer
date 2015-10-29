use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
our($maf,$everything);
my $dirname = dirname(__FILE__);
GetOptions('maf=s'=>\$maf, 
           'everything'=>\$everything);
defined $maf or $maf =0.02;
@ARGV == 2 or die("The arguments should be <phenolyzer_gene_list> <wANNOVAR query.exome/genome_summary.txt>");
open(PHENOLYZER,$ARGV[0]) or die("ERROR:Can't open $ARGV[0]!");
open(ANNOVAR,  $ARGV[1]) or die("ERROR:Can't open $ARGV[1]!");
*OUTPUT = *STDOUT;
my %phenolyzer = ();
my %index_hash = ();
my %output = ();
my %disease_hash = processDisease();
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
   $head =~ s/[\n\r]+//g;
my @words = split("\t", $head);
my (@kg_index, @exac_index, @esp6500_index);
my (@kg_headers, @exac_headers, @esp6500_headers);
for(my $i=0;$i<@words;$i++)
{
      $index_hash{"gene"} = $i if($words[$i] =~ /^gene\b/i);  
      $index_hash{"RadialSVM"} = $i if($words[$i]=~/(^RadialSVM_score$)|(^ljb\d+_metasvm$)/i);
      if($words[$i]=~/^1000g/i){
      	   push @kg_index,$i;
      	   push @kg_headers,$words[$i];
      }
      if($words[$i]=~/^exac/i){
      	   push @exac_index,$i;
      	   push @exac_headers,$words[$i];   
      }
      if($words[$i]=~/^esp6500/i){
      	   push @esp6500_index,$i; 
      	   push @esp6500_headers,$words[$i];
      }
}
my $min_svm = "+inf";
my $max_svm = "-inf";

while(<ANNOVAR>)
{
     s/[\r\n]+//g;
     next if($_!~/\b(exonic|splicing)\b/);
     next if(/\bsynonymous SNV\b/);
     my @words = split("\t");
     my $gene = uc $words[$index_hash{"gene"}];
     my $RadialSVM;
     my $next = 0;
     for ((@kg_index, @esp6500_index, @exac_index)){
     	$next =1 unless ($words[$_]=~/^\W*$/ or $words[$_]<=$maf);
     }
     next if ($next);
        if($words[$index_hash{"RadialSVM"}]){
        	 $RadialSVM = $words[$index_hash{"RadialSVM"}];
        }
        else{
        	$RadialSVM = 0.0001;
        }
        
     my $pheno = $phenolyzer{$gene};
     $pheno = 0 if( not $pheno);
     my $variant = join("\t", (@words[0..4],$gene));
     my @RadialSVMs = split(",", $RadialSVM);
     $RadialSVM = $RadialSVMs[0];
     $output{$variant}{"RadialSVM"} = $RadialSVM if($RadialSVM);
     $output{$variant}{"Phenolyzer"} = $pheno;
      $output{$variant}{'disease'} = "\t";
     $output{$variant}{'disease'} = $disease_hash{$gene} if defined $disease_hash{$gene};
     $output{$variant}{'MAF'} = [@words[@kg_index],@words[@exac_index],@words[@esp6500_index]];
     if($everything){
     	$output{$variant}{'everything'}=join("\t", @words);
     }
     if($RadialSVM and $RadialSVM ne "."){
     $min_svm = $RadialSVM if($RadialSVM < $min_svm); 
     $max_svm = $RadialSVM if($RadialSVM > $max_svm); 
     } 
}
$min_svm-=0.0001;
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
if(not $everything)
{
print OUTPUT join("\t", (qw/Chr Start End Ref Alt Gene RadialSVM_Score Phenolyzer_Score Integrated_Score/,
                         @kg_headers, @exac_headers, @esp6500_headers))."\n"; 
}
else{
	print OUTPUT join("\t", (qw/Chr Start End Ref Alt Gene RadialSVM_Score Phenolyzer_Score Integrated_Score
	                        Disease DiseaseID/,
                           $head))."\n"; 
}
for (sort { $output{$b}{"Score"} <=>  $output{$a}{"Score"}  }keys %output)
{
	   next if(not defined $output{$_}{"Score"});
	   if(not $everything)
	   {
       print OUTPUT join("\t", ($_, $output{$_}{"RadialSVM"}, $output{$_}{"Phenolyzer"}, $output{$_}{"Score"}
                                  , @{$output{$_}{'MAF'}}))."\n";
	   }
	   else{
	   	print OUTPUT join("\t", ($_, $output{$_}{"RadialSVM"}, $output{$_}{"Phenolyzer"}, $output{$_}{"Score"},
                                 $output{$_}{"disease"} , $output{$_}{"everything"}))."\n";
	   }
}
sub processDisease{
	open(DISEASE,"$dirname/lib/compiled_database/DB_COMPILED_GENE_DISEASE_SCORE") or die;
	my %output = ();
	for my $line (<DISEASE>){
		next if($line=~/^GENE\t/);
		$line =~s/[\r\n]+//g;
		my ($gene,$disease,$id,$score,$source) = split("\t",$line);
	    next if($score<0.75 or $source eq "GWAS");
	    my @genes = split(",", $gene);
	    if (not defined $output{$genes[0]}) {
	    	$output{$genes[0]} = "$disease\t$id" ;  }
	    else{
	    	my @entry = split('\t',$output{$genes[0]});
	        $entry[0] .= "|$disease";
	        $entry[1] .= "|$id";
	    	$output{$genes[0]} = join("\t",@entry);  }
	} 
	return %output;
}

=head1 SYNOPSIS

 The arguments should be <phenolyzer_gene_list> <wANNOVAR query.exome/genome_summary.txt> [arguments]
     arguments:
           -maf           The minor allele frequency for filtering 
           -everything    Print every annotation of the variant
=cut


