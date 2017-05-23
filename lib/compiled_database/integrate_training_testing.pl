use strict;
open (DIS_GEN_NET, "../DisGenNet_Curated.txt");
open (GAD, "../GAD_gene_disease.txt");
open (TRAINING, ">DB_TRAINING_GENE_DISEASE_SCORE");
open (TESTING,  ">DB_TESTING_GENE_DISEASE_SCORE");
open (TESTING_INDEX, ">TESING_INDEX");
my %repeat_check = ();
my $i-0;
for my $line (<DIS_GEN_NET>){
	if($i==0) {$i++; next;}
	chomp($line);
	my @words = split("\t", $line);	
	my ($gene, $disease, $score, $disease_id) = @words[6,9,1,8];
	$disease =~s/"//g;
	my $repeat_line = join ("\t", ($gene, $disease_id, $score, "DISGENET"));
    my $output = join ("\t", ($gene, $disease, $disease_id, $score, "DISGENET"));
    print TRAINING $output."\n" if ($score and not $repeat_check{$repeat_line});
    $repeat_check{$repeat_line} = 1;
}
$i=0;
my %conflict_check = ();
for my $line (<GAD>){
    if($i<=2) {$i++; next;}	
	chomp($line);
	my @words = split("\t", $line);
	my ($gene, $disease, $if_association ) = @words[8,5,1];
	my @diseases = split('\|',$disease);
	$disease = $diseases[0];
	$disease =~s/"|('s?)//g;
	$disease =~s/\bII\b/2/g;
	$disease =~s/\bIII\b/3/g;
	$disease =~s/\bI\b/1/g;
	if($if_association and $gene and $disease!~/^\W*$/)
	{
		
		next if($if_association !~ /^Y|N$/);
		my $gene_disease = $gene."\t".$disease;
        $conflict_check{$gene_disease} = $if_association  if (not $conflict_check{$gene_disease});
        $conflict_check{$gene_disease}.= $if_association  if ( $conflict_check{$gene_disease} ne $if_association); 
       		
	};
 }
 my @output;
 my %disease_check = ();
 $i = 0;
 for my $each (keys %conflict_check){
 	if(length($conflict_check{$each})==1)
 	{
 	 push @output,$each."\t".$conflict_check{$each}; 
     my ($gene, $disease) = split("\t", $each);
     $i+=1 and $disease_check{$disease}{"index"} = sprintf('%04s',$i) if(not $disease_check{$disease}{"index"});
     if( not $disease_check{$disease}{"count"})
     {
     $disease_check{$disease}{"count"} = 1; 
     }
 	 else{
 	 $disease_check{$disease}{"count"} +=1;	
        }
 	}
 }
@output = sort {my @words1 = split("\t",$a);
 	            my @words2 = split("\t",$b);
 	             $words1[1].$words1[2] cmp $words2[1].$words2[2]; }  @output;
 	 print TESTING $_."\n" for (@output);
@output = sort { $disease_check{$b}{"count"} <=> $disease_check{$a}{"count"};  
	                   }  keys %disease_check;
 	 print TESTING_INDEX $_."\t".$disease_check{$_}{"count"}."\t".$disease_check{$_}{"index"}."\n" for (@output);
 	 
 	 
 	 
 	 
 	 
 	 
 	 
 	 