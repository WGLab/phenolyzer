use strict;
open(HTRI, "../HTRIdb_data.txt");
open(OUTPUT, ">DB_HTRI_TRANSCRIPTION_INTERACTION");
my $i=0;
my %output = ();
my %evidence_out = ();
my $max_evidence_num = 0.0;
print OUTPUT join("\t", qw/TF TG EVIDENCE PUBMED SCORE/)."\n";

for my $line (<HTRI>)
{
	if($i==0) {$i++; next; }
	chomp($line);
	#print $line."\n";
	my @words = split("\t+", $line);
	my ($tf, $tg, $evidence, $pubmed) = @words[2,4,5,6];
	my $repeat_key = join("\t",($tf,$tg) );
	$output{$repeat_key}{"evidence"} .= $evidence.";";
    $output{$repeat_key}{"pubmed"} .= $pubmed.";";
    $output{$repeat_key}{"evidence_num"} = 1.0 if(not $output{$repeat_key});
    $output{$repeat_key}{"evidence_num"} += 1.0 if($output{$repeat_key});    
    $max_evidence_num = $output{$repeat_key}{"evidence_num"} if(  $max_evidence_num < $output{$repeat_key}{"evidence_num"});
}

for (keys %output){
	my ($tf_tg, $evidence, $pubmed, $evidence_num) = ($_, $output{$_}{"evidence"}, $output{$_}{"pubmed"}, $output{$_}{"evidence_num"}); 
	    $evidence =~s/;$//;
	    $pubmed   =~s/;$//;
    my $score = $evidence_num/$max_evidence_num;
    print OUTPUT join("\t", ($_, $evidence, $pubmed, $score))."\n";	
}
 
