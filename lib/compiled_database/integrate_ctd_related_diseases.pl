use strict;
open(CTD,"DB_COMPILED_CTD_DISEASES");
open(DISEASE,"DB_COMPILED_DISEASE_COUNT");
open(OUT,">DB_COMPILED_CTD_DISEASES_USED_TEMP");
my @disease_ctd=<CTD>;
my @disease=<DISEASE>;
my %diseases=();
my %disease_orig=();
my $i=0;
for my $line (@disease){
	if($i==0){$i++;next;}
    chomp($line);   
	my @words=split("\t",$line);            #build a hash for DISEASE
	my $name=$words[0];
	my $key = lc $name;
	$key=~s/\W+/_/;
	print $name."\n";
	$disease_orig{$key}=$name;
}
$i=0;
for  my $line (@disease_ctd){
	if($i==0){$i++;next;}
	chomp($line);
	my @words=split("\t",$line);
	my @synonyms=split('\|',$words[1]);
	my $name = $words[0];
	my $key = lc $name;
	$key=~s/\W+/_/;
	if(defined $disease_orig{$key})
	{
		$diseases{$line}=1 
		unless defined $diseases{$line};
	}
		else{
			for (@synonyms){
				 s/\W+/_/;   
				if(defined $disease_orig{$_}){
					$diseases{$line}=1 
			#	and print STDERR  "NOTICE:".keys(%diseases)." items name found!! $_ \n"
			     unless defined $diseases{$line};
	             }
          	}
       }
	}

$i=0;
for my $line(@disease_ctd) {
	if($i==0){$i++;next;}
	 chomp($line);
	 my @words=split("\t",$line);   #$words[2] is tree number, then tell if the previous found tree is its children
	 my $father = $words[2];
	 for (keys %diseases){
	 	 chomp();
	 	 my @words2=split("\t");
	 	 my $children = $words2[2];
	 	if($children =~/(^|\|)  $father  [^\|\/]/ix)  #$words[2] should be father, tell if $words2[2] is children or not
	    {
          $diseases{$line}=1 
          unless defined $diseases{$line};
	    }
      }
}

print OUT join("\t",qw/DISEASE SYNONYMS TREE_NUM PARENT_TREE_NUM/)."\n";
print OUT $_."\n" for (keys %diseases);
system ("mv DB_COMPILED_CTD_DISEASES_USED_TEMP DB_COMPILED_CTD_DISEASES_USED");




