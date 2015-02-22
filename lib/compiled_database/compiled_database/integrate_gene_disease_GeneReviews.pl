use strict;
open (GENE,"../GeneReview_NBKid_shortname_genesymbol.txt");
open (OMIM,"../GeneReview_NBKid_shortname_OMIM.txt");
open (MORBID_MAP,"../omim_morbidmap_dev");
open (TEMP,">DB_COMPILED_GENEREVIEWS_TEMP");
open (OUTPUT,">DB_COMPILED_GENEREVIEWS");
open (OMIM_ID_DISEASE,">DB_COMPILED_OMIM_ID_DISEASE");
open (TEMP2,">temp");
# Keys shortname, values is a hash (gene=>list of related genes, omim=>list of disease names)
my %shortname_hash=();           
my %omim_hash=();
for my $line (<MORBID_MAP>){
	chomp($line);
	my @words = split ('\|',$line);
	
	$words[0] =~ /^  [\[\?\{]*  (.*?) [\]\}]* (?:,\s(\d{3,}))?\s\(\d\)$/x;
	push @{$omim_hash{$2}},$1 if($2);
	push @{$omim_hash{$words[1]}},$1 if(not $2);
}

print OMIM_ID_DISEASE join("\t",qw/OMIM_ID DISEASE/)."\n";
for my $omim_id(keys %omim_hash) {
	my @diseases = map {  $_=GetRidOfSusceptibility($_); $_; }  @{$omim_hash{$omim_id}};
	print "@diseases"."\n";
	my $disease = join(";",@diseases);
	   $disease = GetRidOfAnnotations($disease);
	 	print TEMP2 $omim_id."\t".$disease."\n";
	@{$omim_hash{$omim_id}} = @diseases; 	
}
system ("sort temp|uniq >> DB_COMPILED_OMIM_ID_DISEASE");
system ("rm temp");
my $i=0;
for my $line(<GENE>){
	if($i==0){$i++;next;}
	chomp($line);
	my @words=split("\t",$line);
	push @{$shortname_hash{$words[1]}{gene}},$words[2];
}
$i=0;
for my $line(<OMIM>){
	if($i==0){$i++;next;}
	chomp($line);
	my @words=split("\t",$line);
	if (defined $omim_hash{$words[2]}){	
		my $disease = ${$omim_hash{$words[2]}}[0];
		push @{$shortname_hash{$words[1]}{omim}},$disease."\t".$words[2];
	}
}
print OUTPUT join("\t",qw/GENE DISEASE OMIM_NUMBER SCORE/)."\n";
for my $name (keys %shortname_hash){
	  for my $gene(@{$shortname_hash{$name}{gene}}){
	  	for my $disease_omim(@{$shortname_hash{$name}{omim}}){
	  		my $num = @{$shortname_hash{$name}{omim}};
	  		$disease_omim = GetRidOfAnnotations($disease_omim);
	  		print TEMP join("\t",($gene,$disease_omim, 1.0/$num))."\n";
      }
} 
}
sub GetRidOfAnnotations{
	@_==1 or die;
	my $disease = $_[0];
	$disease =~s/\(.*?\)//g;
	return $disease;
}
sub GetRidOfSusceptibility{
	@_==1 or die "input illegal!!";
	$_[0] =~ s/^(.*?)\W*susc?eptibi?lity( to)?,?.*/$1/i;
	$_[0] =~ s/autism \d+/autism/gi;
	return $_[0];
}

system ("sort DB_COMPILED_GENEREVIEWS_TEMP|uniq >>DB_COMPILED_GENEREVIEWS");
system ("rm DB_COMPILED_GENEREVIEWS_TEMP") if (-f "DB_COMPILED_GENEREVIEWS_TEMP");
