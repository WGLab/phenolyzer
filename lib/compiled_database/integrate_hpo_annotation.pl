use strict;
open(OUTPUT,">TEMP") or die;
open(OUTPUT1,">DB_HPO_ANNOTATION") or die;

print OUTPUT1 join("\t",("HPO_ID","SOURCE","HPO_DISEASE_NAME","FREQUENCY"))."\n";
processing("../hpo_annotation_1.txt");
processing("../hpo_annotation_2.txt");
processing("../hpo_annotation_hpoteam.txt");

sub processing
{
my $filename = $_[0];
open(FILE,$filename) or die;
while(<FILE>){
	chomp;
	my @line=split("\t");
	$line[5]=join (":",($line[0],$line[1]));
	$line[2]=~s/^"?(?:\W?\d{3,}\s+)?(.+\w)"?$/\1/;
	$line[2]=getRidOfSusceptibility($line[2]);
	$line[2]=~s/(^|;)\W*(.*?)\W*/$1$2$3/g;

	$line[4]=~s/^HP:(\d+)$/\1/;
	
	print OUTPUT join("\t",($line[4],$line[5],$line[2],$line[8]))."\n";
}
}
system("sort TEMP|uniq >> DB_HPO_ANNOTATION");
system("rm TEMP") if (-f "TEMP");

sub getRidOfSusceptibility{
	@_==1 or die "input illegal!!";
	$_[0] =~ s/^(.*?)\W*susc?eptibi?lity( to)?,?\s*\d*/$1/i;
	return $_[0];
	
}
