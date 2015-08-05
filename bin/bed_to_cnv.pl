use strict;
use warnings;
while(<>){
	s/[\r\n]+//g;
	my @words = split("\t");
	next if @words <=3;
	my ($chr, $start, $end) = @words[0..2];
	print "$chr:$start-$end\n";
}