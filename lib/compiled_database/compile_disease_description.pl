use strict;
open (OMIM_DISEASE, "../omim_disease.txt") or die;
open (OUTPUT, ">DB_OMIM_DESCRIPTION") or die;
open(OMIM, "DB_COMPILED_OMIM_ID_DISEASE") or die;

my %omim_id_disease;
for my $line (<OMIM>)
{
    next if($line=~/^OMIM_ID/);	
	chomp($line);
	my ($id, $disease) = split("\t", $line);
	$omim_id_disease{$id} = $disease;
}

my $read_description = 0;
my $read_id = 0;
my %omim_description;
my $omim_id;
my %output;

for my $line (<OMIM_DISEASE>)
{
        chomp($line);	
        next if(not $line);   
	 
	    if($line eq '*FIELD* NO')
	    {
	     $read_id = 1;	
	     $read_description = 0;
	     next;
	    }
	    if($line eq '*FIELD* CS')
	    {
	     $read_description = 1;
	     $read_id = 0;
	     next;
	    }
 	    if($line =~ /^\*/ and $line ne '*FIELD* CS')
        {
          $read_description = 0;	
          $read_id = 0;
          next;
        }
        if($read_id)
        {
           $omim_id = $line;	
           #print $line."\n";
        }
        if($read_description)
        {
        	next if($line=~/(^.*?:$)|(\[.*?\];?)/);
        	$line =~ s/\(.*?\)//g;
        	$line =~ s/^\W+(.*?)\W+/$1/;
        	$line =~s/;//g;
        	$omim_description{$omim_id}.=$line.";";
        	#print $line."\n";
        }
        
}
print OUTPUT join("\t", qw/OMIM_ID DISEASE DESCRIPTION/)."\n";
for my $omim_id (sort keys %omim_description)
{
	my @descriptions = split(";", $omim_description{$omim_id});
    pop(@descriptions);
    next if (not $omim_id_disease{$omim_id});
    print OUTPUT join("\t", ($omim_id, $omim_id_disease{$omim_id}, $_))."\n" for @descriptions;
}






