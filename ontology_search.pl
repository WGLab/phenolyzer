use Bio::OntologyIO;
use strict;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			 '1.00';
our $LAST_CHANGED_DATE = '$LastChangedDate: 2014-04-022 (22nd, April, 2014) $';
our ($help, $man, $is_phenotype, $ontology_file, $format, $if_exact_match);
GetOptions('help|h'=>\$help, 'man|m' => \$man, 'phenotype|p' => \$is_phenotype,
              'ontology_file|o=s'=>\$ontology_file, 'format|f=s' => \$format,
              'exact'=>\$if_exact_match);
   $help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
   $man and pod2usage  (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
   @ARGV or pod2usage      ("       ERROR: Please enter disease names!!");
   @ARGV == 1 or pod2usage ("       ERROR: too many input arguments");
   $ontology_file||="./lib/hpo.obo"  if ($is_phenotype);             #Phenotype input  
   $ontology_file||="./lib/doid.obo" if (not $is_phenotype);         #Disease input
   $format ||= "name";
   my $parser = Bio::OntologyIO->new ( -format  =>  "obo",
                                       -file    =>  $ontology_file   );
   #Enter the input term                                    
   my $input_term = $ARGV[0];
   my $query_term = $input_term;
      $query_term =~ s/\bs\b//g;
      $query_term =~ s/\W+/ /g;  
   my $ont = $parser->next_ontology();
   my $is_a = Bio::Ontology::RelationshipType->get_instance("IS_A");
   my @roots = $ont->get_root_terms();
   my @all_terms = $ont->get_descendant_terms($roots[0], $is_a);
   unshift @all_terms, @roots;
   my %found_terms = (); 
   for my $term ( @all_terms )
   {
   	   my $id = $term -> identifier();
   	   next if $found_terms{$id};
   	   my $name = $term-> name();
   	   my @synonyms = $term -> get_synonyms();
       my @search_list = ($name, @synonyms);	   
   	   for my $each_name (@search_list)
   	   {
   	   	  my $each_name_change = $each_name;
   	   	     $each_name_change =~ s/\bs\b//g;
   	   	     $each_name_change =~ s/\W+/ /g; 
   	   	  if( $each_name_change =~ /\b$query_term\b/i )
   	   	  {
   	   	  	next if($if_exact_match and $each_name !~ /^input_term$/i);
   	   	  	$found_terms{$id} = $term;
   	   	    last;
   	   	  }
   	   }
       if ($found_terms{$id})                 #If the term is matched in this turn, then find its descendants
   	   {
   	   	 $format=~  /\bid\b/ and print $id."\n";
   	   	 
   	   	 if($format=~/\bname\b/)
   	   	 { 
   	   	 	for (@search_list)
   	   	    {
   	   	    s/^(.*?)\s[A-Z]{4,}.*$/\1/;
   	   	    print $_."\n"  ;
   	   	    }
   	   	 }
   	   	 
   	     my @descendant_terms = $ont->get_descendant_terms($found_terms{$id}, $is_a); 
   	   	 for my $descendant (@descendant_terms)
   	   	 {
   	   	 	my $descendant_id = $descendant -> identifier();
   	   	 	next if ($found_terms{$descendant_id});
   	   	 	$name = $descendant -> name();
   	   	 	@synonyms = $descendant -> get_synonyms();
   	   	 	@search_list = ($name, @synonyms);
   	   	 	$found_terms{$descendant_id} = $descendant;
   	   	 	$format=~  /\bid\b/ and print $descendant_id."\n";
   	   	 	if($format=~/\bname\b/)
   	   	    { 
   	   	 	   for (@search_list)
   	   	       {
   	   	       s/^(.*?)\s[A-Z]{4,}.*$/\1/;
   	   	       print $_."\n"  ;
   	   	       }
   	   	    }
         }
   	   }
   	   
   }
 

=head1 SYNOPSIS

  ontology_search.pl [arguments] <disease_term or phenotype_term>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -p, --phenotype                 the term is a phenotype term
        -o, --ontology_file <string>    the ontology file (default:disease_ontology,
                                        ./lib/hpo.obo for phenotype and ./lib/doid.obo for disease)
        -f, --format <string>           the output format, 'name' or 'id' (default:name)                     
 Function:       
          automatically match the term into a disease or phenotype term in the HPO or DO Database, then found 
          their descendants and then print all the names and synonyms, and ids. 
 Example:
         perl ontology_search.pl alzheimer -format name,id
  