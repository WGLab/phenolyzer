use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
open(CLINVAR,"../clinVar_genemap_20150515.txt") or die;
open(OMIM1,"../omim_genemap_20151119.txt") or die;
open(OMIM2,"../omim_genemap2_20151119.txt") or die;
open(ORPHANET,"../orphanet_genemap.txt") or die;
open(GWAS, "../gwas_catalog_dev") or die;                         
open(CLINVAR_OUT,">DB_CLINVAR_GENE_DISEASE") or die;
open(OMIM_OUT,">DB_OMIM_GENE_DISEASE_TEMP") or die;
open(OMIM_OUT_FINAL,">DB_OMIM_GENE_DISEASE") or die;
open(ORPHANET_OUT,">DB_ORPHANET_GENE_DISEASE_TEMP") or die;
open(ORPHANET_OUT_FINAL,">DB_ORPHANET_GENE_DISEASE") or die;
open(GENE_ID, "DB_HUMAN_GENE_ID") or die;


my $i=0;
my %gene_hash;
for my $line (<GENE_ID>)
{
	if($i==0) { $i++;next;  }
	chomp($line);
	my ($id, $gene, $synonyms) = split("\t", $line);
	if($synonyms eq "-") {$synonyms = "";
		$gene_hash{$gene} = $gene;
	}
	else {
		$synonyms =~s/\|/,/g;
		$gene_hash{$gene} = "$gene,$synonyms";
	}
}
seek GENE_ID, 0,0;
for my $line (<GENE_ID>)
{
	if($i==0) { $i++;next;  }
	chomp($line);
	my ($id, $gene, $synonyms) = split("\t", $line);
	next if(not $synonyms);
	if($synonyms eq "-") {
		next;
	}
	else {
		my @synonyms = split('\|', $synonyms);
		$synonyms = join(",", @synonyms); 
		$gene_hash{$_} = "$gene,$synonyms" for @synonyms;
	}
}


$i=0;
# PUBMEDID	Disease/Trait	Reported Gene(s)	p-Value
open(GWAS_OUT, ">DB_GWAS_GENE_DISEASE");
my %gene_disease_info;
print GWAS_OUT join("\t", qw/GENE DISEASE PUBMED_NUMBER RAW_SCORE/)."\n";
for my $line(<GWAS>){        #process gwas_catalog
       if($i==0) { $i++;next; }
       chomp($line);
       next if($line=~/^\W*$/);
       my @words = split("\t", $line);
       my @gene =  split(",", $words[2]); 
       next if(not looks_like_number($words[3]));
       next if ( $gene[0] eq "NR" or $words[3] >= 1);
       for my $each_gene(@gene){
      	             $each_gene=~s/^\s*(.*?)\s*$/$1/;
      	             $words[1] = GetRidOfAnnotations($words[1]);
                     my $pair = $each_gene."\t".$words[1];
                     if($gene_disease_info{$pair})
                     {
                     	$gene_disease_info{$pair}[0] .= " $words[0]"; 
                     	$gene_disease_info{$pair}[1] = 1-$words[3]
                     	 if( (1-$words[3]) > $gene_disease_info{$pair}[1] );
                     }
                     else{
                     	$gene_disease_info{$pair}[0] = $words[0]; 
                     	$gene_disease_info{$pair}[1] = 1-$words[3];
                     }
                 }
}
print GWAS_OUT $_."\t".$gene_disease_info{$_}[0]."\t".$gene_disease_info{$_}[1]."\n"
for (keys %gene_disease_info);
my %clinvar_count =();
#GeneID	GeneSymbol	ConceptID	DiseaseName	SourceName	SourceID	DiseaseMIM	LastUpdated
print CLINVAR_OUT join("\t",qw/GENE DISEASE MIM_NUMBER SOURCE_COUNT/)."\n";
for my $line(<CLINVAR>){    #process clinVar_genemap.txt
	next if($line=~/^#/);
	chomp($line);
	 my @words=split("\t",$line);
     my($gene, $disease, $mim_num) = ($words[1], $words[3], $words[6]); 
     $disease = GetRidOfAnnotations($disease); 
	 $disease=~s/^\W*(.*?)\W*$/$1/;
	 my $gene_disease_mim = join("\t", ($gene, $disease, $mim_num));
	 $clinvar_count{$gene_disease_mim} += 1;
	}
 print CLINVAR_OUT $_."\t".$clinvar_count{$_}."\n" 
 for keys %clinvar_count;
	
print OMIM_OUT_FINAL join("\t",qw/GENE DISEASE MIM_NUMBER SOURCE_CODE LINGKAGE_INFO/)."\n";
for my $line(<OMIM1>)      #process omim_genemap.txt
{
	chomp($line);
	my @words=split('\|',$line);
	if ($words[13]=~/^\s*$/) {next;}
	my $disease= "$words[13] $words[14] $words[15]";
	my @disease_info=split(";",$disease);
	for (@disease_info){
		if(/(.*?)(?:,\s(\d{3,}))?\s\( (\d) \)/x) {
			my ($disease_name,$mim,$linkage)=($1,$2,$3);
			next if($disease_name =~ /^\?/);
			$disease_name=~s/[\[\]\{\}\?]//g;
			$disease_name=GetRidOfAnnotations($disease_name);
			if(defined $mim) {print OMIM_OUT join ("\t",($words[5],$disease_name,$mim,$words[6],$linkage))."\n";} 
		    else {print OMIM_OUT join ("\t",($words[5],$disease_name,$words[9],$words[6],$linkage))."\n";} 		
		}
	}
}

for my $line(<OMIM2>)      #process omim_genemap.txt
{
	
	chomp($line);
	if ($line=~/^\s*$/) {next;}
	my @words=split('\|',$line);
	if(@words< 12) {next;}
	if ($words[11]=~/^\s*$/) {next;}
	my @disease_info=split(";",$words[11]);
	for (@disease_info){
	
		if(/(.*?)(?:,\s(\d{3,}))?\s\( (\d) \)/x) {
			my ($disease_name,$mim,$linkage)=($1,$2,$3);
			next if($disease_name =~ /^\?/);
			$disease_name=~s/[\[\]\{\}\?]//g;
			$disease_name= GetRidOfAnnotations($disease_name);
			if(defined $mim) {print OMIM_OUT join ("\t",($words[5],$disease_name,$mim,$words[6],$linkage))."\n";} 
		    else {print OMIM_OUT join ("\t",($words[5],$disease_name,$words[8],$words[6],$linkage))."\n";}		
		}
	}
}
system("sort DB_OMIM_GENE_DISEASE_TEMP|uniq >> DB_OMIM_GENE_DISEASE");

print ORPHANET_OUT_FINAL join("\t",qw/GENE DISEASE ORPHANET_NUMBER SOURCE_COUNT LINKAGE_INFO/)."\n";
for my $line(<ORPHANET>){
	chomp($line);
	my @words=split("\t",$line);
	$words[6]= GetRidOfAnnotations($words[6]);
	$words[6]=~s/^\W*(.*?)\W*$/$1/;
	
	if($words[17]=~/^\d+$/){
		print ORPHANET_OUT join("\t",($words[13],$words[6],$words[5],$words[17],$words[22]))."\n";
		}
}
sub GetRidOfAnnotations{
	@_==1 or die;
	my $disease = $_[0];
	$disease =~s/\(.*?\)//g;
	return $disease;
}

system("sort DB_ORPHANET_GENE_DISEASE_TEMP|uniq>>DB_ORPHANET_GENE_DISEASE");
system("rm DB_ORPHANET_GENE_DISEASE_TEMP DB_OMIM_GENE_DISEASE_TEMP");



