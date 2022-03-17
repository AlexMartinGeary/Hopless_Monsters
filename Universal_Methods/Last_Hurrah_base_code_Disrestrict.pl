use strict;
use warnings;
use Cwd qw(getcwd);

my $row;
my %gnames;
my %gnames2;
my %age;
my %ohnolog;
my %paralog;
my %haps;
my %gp;
my %ppi;
my %tppi;
my %omim;
my %fam;
my @nil;


open my $fh0, "<", "/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/Ortho_Ages/WG_Ortho_ages.csv" or die;
while ($row = <$fh0>){
	foreach ($row){
		chomp ($row);
		if ($row =~ /^Gene_id/){
			next;
			}
		elsif ($row =~ /^([^,]*),([^,]*),([^.]*)/){
			$gnames{$2}{$1} = ""; #Name, ENS id
			$gnames2{$1}{$2} = ""; #ENS id, Name
			$age{$1}{$3} = ""; #ENS id, age
			}
		}
	}
close $fh0 or die;

open my $fh1, "<", "/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/Ohnologs/hsapiens.Pairs.Strict.2R.txt" or die;
while ($row = <$fh1>){
	foreach ($row){
		chomp ($row);
		if ($row =~ /^Gene_id/){
			next;
			}
		elsif ($row =~ /^([^\t]*)\t([^\t]*)\t/){
			$ohnolog{$2}{$1} = ""; #ENS partner 1, ENS partner 2
			$ohnolog{$1}{$2} = ""; #Ensure all are caught
			}
		}
	}
close $fh1 or die;

open my $fh1a, "<", "/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/Ohnologs/hsapiens.Pairs.Intermediate.2R.txt" or die;
while ($row = <$fh1a>){
	foreach ($row){
		chomp ($row);
		if ($row =~ /^Gene_id/){
			next;
			}
		elsif ($row =~ /^([^\t]*)\t([^\t]*)\t/){
			$ohnolog{$2}{$1} = ""; #ENS partner 1, ENS partner 2
			$ohnolog{$1}{$2} = ""; #Ensure all are caught
			}
		}
	}
close $fh1a or die;

my @files3 = glob("/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/Ensembl_Paralogs/*.txt");
foreach my $file (@files3) {
	open my $fh2, "<", $file or die; 
	while ($row = <$fh2>) { 
   		foreach ($row) {
       		chomp($row);
       		if ($row =~ /^Gene/){
       			next;
       			}
       		elsif($row =~ /^([^\t]*)\t([^\t]*)\t/){
				$paralog{$1} = ""; #ENS paralog gene 
				$paralog{$2} = ""; #ENS paralog gene 
				}
			}
		}
	}

open my $fh3, "<", "/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/Haploinsufficiency/HI_Predictions_Version3.bed" or die;
while ($row = <$fh3>){
	foreach ($row){
		chomp ($row);
		if ($row =~ /^track/){
			next;
			}
		elsif ($row =~ /^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t/){
			if ($4 =~ /^([A-Z0-9]*)\|([0-9.]*)\|([0-9.%]*)$/){
				$haps{$1}{$3} = ""; #gene name, Hap rank
				}
			}
		}
	}
close $fh3 or die;

open my $fh4, "<", "/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/PPI/Gene_Protein.txt" or die;
while ($row = <$fh4>){
	foreach ($row){
		chomp ($row);
		if ($row =~ /^Gene stable/){
			next;
			}
		elsif ($row =~ /^([^\t]*)\t([^.]*)/){
			$gp{$2}{$1} = ""; #ENS prot name, ENS gene name
			}
		}
	}
close $fh4 or die;

my $phencount = 0;	
my $phencount2 = 0;	
open my $in6, "<", "/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/Disease/Genemap2.txt" or die; 
while ($row = <$in6>) {
    foreach ($row) {
        chomp($row);
        if ($row =~ /^\#/){
        	next;
        	} 
        my $type = "Unknown";
        my $name = ""; 
        my $phenotype = ""; 
        if ($row =~ /ENSG/){ 
        	if ($row =~ /^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/){
        		$name = $11;
        		$phenotype = $13;
        		if ($phenotype =~ /^$/){
					$phenotype = "NONE";
					}
        		}
        	}
        if ($row !~ /ENSG/){		  
			if ($row =~ /^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/){
				my $hgnc = $7;
				if ($hgnc =~ /(A-Z)/i){
					$name = $gnames{$hgnc};
					$phenotype = $13;
					if ($phenotype =~ /^$/){
						$phenotype = "NONE";
						}
					}
				}
			}
		my $phen2 = "";
		if ($phenotype !~ /NONE/){
			my @words = split /;/, $phenotype; 
			foreach my $word (@words) {		#If the phenotype does not indicate a "non-disease" 
				if ($word !~ /\[/) { 
					$phen2 = "$phen2".";"."$word";
					}
				}
			if ($phen2 !~ /^$/){ 
				if ($row =~ /recessive/i){
        			$type = "Recessive";
        			}
        		if ($row =~ /dominant/i){
        			$type = "Dominant";
        			}
        		if (($row =~ /recessive/i) && ($row =~ /dominant/i)){
        			$type = "Both";
        			}
       			$omim{$name}{$type}= "";
       			$phencount2++;
       			}
       		else {
       			$phencount++;
       			}
       		}
		}
	}
close $in6 or die;
print "$phencount\t$phencount2\n";

open my $fh7, "<", "/Users/mqbpqam6/Desktop/PhD/Last_Hurrah_Data_Gen/Last_Hurrah_AGES.csv" or die;
while ($row = <$fh7>){
	foreach ($row){
		chomp ($row);
		if ($row =~ /^Gene/){
			next;
			}
		elsif ($row =~ /^([^,]*),([^.]*)/){
			$fam{$1}{$2} = ""; #ENS gene name, family root age, family identifier
			}
		}
	}
close $fh7 or die;


open my $out, ">", "/Users/mqbpqam6/Desktop/Last_Hurrah_Dataset_Disrestrict.csv" or die;
print $out "Gene_id,Gene_Name,Ortho_Age,Paralog_status,Haplosufficiency_Rank,Disease_Association,Family_Root,Family_ID\n";
foreach my $g (sort keys %gnames2){
	print $out "$g,";
	my $hgnc;
	foreach my $n (sort keys %{$gnames2{$g}}){
		print $out "$n,"; 
		$hgnc = $n;
		}
	foreach my $a (sort keys %{$age{$g}}){
		print $out "$a,";
		}
	if (exists($ohnolog{$g})){
		print $out "Ohnolog,";
		}
	elsif (exists($paralog{$g})){
		print $out "SSD,";
		}
	else {
		print $out "Singleton,";
		}
	if (exists($haps{$hgnc})){
		foreach my $h (sort keys %{$haps{$hgnc}}){
			$h =~ s/\%//g;
			print $out "$h,";
			}
		}
	else {
		print $out "NA,";
		}
	if (exists($omim{$g})){	
		my $dc = 0;
		foreach my $d (sort keys %{$omim{$g}}){
			if ($d =~ /Dominant/){
				$dc = 1;
				}
			if ($d =~ /Recessive/){
				$dc = 2;
				}
			if ($d =~ /Both/){
				$dc = 3;
				}
			if ($d =~ /Unknown/){
				$dc = 4;
				}
			}
		if ($dc == 1){
			print $out "Dominant,";
			}
		elsif ($dc == 2){
			print $out "Recessive,";
			}
		elsif ($dc == 3){
			print $out "Both,";
			}
		elsif ($dc == 4){
			print $out "Unknown,";
			}
		}
	else {
		print $out "None,";
		}
	if (exists($fam{$g})){
		foreach my $d (sort keys %{$fam{$g}}){
			print $out "$d"
			}
		}
	else {
		print $out "Singleton,0";
		}	
	print $out "\n";
	}

close $out or die;

exit;	
		
	
	
