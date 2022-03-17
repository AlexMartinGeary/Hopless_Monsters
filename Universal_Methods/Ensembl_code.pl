use strict;
use warnings; 
use Bio::EnsEMBL::Registry;
use 5.010;
my $row; 
my $host   = 'ensembldb.ensembl.org';
my $user   = 'anonymous';
my $port   = 5306;
my $dbname = 'ensembl_compara_97';

Bio::EnsEMBL::Registry->load_registry_from_db(
  -host    => $host,
  -user    => $user,
  -dbname  => $dbname,
  -verbose => '1');

open my $fh, "<", "TGFb_sups.csv" or die;
my @tgfb_sups;
while ($row = <$fh>){
	foreach ($row){
		chomp ($row);
		push @tgfb_sups, $row;
		}
	}

#my $reg = "Bio::EnsEMBL::Registry";
#$reg->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org');

open my $out, ">", "Ortho_out.csv" or die;
print $out "Gene_id,Taxon_between,Type,Taxonomic_level,is_high_confidence,is_tree_compliant,GOC_score,WGA_coverage\n";
#foreach gene in list
foreach my $g1 (@tgfb_sups){ # first you have to get a GeneMember object. In case of homology is a gene, in
	print "$g1\n"; 
	my $gene_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'GeneMember');
	my $gene_member = $gene_member_adaptor->fetch_by_stable_id($g1); # then you get the homologies where the member is involved
	my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');
	my $homologies = $homology_adaptor->fetch_all_by_Member($gene_member); # That will return a reference to an array with all homologies (orthologues in other species and paralogues in the same one)
	foreach my $homology (@{$homologies}) { # Then for each homology, you can get all the Members implicated 
		my $type = $homology->description;
		my $goc = $homology->goc_score;
		unless (length $goc // ''){ 
			$goc = 0;
			}
		my $wga = $homology->wga_coverage;
		unless (length $wga // ''){ 
			$wga = 0;
			}
		my $conf = $homology->is_high_confidence;
		unless (length $conf // ''){ 
			$conf = 0;
			}
		my $tree = $homology->is_tree_compliant;
		unless (length $tree // ''){ 
			$tree = 0;
			}
		if (($conf == 1) and ($tree == 1)){
			if (($type =~ /ortholog_one2one/) or ($type =~ /ortholog_one2many/)){
  				print $out $gene_member->stable_id,",", $homology->method_link_species_set->name ,",", $homology->description,",", $homology->taxonomy_level ,",", $homology->is_high_confidence ,",", $homology->is_tree_compliant ,",", $goc ,",", $wga ,"\n";
				}
			}
		}
	}



