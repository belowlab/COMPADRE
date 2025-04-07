package PADRE;

## copied from ersa_plus_PRIMUS.pl on Sept 16 2014
## converted into a .pm file

use strict;
use node_v7;
use Statistics::Distributions;

my $min_degree_to_consider = 4; ## THIS NEEDS TO ADJUST DEPENDING ON THE THRESHOLD USED IN RECONSTRUCTION
my $max_degree_for_related = 9; ## 10th degree relatives are more likely to not share any segments IBD, then they are to share a segment.
my $cap_degree = 9; ## All relationships above this degree, will be assigned the likelihood based on this degree: see below.
my $allow_direct_descendant_comparisons = 0;
my $AIC_correction = 2;
my $multuple_test_correct = 0;
my $p_val_cutoff = 0.05;
my $verbose = 1;
my $USE_PEDIGREE_LIKELIHOODS = 1;
my $line_position_for_first_pedigree = 9; ## the line position in the network summary file for the first pedigree

my $num_founder_founder_comparisons = 0;
my $num_unrelated_founder_comparisons = 0;

my $GLOBAL_COUNTER = 0;


sub run_PADRE_project_summary
{
	my $summary_file = shift;
	my $ersa_likelihoods_file = shift;
	my $ersa_results_file = shift;
	my $max_degree_relationship_used_to_reconstruct_pedigrees = shift;
	my $output_dir = shift;
	my $multiple_test_correct = shift;

	#$AIC_correction = $user_defined_AIC_correction if $user_defined_AIC_correction ne "";

	my %relatedness_network;

	system("mkdir $output_dir") if $output_dir ne "" && !-d $output_dir;

	## Update the min degree to consider based on what degree was used to reconstruct
	$min_degree_to_consider = $max_degree_relationship_used_to_reconstruct_pedigrees+1 if $max_degree_relationship_used_to_reconstruct_pedigrees ne "";
	


	$summary_file =~ /(.*)\/Summary_(.+).txt/;
	my $project_dir = $1;
	my $project_name = $2;
	#print "project_dir: $project_dir\n";
	my $results_file = "$project_dir/PADRE_network_connections.txt";
	$results_file = "$output_dir/PADRE_network_connections.txt" if $output_dir ne "";
	my $all_degrees_file = "$project_dir/PADRE_sample_relatedness.txt";
	$all_degrees_file = "$output_dir/PADRE_sample_relatedness.txt" if $output_dir ne "";
	my $dot_file = "$project_dir/PADRE_network_connections.dot";
	$dot_file = "$output_dir/PADRE_network_connections.dot" if $output_dir ne "";
	
	my $unrelated_samples_file = "$project_dir/$project_name\_unrelated_samples.txt";

	print "\nRUNNING PADRE\n";
	print "Project summary file: $summary_file\n";
	print "ERSA relationship likelihoods: $ersa_likelihoods_file\n";
	print "ERSA max likelihood relationship results: $ersa_results_file\n";
	print "Max relationship degree use in reconstruction: $max_degree_relationship_used_to_reconstruct_pedigrees\n";
	print "min degree to consider: $min_degree_to_consider\n";
	print "max degree for related: $max_degree_for_related\n";
	print "Results: $results_file\n";
	print "Project name: $project_name\n\n";

	## Adjust p_val_cutoff for multiple testing if user specifies to do so
	$p_val_cutoff = get_multiple_test_corrected_p_val($summary_file,$unrelated_samples_file,$project_dir) if $multiple_test_correct ne 0;
	
	## Load likelihoods
	my $likelihoods = load_ersa_likelihoods($ersa_likelihoods_file,$ersa_results_file);

	## Load network data from summary file
	my %network_data;
	open(IN,$summary_file);
	my ($junk,$family) = split(/\s+/,<IN>);
	my ($junk,$num_networks) = split(/\s+/,<IN>);
	<IN>;
	my $header = split(/\s+/,<IN>);
	while(my $line = <IN>)
	{
		chomp($line);
		my ($net_num,$net_dir,$num_possible,$num_unflagged,$num_samples,$num_at_max,$sample_IDs,$error_message) = split(/\s+/,$line);
		$network_data{$net_num}{'dir'} = $net_dir;
		my @ids = (split(/,/,$sample_IDs));
		for(my $i = 0; $i < @ids; $i++){$ids[$i] =~ s/.+__//;}
		my %IDs = map { $_ => 1 } @ids;
		$network_data{$net_num}{'samples'} = \%IDs;
		#foreach my $id (keys %{$network_data{$net_num}{'samples'} })
		#{
			#print "id $id\n";
		#}
	}
	close(IN);

	open(PROJECT_LEVEL,">$results_file");
	open(ALL,">$all_degrees_file");
	print PROJECT_LEVEL "Network\tPedigree\tFounder\tNetwork\tPedigree\tFounder\tdegree_related\tp_value\n";
	print ALL "Network\tPedigree\tID\tNetwork\tPedigree\tID\tdegree_related\tp_value\n";
	
	## Write out within network degree relationships based on the highest ranked pedigree
	for(my $i = 1; $i <= $num_networks; $i++)
	{
		my $network1_summary_file = "$project_dir/$network_data{$i}{'dir'}/Summary_$network_data{$i}{'dir'}.txt";

		my %degrees;
		get_degrees_within_network($network1_summary_file,\%degrees);
		
		foreach my $id1 (sort keys %degrees)
		{
			next if $id1 =~ /Missing/i;
			foreach my $id2 (sort keys %{$degrees{$id1} })
			{
				last if $id1 eq $id2;
				next if $id2 =~ /Missing/i;
				$degrees{$id1}{$id2} = "UN" if !exists $degrees{$id1}{$id2};
				print ALL "$network_data{$i}{'dir'}\t1\t$id1\t$network_data{$i}{'dir'}\t1\t$id2\t$degrees{$id1}{$id2}\tPRIMUS\n";
			}
		}
	}

	## Compare each pair of networks
	for(my $i = 1; $i <= $num_networks-1; $i++)
	{
		my $network1_summary_file = "$project_dir/$network_data{$i}{'dir'}/Summary_$network_data{$i}{'dir'}.txt";

		#print "net1_sum: $network1_summary_file\n";
		for(my $j = $i + 1; $j <= $num_networks; $j++)
		{
			my $network2_summary_file = "$project_dir/$network_data{$j}{'dir'}/Summary_$network_data{$j}{'dir'}.txt";
			
			## load likelihoods for just these two networks
			#my $likelihoods = load_ersa_likelihoods($ersa_likelihoods_file,$ersa_results_file,$network_data{$j}{'samples'},$network_data{$j}{'samples'});
			
			## Get min AICcombined for these two networks
			my ($min_AICcombined,$min_pair_ref,$min_AICcombined_degrees_ref,$min_ped1,$min_ped2,$fam_file1,$fam_file2) = get_min_AICcombined_score_two_networks($network1_summary_file,$network2_summary_file,$likelihoods);

			#print "$min_AICcombined,$min_pair_ref,$min_AICcombined_degrees_ref,$min_ped1,$min_ped2,$fam_file1,$fam_file2\n";

			## Check that the connection is supported by significant ERSA results
			my $is_significant = is_network_to_network_connection_significant($min_pair_ref,$likelihoods,$fam_file1,$fam_file2);
			if ($is_significant eq "nonsig")
			{
				$$min_pair_ref[0] = "all";
				$$min_pair_ref[1] = "all";
				$$min_pair_ref[2] = "UN";
			}
      if($is_significant == 0)
      {
        $is_significant = "<10^-80";
      }

			## Write results to network file
			print PROJECT_LEVEL "$network_data{$i}{'dir'}\t$min_ped1\t$$min_pair_ref[0]\t$network_data{$j}{'dir'}\t$min_ped2\t$$min_pair_ref[1]\t$$min_pair_ref[2]\t$is_significant\n";
			
			## Store results in relatedness_network to make the .dot file
			$relatedness_network{"Network$i"}{"Network$j"}=$$min_pair_ref[2];

			## Write results to sample file
			foreach my $id1 (keys %{$min_AICcombined_degrees_ref})
			{
				next if $id1 =~ /Missing/i;
				next if !exists $network_data{$i}{'samples'}{$id1};
				foreach my $id2 (keys %{$$min_AICcombined_degrees_ref{$id1} })
				{
					next if !exists $network_data{$j}{'samples'}{$id2};
					next if $id2 =~ /Missing/i;
					print ALL "$network_data{$i}{'dir'}\t$min_ped1\t$id1\t$network_data{$j}{'dir'}\t$min_ped2\t$id2\t$$min_AICcombined_degrees_ref{$id1}{$id2}\t$is_significant\n";
				}

			}
		}
	}

	## Compare each unrelated sample to each network
	open(UNRELATEDS,"$unrelated_samples_file") or die "can't open $unrelated_samples_file; $!\n";
	my $header = <UNRELATEDS>;
	my @lines = <UNRELATEDS>;
	my @unrelateds;
	close(UNRELATEDS);
	foreach my $line (@lines)
	{
		chomp($line);
		my ($fid, $iid) = split(/\s+/,$line);
		push(@unrelateds,$iid);
		for(my $i = 1; $i <= $num_networks-1; $i++)
		{
			my $network1_summary_file = "$project_dir/$network_data{$i}{'dir'}/Summary_$network_data{$i}{'dir'}.txt";
			
			## Get min AICcombined for these two networks
			my ($min_AICcombined,$min_pair_ref,$min_AICcombined_degrees_ref,$min_ped1,$min_ped2,$fam_file1,$fam_file2) = get_min_AICcombined_score_unrelated_sample_and_a_network($iid,$network1_summary_file,$likelihoods,$project_dir);	
			
			## Check that the connection is supported by significant ERSA results
			my $is_significant = is_network_to_network_connection_significant($min_pair_ref,$likelihoods,$fam_file1,$fam_file2);
			if ($is_significant eq "nonsig")
			{
				$$min_pair_ref[0] = "all";
				$$min_pair_ref[1] = "all";
				$$min_pair_ref[2] = "UN";
			}
      if($is_significant == 0)
      {
        $is_significant = "<10^-80";
      }

			## Write results to network file
			print PROJECT_LEVEL "$network_data{$i}{'dir'}\t$min_ped1\t$$min_pair_ref[0]\tNA\t$iid\t$$min_pair_ref[1]\t$$min_pair_ref[2]\t$is_significant\n";
			## Store results in relatedness_network to make the .dot file
			$relatedness_network{"Network$i"}{$iid}=$$min_pair_ref[2];

			## Write results to sample file
			foreach my $id1 (keys %{$min_AICcombined_degrees_ref})
			{
				next if $id1 =~ /Missing/i;
				next if !exists $network_data{$i}{'samples'}{$id1};
				print ALL "$network_data{$i}{'dir'}\t$min_ped1\t$id1\tunrelated_samples\tNA\t$iid\t$$min_AICcombined_degrees_ref{$id1}{$iid}\t$is_significant\n";
			}
		}
	}
	close(PROJECT_LEVEL);

	## Compare unrelateds to each other
	foreach my $id1 (@unrelateds)
	{
		next if $id1 =~ /Missing/i;
		foreach my $id2 (@unrelateds)
		{
			last if $id1 eq $id2;
			next if $id2 =~ /Missing/i;
			next if $$likelihoods{$id1}{$id2}{'best'} eq "UN";
			
			my $is_significant = is_sample_to_sample_connection_significant($id1,$id2,$$likelihoods{$id1}{$id2}{'best'},$likelihoods);
			if ($is_significant != -1)
			{
				print ALL "unrelated_samples\tNA\t$id1\tunrelated_samples\tNA\t$id2\t$$likelihoods{$id1}{$id2}{'best'}\t$is_significant\n";
			}
			
			## Store results in relatedness_network to make the .dot file
			#$relatedness_network{$id1}{$id2}=$$likelihoods{$id1}{$id2}{'best'};
		}
	}
	close(ALL);

	## Prune out the network pairs below significance cutoff
	

	## Make the .dot file for the networks and unrelated samples
	write_out_dot_file($dot_file,\%relatedness_network);

  #print "GLOBAL COUNTER: $GLOBAL_COUNTER\n";

	return $results_file;
}

##############################################################################################
##############################################################################################
##############################################################################################
sub is_network_to_network_connection_significant
{
	my ($min_pair_ref,$likelihoods,$fam_file1,$fam_file2) = @_;
	my ($id1,$id2,$degree) = @$min_pair_ref;

	return 0 if $degree eq "UN";

	my $adjusted_degree = $degree;

	my @network1_individuals_to_compare;
	my @network2_individuals_to_compare;

	## Need to compare non-missing to non-missing; so if founder is Missing, then use the non-missing children of the founder
	## If the founder doesn't have any non-missing children, then keep going deeper until you find one.
	if($id1 =~ /Missing/)
	{
        #print "fam1: $fam_file1\n";
		my $network1 = build_network_from_ped_file($fam_file1);
		my ($arr_ref,$depth) = get_non_missing_decendants($id1,$network1,0);
		@network1_individuals_to_compare = @{$arr_ref};
		$adjusted_degree += $depth;
        #print "$id1\'s non_missing decendants: @network1_individuals_to_compare ($depth)\n";
	}
	else
	{
		push(@network1_individuals_to_compare,$id1);
	}
	if($id2 =~ /Missing/)
	{
        #print "fam2: $fam_file2\n";
		my $network2 = build_network_from_ped_file($fam_file2);
		my ($arr_ref,$depth) = get_non_missing_decendants($id2,$network2,0);
		@network2_individuals_to_compare = @{$arr_ref};
		$adjusted_degree += $depth;
        #print "$id2\'s non_missing decendants: @network2_individuals_to_compare ($depth)\n";
	}
	else
	{
		push(@network2_individuals_to_compare,$id2);
	}
	
    #print "net1 compare: @network1_individuals_to_compare\n";
    #print "net2 compare: @network2_individuals_to_compare\n";

	## Sum the chi-squared statistics and count the number of comparisons to see the degree of freedom correctly
	my $chi_square_sum = 0;
	my $num_comparisons = 0;
	foreach my $id1 (@network1_individuals_to_compare)
	{
		foreach my $id2 (@network2_individuals_to_compare)
		{
			my $lnl_related = $$likelihoods{$id1}{$id2}{$adjusted_degree} + $AIC_correction;
			my $lnl_unrelated = $$likelihoods{$id1}{$id2}{"UN"};
			my $chi_square_val = -2*($lnl_unrelated-$lnl_related);
			$chi_square_sum += $chi_square_val;
			$num_comparisons++;
			#my $is_significant = is_sample_to_sample_connection_significant($id1,$id2,$adjusted_degree,$likelihoods);
			
			#print "$id1 <-> $id2 = $is_significant\n";
			#if($is_significant eq 1)
			#{
			#	return 1;
			#}
		}
	}
	die "no valid comparisons ($fam_file1,$fam_file2)\n" if $num_comparisons == 0;

	## Determine significance
	my $degrees_of_freedom = 2 * $num_comparisons;
	my $chi_square_prob = Statistics::Distributions::chisqrprob($degrees_of_freedom,$chi_square_sum);
    #print "chi_square_sum: $chi_square_sum\n";
    #print "degrees of freedom: $degrees_of_freedom\n";
    #print "chi_square_prob: $chi_square_prob\n";
	if($chi_square_prob < $p_val_cutoff)
	{
		return $chi_square_prob;
	}
	else
	{
    #print "HHHHHHHHHHHHHEEEEEEEEEEEEEERRRRRRRRRRRRRRRREEEEEEEEEEEEEEEEEEEEEEEE\n";
		return "nonsig";
	}
}


sub get_non_missing_decendants ## bredth first search, better
{
	my ($parent,$network_ref,$depth) = @_;

	#print "\nParent: $parent\n";
	my @non_missing_individuals;
	my %id_depth;
	my %depth_ids;
	push(@{$depth_ids{0} },$parent);
	
	#print "depth 0: @{$depth_ids{0} }\n";

	my $search_depth = 0;
	my $continue = 1;
	while($continue)
	{
		$search_depth++;
		foreach my $parent (@{$depth_ids{($search_depth-1)} })
		{
			my @children = $$network_ref{$parent}->children();
			#print "$parent: @children\n";
			foreach my $child (@children)
			{
				push(@{$depth_ids{$search_depth} },$child);
				next if $child =~ /Missing/;
				push(@non_missing_individuals,$child);
			}
			if(@non_missing_individuals > 0)
			{
				#print "HERE!\n";
				$continue = 0;
			}
			#print "depth $search_depth: @{$depth_ids{$search_depth} }\n";
		}
	}
	$depth = $search_depth;
	return (\@non_missing_individuals,$depth);
}

sub get_non_missing_decendants_OLD ## depth first search... not ideal
{
	my ($parent,$network_ref,$depth) = @_;
    #print "\nParent: $parent\n";
	my @children = $$network_ref{$parent}->children();
	my @non_missing_individuals;

	if(@children < 1)
	{
		die "$parent has no children\n";
	}

	$depth++;

	## Get all non-missing children
	foreach my $child (@children)
	{
        #print "$parent\'s child: $child\n";
		next if $child =~ /Missing/;
		push(@non_missing_individuals,$child);
	}

	## Decend if no non-missing children
	if(@non_missing_individuals < 1)
	{
        #print "HERE!\n";
		foreach my $child (@children)
		{
            #print "$parent\'s child: $child\n";
			my ($non_missing_individuals_ref,$depth_temp) = get_non_missing_decendants($child,$network_ref,$depth);
			$depth = $depth_temp;
			@non_missing_individuals = @{$non_missing_individuals_ref};
		}
	}
    #print "depth: $depth\n";
	return (\@non_missing_individuals,$depth);
}

sub is_sample_to_sample_connection_significant
{
	my ($id1,$id2,$degree,$likelihoods) = @_;
	
	return 0 if $degree eq "UN";
	
	my $lnl_related = $$likelihoods{$id1}{$id2}{$degree} + $AIC_correction;
	my $lnl_unrelated = $$likelihoods{$id1}{$id2}{"UN"};
	my $chi_square_val = -2*($lnl_unrelated-$lnl_related);
	my $chi_square_prob = Statistics::Distributions::chisqrprob(2,$chi_square_val);
	
    #print "$id1:$id2 = $degree\n";
    #print "$lnl_unrelated-$lnl_related\n";
    #print "chi_square_val: $chi_square_val\n";
    #print "chi_square_prob: $chi_square_prob\n";
	if($chi_square_prob < $p_val_cutoff)
	{
		return $chi_square_prob;
	}
	else
	{
        #print "HHHHHHHHHHHHHEEEEEEEEEEEEEERRRRRRRRRRRRRRRREEEEEEEEEEEEEEEEEEEEEEEE\n";
		return -1;
	}
}

sub get_multiple_test_corrected_p_val
{
	my $summary_file = shift;
	my $unrelateds_file = shift;
	my $project_dir = shift;
	
	## Get number unrelated samples
	my $num_unrelated = 0;
	open(MULT,$unrelateds_file);
	<MULT>; #header
	while(<MULT>)
	{
		$num_unrelated++;
	}
	close(MULT);

	## Get num_founders
	my %founder_data;
	open(MULT,$summary_file) or die "Can't open $summary_file: $!\n";
	my ($junk,$family) = split(/\s+/,<MULT>);
	my ($junk,$num_networks) = split(/\s+/,<MULT>);
	<MULT>;
	my $header = split(/\s+/,<MULT>);
	while(my $line = <MULT>)
	{
		#print "summary_file line: $line";
		chomp($line);
		last if $line eq "";
		my ($net_num,$net_dir,$num_possible,$num_unflagged,$num_samples,$num_at_min,$sample_IDs,$error_message) = split(/\s+/,$line);
		my $network_summary_file = "$project_dir/$net_dir/Summary_$net_dir.txt";
		get_founders_within_network($network_summary_file,\%founder_data,$net_num);
	}
	close(MULT);
	
	## Get number of pairwise comparisons between founders
	my $num_networks = keys %founder_data;

	#print "num_networks: $num_networks\n";
	my $num_pairwise_comparisons = 1;

	for(my $net1 = 1; $net1 <= $num_networks-1; $net1++)
	{
		#print "net1 $net1\n";
		my $num_founders1 = keys %{ $founder_data{$net1} };
		for(my $net2 = $net1+1; $net2 <= $num_networks; $net2++)
		{
			#print "net2 $net2\n";
			my $num_founders2 = keys %{ $founder_data{$net2} };
			my $pairs = $num_founders1 * $num_founders2;
			#print " F1*F2:$num_founders1 * $num_founders2\n";
			$num_pairwise_comparisons += $pairs;
		}
	}
	#print "num pairwise comparisons: $num_pairwise_comparisons\n";

	## get number of pairwise comparison between founders and unrelateds
	for(my $net1 = 1; $net1 < $num_networks; $net1++)
	{
		my $num_founders1 = keys %{ $founder_data{$net1} };
		my $pairs = $num_founders1 * $num_unrelated;
		$num_pairwise_comparisons += $pairs;
	}
	#print "num pairwise comparisons: $num_pairwise_comparisons\n";

	## get number of unrelated to unrelateds
	my $pairs = $num_unrelated * $num_unrelated;
	$num_pairwise_comparisons += $pairs;
	
  print "num pairwise comparisons: $num_pairwise_comparisons\n";

	### get corrected p-value
	my $p_val = 0.05/$num_pairwise_comparisons;
    #print "p_val = $p_val\n";
	#my $percentile = 1-$p_val;
	#my $chis=Statistics::Distributions::chisqrdistr (2,$p_val);
	#my $penalty = $chis/2;
	#print "Chi-squared-crit (2 degrees of freedom, $percentile"."th percentile = $p_val level) = $chis\n";	
	#print "AIC penalty value: $penalty\n";
	return $p_val;
}

sub get_founders_within_network
{
	my $network_summary = shift;
	my $founder_ref = shift;
	my $net_num = shift;
	
	#print "net num: $net_num\n";

	$network_summary =~ /(.*)\/Summary_/;
	my $network_dir = $1;

	## LOAD network data
	open(NET,$network_summary) or die "can't open $network_summary: $!\n";
	my @lines = <NET>;
	close(NET);
	#print "@lines1";
	my ($junk1,$junk2, $network_name) = split(/\s+/,@lines[0]);
	my ($junk1,$junk2, $num_possible) = split(/\s+/,@lines[1]);

	my %ped_lnl_data;

	#print "net1_name: $network1_name\n";
	#print "net1_possible: $num_possible1\n";
	for(my $i = $line_position_for_first_pedigree; $i < $line_position_for_first_pedigree + $num_possible; $i++)
	{
		my($ped_num, $num_dummies, $num_generations, $rank, $pedigree_lnl, $age_flags,@rest) = split(/\s+/,@lines[$i]); 
		$ped_lnl_data{$ped_num} = $pedigree_lnl;
	}

	## Compare each pair of pedigrees from the two networks
	foreach my $ped (sort {$a<=>$b} keys %ped_lnl_data)
	{
		my $fam_file = "$network_dir/$network_name\_$ped.fam";
		
		my $network = build_network_from_ped_file($fam_file);
		foreach my $id (sort {$a cmp $b} keys %$network)
		{
			## Check that $id1 is a founder, regardless of whether it is missing
			my $id_is_founder = is_founder($id,$network);

			$$founder_ref{$net_num}{$id}++ if $id_is_founder;
		}
	}
}

sub get_min_AICcombined_score_unrelated_sample_and_a_network
{
	my $unrelated_sample = shift;
	my $network_summary1 = shift;
	my $likelihoods = shift;
	my $project_dir = shift;
	
	$network_summary1 =~ /(.*)\/Summary_/;
	my $network1_dir = $1;
	
	#print "Unrelated sample: $unrelated_sample\n";
	#print "net_dir1: $network_dir1\n";

	## Make the .fam file for the unrelated sample
	my $unrelated_sample_fam_file = "$project_dir/$unrelated_sample.fam";
	open(OUT,">$unrelated_sample_fam_file");
	print OUT "0 $unrelated_sample 0 0";
	close(OUT);

	my %ped_lnl_data1;

	## LOAD network data
	open(NET1,$network_summary1) or die "can't open $network_summary1: $!\n";
	my @lines1 = <NET1>;
	close(NET1);
	#print "@lines1";
	my ($junk1,$junk2, $network1_name) = split(/\s+/,@lines1[0]);
	my ($junk1,$junk2, $num_possible1) = split(/\s+/,@lines1[1]);

	#print "net1_name: $network1_name\n";
	#print "net1_possible: $num_possible1\n";
	for(my $i = $line_position_for_first_pedigree; $i < $line_position_for_first_pedigree + $num_possible1; $i++)
	{
		my($ped_num, $num_dummies, $num_generations, $rank, $pedigree_lnl, $age_flags,@rest) = split(/\s+/,@lines1[$i]); 
		$ped_lnl_data1{$ped_num} = $pedigree_lnl;
	}


	## Compare each pair of pedigrees from the two networks
	my $min_AICcombined = ""; 
	my $min_pair_ref;
	my $min_AICcombined_degrees_ref;
	my $min_ped;
	foreach my $ped1 (sort {$a<=>$b} keys %ped_lnl_data1)
	{
		my $fam_file1 = "$network1_dir/$network1_name\_$ped1.fam";
		my $fam_file2 = "$unrelated_sample_fam_file";
		
		## Get lnl for the ERSA relationhips
		my ($ersa_min_AICcombined,$pair_ref,$AIC_combined_degrees_ref) = get_min_AICcombined_score_two_pedigrees($fam_file1,$fam_file2,$likelihoods);
		
		## Combine that with the pedigree lnl caculated in the PRIMUS reconstructions
		my $possible_min_AICcombined = $ersa_min_AICcombined + $ped_lnl_data1{$ped1};
		#print "AIC_combined: $ersa_min_AICcombined + $ped_lnl_data1{$ped1} = $AIC_combined\n";

		if($min_AICcombined eq "")
		{
			$min_AICcombined = $possible_min_AICcombined;
			$min_pair_ref = $pair_ref;
			$min_AICcombined_degrees_ref = $AIC_combined_degrees_ref;
			$min_ped = $ped1;
		}
		elsif($possible_min_AICcombined < $min_AICcombined)
		{
			$min_AICcombined = $possible_min_AICcombined;
			$min_pair_ref = $pair_ref;
			$min_AICcombined_degrees_ref = $AIC_combined_degrees_ref;
			$min_ped = $ped1;
		}
	}
	#system("rm $unrelated_sample_fam_file");
	my $fam_file1 = "$network1_dir/$network1_name\_$min_ped.fam";
	my $fam_file2 = "$unrelated_sample_fam_file";
	return ($min_AICcombined,$min_pair_ref,$min_AICcombined_degrees_ref,$min_ped,"unrelateds",$fam_file1,$fam_file2);
}


sub get_min_AICcombined_score_two_networks
{
	my $network_summary1 = shift;
	my $network_summary2 = shift;
	my $likelihoods = shift;
	
	$network_summary1 =~ /(.*)\/Summary_/;
	my $network1_dir = $1;
	$network_summary2 =~ /(.*)\/Summary_/;
	my $network2_dir = $1;
	
	#print "net1_dir: $network1_dir\n";
	#print "net2_dir: $network2_dir\n";

	my %ped_lnl_data1;
	my %ped_lnl_data2;
	my $line_position_for_first_pedigree = 9;

	## LOAD network data
	open(NET1,$network_summary1) or die "can't open $network_summary1: $!\n";
	my @lines1 = <NET1>;
	close(NET1);
	#print "@lines1";
	my ($junk1,$junk2, $network1_name) = split(/\s+/,@lines1[0]);
	my ($junk1,$junk2, $num_possible1) = split(/\s+/,@lines1[1]);

	#print "net1_name: $network1_name\n";
	#print "net1_possible: $num_possible1\n";
	for(my $i = $line_position_for_first_pedigree; $i < $line_position_for_first_pedigree + $num_possible1; $i++)
	{
		my($ped_num, $num_dummies, $num_generations, $rank, $pedigree_lnl, $age_flags,@rest) = split(/\s+/,@lines1[$i]); ## UPDATE THIS INCLUDE PEDIGREE LIKELIHOOD ONCE THAT MAKES IT INTO FILE
		$ped_lnl_data1{$ped_num} = $pedigree_lnl; ## UPDATE THIS TO PEDIGREE LIKELIHOOD ONCE THAT MAKES IT INTO FILE
	}

	open(NET2,$network_summary2);
	my @lines2 = <NET2>;
	close(NET2);
	my ($junk1,$junk2, $network2_name) = split(/\s+/,@lines2[0]);
	my ($junk1,$junk2, $num_possible2) = split(/\s+/,@lines2[1]);
	for(my $i = $line_position_for_first_pedigree; $i < $line_position_for_first_pedigree + $num_possible2; $i++)
	{
		my($ped_num, $num_dummies, $num_generations, $rank, $pedigree_lnl, $age_flags,@rest) = split(/\s+/,@lines2[$i]);
		#print "ped $ped_num $rank, $pedigree_lnl\n";
		$ped_lnl_data2{$ped_num} = $pedigree_lnl; 
	}


	## Compare each pair of pedigrees from the two networks
	my $min_AICcombined = ""; 
	my $min_pair_ref;
	my $min_AICcombined_degrees_ref;
	my $min_ped1;
	my $min_ped2;
	foreach my $ped1 (sort {$a<=>$b} keys %ped_lnl_data1)
	{
		my $fam_file1 = "$network1_dir/$network1_name\_$ped1.fam";
		foreach my $ped2 (sort {$a<=>$b} keys %ped_lnl_data2)
		{
			my $fam_file2 = "$network2_dir/$network2_name\_$ped2.fam";
			
			## Get lnl for the ERSA relationhips
			my ($ersa_min_AICcombined,$pair_ref,$AIC_combined_degrees_ref) = get_min_AICcombined_score_two_pedigrees($fam_file1,$fam_file2,$likelihoods);
		
			## Combine that with the pedigree lnl caculated in the PRIMUS reconstructions
			
			
			my $possible_min_AICcombined = $ersa_min_AICcombined;
			$possible_min_AICcombined = $possible_min_AICcombined + $ped_lnl_data1{$ped1} + $ped_lnl_data2{$ped2} if $USE_PEDIGREE_LIKELIHOODS;
			
			#print "AIC_combined: $ersa_min_AICcombined + $ped_lnl_data1{$ped1} + $ped_lnl_data2{$ped2} = $AIC_combined\n";
			#exit;
			if($min_AICcombined eq "")
			{
				$min_AICcombined = $possible_min_AICcombined;
				$min_pair_ref = $pair_ref;
				$min_AICcombined_degrees_ref = $AIC_combined_degrees_ref;
				$min_ped1 = $ped1;
				$min_ped2 = $ped2;
			}
			elsif($possible_min_AICcombined < $min_AICcombined)
			{
				$min_AICcombined = $possible_min_AICcombined;
				$min_pair_ref = $pair_ref;
				$min_AICcombined_degrees_ref = $AIC_combined_degrees_ref;
				$min_ped1 = $ped1;
				$min_ped2 = $ped2;
			}
		}
	}

	my $fam_file1 = "$network1_dir/$network1_name\_$min_ped1.fam";
	my $fam_file2 = "$network2_dir/$network2_name\_$min_ped2.fam";
	return ($min_AICcombined,$min_pair_ref,$min_AICcombined_degrees_ref,$min_ped1,$min_ped2,$fam_file1,$fam_file2);
}


## Cycle though each pair of individuals from the two pedigrees to find the min AIC_combined log likelihood
sub get_min_AICcombined_score_two_pedigrees
{
	my $fam_file1 = shift;
	my $fam_file2 = shift;
	my $likelihoods = shift;
	
	my %AIC_combineds;
	my $min_AICcombined;
	my @min_pair;

	#print "$fam_file1 and $fam_file2\n";
	my $network1 = build_network_from_ped_file($fam_file1);
	my $network2 = build_network_from_ped_file($fam_file2);
	
	## Get the unrelated version first
	$min_AICcombined = get_AIC_combined("all","all","UN",$network1,$network2,$likelihoods);
	@min_pair = ("all","all","UN");
	$AIC_combineds{"all"}{"all"}{"UN"} = $min_AICcombined;
			
	foreach my $id1 (sort {$a cmp $b} keys %$network1)
	{
		## Check that $id1 is a founder, regardless of whether it is missing
		my $id1_is_founder = is_founder($id1,$network1);
		my $id1_is_leaf = is_leaf($id1,$network1);
		if(!$id1_is_founder && !$id1_is_leaf){next}
		foreach my $id2 (sort {$a cmp $b} keys %$network2)
		{

            #print "id1:id2 = $id1:$id2\n";
			my $id2_is_founder = is_founder($id2,$network2);
			my $id2_is_leaf = is_leaf($id2,$network2);
			## Check that $id2 is a founder or leaf, regardless of whether it is missing
			if(!$id2_is_founder && !$id2_is_leaf){next}
			
			## Check that id1 or id2 founders
			if(!$id1_is_founder && !$id2_is_founder){next}

			## If the run is not allowing direction decendant (leaf to founder) comparisons, then not allow it.
			if(!$allow_direct_descendant_comparisons && (!$id1_is_founder || !$id2_is_founder)){next}

			## Adjust the degree to consider based on depth to nearest non-missing decedant; the adjustment will not exceed 6 due to the depth of possible reconstructions
			my $adj_min_degree_to_consider = $min_degree_to_consider;
			my $adj_max_degree_to_consider = $max_degree_for_related;
			if($id1 =~ /Missing/i)
			{
				my ($arr_ref,$depth) = get_non_missing_decendants($id1,$network1,0);
				$adj_min_degree_to_consider = $adj_min_degree_to_consider - $depth;
				$adj_max_degree_to_consider = $adj_max_degree_to_consider - $depth;
			}
			if($id2 =~ /Missing/i)
			{
				my ($arr_ref,$depth) = get_non_missing_decendants($id2,$network2,0);
				$adj_min_degree_to_consider = $adj_min_degree_to_consider - $depth;
				$adj_max_degree_to_consider = $adj_max_degree_to_consider - $depth;
			}
			
			$adj_min_degree_to_consider = 1 if $adj_min_degree_to_consider < 1;
            #print "adj_min_degree_to_consider $adj_min_degree_to_consider\n";
            #print "adj_max_degree_to_consider $adj_max_degree_to_consider\n";
			

			for(my $degree = $adj_min_degree_to_consider; $degree <= $adj_max_degree_to_consider; $degree++)
			{
				my $lnl = $$likelihoods{$id1}{$id2}{$degree};
				#print "before getting AIC_combined score: $id1 <-> $id2 $degree = $lnl\n";
				
				my $AIC_combined = get_AIC_combined($id1,$id2,$degree,$network1,$network2,$likelihoods);
                #die "ERROR: Composite likelihood is 0; $id1,$id2,$degree\n" if $AIC_combined == 0;
                if($AIC_combined == 0)
                {
                    #print "WARNING: Composite likelihood is 0; $id1,$id2,$degree\n";
                    $AIC_combined = 0.00000001;
                }
				
                #print "Composite $id1,$id2,$degree: $AIC_combined\n" if $id1 eq "A2" and $id2 eq "B2";
				$AIC_combineds{$id1}{$id2}{$degree} = $AIC_combined;
				if($AIC_combined < $min_AICcombined)
				{
					$min_AICcombined = $AIC_combined;
					@min_pair = ($id1,$id2,$degree);
				}
			}
		}
	}
    #print "min Composite = $min_AICcombined (@min_pair)\n";
	my $min_AICcombined_degrees_ref = get_degrees($min_pair[0],$min_pair[1],$min_pair[2],$network1,$network2);
	return ($min_AICcombined,\@min_pair,$min_AICcombined_degrees_ref);
}

sub get_AIC_combined
{
	my $id1_founder = shift;
	my $id2_founder = shift;
	my $degree = shift;
	my $network1 = shift;
	my $network2 = shift;
	my $likelihoods = shift;

	#print "$id1_founder,$id2_founder,$degree\n";

	## Given the input degree, get the degree relatedness among all other samples
	my $degrees = get_degrees($id1_founder,$id2_founder,$degree,$network1,$network2);
	
	## Print degrees
	foreach my $id1 (sort {$b cmp $a} keys %$network1)
	{
		foreach my $id2 (sort {$b cmp $a} keys %$network2)
		{
			#print "Degree $id1 <-> $id2 = $$degrees{$id1}{$id2}\n";
			for(1..9)
			{
				#print "lnl $id1 <-> $id2 degree $_ = $$likelihoods{$id1}{$id2}{$_}\n" if $id1_founder eq "A2" and $id2_founder eq "B2";
			}
		}
	}

	
	## Add up the likelihoods for all samples given the degree of relatedness
	my $AIC_combined = 0;
	foreach my $id1 (sort {$b cmp $a} keys %$network1)
	{
		next if $id1 =~ /Missing/i;
		foreach my $id2 (sort {$b cmp $a} keys %$network2)
		{
			next if $id2 =~ /Missing/i;
      $GLOBAL_COUNTER++;
			my $lnl;
      my $AIC;
			my $degree = $$degrees{$id1}{$id2};
			if($degree eq "")
			{
				$degree = "UN";
				$lnl = $$likelihoods{$id1}{$id2}{'UN'};
			}
			elsif($degree > $cap_degree)
			{
				## To a first approximation, given that two individuals A and B are genetic 9th degree relatives, the unconditional likelihood of the 10th degree relationship for individual A and the offspring of individual B is as follows.  With probability one half, the segment is passed on and the likelihood is equal to the 9th degree relationship likelihood.  Otherwise, the likelihood is equal to the unrelated likelihood.  This leads to the following formula for 10+ degrees of relationship:
				## ln((0.5)^(degree_of_relationship-9)*exp(9th_degree_log_likelihood)+(1-(0.5)^(degree_of_relationship-9))*exp(unrelated_log_likelihood))
				my $is_inherited  = ((0.5)**($degree-$cap_degree))*exp($$likelihoods{$id1}{$id2}{$cap_degree});
				my $not_inherited = (1-(0.5)**($degree-$cap_degree))*exp($$likelihoods{$id1}{$id2}{'UN'});
				$lnl = log($is_inherited+$not_inherited);
			}
			else
			{
				$lnl = $$likelihoods{$id1}{$id2}{$degree};
			}
			#print "lnl $id1 <-> $id2 $degree = $lnl\n" if $id1_founder eq "A2" and $id2_founder eq "B2";
			$AIC = 4 - 2 * $lnl;
      $AIC_combined += $AIC;
			#print "$AIC_combined\n" if $id1_founder eq "A2" and $id2_founder eq "B2";
		}
	}
  #print "$id1_founder<->$id2_founder($degree) = $AIC_combined\n" if $verbose > 0;
	return $AIC_combined;
}

sub get_degrees
{
	my $id1 = shift;
	my $id2 = shift;
	my $degree = shift;
	my $network1 = shift;
	my $network2 = shift;

	my %all_degrees;
	
	## If degree is unrelated then return the empty hash
	if($degree eq "UN"){return \%all_degrees}
	
	my @ids_to_visit_net1;
	my @ids_to_visit_net2;
	
	## Prime the arrays with $id1 and set to degree provided
	$all_degrees{$id1}{$id2} = $degree;
	$all_degrees{$id2}{$id1} = $degree;
	push(@ids_to_visit_net1,$id1);
	
	## Proceed to each individual in this ids_to_visit_net1 array and add their relevant relatives to the arrays and set their degree
	my %visited_net1;
	while(@ids_to_visit_net1)
	{
		my $net1_id = shift(@ids_to_visit_net1);
		if(exists $visited_net1{$net1_id}){next}
		$visited_net1{$net1_id} = 1;
		#print "$net1_id\n";
		
		## Set degree of sibs to id2 and push sib onto array to visit
		my @sibs = $$network1{$net1_id}->get_full_sibs($network1);
		foreach my $sib (@sibs)
		{
			if(exists $visited_net1{$sib}){next}
			$all_degrees{$sib}{$id2} = $all_degrees{$net1_id}{$id2};
			$all_degrees{$id2}{$sib} = $all_degrees{$net1_id}{$id2};
			push(@ids_to_visit_net1,$sib);
		}
		
		## Set degree of children and push onto array to visit
		my @children = $$network1{$net1_id}->children();
		foreach my $child (@children)
		{
			if(exists $visited_net1{$child}){next}
			$all_degrees{$child}{$id2} = $all_degrees{$net1_id}{$id2} + 1;
			$all_degrees{$id2}{$child} = $all_degrees{$net1_id}{$id2} + 1;
			push(@ids_to_visit_net1,$child);
		}
		
		## Set degree of parents and push onto array to visit, if neither parents has been seen before, meaning we are dealing with a leaf to founder comparison; the order that they will be pulled off of the array will matter. We want to pull them off working up the tree, so the sibs get the correct degree of relationships. If we move to the top of the tree and work back down, then the sibs will be one degree more distant than they should be.
		my @parents = $$network1{$net1_id}->parents();
		if(@parents > 0 && !exists $visited_net1{$parents[0]} && !exists $visited_net1{$parents[1]})
		{
			foreach my $parent (@parents)
			{
				if(exists $visited_net1{$parent}){next}
				$all_degrees{$parent}{$id2} = $all_degrees{$net1_id}{$id2} + 1;
				$all_degrees{$id2}{$parent} = $all_degrees{$net1_id}{$id2} + 1;
				push(@ids_to_visit_net1,$parent);
			}
		}
		
		push(@ids_to_visit_net2,$id2);
		my %visited_net2;
		## Cycle through each individual in net2
		while(@ids_to_visit_net2)
		{
			my $net2_id = shift(@ids_to_visit_net2);
			if(exists $visited_net2{$net2_id}){next}
			$visited_net2{$net2_id} = 1;
			#print "$net2_id\n";
			
			## Set degree of sibs and push onto array to visit
			my @sibs = $$network2{$net2_id}->get_full_sibs($network2);
			foreach my $sib (@sibs)
			{
				if(exists $visited_net2{$sib}){next}
				$all_degrees{$net1_id}{$sib} = $all_degrees{$net1_id}{$net2_id};
				$all_degrees{$sib}{$net1_id} = $all_degrees{$net1_id}{$net2_id};
				push(@ids_to_visit_net2,$sib);
			}
			
			## Set degree of children and push onto array to visit
			my @children = $$network2{$net2_id}->children();
			foreach my $child (@children)
			{
				if(exists $visited_net2{$child}){next}
				$all_degrees{$net1_id}{$child} = $all_degrees{$net1_id}{$net2_id} + 1;
				$all_degrees{$child}{$net1_id} = $all_degrees{$net1_id}{$net2_id} + 1;
				push(@ids_to_visit_net2,$child);
			}

			## Set degree of parents and push onto array to visit, if neither parents has been seen before, meaning we are dealing with a leaf to founder comparison
			my @parents = $$network2{$net2_id}->parents();
			if(@parents > 0 && !exists $visited_net2{$parents[0]} && !exists $visited_net2{$parents[1]})
			{
				foreach my $parent (@parents)
				{
					if(exists $visited_net2{$parent}){next}
					$all_degrees{$net1_id}{$parent} = $all_degrees{$net1_id}{$net2_id} + 1;
					$all_degrees{$parent}{$net1_id} = $all_degrees{$net1_id}{$net2_id} + 1;
					push(@ids_to_visit_net2,$parent);
				}
			}
		}
	}
	return \%all_degrees;
}

sub get_degrees_within_network
{
	my $network_summary = shift;
	my $all_degrees = shift;

	my $line_position_for_first_pedigree = 9;
	$network_summary =~ /(.*)\/Summary_/;
	my $network1_dir = $1;

	## LOAD network data
	open(NET1,$network_summary) or die "can't open $network_summary: $!\n";
	my @lines1 = <NET1>;
	close(NET1);
	my ($junk1,$junk2, $network1_name) = split(/\s+/,@lines1[0]);
	my $fam_file = "$network1_dir/$network1_name\_1.fam";

	my $network_ref = build_network_from_ped_file($fam_file);
	
	foreach my $ref_id (keys %{$network_ref} )
	{
		#print "ref_id: $ref_id\n";
		## PRIME the array
		$$all_degrees{$ref_id}{$ref_id} = 0;
		my @ids_to_visit;
		push(@ids_to_visit,$ref_id);
		
		## Proceed to each individual in this ids_to_visit_net1 array and add their relevant relatives to the arrays and set their degree
		my %visited;
		while(@ids_to_visit)
		{
			my $id2 = shift(@ids_to_visit);

			if(exists $visited{$id2}){next}
			$visited{$id2} = 1;
			#print "$ref_id:$id2 = $$all_degrees{$id2}{$ref_id}\n";
			## Set degree of sibs to id2 and push sib onto array to visit
			my @sibs = $$network_ref{$id2}->get_full_sibs($network_ref);
			foreach my $sib (@sibs)
			{
				if(exists $visited{$sib}){next}
				push(@ids_to_visit,$sib);
				if(exists $$all_degrees{$sib}{$ref_id} && $$all_degrees{$id2}{$ref_id} + 1 >= $$all_degrees{$sib}{$ref_id}){next}
				#print "sib: $sib\n";
				$$all_degrees{$sib}{$ref_id} = $$all_degrees{$id2}{$ref_id} + 1;
				$$all_degrees{$ref_id}{$sib} = $$all_degrees{$id2}{$ref_id} + 1;
			}
			
			## Set degree of children and push onto array to visit
			my @children = $$network_ref{$id2}->children();
			foreach my $child (@children)
			{
				if(exists $visited{$child}){next}
				push(@ids_to_visit,$child);
				if(exists $$all_degrees{$child}{$ref_id} && $$all_degrees{$id2}{$ref_id} + 1 >= $$all_degrees{$child}{$ref_id}){next}
				$$all_degrees{$child}{$ref_id} = $$all_degrees{$id2}{$ref_id} + 1;
				$$all_degrees{$ref_id}{$child} = $$all_degrees{$id2}{$ref_id} + 1;
			}
			
			## Set degree of parents and push onto array to visit, if neither parents has been seen before, meaning we are dealing with a leaf to founder comparison; the order that they will be pulled off of the array will matter. We want to pull them off working up the tree, so the sibs get the correct degree of relationships. If we move to the top of the tree and work back down, then the sibs will be one degree more distant than they should be.
			my @parents = $$network_ref{$id2}->parents();
			if(@parents > 0 && !exists $visited{$parents[0]} && !exists $visited{$parents[1]})
			{
				foreach my $parent (@parents)
				{
					if(exists $visited{$parent}){next}
					push(@ids_to_visit,$parent);
					if(exists $$all_degrees{$parent}{$ref_id} && $$all_degrees{$id2}{$ref_id} + 1 >= $$all_degrees{$parent}{$ref_id}){next}
					$$all_degrees{$parent}{$ref_id} = $$all_degrees{$id2}{$ref_id} + 1;
					$$all_degrees{$ref_id}{$parent} = $$all_degrees{$id2}{$ref_id} + 1;
				}
			}
		}
	}
	return $all_degrees;
}


sub is_founder
{
	my $id = shift;
	my $net = shift;
	#print "id $id\n";
		
	my @parents = $$net{$id}->parents();
	if(@parents > 0){return 0}
	else{return 1}
}

sub is_leaf
{
	my $id = shift;
	my $net = shift;
	#print "id $id\n";
		
	my @children = $$net{$id}->children();
	if(@children > 0){return 0}
	else{return 1}
}

sub load_ersa_likelihoods
{
	print "Loading ERSA likelihoods\n";
	
	my $lnl_file = shift;
	my $results_file = shift;
	my $network1 = shift;
	my $network2 = shift;

	my %temp_likelihoods;
	my %likelihoods;
	my %pairs;
	my %individuals;

	#open(OUT, ">$file\_networks1_and_170_likelihoods");
	open(IN,$lnl_file);
	<IN>;
	while(my $line = <IN>)
	{
		chomp($line);
		#print "$line\n";
		my ($id1,$id2,$num_shared,$degree,$maxlnl) = split(/\s+/,$line);
		$pairs{$id1}{$id2}=1;
		$individuals{$id1}=1;
		$individuals{$id2}=1;
		if($network1 eq "" || $network2 eq "")
		{
			$temp_likelihoods{$id1}{$id2}{$degree}{$num_shared} = $maxlnl;
			$temp_likelihoods{$id2}{$id1}{$degree}{$num_shared} = $maxlnl;
		}
		elsif((exists $$network1{$id1} && exists $$network2{$id2}) || (exists $$network1{$id2} && exists $$network2{$id1}))
		{
			#print OUT "$line\n";
			#print "$id1 <-> $id2 = $degree\n";
			$temp_likelihoods{$id1}{$id2}{$degree}{$num_shared} = $maxlnl;
			$temp_likelihoods{$id2}{$id1}{$degree}{$num_shared} = $maxlnl;
		}
	}
	close(IN);
	#close(OUT);

	## Load unrelated lnl for each pair from results file
	open(IN,$results_file);
	while(my $line = <IN>)
	{
		next if $line =~ /^#/ || $line =~ /^individual_1/;
		chomp($line);
		my @temp = split(/\s+/,$line);
		my $id1 = @temp[0];
		my $id2 = @temp[1];
		my $lnl_UN = @temp[-1];
		my $best = @temp[3];
		$best =~ s/no_sig_rel/UN/;
		$likelihoods{$id2}{$id1}{'best'} = $best;
		$likelihoods{$id1}{$id2}{'best'} = $best;
		$likelihoods{$id2}{$id1}{'UN'} = $lnl_UN;
		$likelihoods{$id1}{$id2}{'UN'} = $lnl_UN;
	}
	close(IN);

	foreach my $id1 (keys %temp_likelihoods)
	{
		foreach my $id2 (keys %{$temp_likelihoods{$id1} })
		{
			## Get the maximum between one parent in common and two parents in common
			#for(my $degree = $min_degree_to_consider; $degree <= $max_degree_for_related; $degree++)
			for(my $degree = $min_degree_to_consider; $degree <= 40; $degree++)
			{
				# If missing from the likelihood file, set to 0 for each degree
				if(!exists $temp_likelihoods{$id1} || !exists $temp_likelihoods{$id1}{$id2})
				{
					$likelihoods{$id1}{$id2}{$degree} = 0;
					$likelihoods{$id2}{$id1}{$degree} = 0;
					next;
				}

				my $maxlnl_1 = $temp_likelihoods{$id1}{$id2}{$degree}{1};
				my $maxlnl_2 = $temp_likelihoods{$id1}{$id2}{$degree}{2};
				#print "$id1 <-> $id2 1 $degree = $maxlnl_1\n";
				#print "$id1 <-> $id2 2 $degree = $maxlnl_2\n";
				if($maxlnl_1 > $maxlnl_2 && $maxlnl_1 ne "")
				{
					
					my $adjusted_maxlnl_1  = $maxlnl_1-$AIC_correction;
					$likelihoods{$id1}{$id2}{$degree} = $adjusted_maxlnl_1;
					$likelihoods{$id2}{$id1}{$degree} = $adjusted_maxlnl_1;
				}
				else
				{
					my $adjusted_maxlnl_2  = $maxlnl_2-$AIC_correction;
					$likelihoods{$id1}{$id2}{$degree} = $adjusted_maxlnl_2;
					$likelihoods{$id2}{$id1}{$degree} = $adjusted_maxlnl_2;
				}
			}
			
			next if exists $likelihoods{$id2}{$id1}{'UN'};
			
			## get the maximum log likelihood among all the unrelated relationships
			# If missing from the likelihood file, set to 0
			if(!exists $temp_likelihoods{$id1} || !exists $temp_likelihoods{$id1}{$id2})
			{
				$likelihoods{$id1}{$id2}{"UN"} = $likelihoods{$id1}{$id2}{40};
				$likelihoods{$id2}{$id1}{"UN"} = $likelihoods{$id1}{$id2}{40};
				next;
			}
		}
	}
	print "done\n";
	return \%likelihoods;
}


sub build_network_from_ped_file
{
	my $fam_file = shift;
	#print "Building network for $ped_file\n";
	
	## Read in file
	open(IN,$fam_file) or die "Can't open $fam_file: $!\n";
	my %all_nodes_network;
	my $network_ref = \%all_nodes_network;
	my @network_refs;
	
	## Build pedigree
	while(my $line = <IN>)
	{
		chomp($line);
		my ($FID,$IID,$PID,$MID,$SEX,$PHENOTYPE) = split(/\s+/,$line);
		my $child = $IID;
		my $dad = $PID;
		my $mom = $MID;
		
		## If the name is seperated with a '__' from PRIMUS output then change it to just the part after the '__'
		if($IID =~ /__/){$child = (split(/__/,$IID))[1];} # "$FID\__$IID";}
		if($MID =~ /__/){$mom =(split(/__/,$MID))[1];} # "$FID\__$MID";}
		if($PID =~ /__/){$dad =(split(/__/,$PID))[1];} # "$FID\__$MID";}
		
		if(!exists $$network_ref{$child})
		{
			my $node = new node_v7($child);
			$$network_ref{$child} = $node;
		}
		if(!exists $$network_ref{$dad} && $PID ne 0)
		{
			my $node = new node_v7($dad);
			$$network_ref{$dad} = $node;
		}
		if(!exists $$network_ref{$mom} && $MID ne 0)
		{
			my $node = new node_v7($mom);
			$$network_ref{$mom} = $node;
		}
		if($PID ne 0)
		{
			$$network_ref{$child}->add_parent($dad);
			$$network_ref{$dad}->add_child($child);
		}
		
		if($MID ne 0)
		{
			$$network_ref{$child}->add_parent($mom);
			$$network_ref{$mom}->add_child($child);
		}
	}
	close(IN);


	## Break network into individual pedigrees and
	## Build relationships
	my @keys = keys %$network_ref;
	#print "# keys = " . @keys . "\n";
	while (@keys > 0)
	{
		my $node_name = @keys[0];
		
		my %names = $$network_ref{$node_name}->get_subpedigree_names($network_ref);
		my %pedigree;
		
		foreach my $node_name (keys %names)
		{
			$pedigree{$node_name} = $$network_ref{$node_name};
			delete $$network_ref{$node_name};
		}
		push(@network_refs,\%pedigree);
		@keys = keys %$network_ref;
		#print "# keys = " . @keys . "\n";
		#print "# network_refs = " . @network_refs . "\n";
	}
	

	
	## Write out relationships
	foreach my $network_ref (@network_refs)
	{
		foreach my $node_name (keys %$network_ref)
		{
			$$network_ref{$node_name}->make_relative_network_from_pedigree($network_ref);
			my %rels = $$network_ref{$node_name}->relatives();
			foreach(keys %rels)
			{
				#print "$node_name -> $_  = @{$rels{$_} }\n";
			}
		}
	}
	
	if(@network_refs > 1){die "should not be more than one network: $fam_file\n";}
	return @network_refs[0];	
}

### Write the .dot file that can be read into a graph visualization program like Graphviz
sub write_out_dot_file
{
	my $file = shift;
	my $networks_ref = shift;

	open(GRAPH_OUT,">$file");
	print GRAPH_OUT "graph Network_Relationships{\n";
	print GRAPH_OUT "\tnode [shape=circle];\n\n";
	
	my %added;

	## Write out all connections
	foreach my $id1 (keys %{$networks_ref})
	{
		#print "id1: $id1\n";
		foreach my $id2 (keys %{$$networks_ref{$id1} })
		{
			#print "id2: $id2\n";
			next if exists $added{$id1}{$id2};
			$added{$id1}{$id2} = 1;
			$added{$id2}{$id1} = 1;
			my $degree = $$networks_ref{$id1}{$id2};
			
			if($degree ne "UN")
			{
				my $color = "black";
				print GRAPH_OUT "\t\"$id1\" -- \"$id2\" [color=$color;label=\"$degree\"];\n";
			}
		}
	}
	print GRAPH_OUT "}";
	close(GRAPH_OUT);
}



return 1;
