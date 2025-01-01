package PRIMUS::reconstruct_pedigree_v7;

#### Version 4 is attempting to resolve the infinite run times by splitting the likelihood vectors on the fly, rather than up front.
#### Version 5 is attempting to make the flow of Phase1 match the flow of phases 2 and 3.
#### Version 6 is attempting to leverage sibling blocks to speed up reconstruction; and fixed the sib-mating check; and error catching during reconstruction (note: this may mask bugs in the code)
#### Version 7 is correcting several bugs, cleaning up the code, and minor fixes in preparation for the initial public release

use strict;
use Getopt::Long qw(GetOptionsFromArray);
use PRIMUS::node_v7;
use PRIMUS::get_age_flags;
use PRIMUS::predict_relationships_2D qw(get_relationship_likelihood_vectors);
use PRIMUS::compare_fam_files;
use List::Util qw(sum);

## Code to build a network and print out the pairwise relationships from a ped/fam file
#my @network_refs_from_fam = build_network_from_fam_file("../data/PCfamilies.fam");
#write_out_relationships("../data/PCrelationships_test",@network_refs_from_fam);
#exit;

## Settings
my $verbose = 1;
my $lib_dir;
my $bin_dir;
my $ersa_data;

## Necessary variables
my $MAX_MATING_GENERATION_GAP = 0;
my $MAX_RUNTIME = 129600; # 43200 = 12 hours; 129,600 = 36 hours
my $MAX_NETWORKS_TO_RESOLVE = 100000;
my $MAX_GENERATIONS = 5;
my $affected_status_value = 2;
my $MIN_LIKELIHOOD = 0.1; ## should be .1, adjusted in OG primus paper to 0.375 
my $NO_MITO = 0;
my $NO_Y = 0;
my $USE_NO_MATCH_MITO = 0; ## Will be updated to 1 if there is a mito file passed in
my $USE_MATCH_MITO = 0; ## Unless you have extremely accurate and high precission MT relationships, this option should remain off.
my $USE_NO_MATCH_Y = 0; ## Will be updated to 1 if there is a mito file passed in
my $USE_MATCH_Y = 0; ## Unless you have extremely accurate and high precission MT relationships, this option should remain off.
my @likelihood_names = qw(PC FS HAG CGH DR UN);
my $cranefoot_binary = "cranefoot";
if($^O =~ /darwin/i){$cranefoot_binary = "cranefoot_mac"}
my $allow_half_sib_dummy_mating = 0;
my $allow_full_sib_dummy_mating = 0;
my $USE_SIBLING_BLOCKS = 1; ## Setting this to 1 will speed up the runtime, but will reduce the number of possible pedigrees that are inbred, but it will still produce all possible outbred pedigrees.
my $dummy_ctr = 1;
my $LOG;
my %age;
my $enforce_age_filtering = 0;
my $degree_rel_cutoff = 3;

open($LOG,"");

sub reconstruct_pedigree
{
	my @commands = @_;
	$dummy_ctr = 1;
	my $network_name;
	my $IBD_file_ref;
	my $MITO_file_ref;
	my $Y_file_ref;
	
	my %affected_status;
	my %gender;
	
	my $output_directory;
	my $network_num = 1;
	my @sample_names;
	my @networks;
	my @all_possible_networks;
	my %scores;
	my %num_dummies;
	my %num_generations;
	my %age_flags;

	GetOptionsFromArray(
		\@_,
		# Diagnostic options
		'verbose=i' => \$verbose,
		
		# Settings
		"network=s" => \$network_name, 
		"int_likelihood_cutoff=f" => \$MIN_LIKELIHOOD, 
		"output_dir=s" => \$output_directory,
		"ibd_estimates=s" => \$IBD_file_ref,
		"FILE=s"=> sub{$$IBD_file_ref{'FILE'}=$_[1]},
		"FID1=s"=> sub{$$IBD_file_ref{'FID1'}=$_[1]},
		"IID1=s"=> sub{$$IBD_file_ref{'IID1'}=$_[1]},
		"FID2=s"=> sub{$$IBD_file_ref{'FID2'}=$_[1]},
		"IID2=s"=> sub{$$IBD_file_ref{'IID2'}=$_[1]},
		"IBD0=s"=> sub{$$IBD_file_ref{'IBD0'}=$_[1]},
		"IBD1=s"=> sub{$$IBD_file_ref{'IBD1'}=$_[1]},
		"IBD2=s"=> sub{$$IBD_file_ref{'IBD2'}=$_[1]},
		
		"mito_hash=s" => \$MITO_file_ref,
		"MITO_FILE=s"=> sub{$$MITO_file_ref{'FILE'}=$_[1]},
		"MITO_FID1=s"=> sub{$$MITO_file_ref{'FID1'}=$_[1]},
		"MITO_IID1=s"=> sub{$$MITO_file_ref{'IID1'}=$_[1]},
		"MITO_FID2=s"=> sub{$$MITO_file_ref{'FID2'}=$_[1]},
		"MITO_IID2=s"=> sub{$$MITO_file_ref{'IID2'}=$_[1]},
		"MITO_MATCH=s"=> sub{$$MITO_file_ref{'MATCH'}=$_[1]},
		"MITO_MATCH_VAL=s"=> sub{$$MITO_file_ref{'MATCH_VAL'}=$_[1]},
		
		"y_hash=s" => \$Y_file_ref,
		"Y_FILE=s"=> sub{$$Y_file_ref{'FILE'}=$_[1]},
		"Y_FID1=s"=> sub{$$Y_file_ref{'FID1'}=$_[1]},
		"Y_IID1=s"=> sub{$$Y_file_ref{'IID1'}=$_[1]},
		"Y_FID2=s"=> sub{$$Y_file_ref{'FID2'}=$_[1]},
		"Y_IID2=s"=> sub{$$Y_file_ref{'IID2'}=$_[1]},
		"Y_MATCH=s"=> sub{$$Y_file_ref{'MATCH'}=$_[1]},
		"Y_MATCH_VAL=s"=> sub{$$Y_file_ref{'MATCH_VAL'}=$_[1]},
		
		"no_mito=i" => \$NO_MITO,
		"no_y=i" => \$NO_Y,
		"use_mito_match=i" => \$USE_MATCH_MITO,
		"use_y_match=i" => \$USE_MATCH_Y,


		"lib=s"=>\$lib_dir,
		"bin=s"=>\$bin_dir,
		"log_file_handle=s"=>\$LOG,
		"sex_file=s"=> sub{
			open(IN,$_[1]);
			while(my $line = <IN>){
				chomp($line);
				my @temp = split(/\s+/,$line);
				$gender{$temp[0]}= $temp[1];
			} 
		},
		"age_file=s"=> sub{
			open(IN,$_[1]);
			while(my $line = <IN>){
				chomp($line);
				my @temp = split(/\s+/,$line);
				$age{$temp[0]}= $temp[1];
			} 
		},
		"affection_file=s"=> sub{
			open(IN,$_[1]);
			while(my $line = <IN>){
				chomp($line);
				my @temp = split(/\s+/,$line);
				$affected_status{$temp[0]}= $temp[1];
			} 
		},
		"PI_HAT|RELATEDNESS=s"=> sub{$$IBD_file_ref{'PI_HAT'}=$_[1]},
		"degree_rel_cutoff=i"=> \$degree_rel_cutoff,
		"sex_ref=s" => sub{%gender = %{ $_[1] } },
		"affection_ref=s" => sub{ %affected_status = %{ $_[1] } },
		"affection_status_value=s"=> \$affected_status_value,
		"age_ref=s" => sub{%age = %{ $_[1] } },
		"network_num=i" => \$network_num,
		"max_gen=i" => \$MAX_GENERATIONS,

	) or die "Failed to parse options for Pedigree Reconstruction\n";
	
	open($LOG,">$output_directory/$network_name.log") if($LOG eq "");

	PRIMUS::node_v7::set_min_likelihood($MIN_LIKELIHOOD);
	PRIMUS::node_v7::set_verbose($verbose);
	
	my $mito_ref = load_mito_data($MITO_file_ref) if exists $$MITO_file_ref{'FILE'};
	my $y_ref = load_y_data($Y_file_ref) if exists $$Y_file_ref{'FILE'};
	
	$age_flags{"age_file"}=1 if(keys %age > 0);
	print "RECONSTRUCTING $network_name\n" if($verbose > 0);
	print "Output directory: $output_directory\n" if($verbose > 0);
	print "Use mito non-match: $USE_NO_MATCH_MITO\n" if($verbose > 0);
	print "Use mito match: $USE_MATCH_MITO\n" if($verbose > 0);
	print "Use Y non-match: $USE_NO_MATCH_Y\n" if($verbose > 0);
	print "Use Y match: $USE_MATCH_Y\n" if($verbose > 0);
	print "min likelihood: $MIN_LIKELIHOOD\n" if($verbose > 1);
	print "MAX GENERATIONS: $MAX_GENERATIONS\n" if($verbose > 1);
	print $LOG "RECONSTRUCTING $network_name\n" if($verbose > 0);
	print $LOG "Output directory: $output_directory\n" if($verbose > 0);
	print $LOG "min likelihood: $MIN_LIKELIHOOD\n" if($verbose > 1);
	print $LOG "MAX GENERATIONS: $MAX_GENERATIONS\n" if($verbose > 1);
	
	## Load network data
	## $relationships_ref = a reference to a hash of hashes of arrays; format is $relationship{$id1}{$id2} = [likelihood of relationships]
	my ($relationships_ref,$raw_relationship_densities_ref,$total_possibilities, @fails);
	eval
	{
		print "\nmin_likelihood: $MIN_LIKELIHOOD\n\n";
		print $LOG "min_likelihood: $MIN_LIKELIHOOD\n" if $verbose > 1;

		###################################
		# This is where the relationship vector data gets read in from predict_relationships_2D.pm

		print "\nIBD file: $MIN_LIKELIHOOD\n\n";
		print $LOG "min_likelihood: $MIN_LIKELIHOOD\n" if $verbose > 1;
		
		($relationships_ref,$raw_relationship_densities_ref,$total_possibilities, @fails) = PRIMUS::predict_relationships_2D::get_relationship_likelihood_vectors($IBD_file_ref,$MIN_LIKELIHOOD,$verbose,$lib_dir,$output_directory);

		# can we just run the LOAD function here instead? why do we need to do this twice 
	
		###################################
	};
	if($@)
	{	
		## If it crashed, I need to get all the sample IDs and the write the summary file.
		my %samples;
		my $IBD_file = $$IBD_file_ref{'FILE'};
		open(IN,$IBD_file) or die "Can't open $IBD_file; $!";
		my $header = <IN>;
		while(my $line = <IN>)
		{
			$line =~ s/^\s+//;
			chomp($line);
			my @temp = split(/\s+/,$line);
			my $FID1 = @temp[$$IBD_file_ref{'FID1'}-1]; 
			my $IID1 = @temp[$$IBD_file_ref{'IID1'}-1];
			my $FID2 = @temp[$$IBD_file_ref{'FID2'}-1];
			my $IID2 = @temp[$$IBD_file_ref{'IID2'}-1];

			my $name1 = "$IID1";
			my $name2 = "$IID2";
			$samples{$name1}=1;
			$samples{$name2}=1;
		}
		my @sample_names = keys %samples;

		print "FAILED RELATIONSHIP PREDICTION: $@\n";
		print $LOG "FAILED RELATIONSHIP PREDICTION: $@\n";
		chomp($@);
		write_summary_file(\%scores,\%num_dummies,\%num_generations,$output_directory,$network_name,\@networks,\@sample_names,"",$@);
		return $total_possibilities;
	};


	my %relationships = %{$relationships_ref};
	my $network_ref = load_network(\%relationships,$mito_ref,$y_ref,\%gender);
	foreach(keys %$network_ref)
	{
		push(@sample_names,$_);
	}
	

	## Reconstruct and catch errors so it can still write out the summary file
	eval
	{
		@networks = reconstruct_network($network_ref);
		1;
	};
	if($@)
	{
		print "FAILED RECONSTRUCTION: $@\n";
		print $LOG "FAILED RECONSTRUCTION: $@\n";
		write_summary_file(\%scores,\%num_dummies,\%num_generations,$output_directory,$network_name,\@networks,\@sample_names,"",$@);
		return $total_possibilities;
	};
		
	
	## Test if timed out
	if (@networks[0] == -100)
	{
		print "TIMED OUT\n";
		print $LOG "TIMED OUT\n";
		write_summary_file(\%scores,\%num_dummies,\%num_generations,$output_directory,$network_name,\@networks,\@sample_names,"","TIMED OUT");
		return $total_possibilities;
	}
	
	## Test if too many possible pedigrees
	if (@networks[0] == -101)
	{
		print "TOO many possible pedigrees\n";
		print $LOG "TOO many possible pedigrees\n";
		write_summary_file(\%scores,\%num_dummies,\%num_generations,$output_directory,$network_name,\@networks,\@sample_names,"","Too many possible pedigrees");
		return $total_possibilities;
	}
	
	## Rank
	my @ranked_networks = rank_networks($raw_relationship_densities_ref,\%scores,\%num_dummies,\%num_generations,\%age,\%age_flags,@networks);

	if($verbose > 0)
	{
		print "networks pre-prune: ".@ranked_networks."\n";
		print $LOG "networks pre-prune: ".@ranked_networks."\n";
	}

	## Prune
	@ranked_networks = post_prune_networks(@ranked_networks);
	if($verbose > 0)
	{
		print "networks post-prune: ".@ranked_networks."\n";
		print $LOG "networks post-prune: ".@ranked_networks."\n";
	}

	## Write out pedigree files
	write_summary_file(\%scores,\%num_dummies,\%num_generations,$output_directory,$network_name,\@ranked_networks,\@sample_names,\%age_flags,"none");

	my $possible_pedigree_ctr = 1;	
	foreach my $network_ref (@ranked_networks)
	{
		my $network_name = "$network_name\_$possible_pedigree_ctr";
		print "Writing .fam file for $network_name\n" if $verbose > 0;
		print $LOG "Writing .fam file for $network_name\n" if $verbose > 0;
		eval
		{
			write_fam_file($network_ref,$network_name,\%affected_status,$output_directory,$network_num);
			$possible_pedigree_ctr++;
			write_cranefoot_files($network_ref,"$network_name",\%affected_status,$output_directory,$network_num,\%age);
			#exit;
		};
		if($@)
		{
			warn "$network_name failed to write\n";
		};
	}

	push(@all_possible_networks,@ranked_networks);
	
	## Reduce min likelihood and run again if no networks reconstructed

	if(@all_possible_networks eq 0 && $MIN_LIKELIHOOD > 0.01) ## SHOULD BE 0.01
	{
		if($MIN_LIKELIHOOD <= .1)
		{
			$MIN_LIKELIHOOD = $MIN_LIKELIHOOD/3;
		}
		else
		{
			$MIN_LIKELIHOOD = $MIN_LIKELIHOOD - 0.1;
		}

		for(my $i = 0; $i < @commands; $i++)
		{
			if(@commands[$i] eq "--int_likelihood_cutoff")
			{
				@commands[$i+1] = $MIN_LIKELIHOOD;
				last;
			}
		}
		($total_possibilities, @all_possible_networks) = reconstruct_pedigree(@commands); 
	}


	$MIN_LIKELIHOOD = 0.3;
	return ($total_possibilities, @all_possible_networks);
}

####################################################################################
####################################################################################
####################################################################################
sub load_y_data
{
	return if $NO_Y == 1;
	my $y_ref = shift;
	my %y = %{$y_ref};
	my $file = $y{'FILE'};
	open(IN,$file) or warn "Can't open y file: $file; $!" ;
	if(tell(IN) == -1)
	{
		return "";
	}
	
	$USE_NO_MATCH_Y = 1;

	my $FID1_col = $y{'FID1'};
	my $IID1_col = $y{'IID1'};
	my $FID2_col = $y{'FID2'};
	my $IID2_col = $y{'IID2'};
	my $match_col = $y{'MATCH'};

	my %match_data;
	#print "file: $file\n";
	#print "$FID1_col :: $IID1_col\n";
	while(my $line = <IN>)
	{
		$line =~ s/^\s+//;
		my @temp = split(/\s+/,$line);
		my $FID1 = @temp[$FID1_col-1];
		my $IID1 = @temp[$IID1_col-1];
		my $FID2 = @temp[$FID2_col-1];
		my $IID2 = @temp[$IID2_col-1];
		my $val = @temp[$match_col-1];
		
		my $match = 0;
		$match = 1 if $val == $y{'MATCH_VAL'};
		$match_data{"$FID1\__$IID1"}{"$FID2\__$IID2"} = $match;
		$match_data{"$FID2\__$IID2"}{"$FID1\__$IID1"} = $match;
		#print "$FID2\__$IID2 <-> $FID1\__$IID1 = ". $match_data{"$FID2\__$IID2"}{"$FID1\__$IID1"} ."\n";
	}
	return \%match_data;
}

sub load_mito_data
{
	return if $NO_MITO == 1;
	my $mito_ref = shift;
	my %mito = %{$mito_ref};
	my $file = $mito{'FILE'};
	open(IN,$file) or warn "Can't open mito file: $file; $!" ;
	if(tell(IN) == -1)
	{
		return "";
	}
	
	print "NO_MITO: $NO_MITO\n" if $verbose > 1;

	$USE_NO_MATCH_MITO = 1;
	print "USE_NO_MATCH_MITO: $USE_NO_MATCH_MITO\n" if $verbose > 1;
	print "USE_MATCH_MITO: $USE_MATCH_MITO\n" if $verbose > 1;

	my $FID1_col = $mito{'FID1'};
	my $IID1_col = $mito{'IID1'};
	my $FID2_col = $mito{'FID2'};
	my $IID2_col = $mito{'IID2'};
	my $match_col = $mito{'MATCH'};

	my %match_data;
	#print "file: $file\n";
	#print "$FID1_col :: $IID1_col\n";
	while(my $line = <IN>)
	{
		$line =~ s/^\s+//;
		my @temp = split(/\s+/,$line);
		my $FID1 = @temp[$FID1_col-1];
		my $IID1 = @temp[$IID1_col-1];
		my $FID2 = @temp[$FID2_col-1];
		my $IID2 = @temp[$IID2_col-1];
		my $val = @temp[$match_col-1];
		
		my $match = 0;
		$match = 1 if $val == $mito{'MATCH_VAL'};
		$match_data{"$FID1\__$IID1"}{"$FID2\__$IID2"} = $match;
		$match_data{"$FID2\__$IID2"}{"$FID1\__$IID1"} = $match;
		#print "$FID2\__$IID2 <-> $FID1\__$IID1 = ". $match_data{"$FID2\__$IID2"}{"$FID1\__$IID1"} ."\n";
	}
	return \%match_data;
}

sub write_summary_file
{
	print "Writing summary file\n" if($verbose > 0);
	print $LOG "Writing summary file\n" if($verbose > 0);
	
	my $scores = shift;
	my $num_dummies = shift;
	my $num_generations = shift;
	my $output_directory = shift;
	my $network_name = shift;
	my $ranked_networks = shift;
	my $non_dummy_samples = shift;
	my $age_flags = shift;
	my $error_message = shift;
	
	my $min = 100000000000;
	my $max = -1000000000000;
	my $mean = "NA";
	my $median = "NA";
	my $num_pedigrees = @$ranked_networks;
	my $num_pedigrees_unflagged = $num_pedigrees;
	my $num_samples = 0;
	my $num_at_max_score = 0;

	my @lines;
	my @scores_arr;
	my $sum;
	my @non_dummies = @$non_dummy_samples;
	my $ctr = 0;
	my $scoring_rank = 0;
	my $ranks_score = 10000000000; 
	my $network_ref = @$ranked_networks[0];
	my $num_samples = @non_dummies;
	foreach my $name (@$ranked_networks)
	{
		if($$scores{$name} > $max){ $max=$$scores{$name} }
		if($$scores{$name} < $min){ $min=$$scores{$name} }
		$sum += $$scores{$name};
		push(@scores_arr,$$scores{$name});
		## Since the score are ranked, once you find a score lower than max, then you are out of the top scoring pedigrees
		if($$scores{$name} < $max && $num_at_max_score eq 0)
		{
			$num_at_max_score = $ctr;
		}
		$ctr++;
		## everytime the score drops, update the rank based on score of the pedigree to be the number of pedigrees seen ($ctr)
		if($$scores{$name} < $ranks_score)
		{
			$ranks_score = $$scores{$name};
			$scoring_rank = $ctr;
		}
		#### Make the summary line for the file
		my $line = "$ctr\t$$num_dummies{$name}\t$$num_generations{$name}\t$scoring_rank\t$$scores{$name}";
		## Look for age flagged pedigrees
		if(keys %{$$age_flags{$name} } > 0)
		{
			$num_pedigrees_unflagged--;
		}
		## add the reason(s) it was flagged
		foreach my $id1 (keys %{$$age_flags{$name} })
		{
			foreach my $id2 (keys %{$$age_flags{$name}{$id1} } )
			{
				$line .= "$id1<->$id2=$$age_flags{$name}{$id1}{$id2};";
			}
		}
		push(@lines,$line);
	}
	## Calculate statistics for the scores, but these aren't output in this version
	if($num_at_max_score eq 0){$num_at_max_score = $ctr}
	if($num_pedigrees % 2 eq 0)
	{
		my $pos = $num_pedigrees/2 - 1;
		my $pos2 = $num_pedigrees/2;
		$median = (@scores_arr[$pos] + @scores_arr[$pos2])/2;
	}
	else
	{
		my $pos = $num_pedigrees/2;
		$median = (@scores_arr[$pos] + @scores_arr[$pos])/2;
	}
	if($ctr ne 0)
	{	
		$mean = sprintf("%.4f", $sum/$ctr);
		$median = sprintf("%.4f", $median);
		$max = sprintf("%.4f", $max);
		$min = sprintf("%.4f", $min);
	}
	else
	{
		$mean="NA";
		$median="NA";
		$max="NA";
		$min="NA";
	}
	
	open(OUT,">$output_directory/Summary_$network_name.txt");
	print OUT "## Network_name: $network_name\n";
	print OUT "## Num_pedigrees: $num_pedigrees\n";
	print OUT "## Num_pedigrees_unflagged: $num_pedigrees_unflagged\n";
	print OUT "## Num_samples: $num_samples\n";
	print OUT "## Num_at_max_score: $num_at_max_score\n";
	print OUT "## Sample_IIDs: " . join(',',@non_dummies) . "\n";
	print OUT "## Error_message: $error_message\n";
	print OUT "\n";
	print OUT "Pedigree_num\tnum_dummies\tnum_generations\trank\tpedigree_composite_lnl\tAge_flags\n";

	foreach my $line (@lines)
	{
		print OUT "$line\n";
	}
	close(OUT);
}

sub rank_networks
{
	my $vectors = shift;
	my $likelihoods_ref = shift;
	my $num_dummies_ref = shift;
	my $num_generations_ref = shift;
	my $ages_ref = shift;
	my $age_flags_ref = shift;
	my @networks = @_;
	my %num_flags;
	my %network_hash;

	## Foreach network
	foreach my $network_ref (@networks)
	{
		$num_flags{$network_ref} = 0;
		$network_hash{$network_ref} = $network_ref;
		my $ctr = 0;
		my $dummy_ctr = 0;
		my $network_score = 0;

		$$num_generations_ref{$network_ref} = get_num_generations($network_ref);

		## Foreach node in network
		foreach my $node_name (sort {$a cmp $b} keys %$network_ref)
		{
			if($node_name =~ /Dummy/i){$dummy_ctr++;next;};
			## Foreach other node in network
			foreach my $rel (sort {$a cmp $b} keys %$network_ref)
			{
				if($node_name eq $rel){last;}
				if($node_name =~ /Dummy/i || $rel =~ /Dummy/i){next;}
				
				## Get the passed in likelihood vector
				my @likelihood_vector;
				if(exists $$vectors{$node_name}{$rel})
				{
					@likelihood_vector = @{ $$vectors{$node_name}{$rel} };
				}
				elsif(exists $$vectors{$rel}{$node_name})
				{
					@likelihood_vector = @{ $$vectors{$rel}{$node_name} };
				}
				else
				{
					## A likelihood vector is not a available
					#print "$node_name and $rel don't have a likelihood vector\n";
				}
				
				my $score = 0.0001; ## non-zero just in class no relations class had a likelihood above the MIN_LIKELIHOOD
				for(my $i = 0; $i < @likelihood_vector; $i++)
				{
					if(@likelihood_vector[$i] >= $MIN_LIKELIHOOD)
					{
						my $val = $$network_ref{$node_name}->does_relationship_exist_in_pedigree($network_ref,$rel,@likelihood_names[$i] , 10);
						if($val ne 0)
						{
							if($score != 0.0001){die "score is not 0.0001 ($score), meaning that there are two possible relationships between two people\n"}
							$score = @likelihood_vector[$i];
							last; ## this is currently a problem be an problem for consanginous or complex pedigrees where two people have more than one way they are connected to each other
						}
					}
				}
				$network_score = $network_score + log($score); 
				#print "network_score= $network_score\n";
				$ctr++;
			}
		}
		if($ctr != 0)
		{
			$$likelihoods_ref{$network_ref} = $network_score;
		}
		my $num_unnecessary_dummies = get_num_unnecessary_dummies($network_ref);
		$$num_dummies_ref{$network_ref} = $dummy_ctr - $num_unnecessary_dummies;
		
		## Test ages to rank by age flags
		if($$age_flags_ref{"age_file"} ne "")
		{
			PRIMUS::get_age_flags::get_age_flags_in_network($network_ref,$ages_ref,$age_flags_ref);
			#print "Flags:\n";
			foreach my $id1 (keys %{ $$age_flags_ref{$network_ref} })
			{
				foreach my $id2 (keys %{ $$age_flags_ref{$network_ref}{$id1} })
				{
					#print "$id1 -> $id2 = $$age_flags_ref{$network_ref}{$id1}{$id2}\n";
					$num_flags{$network_ref}++;
				}
			}
			#print "Done with age flags\n";
			#exit;
		}
		else
		{
			#print "no age file\n";
			#exit;
		}
	}

	my @sorted_networks;
	foreach my $network (reverse sort { $num_flags{$b} <=> $num_flags{$a} || $$likelihoods_ref{$a} <=> $$likelihoods_ref{$b} || $$num_generations_ref{$b} <=> $$num_generations_ref{$a} || $$num_dummies_ref{$b} <=> $$num_dummies_ref{$a} } keys %$likelihoods_ref)
	{
		push(@sorted_networks,$network_hash{$network});
	}
	return @sorted_networks;
}

sub get_num_generations
{
	my $network_ref = shift;
	my @names = keys %$network_ref;
	my $self_name = @names[0];
	my $ctr = 1;
	while($self_name =~ /Dummy/)
	{
		$self_name = @names[$ctr];
		$ctr++;
	}
	my @nodes_to_explore;
	
	my $min_gen = 0;
	my $max_gen = 0;

	my %relatives;
	$relatives{$self_name} = 0;
	push(@nodes_to_explore,$self_name);
	
	## Set generation values
	while(my $name = shift(@nodes_to_explore) )
	{
		my @parents = $$network_ref{$name}->parents();
		foreach my $parent (@parents)
		{
			if(exists $relatives{$parent})
			{
				next;
			}
			push(@nodes_to_explore,$parent);
			$relatives{$parent} = $relatives{$name} + 1;
			if($parent =~ /Dummy/i){next;}
			if($parent !~ /Dummy/i && $relatives{$parent} > $max_gen){$max_gen = $relatives{$parent}; }
			if($parent !~ /Dummy/i && $relatives{$parent} < $min_gen){$min_gen = $relatives{$parent}; }
		}
		
		my @children = $$network_ref{$name}->children();
		foreach my $child(@children)
		{
			if(exists $relatives{$child})
			{
				next;
			}
			push(@nodes_to_explore,$child);
			$relatives{$child} = $relatives{$name} - 1;
			if($child =~ /Dummy/i){next;}
			if($child !~ /Dummy/i && $relatives{$child} < $min_gen){$min_gen = $relatives{$child}; }
			if($child !~ /Dummy/i && $relatives{$child} > $max_gen){$max_gen = $relatives{$child}; }
		}
	}

	## Get the max and min generations for the total generations
	my $total_gen = $max_gen + abs($min_gen) + 1; ## + 1 for the founder generation
	return $total_gen;
}

sub get_num_unnecessary_dummies
{
	my $network_ref = shift;
	my %unnecessary_dummies;
	
	## Remove unnecessary dummy parents;
	## For each node, remove the parents if:
	## 1. both parents are dummy parents
	## 2. Both dummy parents don't have parents
	## 3. Both parents have node as only child
	foreach my $node_name (keys %$network_ref)
	{
		if(!exists $$network_ref{$node_name}){next;}
		my @parents = $$network_ref{$node_name}->parents();
		if(@parents eq 0 || @parents[0] !~ /Dummy/ || @parents[1] !~ /Dummy/){next;}

		my @parents0 = $$network_ref{@parents[0]}->parents();
		my @parents1 = $$network_ref{@parents[1]}->parents();
		if(@parents0 > 0 || @parents1 > 0){next;}
		
		my @children0 = $$network_ref{@parents[0]}->children();
		my @children1 = $$network_ref{@parents[1]}->children();
		if(@children0 > 1 || @children1 > 1){next;}
		
		## Passed all criteria, label as unnecessary dummy parents
		$unnecessary_dummies{@parents[0]} = 1;
		$unnecessary_dummies{@parents[1]} = 1;
	}
	my @temp = keys %unnecessary_dummies;
	my $num = @temp;
	return $num;
}

sub post_prune_networks
{
	my @networks = @_;
	my @gender_pass_networks;
	my @unique_networks;
	foreach my $network_ref (@networks)
	{
		## Remove unnecessary dummy parents;
		## For each node, remove the parents if:
		## 1. both parents are dummy parents
		## 2. Both dummy parents don't have parents
		## 3. Both parents have node as only child
		foreach my $node_name (keys %$network_ref)
		{
			if(!exists $$network_ref{$node_name}){next;}
			my @parents = $$network_ref{$node_name}->parents();
			if(@parents eq 0 || @parents[0] !~ /Dummy/ || @parents[1] !~ /Dummy/){next;}

			my @parents0 = $$network_ref{@parents[0]}->parents();
			my @parents1 = $$network_ref{@parents[1]}->parents();
			if(@parents0 > 0 || @parents1 > 0){next;}
			
			my @children0 = $$network_ref{@parents[0]}->children();
			my @children1 = $$network_ref{@parents[1]}->children();
			if(@children0 > 1 || @children1 > 1){next;}
			
			## Passed all criteria, remove dummy parents
			$$network_ref{$node_name}->delete_parents();
			delete $$network_ref{@parents[0]};
			delete $$network_ref{@parents[1]};
		}
	}
	
	## Do the gender checks
	foreach my $network_ref (@networks)
	{
		my $genders_pass = get_genders($network_ref);
		if($genders_pass == 0){next;}
		push(@gender_pass_networks,$network_ref);
	}
	
	my @unique_networks;
	my %fams;
	## This next section is trying to remove duplicate pedigrees but may not be necessary if I don't make them during the program
	foreach my $network_ref (@gender_pass_networks)
	{
		## Do the gender checks
		my $genders_pass = get_genders($network_ref);
		if($genders_pass == 0){next;}
		
		if(!exists $fams{$network_ref})
		{
			$fams{$network_ref} = get_fam($network_ref,$genders_pass);
		}

		## Store unique networks in array of network_refs
		## To check if fam files are the same. 
		## 1. Then compare each network to every other network already added to the unique least
		## 2. If not already in the unique list, add it.
		
		## Check that the network ref does not match any of the unique networks
		my $unique = 1;
		foreach my $u_network_ref (@unique_networks)
		{
			if(PRIMUS::compare_fam_files::are_fams_same($fams{$u_network_ref},$fams{$network_ref}))
			{
				$unique = 0;
				last;
			}
		}
		
		if($unique == 1)
		{
			push(@unique_networks,$network_ref);
		}
	}
	return @unique_networks;
}

## get the data in the .fam format
sub get_fam
{
	my $network_ref = shift;
	my $gender_ref = shift;
	my %fam;

	foreach my $ID (sort keys %$network_ref)
	{
		my $FID;
		my $node_name = $ID;
		my @parent_IDs = $$network_ref{$ID}->parents();
		my @parents; ## without the FID
		foreach(@parent_IDs)
		{
			push(@parents,$_);
		}
		
		my $g = $$network_ref{$ID}->sex();

		my $a = 0;
		if(@parents eq 0){@parents = (0,0);}
		if(@parents eq 1){push @parents, 0;}
		if(@parents ne 2)
		{
			die "ERROR!!! Incorrect number of parents ". @parents .": @parents\n";
		}
		my $pat;
		my $mat;
		if(exists $$gender_ref{@parent_IDs[0]} && $$gender_ref{@parent_IDs[0]} == $$gender_ref{@parent_IDs[1]} && $$gender_ref{@parent_IDs[0]} ne 0)
		{
			if($verbose > 1)
			{
				warn "WARNING3!!! Parents @parents[0] and @parents[1] are the same gender ". $$gender_ref{@parents[1]}."!!!\n";
			}
		}
		elsif($$gender_ref{@parent_IDs[0]} == 1)
		{
			$pat = @parents[0];
			$mat = @parents[1]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 2)
		{
			$pat = @parents[1];
			$mat = @parents[0]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 0 && $$gender_ref{@parent_IDs[1]} == 2)
		{
			$pat = @parents[0]; 
			$mat = @parents[1]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 0 && $$gender_ref{@parent_IDs[1]} == 1)
		{
			$pat = @parents[1]; 
			$mat = @parents[0]; 
		}
		else
		{
			$pat = @parents[1]; 
			$mat = @parents[0];
		}
		$pat =~ s/Dummy/Missing/g;
		$mat =~ s/Dummy/Missing/g;
		$node_name =~ s/Dummy/Missing/g;
		$$gender_ref{$pat} = 1;
		$$gender_ref{$mat} = 2;
		$fam{$node_name}{"PID"}=$pat;
		$fam{$node_name}{"MID"}=$mat;
	}
	return \%fam;
}

sub reconstruct_network
{
	#### Phase 1: resolve all possible relationships that are 1st degree
	my @networks_unresolved = @_;
	my $start_time = time;
	my @networks_resolved0;

	print "Entering Resolve PC trios. # of possible pedigrees: ". @networks_unresolved ."\n" if $verbose > 0;
	
	while(@networks_unresolved > 0)
	{
		print "# of networks to resolve: " . @networks_unresolved . "\n" if($verbose > 1);
		print $LOG "# of networks to resolve: " . @networks_unresolved . "\n" if($verbose > 1);
		my $network_ref = shift(@networks_unresolved);
		
    my ($resolved, @temp_networks) = resolve_PC_trio($network_ref);
    
    print "# temp_networks: " . @temp_networks . "\n" if($verbose > 1);
    push(@networks_unresolved,@temp_networks);

		if($resolved eq 1)
		{
      print "RESOLVED\n" if($verbose > 1);
			push(@networks_resolved0,$network_ref);
		}
		if(@networks_unresolved > $MAX_NETWORKS_TO_RESOLVE)
		{
			die "> $MAX_NETWORKS_TO_RESOLVE networks to resolve (Phase 1)";
		}
	}

	foreach my $network_ref (@networks_resolved0)
	{
		## Make sure everyone has dummy parents
		foreach my $node_name (keys %$network_ref)
		{
			#print "Node: $node_name\n";
			my @parents = $$network_ref{$node_name}->parents();
			while(@parents < 2)
			{
				my $dummy1 = add_dummy($network_ref);
				$$network_ref{$node_name}->add_parent($dummy1);
				$$network_ref{$dummy1}->add_child($node_name);
				@parents = $$network_ref{$node_name}->parents();
			}
		}
	}

	###########################################
	# at this point it would be nice if it could construct a pedigree straight up 
	
	
	my @networks_resolved1;
	print "Entering Phase 1. # of possible pedigrees: ". @networks_resolved0 ."\n" if $verbose > 0;
	print $LOG "Entering Phase 1. # of possible pedigrees: ". @networks_resolved0 ."\n" if $verbose > 0;


	while(@networks_resolved0 > 0)
	{
		print "# of networks to resolve: " . @networks_resolved0 . "\n" if($verbose > 1);
		print $LOG "# of networks to resolve: " . @networks_resolved0 . "\n" if($verbose > 1);
		my $network_ref = shift(@networks_resolved0);
		my ($continue, @temp_networks) = phase_1($network_ref);

		if($continue eq 1)
		{
			push(@networks_resolved0,@temp_networks);
		}
		else
		{
			## Test all the networks that came out of Phase 1 to make sure they all fit the phase 2 requirements
			foreach(@temp_networks)
			{
				if(0 ne check_network($_,2)){next;}
				push(@networks_resolved1,$_);
			}
		}
		if(@networks_unresolved > $MAX_NETWORKS_TO_RESOLVE)
		{
			die "> $MAX_NETWORKS_TO_RESOLVE networks to resolve (Phase 1)";
		}
	}
	my $curr_time = time;
	my $timer = $curr_time - $start_time;
	print "phase 1 complete (t=$timer\s)\n" if $verbose > 1;
	print $LOG "phase 1 complete (t=$timer\s)\n" if $verbose > 1;

  if($degree_rel_cutoff < 2)
  {
    return (@networks_resolved1);
  }


	##################################################################################################
	#### Phase 2: resolve all possible relationships that are still not resolved and check existing relationships by cycling through HS/GG/AV (HAG) relationship
	
	my @networks_resolved2;
	print "Entering Phase 2. # of possible pedigrees: ". @networks_resolved1 ."\n" if($verbose > 0);
	print $LOG "Entering Phase 2. # of possible pedigrees: ". @networks_resolved1 ."\n" if($verbose > 0);
	while(@networks_resolved1 > 0)
	{
		print "# of networks to resolve: " . @networks_resolved1 . "\n" if($verbose > 1);
		print $LOG "# of networks to resolve: " . @networks_resolved1 . "\n" if($verbose > 1);
		my $network_ref = shift(@networks_resolved1);

		## Check to see if the phase_2 reconstruction of the network failed due to an error, if so catch/report the error and reject the network.
		my ($continue, @temp_networks);
		eval
		{
			($continue, @temp_networks) = phase_2($network_ref);
			1;
		};
		if($@)
		{
			$continue = 1;
			@temp_networks = ();
			warn "ERROR IN PHASE2 RECONSTRUCTION: $@\n";
		};
		if($continue eq 1)
		{
			push(@networks_resolved1,@temp_networks);
		}
		else
		{
			## Test all the networks that came out of Phase 2 to make sure they all fit the phase 3 requirements
			foreach(@temp_networks)
			{
				eval
				{
					if(0 ne check_network($_,3)){next;}
					push(@networks_resolved2,$_);
					1;
				};
				if($@)
				{
					warn "ERROR IN check_network before Phase 3: $@\n";
				};
			}
		}
		my $curr_time = time;
		my $timer = $curr_time - $start_time;
		if($timer > $MAX_RUNTIME){die "TIMED OUT: $timer\s > $MAX_RUNTIME\s"}
		if(@networks_resolved1 > $MAX_NETWORKS_TO_RESOLVE)
		{
			die "> $MAX_NETWORKS_TO_RESOLVE networks to resolve (Phase 2)";
		}
	}
	my $curr_time = time;
	my $timer = $curr_time - $start_time;
	print "Phase 2 complete (t=$timer\s)\n" if $verbose > 1;
	print $LOG "Phase 2 complete (t=$timer\s)\n" if $verbose > 1;

  if($degree_rel_cutoff < 3)
  {
    return (@networks_resolved2);
  }

	##################################################################################################
	#### Phase 3: combine family networks only conencted by CGH relationship (cousin, great grandparent, half avuncular, and great avuncular)
	print "Entering Phase 3. # of possible pedigrees: ". @networks_resolved2 ."\n" if($verbose > 0);
	print $LOG "Entering Phase 3. # of possible pedigrees: ". @networks_resolved2 ."\n" if($verbose > 0);
	my @networks_resolved3;
	
	while(@networks_resolved2 > 0)
	{
		print "# of networks to resolve: " . @networks_resolved2 . "\n" if($verbose > 1);
		print $LOG "# of networks to resolve: " . @networks_resolved2 . "\n" if($verbose > 1);
		my $network_ref = shift(@networks_resolved2);

		## Check to see if the phase_3 reconstruction of the network failed due to an error, if so catch/report the error and reject the network.
		my ($continue, @temp_networks);
		eval
		{
			($continue, @temp_networks) = phase_3($network_ref);
			1;
		};
		if($@)
		{
			$continue = 1;
			@temp_networks = ();
			warn "ERROR IN PHASE3 RECONSTRUCTION: $@\n";
		};

		if($continue eq 1)
		{
			push(@networks_resolved2,@temp_networks);
		}
		else
		{
			push(@networks_resolved3,@temp_networks);
		}
		my $curr_time = time;
		my $timer = $curr_time - $start_time;

		if($timer > $MAX_RUNTIME){die "TIMED OUT: $timer\s > $MAX_RUNTIME\s"}
		if(@networks_resolved2 > $MAX_NETWORKS_TO_RESOLVE)
		{
			die "More than $MAX_NETWORKS_TO_RESOLVE networks to resolve (Phase 3)";
		}
	}
	
	my $curr_time = time;
	my $timer = $curr_time - $start_time;
	print "Phase 3 complete (t=$timer\s)\n" if $verbose > 1;
	print $LOG "Phase 3 complete (t=$timer\s)\n" if $verbose > 1;
	
	##################################################################################################
	#### TESTING PHASE
	my @networks_resolved;
	foreach my $network_ref (@networks_resolved3)
	{
		eval
		{
			if(0 eq check_network($network_ref,4))
			{
				push(@networks_resolved, $network_ref);
			}
			1;
		};
		if($@)
		{
			warn "ERROR IN check_network after Phase 3: $@\n";
		};
	}
	return (@networks_resolved);
}

##################################
##################################
##################################
## Each section of this phase could be written in a double for loop
sub phase_1
{
	my $network_ref = shift;
	foreach my $node_name (keys %$network_ref)
	{
		$$network_ref{$node_name}->update_unres_first_degree($network_ref);
		while(my $first_degree_rel = $$network_ref{$node_name}->get_unres_first_degree($network_ref))
		{
			my @return_vals = ();
		  my $hag_network_ref = separate_first_degree_and_HAG($network_ref,$node_name,$first_degree_rel);
		   	
			# Split 1st degree from the rest of the likelihood vector and need to push rest back on to be resolved
		  push(@return_vals,$hag_network_ref) if $hag_network_ref ne "";
			
			## Check which first degree relationship they have
			my @possible_rels = @{$$network_ref{$node_name}->possible_relationships($first_degree_rel) };
			if(@possible_rels > 1)
			{
				#print "Splitting PC and FS\n";
				my $fs_network_ref = separate_PC_and_FS($network_ref,$node_name,$first_degree_rel);
				if($fs_network_ref ne "")
				{
					push(@return_vals,$fs_network_ref) if $fs_network_ref ne "";
				}
				@possible_rels = @{$$network_ref{$node_name}->possible_relationships($first_degree_rel) };
			}
			my $relationship = @possible_rels[0];

			## Get parent info for self and rel. If no parents, add dummy parents
			my @parents = $$network_ref{$node_name}->parents();
			my @rel_parents = $$network_ref{$first_degree_rel}->parents();
			if(@parents eq 0)
			{
				my $dummy1 = add_dummy($network_ref);
				my $dummy2 = add_dummy($network_ref);
				$$network_ref{$node_name}->add_parent($dummy1);
				$$network_ref{$node_name}->add_parent($dummy2);
				$$network_ref{$dummy1}->add_child($node_name);
				$$network_ref{$dummy2}->add_child($node_name);
				@parents = $$network_ref{$node_name}->parents();
			}
			if(@rel_parents eq 0)
			{
				my $dummy1 = add_dummy($network_ref);
				my $dummy2 = add_dummy($network_ref);
				$$network_ref{$first_degree_rel}->add_parent($dummy1);
				$$network_ref{$first_degree_rel}->add_parent($dummy2);
				$$network_ref{$dummy1}->add_child($first_degree_rel);
				$$network_ref{$dummy2}->add_child($first_degree_rel);
				@rel_parents = $$network_ref{$first_degree_rel}->parents();
			}

			#### Resolve the relationship
			if($relationship eq "PC")
			{
				## Set rel as parent of node_name
				if($$network_ref{$node_name}->num_real_parents() < 2)
				{
					my $new_network_ref = new_network($network_ref);
					my @temp = add_parent($new_network_ref,$first_degree_rel,$node_name);
					if(@temp > 1){die "More pedigrees than I bargained for\n";}
					$new_network_ref = @temp[0];
					push(@return_vals,$new_network_ref) if $new_network_ref ne "";
				}


				## Set rel as child of self
				if($$network_ref{$first_degree_rel}->num_real_parents() < 2)
				{
					my $new_network_ref = new_network($network_ref);
					my @temp = add_parent($new_network_ref,$node_name,$first_degree_rel);
					if(@temp > 1){die "More pedigrees than I bargained for\n";}
					$new_network_ref = @temp[0];
					push(@return_vals,$new_network_ref) if $new_network_ref ne "";
				}
			}
			elsif($relationship eq "FS")
			{
				## Merge parents of self and rel
				my $new_network_ref = new_network($network_ref);
				my @network_refs = merge_FS_parents($new_network_ref,$node_name,$first_degree_rel);
				my $ctr = 1;
				foreach(@network_refs)
				{
					$ctr++;
				}
				if(@network_refs[0] ne -1)
				{	
					push(@return_vals,@network_refs)
				}
			}
			else
			{
				die "$node_name and $first_degree_rel have a relationship other than PC and FS in phase 1: $relationship\n";
			}

			my @networks_to_return;
			foreach(@return_vals)
			{
				my $failed = check_network($_,1);
				if($failed == 0){push(@networks_to_return,$_)}
			}
			if(@networks_to_return > 0)
			{
				return (1,@networks_to_return);
			}
		}
	}
	return (0,$network_ref);
}

sub do_sib_blocks_match
{
	my $network_ref = shift;
	my $sb1_ref = shift;
	my $sb2_ref = shift;
	my $degree_relatedness = shift;

	foreach my $node1 (@$sb1_ref)
	{
		if($node1 =~ /Dummy/i || $node1 =~ /Missing/i){next}
		my %hag1 = $$network_ref{$node1}->hag();
		my %cgh1 = $$network_ref{$node1}->cgh();
		
		foreach my $node2 (@$sb2_ref)
		{
			if($node2 =~ /Dummy/i || $node2 =~ /Missing/i){next}
			my %hag2 = $$network_ref{$node2}->hag();
			my %cgh2 = $$network_ref{$node2}->cgh();
			if($degree_relatedness eq 2 && !exists $hag1{$node2})
			{
				return 0;
			}
			if($degree_relatedness eq 3 && !exists $cgh1{$node2})
			{
				return 0;
			}
		}
	}
	return 1;
}

sub is_grandparent_to_sib_block
{
	my $network_ref = shift;
	my $node = shift;
	my $sb1_ref = shift;
	my $sb2_ref = shift;
	my $degree_relatedness = shift;

	my %hag = $$network_ref{$node}->hag();
	my %cgh = $$network_ref{$node}->cgh();

	my @hag = keys %hag;
	my @cgh = keys %cgh;
	
	## Check that node is the degree_relatedness to everyone in sb2
	foreach my $node2 (@$sb2_ref)
	{
		if($node2 =~ /Dummy/i || $node2 =~ /Missing/i){next}
		if($degree_relatedness eq 2 && !exists $hag{$node2})
		{
			return 0;
		}
		if($degree_relatedness eq 3 && !exists $cgh{$node2})
		{
			return 0;
		}
	}

	## Check that everyone except $node in sb1 can be degree_relatedness + 1 to everyone in sb2 
	foreach my $node1 (@$sb1_ref)
	{
		if($node1 eq $node){next}
		if($node1 =~ /Dummy/i || $node1 =~ /Missing/i){next}

		my %hag1 = $$network_ref{$node1}->hag();
		my %cgh1 = $$network_ref{$node1}->cgh();
		
		foreach my $node2 (@$sb2_ref)
		{
			if($node2 =~ /Dummy/i || $node2 =~ /Missing/i){next}
			if($degree_relatedness+1 eq 2 && !exists $hag1{$node2})
			{
				return 0;
			}
			if($degree_relatedness+1 eq 3 && !exists $cgh1{$node2})
			{
				return 0;
			}
			if($degree_relatedness+1 eq 4)
			{
				my $possible_relationship_ref = $$network_ref{$node1}->possible_relationships($node2);
				if(!grep($_ eq "DR",@$possible_relationship_ref) && !grep($_ eq "UN",@$possible_relationship_ref))
				{
					return 0;
				}
			}
		}
	}
	return 1;
}


sub phase_3
{
	my $network_ref = shift;
	my $temp = new_network($network_ref);
	foreach my $node_name (keys %$network_ref)
	{
		my %CGH_rels = $$network_ref{$node_name}->cgh();
		$$network_ref{$node_name}->update_unres_cgh($network_ref);
		while(my $cgh_rel = $$network_ref{$node_name}->get_unres_cgh($network_ref))
		{
			my @networks_to_return;
			my @return_vals = ();
		   	my $rest_network_ref = seperate_CGH_and_REST($network_ref,$node_name,$cgh_rel);
		   	
		   	if($rest_network_ref ne "")
		   	{
				push(@networks_to_return,$rest_network_ref);
		   	}
			
			my @parents = $$network_ref{$node_name}->parents();
			my @grandparents0 = $$network_ref{$parents[0]}->parents();
			my @grandparents1 = $$network_ref{$parents[1]}->parents();

			my @rel_parents = $$network_ref{$cgh_rel}->parents();
			my @rel_grandparents0 = $$network_ref{$rel_parents[0]}->parents();
			my @rel_grandparents1 = $$network_ref{$rel_parents[1]}->parents();

			## Leverage the sibling blocks to narrow down the possible pedigrees
			my @sibs = $$network_ref{$node_name}->get_full_sibs($network_ref);
			my @rel_sibs = $$network_ref{$cgh_rel}->get_full_sibs($network_ref);
			my $sib_blocks_all_cghs = 1;
			my $self_is_GGP = 1;
			my $rel_is_GGP = 1;

			push(@sibs,$node_name);
			push(@rel_sibs,$cgh_rel);
			
			if($USE_SIBLING_BLOCKS)
			{
				$sib_blocks_all_cghs = do_sib_blocks_match($network_ref,\@sibs,\@rel_sibs,3);
				$self_is_GGP = is_grandparent_to_sib_block($network_ref,$node_name,\@sibs,\@rel_sibs,3);
				$rel_is_GGP = is_grandparent_to_sib_block($network_ref,$cgh_rel,\@rel_sibs,\@sibs,3);
			}

			#### 3A: First cousin
			## rel and self must share a pair of grandparents to be first cousins
			## Both of them must have 1 dummy parent, because if either have both parents, then this relationship should have already been formed in phase 2
			## Additionally, the pair of grandparents in common must also be dummies. So were are looking for a dummy with two dummy parents 
			## or a dummy without parents
			if($$network_ref{$node_name}->num_real_parents() < 2 && $$network_ref{$cgh_rel}->num_real_parents() < 2 && $sib_blocks_all_cghs)
			{
				my $just_one_parent = expand_on_just_one_parent($network_ref,@parents);
				
				foreach my $parent (@parents)
				{
					if($parent !~ /Dummy/i){next;}
					if($$network_ref{$parent}->num_real_parents() > 0){next;}
					
					my $just_one_rel_parent = expand_on_just_one_parent($network_ref,@rel_parents);
					
					foreach my $rel_parent (@rel_parents)
					{
						if($rel_parent !~ /Dummy/i){next;}
						if($$network_ref{$rel_parent}->num_real_parents() > 0){next;}

						my $new_network_ref = new_network($network_ref);
						my @network_refs = merge_dummy_parents($new_network_ref, $parent,$rel_parent);
						if(@network_refs[0] ne -1)
						{	
							push(@return_vals,@network_refs);
						}
						if($just_one_rel_parent == 1){last;}
					}
					if($just_one_parent == 1){last;}
				}			
			}

			
			#### 3B: Grand-avuncular
			## grand nephew needs to have a dummy parent and that dummy parent must have a dummy parent
			## and that dummy grandparent will then be the sibling of the grand-uncle
			if($sib_blocks_all_cghs)
			{
				my $new_network_ref = new_network($network_ref);
				my @network_refs = merge_grand_avuncular($new_network_ref,$node_name,$cgh_rel);
				if(@network_refs[0] ne -1)
				{	
					push(@return_vals,@network_refs);
				}
				
				my $new_network_ref = new_network($network_ref);
				my @network_refs = merge_grand_avuncular($new_network_ref,$cgh_rel,$node_name);
				if(@network_refs[0] ne -1)
				{	
					push(@return_vals,@network_refs);
				}
			}
			

			#### 3C: Great-grandparental
			## grand child needs to have a dummy parent and that dummy parent must have a dummy parent
			## and that dummy grandparent will then be the child of the great-grand parent
			if($rel_is_GGP)
			{
				my $new_network_ref = new_network($network_ref);
				my @network_refs = merge_great_grandparent($new_network_ref,$node_name,$cgh_rel);
				if(@network_refs[0] ne -1)
				{	
					push(@return_vals,@network_refs);
				}
			}
			if($self_is_GGP)
			{
				my $new_network_ref = new_network($network_ref);
				my @network_refs = merge_great_grandparent($new_network_ref,$cgh_rel,$node_name);
				if(@network_refs[0] ne -1)
				{	
					push(@return_vals,@network_refs);
				}			
			}

			#### 3D: Half-avuncular
			## Half uncle must have at least 1 dummy parent (the relative in common)
			## Half nephew should have a dummy parent who also has a dummy parent
			if($sib_blocks_all_cghs)
			{
				my $new_network_ref = new_network($network_ref);
				my @network_refs = merge_half_avuncular($new_network_ref,$node_name,$cgh_rel);
				if(@network_refs[0] ne -1)
				{	
					push(@return_vals,@network_refs);
				}
				my $new_network_ref = new_network($network_ref);
				my @network_refs = merge_half_avuncular($new_network_ref,$cgh_rel,$node_name);
				if(@network_refs[0] ne -1)
				{	
					push(@return_vals,@network_refs);
				}
			}
			
			## Check networks and return if valid
			foreach(@return_vals)
			{
				if(check_network($_,3) == 0){push(@networks_to_return,$_)}
			}
			if(@networks_to_return > 0)
			{
				return (1,@networks_to_return);
			}
			else
			{
				## Maybe I need to delete the relationship because I couldn't build on it 
				##and if I leave it in then I will get into an infinite loop
				#$$network_ref{$node_name}->delete_cgh($cgh_rel);
				#
				#
				# I guess not since I never coded this and it seems to work just fine
			}
		}
	}
	return (0,$network_ref);
}

sub phase_2
{
	my $network_ref = shift;
	foreach my $node_name (keys %$network_ref)
	{
		$$network_ref{$node_name}->update_unres_hag($network_ref);
		while(my $hag_rel = $$network_ref{$node_name}->get_unres_hag($network_ref))
		{
			my @networks_to_return;
		   	my @return_vals = (); ## These are the networks that need to be tested
		   	my $cgh_network_ref = seperate_HAG_and_CGH($network_ref,$node_name,$hag_rel);
		   	
		   	if($cgh_network_ref ne "")
		   	{
		   		push(@networks_to_return,$cgh_network_ref);
		   	}
		   	
			my @parents = $$network_ref{$node_name}->parents();
			my @rel_parents = $$network_ref{$hag_rel}->parents();

			my $just_one_parent = expand_on_just_one_parent($network_ref,@parents);
			my $just_one_rel_parent = expand_on_just_one_parent($network_ref,@rel_parents);
	
			my @sibs = $$network_ref{$node_name}->get_full_sibs($network_ref);
			my @rel_sibs = $$network_ref{$hag_rel}->get_full_sibs($network_ref);
			push(@sibs,$node_name);
			push(@rel_sibs,$hag_rel);

			## Leverage the sibling blocks to narrow down the possible pedigrees
			my $sib_blocks_all_hags = 1;
			my $self_is_GP = 1;
			my $rel_is_GP = 1;

			if($USE_SIBLING_BLOCKS)
			{
				$sib_blocks_all_hags = do_sib_blocks_match($network_ref,\@sibs,\@rel_sibs,2);
				$self_is_GP = is_grandparent_to_sib_block($network_ref,$node_name,\@sibs,\@rel_sibs,2);
				$rel_is_GP = is_grandparent_to_sib_block($network_ref,$hag_rel,\@rel_sibs,\@sibs,2);
			}

			### A: Half sib relationship #######################################################
			## Can be half sibs if neither have both parents
			if($$network_ref{$node_name}->num_real_parents() < 2 && $$network_ref{$hag_rel}->num_real_parents() < 2 && $sib_blocks_all_hags)
			{
				
				## Parent in common is P1 (parent[0]) and RP1 (rel_parent[0])
				if(@parents[0] =~ /Dummy/i && @rel_parents[0] =~ /Dummy/i)
				{
					my $network_ref11 = new_network($network_ref);
					my @networks = merge_dummies($network_ref11,@parents[0],@rel_parents[0]);
					if(@networks[0] != -1)
					{
						foreach(@networks)
						{
							if(check_network($_,2) == 0){push(@return_vals,$_)}
						}
					}
				}
				
				## Parent in common is P1 (parent[0]) and RP2 (rel_parent[1])
				if(@parents[0] =~ /Dummy/i && @rel_parents[1] =~ /Dummy/i && $just_one_rel_parent == 0)
				{
					my $network_ref12 = new_network($network_ref);
					my @networks = merge_dummies($network_ref12,@parents[0],@rel_parents[1]);
					if(@networks[0] != -1)
					{
						foreach(@networks)
						{
							if(check_network($_,2) == 0){push(@return_vals,$_)}
						}
					}
				}
				
				## Parent in common is P2 (parent[1]) and RP1 (rel_parent[0])
				if(@parents[1] =~ /Dummy/i && @rel_parents[0] =~ /Dummy/i && $just_one_parent == 0)
				{
					my @parent1_parents = $$network_ref{@parents[1]}->parents();
					my @rel_parent0_parents = $$network_ref{@rel_parents[0]}->parents();
					my $network_ref21 = new_network($network_ref);
					my @networks = merge_dummies($network_ref21,@parents[1],@rel_parents[0]);
					if(@networks[0] != -1)
					{
						foreach(@networks)
						{
							if(check_network($_,2) == 0){push(@return_vals,$_)}
						}
					}
				}
								
				## Parent in common is P2 (parent[1]) and RP2 (rel_parent[1])
				if(@parents[1] =~ /Dummy/i && @rel_parents[1] =~ /Dummy/i && $just_one_parent == 0 && $just_one_rel_parent == 0)
				{
					my $network_ref22 = new_network($network_ref);
					my @networks = merge_dummies($network_ref22,@parents[1],@rel_parents[1]);
					if(@networks[0] != -1)
					{
						foreach(@networks)
						{
							if(check_network($_,2) == 0){push(@return_vals,$_)}
						}
					}
				}
			}
			
			
			#print "A done\n";
			### B: avuncular relationships #######################################################
			## Dummy is selfs full sibling and is dummy parent1 of rel_hag = neice or nephew; 
			## If dummy parent1 has a real parent, then it must be shared between it and the uncle
			## because otherwise the PC relationship would have already been made between the real parent and the uncle.
			## Cannot have 2 real parents or this relationship would already have existed
			if(@rel_parents[0] =~ /Dummy/i && $$network_ref{@rel_parents[0]}->num_real_parents() ne 2 && $sib_blocks_all_hags)
			{
				my $new_network_ref = new_network($network_ref);
				my @new_network_ref = merge_avuncular($new_network_ref,$node_name,@rel_parents[0]);
				if(@new_network_ref[0] ne -1)
				{	
					push(@return_vals,@new_network_ref);
				}
			}
			
			## Dummy is selfs full sibling and is dummy parent2 of rel_hag = neice or nephew;
			if(@rel_parents[1] =~ /Dummy/i && $$network_ref{@rel_parents[1]}->num_real_parents() ne 2 && $just_one_rel_parent == 0 && $sib_blocks_all_hags)
			{
				my @rel_parent1_parents = $$network_ref{@rel_parents[1]}->parents();
				my $new_network_ref = new_network($network_ref);
				my @new_network_ref = merge_avuncular($new_network_ref,$node_name,@rel_parents[1]);
				if(@new_network_ref[0] ne -1)
				{	
					push(@return_vals,@new_network_ref);
				}
			}
			
			## Dummy is selfs parent1, and parent1 is full sibling with hag_rel = uncle/aunt
			if(@parents[0] =~ /Dummy/i && $$network_ref{@parents[0]}->num_real_parents() ne 2 && $sib_blocks_all_hags)
			{
				my $new_network_ref = new_network($network_ref);
				my @new_network_ref = merge_avuncular($new_network_ref,$hag_rel,@parents[0]);
				if(@new_network_ref[0] ne -1)
				{	
					push(@return_vals,@new_network_ref);
				}
			}
			
			## Dummy is selfs parent2, and parent2 is full sibling with rel_hag = uncle/aunt
			if(@parents[1] =~ /Dummy/i && $$network_ref{@parents[1]}->num_real_parents() ne 2 && $just_one_parent == 0 && $sib_blocks_all_hags)
			{
				my $new_network_ref = new_network($network_ref);
				my @new_network_ref = merge_avuncular($new_network_ref,$hag_rel,@parents[1]);
				if(@new_network_ref[0] ne -1)
				{	
					push(@return_vals,@new_network_ref);
				}
			}

			### C: Grandparent/grandchild #######################################################
			## Dummy is self parent1 and dummy child of rel_hag; rel_hag = grandparent
			## Dummy parent/child cannot have two real parents because if it does, then it can't 
			## add a third to make the relationship grandparental
			if(@parents[0] =~ /Dummy/i && $$network_ref{@parents[0]}->num_real_parents() ne 2 && $rel_is_GP)
			{
				my @rel_parent1_parents = $$network_ref{@rel_parents[1]}->parents();
				my $new_network_ref = new_network($network_ref);
				my @new_network_ref = merge_grandparental($new_network_ref,$hag_rel,@parents[0]);
				if(@new_network_ref[0] ne -1)
				{	
					push(@return_vals,@new_network_ref);
				}
			}
			
			## Dummy is self parent2 and dummy child of rel_hag; rel_hag = grandparent
			if(@parents[1] =~ /Dummy/i && $$network_ref{@parents[1]}->num_real_parents() ne 2 && $just_one_parent == 0 && $rel_is_GP)
			{
				my $new_network_ref = new_network($network_ref);
				my @new_network_ref = merge_grandparental($new_network_ref,$hag_rel,@parents[1]);
				if(@new_network_ref[0] ne -1)
				{	
					push(@return_vals,@new_network_ref);
				}
			}
			
			## Dummy is child and dummy parent1 of rel_hag; rel_hag = grand child
			if(@rel_parents[0] =~ /Dummy/i && $$network_ref{@rel_parents[0]}->num_real_parents() ne 2 && $self_is_GP)
			{
				my $new_network_ref = new_network($network_ref);
				my @new_network_ref = merge_grandparental($new_network_ref,$node_name,@rel_parents[0]);
				if(@new_network_ref[0] ne -1)
				{	
					push(@return_vals,@new_network_ref);
				}

			}
			
			## Dummy is child and is dummy parent2 of rel_hag; rel_hag = grand child
			if(@rel_parents[1] =~ /Dummy/i && $$network_ref{@rel_parents[1]}->num_real_parents() ne 2 && $just_one_rel_parent == 0 && $self_is_GP)
			{
				my $new_network_ref = new_network($network_ref);
				my @new_network_ref = merge_grandparental($new_network_ref,$node_name,@rel_parents[1]);
				if(@new_network_ref[0] ne -1)
				{	
					push(@return_vals,@new_network_ref);
				}			
			}
			
			foreach(@return_vals)
			{
				if(check_network($_,2) == 0){push(@networks_to_return,$_)}
			}
			if(@networks_to_return > 0)
			{
				return (1,@networks_to_return);
			}
		}
	}
	return (0,$network_ref);
}

## Find B with PC rel A-B where A has both parents
sub phase_1A
{
	my $network_ref = shift;
	foreach my $node_name (keys %$network_ref)
	{
		## B = $node_name
		my @parents = $$network_ref{$node_name}->parents();
		## B should not have any parents yet
		if(@parents == 0)
		{
			my %PC = $$network_ref{$node_name}->pc();
			## Find the PC_rel with two parents
			foreach my $PC_rel (keys %PC)
			{
				my @PC_rel_parents = $$network_ref{$PC_rel}->parents();
				if($$network_ref{$PC_rel}->num_real_parents() == 2 && !grep {$_ eq $node_name} @PC_rel_parents)
				{
					## B is $node_name and A is PC_rel
					## Found B's grandparents
					## Set B as child of A and A as parent of B and make dummy parent for B
					$$network_ref{$node_name}->add_parent($PC_rel);
					my $temp_dummy = add_dummy($network_ref);
					$$network_ref{$node_name}->add_parent($temp_dummy);
					$$network_ref{$temp_dummy}->add_child($node_name);
					$$network_ref{$PC_rel}->add_child($node_name);
				}
			}			
		}
	}
	return (0,$network_ref);
}

## Find node A that has a PC rel to B and B has children. 
## Check that A is HAG to all children of B.
## Make network with A is aprent of B and 
## Make network with B is parent of A
sub phase_1B
{
	my $network_ref = shift;
	my $entered_loop = 0;
	foreach my $node_name (keys %$network_ref)
	{
		my @parents = $$network_ref{$node_name}->parents();
		## B should not have any more than one real parents
		if($$network_ref{$node_name}->num_real_parents() == 2 && !grep ($_ =~ /Dummy/, @parents))
		{
			next;
		}
		my @children = $$network_ref{$node_name}->children();

		## Look for unresolved pc rel
		my %PC_rels = $$network_ref{$node_name}->pc();
		my %HAG_rels = $$network_ref{$node_name}->hag();

		OUTER: foreach my $PC_rel (keys %PC_rels)
		{
			## we are not interested in the possible one parent or known children
			if(grep {$_ eq $PC_rel} @parents){next;}
			if(grep {$_ eq $PC_rel} @children){next;}
			
			## Looking for a PC rel with known children
			my @PC_rel_children = $$network_ref{$PC_rel}->children();
			if(@PC_rel_children eq 0){next;}
			
			## Check that $node_name has a HAG_rel with each of the PC_rel_children
			foreach my $PC_rel_child (@PC_rel_children)
			{
				if(!exists $HAG_rels{$PC_rel_child}){
					next OUTER;
				}
			}
			
			## This is necessary in the case that there is a network that has 
			## this unresolved PC_rel but neither combination fits all relationships
			$entered_loop = 1;
			## FOUND CASE 9/16 = a trio with one sample that could be HS or Grandparent to the child
			my $new_network9 = new_network($network_ref);
			my @temp = add_parent($new_network9,$PC_rel,$node_name);
			if(@temp > 1){die "More pedigrees than I bargained for\n";}
			$new_network9 = @temp[0];

			my $new_network16 = new_network($network_ref);
			my @temp = add_parent($new_network16,$node_name,$PC_rel);
			if(@temp > 1){die "More pedigrees than I bargained for\n";}
			$new_network16 = @temp[0];
			my $valid9 = check_network($new_network9,2);
			my $valid16 = check_network($new_network16,2);
			
			my @return_vals = ();
			if($valid9 eq 0){push(@return_vals,$new_network9)}
			if($valid16 eq 0){push(@return_vals,$new_network16)}
			if(@return_vals > 0)
			{
				return (1,@return_vals);
			}
		}
	}
	if($entered_loop == 1)
	{
		return (0);
	}

	return (0,$network_ref);
}

#### looking for unresolved A-B-C PC rels where A-C is HAG; 
#### A = node_name; B = PC_rel; C = PC_rel_PC_rel
sub phase_1C
{
	my $network_ref = shift;
	foreach my $node_name (keys %$network_ref)
	{
		my $num_real_parents = $$network_ref{$node_name}->num_real_parents();
		## B should not have any more than one real parents
		if($num_real_parents > 1){next;}

		## Look for unresolved pc rel
		my %PC_rels = $$network_ref{$node_name}->pc();
		my %HAG_rels = $$network_ref{$node_name}->hag();
		foreach my $PC_rel (keys %PC_rels)
		{
			my @PC_rel_parents = $$network_ref{$PC_rel}->parents();
			
			## PC_rel should not have parents, because if so, this would be resolved by a different subroutine
			if(@PC_rel_parents >= 1){next;}
			
			## Found B = PC_rel
			## Find  C = PC_rel_PC_rel 
			my %PC_rel_PC_rels = $$network_ref{$PC_rel}->pc();
			
			## Look for unresolved PC rel of PC_rel that has HAG to node_name
			foreach my $PC_rel_PC_rel (keys %PC_rel_PC_rels)
			{
				## We don't want to go back to node_name (A)
				if ($PC_rel_PC_rel eq $node_name){next;}
				
				my @PC_rel_PC_rel_parents = $$network_ref{$PC_rel_PC_rel}->parents();
				## PC_rel_PC_rel should not have parents, because if so, this would be resolved by a different subroutine
				if(@PC_rel_PC_rel_parents >= 1){next;}
				
				## Check if PC_rel_PC_rel (C) has HAG rel to node_name (A)
				if(!exists $HAG_rels{$PC_rel_PC_rel}){next;}
				
				## Found C = PC_rel_PC_rel and CASE 12/18x2 = similar to the 9/16, but you are missing the parent not in common
				my $new_network12 = new_network($network_ref);
				my @temp = add_parent($new_network12,$PC_rel,$node_name);
				if(@temp > 1){die "More pedigrees than I bargained for\n";}
				$new_network12 = @temp[0];
				$$new_network12{$PC_rel}->add_child($node_name);
				my @temp = add_parent($new_network12,$PC_rel,$PC_rel_PC_rel);
				if(@temp > 1){die "More pedigrees than I bargained for\n";}
				$new_network12 = @temp[0];
				$$new_network12{$PC_rel}->add_child($PC_rel_PC_rel);
				my @parents = $$new_network12{$node_name}->parents();
				if(@parents < 2)
				{
					my $temp_dummy12a = add_dummy($new_network12);
					my @temp = add_parent($new_network12,$temp_dummy12a,$node_name);
					if(@temp > 1){die "More pedigrees than I bargained for\n";}
					$new_network12 = @temp[0];
					$$new_network12{$temp_dummy12a}->add_child($node_name);
				
				}
				my @parents = $$new_network12{$PC_rel_PC_rel}->parents();
				if(@parents < 2)
				{
					my $temp_dummy12b = add_dummy($new_network12);
					my @temp = add_parent($new_network12,$temp_dummy12b,$PC_rel_PC_rel);
					if(@temp > 1){die "More pedigrees than I bargained for\n";}
					$new_network12 = @temp[0];
					$$new_network12{$temp_dummy12b}->add_child($PC_rel_PC_rel);
				}

				my $new_network18a = new_network($network_ref);
				my @temp = add_parent($new_network18a,$node_name,$PC_rel);
				if(@temp > 1){die "More pedigrees than I bargained for\n";}
				$new_network18a = @temp[0];
				$$new_network18a{$node_name}->add_child($PC_rel);
				my @temp = add_parent($new_network18a,$PC_rel,$PC_rel_PC_rel);
				if(@temp > 1){die "More pedigrees than I bargained for\n";}
				$new_network18a = @temp[0];
				$$new_network18a{$PC_rel}->add_child($PC_rel_PC_rel);
				my @parents = $$new_network12{$PC_rel}->parents();
				if(@parents < 2)
				{
					my $temp_dummy18a1 = add_dummy($new_network18a);
					my @temp = add_parent($new_network18a,$temp_dummy18a1,$PC_rel);
					if(@temp > 1){die "More pedigrees than I bargained for\n";}
					$new_network18a = @temp[0];
					$$new_network18a{$temp_dummy18a1}->add_child($PC_rel);
				}
				my @parents = $$new_network12{$PC_rel_PC_rel}->parents();
				if(@parents < 2)
				{
					my $temp_dummy18a2 = add_dummy($new_network18a);
					my @temp = add_parent($new_network18a,$temp_dummy18a2,$PC_rel_PC_rel);
					if(@temp > 1){die "More pedigrees than I bargained for\n";}
					$new_network18a = @temp[0];
					$$new_network18a{$temp_dummy18a2}->add_child($PC_rel_PC_rel);
				}

				my $new_network18b = new_network($network_ref);
				$$new_network18b{$PC_rel}->add_child($node_name);
				my @temp = add_parent($new_network18b,$PC_rel,$node_name);
				if(@temp > 1){die "More pedigrees than I bargained for\n";}
				$new_network18b = @temp[0];
				$$new_network18b{$PC_rel_PC_rel}->add_child($PC_rel);
				my @temp = add_parent($new_network18b,$PC_rel_PC_rel,$PC_rel);
				if(@temp > 1){die "More pedigrees than I bargained for\n";}
				$new_network18b = @temp[0];
				my @parents = $$new_network12{$PC_rel}->parents();
				if(@parents < 2)
				{
					my $temp_dummy18b1 = add_dummy($new_network18b);
					my @temp = add_parent($new_network18b,$temp_dummy18b1,$PC_rel);
					if(@temp > 1){die "More pedigrees than I bargained for\n";}
					$new_network18b = @temp[0];
					$$new_network18b{$temp_dummy18b1}->add_child($PC_rel);
				}
				my @parents = $$new_network12{$node_name}->parents();
				if(@parents < 2)
				{
					my $temp_dummy18b2 = add_dummy($new_network18b);
					my @temp = add_parent($new_network18b,$temp_dummy18b2,$node_name);
					if(@temp > 1){die "More pedigrees than I bargained for\n";}
					$new_network18b = @temp[0];
					$$new_network18b{$temp_dummy18b2}->add_child($node_name);
				}
				my $valid12 = check_network($new_network12,2);
				my $valid18a = check_network($new_network18a,2);
				my $valid18b = check_network($new_network18b,2);
				my @return_vals = ();
				if($valid12 eq 0){push(@return_vals,$new_network12)}
				if($valid18a eq 0){push(@return_vals,$new_network18a)}
				if($valid18b eq 0){push(@return_vals,$new_network18b)}
				if(@return_vals > 0)
				{
					return (1,@return_vals);
				}
			}
		}
	}
	return (0,$network_ref);
}


### Subroutine fills in all the missing dummy parents in the pedigree and resolved the directionality of PC relationships
sub phase_1D
{
	my $network_ref = shift;
	foreach my $node_name (keys %$network_ref)
	{
		if($node_name =~ /Dummy/i){next;}
		
		my @parents = $$network_ref{$node_name}->parents();
		my @children = $$network_ref{$node_name}->children();
		my %PC = $$network_ref{$node_name}->pc();
		
		## B should not have any parents yet
		if(@parents == 2){next;}
		if(@parents == 1) ## directionality already established
		{
			## Parent should not have any more than one child
			my @P_children = $$network_ref{@parents[0]}->children();
			if(@P_children > 1){die "ERROR!!! @parents has more than one child: @P_children\n";}

			my $dummy_parent_name = add_dummy($network_ref);
			$$network_ref{$node_name}->add_parent($dummy_parent_name);
			$$network_ref{$dummy_parent_name}->add_child($node_name);
		}
		
		if(@parents == 0 && keys  %PC > 0)
		{
			
			## is a founder
			if(@children > 0)
			{
				## Add dummy parents
				while (@parents < 2)
				{
					my $dummy_parent_name = add_dummy($network_ref);
					$$network_ref{$node_name}->add_parent($dummy_parent_name);
					$$network_ref{$dummy_parent_name}->add_child($node_name);
					@parents = $$network_ref{$node_name}->parents();
				}
			}
			else  ## unestablished directionality
			{
				my @PC_rels = keys %PC;
				if(@PC_rels > 1){print $LOG "ERROR!!! $node_name has more than one PC_rels: @PC_rels\n";}	## longer PC strings would have been resolved in Phase1C and then Phase1A
				if(@PC_rels > 1){print "ERROR!!! $node_name has more than one PC_rels: @PC_rels\n";return -1;}	## longer PC strings would have been resolved in Phase1C and then Phase1A
				

				## set node as child
				my $new_network2 = new_network($network_ref);
				my @networks = add_parent($new_network2,$PC_rels[0],$node_name);
				if(@networks > 1){die "IN 1D; new_network2 SHOULD NOT HAVE RETURNED MORE THAN ONE NETWORK\n";}
				$new_network2 = shift(@networks);
				
				## set node as parent of PC_rel if PC_rel doesn't already have two parents
				if($$network_ref{$PC_rels[0]}->num_real_parents() >= 2)
				{
					return(1,$new_network2);
				}
				my $new_network1 = new_network($network_ref);
				my @networks = add_parent($new_network1,$node_name,$PC_rels[0]);
				if(@networks > 1){die "IN 1D; SHOULD NOT HAVE RETURNED MORE THAN ONE NETWORK\n";}
				$new_network1 = shift(@networks);
				
				return (1,$new_network1,$new_network2);
			}
		}
		else
		{
			while (@parents < 2)
			{
				my $dummy_parent_name = add_dummy($network_ref);
				my @networks = add_parent($network_ref,$dummy_parent_name,$node_name);
				if(@networks > 1){die "IN 1D; new_network2 SHOULD NOT HAVE RETURNED MORE THAN ONE NETWORK\n";}
				$network_ref = shift(@networks);
				@parents = $$network_ref{$node_name}->parents();
			}
		}
	}
	return (0,$network_ref);
}

sub resolve_FS_dummy_parents
{
	my $network_ref = shift;
	my $node_name = shift;
	my @parents = $$network_ref{$node_name}->parents();
	my %FS = $$network_ref{$node_name}->fs();
	my @sib_names = keys %FS;

	## Check that parents match between all sibs
	foreach my $sib_name (@sib_names)
	{
		my @sib_parents = $$network_ref{$sib_name}->parents();
		my $match = do_arrays_match(\@parents,\@sib_parents);
		if($match eq 0)
		{
			print "$node_name parents (@parents) do not match sibling $sib_name parents (@sib_parents)\n" if $verbose > 1;
			print $LOG "$node_name parents (@parents) do not match sibling $sib_name parents (@sib_parents)\n" if $verbose > 1;
			return -1;
		}
	}
	
	## Add dummy parents
	while (@parents < 2)
	{
		my $dummy_parent_name = add_dummy($network_ref);

		$$network_ref{$node_name}->add_parent($dummy_parent_name);
		$$network_ref{$dummy_parent_name}->add_child($node_name);
		foreach my $sib_name (@sib_names)
		{
			$$network_ref{$sib_name}->add_parent($dummy_parent_name);
			$$network_ref{$dummy_parent_name}->add_child($sib_name);
		}
		@parents = $$network_ref{$node_name}->parents();
	}

	## Check that parents match between all sibs
	foreach my $sib_name (@sib_names)
	{
		my @sib_parents = $$network_ref{$sib_name}->parents();
		my $match = do_arrays_match(\@parents,\@sib_parents);
		if($match eq 0)
		{
			print "$node_name parents (@parents) do not match sibling $sib_name parents (@sib_parents)\n" if $verbose > 1;
			print $LOG "$node_name parents (@parents) do not match sibling $sib_name parents (@sib_parents)\n" if $verbose > 1;
			return -1;
		}
	}
}

## Return -1 if the pass in network_ref failed; return 1 if it is resolved; Also return any new networks generated from spliting FS and HAG relationships
sub resolve_PC_trio
{
	my $network_ref = shift;
  my @return_networks = ();
  foreach my $node_name (keys %$network_ref)
  {
    my %PC = $$network_ref{$node_name}->pc();
    next if keys %PC < 2;
	
    my @PC_rels = keys %PC;
    my $relationships_ref = $$network_ref{$node_name}->possible_relationships(@PC_rels[0]);

    ## check if PC_rels are FS_rels, if so, $node_name is parent, and they are children;
    # else if they could be unrelated, then $node_name is child and they are parents
    # else if they are not unrelated or FS, I need TO DO THIS##### 
    # but are not FS something more complicated that I will come to later
    for(my $i = 0; $i < @PC_rels - 1; $i++)
    {
      for(my $j = $i+1; $j < @PC_rels; $j++)
      {
        my $rel1_name = @PC_rels[$i];
        my $rel2_name = @PC_rels[$j];
        my @num_real_parents = $$network_ref{$node_name}->num_real_parents();
    
        my $relationships_ref = $$network_ref{$rel1_name}->possible_relationships($rel2_name);
        #if(grep($_ eq "", @{$relationships_ref}) || (!grep($_ eq "FS", @{$relationships_ref})) || @$relationships_ref eq 0) ##check if parents of $node_name;
        if(grep($_ eq "", @{$relationships_ref}) || (!grep($_ eq "FS", @{$relationships_ref}) && !grep($_ eq "HAG", @{$relationships_ref})) || @$relationships_ref eq 0) ##check if parents of $node_name;
        {		
          #print "$rel1_name : $rel2_name = !FS && !HAG; Parents of $node_name\n";
          if(@num_real_parents > 1)
          {	
            my @parents = $$network_ref{$node_name}->parents();
            print "$node_name already has two parents (@parents) can't add these two\n" if $verbose > 1;
            print $LOG "$node_name already has two parents (@parents) can't add these two\n" if $verbose > 1;
            return (-1,@return_networks);
          }
          $$network_ref{$node_name}->parents($rel1_name,$rel2_name);
          $$network_ref{$rel1_name}->add_child($node_name);
          $$network_ref{$rel2_name}->add_child($node_name);
        }
        elsif(grep($_ eq "FS", @{$relationships_ref})) ##check if children of $node_name
        {
          #
          # Split 1st degree from the rest of the likelihood vector and need to push rest back on to be resolved
          my $hag_network_ref = separate_first_degree_and_HAG($network_ref,$rel1_name,$rel2_name);
          push(@return_networks,$hag_network_ref) if $hag_network_ref ne "";

          #print "$rel1_name : $rel2_name = FS; PARENT = $node_name\n";
          #print "HAG net not empty\n" if $hag_network_ref ne "";
          if($$network_ref{$rel1_name}->num_real_parents() > 1)
          {
            my @parents = $$network_ref{$rel1_name}->parents();
            if(!grep ($_ eq $node_name, @parents ))
            {
              print "$rel1_name already has two parents (@parents) can't add $node_name\n" if $verbose > 1;
              print $LOG "$rel1_name already has two parents (@parents) can't add $node_name\n" if $verbose > 1;
              return (-1,@return_networks);
            }		
          }
          if($$network_ref{$rel2_name}->num_real_parents() > 1)
          {
            my @parents = $$network_ref{$rel2_name}->parents();
            if(!grep ($_ eq $node_name, @parents ))
            {
              print "$rel2_name already has two parents (@parents) can't add $node_name\n" if $verbose > 1;
              print $LOG "$rel2_name already has two parents (@parents) can't add $node_name\n" if $verbose > 1;
              return (-1,@return_networks);
            }
          }
          $$network_ref{$node_name}->add_child($rel1_name);
          $$network_ref{$node_name}->add_child($rel2_name);
          $$network_ref{$rel1_name}->add_parent($node_name);
          
          $$network_ref{$rel2_name}->add_parent($node_name);
        }
        elsif(grep($_ eq "HAG", @{$relationships_ref}))
        {
          if($verbose > 2)
          {
            print "$rel1_name : $rel2_name = HAG; do nothing\n" if $verbose > 1; ## for now. This precents PRIMUS from reconsructing pedigrees where the parents are 2nd degree relatives (HAGs)
            print $LOG "$rel1_name : $rel2_name = HAG; do nothing\n" if $verbose > 1; ## for now. This precents PRIMUS from reconsructing pedigrees where the parents are 2nd degree relatives (HAGs)
          }
           ## This means that this could be a parent or a child. 
        }
        else
        {
          ## Ignore these cases for now...
        }
      }
    }
  }
  return (1,@return_networks);
}

sub load_network
{	
	my $data_ref = shift; ## Relationship hash reference
	my $mito_ref = shift; 
	my $y_ref = shift;
	my $gender_ref = shift;
	my %data = %{$data_ref};

	my %network = ();
	my $network_ref = \%network;
	
  #### Merge MZ twins into the same node by merging their names
  ## Get the MZ twin networks
	my %mz_twin_network_lookup;
  my $mz_network_counter = 1;
  foreach my $self (sort keys %data)
	{
		foreach my $rel (sort keys %{ $data{$self} })
		{
			my @likelihoods = @{ $data{$self}{$rel} };
		  next if(@likelihoods[6] == 0);
      ## Check to see if these two samples are in the same mz_twin_network
      if(exists $mz_twin_network_lookup{$self} && exists $mz_twin_network_lookup{$rel})
      {
        ## If they are NOT in the same network, change everyone in rel's network to self's network; else don't do anything
        my $selfs_network = $mz_twin_network_lookup{$self};
        my $rels_network = $mz_twin_network_lookup{$rel};
        next if $selfs_network == $rels_network;
        foreach my $twin (keys %mz_twin_network_lookup)
        {
          $mz_twin_network_lookup{$twin} = $selfs_network if $mz_twin_network_lookup{$twin} == $rels_network;
        }
      }
      elsif (exists $mz_twin_network_lookup{$self} && ! exists $mz_twin_network_lookup{$rel})
      {
        ## put rel in self's network
        my $selfs_network = $mz_twin_network_lookup{$self};
        $mz_twin_network_lookup{$rel} = $selfs_network;
      }
      elsif (! exists $mz_twin_network_lookup{$self} && exists $mz_twin_network_lookup{$rel})
      {
        ## put self in rel's network
        my $rels_network = $mz_twin_network_lookup{$rel};
        $mz_twin_network_lookup{$self} = $rels_network;
      }
      elsif (! exists $mz_twin_network_lookup{$self} && ! exists $mz_twin_network_lookup{$rel})
      {
        ## add both to network mz_network_counter and then increment mz_network_counter
        $mz_twin_network_lookup{$self} = $mz_network_counter;
        $mz_twin_network_lookup{$rel} = $mz_network_counter;
        $mz_network_counter++;
      }
    }
  }
  ## Should probably Check that each mz_twin_network is a clique
  
  ## Collapse each MZ twin network into a single "person" or node
  for(my $i = 1; $i < $mz_network_counter; $i++)
  {
    print "Collapsing network $i\n";
    # Get the collapsed name
    my @names_to_collapse;
    foreach my $twin (keys %mz_twin_network_lookup)
    {
      push(@names_to_collapse,$twin) if $mz_twin_network_lookup{$twin} == $i;
    }
    my $new_name = join(':',@names_to_collapse);

    # for each collapsed sample, rename it in dataset and remove the old name; this will result in a reduction in size of the data hash table as relationship between collapsed samples and their relatives overwrite each other
    foreach my $old_name (@names_to_collapse)
    {
      # If old name is first name in pair
      $data{$new_name} = $data{$old_name};
      delete $data{$old_name};
      
      # if old name is second name in pair
      foreach my $rel (sort keys %data)
      {
        next if ! exists $data{$rel}{$old_name};
        $data{$rel}{$new_name} = $data{$rel}{$old_name};
        delete $data{$rel}{$old_name};
      }
    }
    delete $data{$new_name}{$new_name};
  }

  #print "AFTER\n";
  foreach my $self (sort keys %data)
	{
		foreach my $rel (sort keys %{ $data{$self} })
		{
			my @likelihoods = @{ $data{$self}{$rel} };
      #print "$self <-> $rel = @likelihoods\n";
    }
  }

  ## Now make the network for PRIMUS
	foreach my $self (sort keys %data)
	{
		if(!exists $$network_ref{$self})
		{
			my $sex = $$gender_ref{$self};
			my $node = new PRIMUS::node_v7($self,$sex);
			$$network_ref{$self} = $node;
		}

		foreach my $rel (sort keys %{ $data{$self} })
		{
			my @likelihoods = @{ $data{$self}{$rel} };
			
			#print "$self <-> $rel\n";

			## Add rel to network
			if(!exists $$network_ref{$rel})
			{
				my $sex = $$gender_ref{$rel};
				my $node = new PRIMUS::node_v7($rel,$sex);
				$$network_ref{$rel} = $node;
			}
			my $mito_match = "";
			my $y_match = "";
			$mito_match = $$mito_ref{$self}{$rel} if exists $$mito_ref{$self}{$rel};
			$y_match = $$y_ref{$self}{$rel} if exists $$y_ref{$self}{$rel};
			#print "$self <-> $rel = $match\n";
			$$network_ref{$self}->add_relative($rel,@likelihoods);
			$$network_ref{$rel}->add_relative($self,@likelihoods);
			$$network_ref{$self}->add_mito_match($rel,$mito_match);
			$$network_ref{$rel}->add_mito_match($self,$mito_match);
			$$network_ref{$self}->add_y_match($rel,$y_match);
			$$network_ref{$rel}->add_y_match($self,$y_match);
		}
	}
	
	## Set the relationship for each network
	foreach my $node_name (keys %$network_ref)
	{			
		my $node = $$network_ref{$node_name};
		my $node_name = $node->name();
		my %relatives = $node->relatives();
		
		## Set relationships
		my %PC;
		my %FS;
		my %HAG;
		my %CGH;
				
		## Get all possible relationships and add them to each of the lists
		foreach my $relative (keys %relatives)
		{
			my $possible_relationships_ref = $node->possible_relationships($relative);
			
			foreach my $relationship (@$possible_relationships_ref)
			{
				if ($relationship =~ /^PC$/i)
				{
					$PC{$relative} = $relatives{$relative};
				}
				elsif ($relationship =~ /^FS$/i)
				{
					$FS{$relative} = $relatives{$relative};
				}
				elsif ($relationship =~ /^HAG$/i)
				{
					$HAG{$relative} = $relatives{$relative};
				}
				elsif ($relationship =~ /^CGH$/i)
				{
					$CGH{$relative} = $relatives{$relative};
				}
			}
		}
		$node->pc(%PC);
		$node->fs(%FS);
		$node->hag(%HAG);
		$node->cgh(%CGH);
	}
	
	foreach(keys %$network_ref)
	{
		my %relatives = $$network_ref{$_}->relatives();
		my @rels = keys %relatives;
	}
	
	return $network_ref;
}

sub add_dummy
{
	my $network_ref = shift;
	my $dummy_name = "Dummy$dummy_ctr";
	$dummy_ctr++;
	
	## Add dummy_parent to network
	my $dummy = new PRIMUS::node_v7($dummy_name);
	$$network_ref{$dummy_name} = $dummy;
	return $dummy_name;
}

sub merge_dummies
{
	my @networks = merge_nodes(@_);
	return @networks;
}

sub merge_nodes
{
	my $old_network_ref = shift;
	my $network_ref = new_network($old_network_ref);
	my $name1 = shift;
	my $name2 = shift;
	my $name;
	my $dummy_name;
	my @networks_to_return;
	
	## If the nodes are the same, then they are already merged
	## This can happen when I am merging parents
	if($name1 eq $name2)
	{
		return $old_network_ref;
	}
	
	if(does_network_pass_generation_test($network_ref) eq 0)
	{
		return;
	}
	
	## Identify which ones are dummies and which aren't. 
	if($name1 =~ /Dummy/i && $name2 =~ /Dummy/i) 
	{
		$name = $name1;
		$dummy_name = $name2;
		
		## Since both are dummies, and one has parents and the other other doesn't, 
		## then I want the one with parents to be $name, and the one without parents to be $dummy_name
		## If they both have parents, then it doesn't matter
		my @temp_parents = $$network_ref{$dummy_name}->parents();
		if(@temp_parents > 0)
		{
			$name = $name2; 
			$dummy_name = $name1;
		}
	}
	elsif($name2 =~ /Dummy/i)
	{
		$name = $name1;
		$dummy_name = $name2;
	}
	elsif($name1 =~ /Dummy/i) 
	{
		$name = $name2;
		$dummy_name = $name1;
	}
	else
	{
		return;
	}
		
	## merge parents
	my @parents1 = $$network_ref{$name}->parents();
	my @parents2 = $$network_ref{$dummy_name}->parents();
		
	## If both nodes have parents, recursively descend combining ancestors
	if(@parents1 > 0 && @parents2 > 0)
	{
		my @networks_to_return;
		
		## This if statement is here to prevent doing a merge on two pairs of parents that share a common individual
		## If this merge is allowed, it would result in a pair of identical parents, which is not possible 
		if(@parents1[0] ne @parents2[1] && @parents1[1] ne @parents2[0])
		{
			my $new_network_ref = new_network($network_ref);
			my @networks = merge_nodes($new_network_ref,@parents1[0],@parents2[0]);
			foreach my $temp_network_ref (@networks)
			{
				my @networks = merge_nodes($temp_network_ref,@parents1[1],@parents2[1]);
				foreach(@networks)
				{
					if($$_{$name}->pass_generation_check($_) eq 1)
					{
						push(@networks_to_return, $_);
					}
				}
			}
		}
		
		## Merge nodes, part 1 complete
		
		## This if statement is here to prevent doing a merge on two pairs of parents that share a common individual
		## If this merge is allowed, it would result in a pair of identical parents, which is not possible 
		if(@parents1[0] ne @parents2[0] && @parents1[1] ne @parents2[1])
		{	
			my $new_network_ref = new_network($network_ref);
			my @networks = merge_nodes($new_network_ref,@parents1[0],@parents2[1]);
			foreach my $temp_network_ref (@networks)
			{
				my @networks = merge_nodes($temp_network_ref,@parents1[1],@parents2[0]);
				foreach(@networks)
				{
					if($$_{$name}->pass_generation_check($_) eq 1)
					{
						push(@networks_to_return, $_);
					}
				}
			}
		}
		## Merge nodes, part 2 complete

		foreach my $network_ref (@networks_to_return)
		{			
			## Merge children
			my @children = $$network_ref{$dummy_name}->children();
			foreach(@children)
			{
				$$network_ref{$name}->add_child($_);
				$$network_ref{$_}->replace_parent($name,$dummy_name);
			}
			
			## Remove dummy_name from the network
			my @dummy_name_parents = $$network_ref{$dummy_name}->parents();
			foreach(@dummy_name_parents)
			{
				$$network_ref{$_}->remove_child($dummy_name);
			}
			delete $$network_ref{$dummy_name};
		}
		return @networks_to_return;
	}
	
	## Merge children
	my @children = $$network_ref{$dummy_name}->children();
	foreach(@children)
	{
		$$network_ref{$name}->add_child($_);
		$$network_ref{$_}->replace_parent($name,$dummy_name);
	}
	
	## Remove dummy_name from the network
	my @dummy_name_parents = $$network_ref{$dummy_name}->parents();
	foreach(@dummy_name_parents)
	{
		$$network_ref{$_}->remove_child($dummy_name);
	}
	delete $$network_ref{$dummy_name};
	return $network_ref;
}

sub merge_dummy_parents
{
	my $old_network_ref = shift;
	my $network_ref = new_network($old_network_ref);
	my $p1 = shift;
	my $p2 = shift;
	my @p1_parents = $$network_ref{$p1}->parents();
	my @p2_parents = $$network_ref{$p2}->parents();
	
	## These two sets of parents should be dummies or not added yet
		
	if(@p1_parents eq 0 && @p2_parents eq 0) ## If the dummy parents don't exists, make them
	{
		my $dummy1 = add_dummy($network_ref);
		my $dummy2 = add_dummy($network_ref);
		$$network_ref{$p1}->add_parent($dummy1);
		$$network_ref{$p2}->add_parent($dummy1);
		$$network_ref{$p1}->add_parent($dummy2);
		$$network_ref{$p2}->add_parent($dummy2);
		
		$$network_ref{$dummy1}->add_child($p1);
		$$network_ref{$dummy2}->add_child($p1);
		$$network_ref{$dummy1}->add_child($p2);
		$$network_ref{$dummy2}->add_child($p2);						

		if(check_network($network_ref,3) == 0){return $network_ref;}
	}
	elsif(@p1_parents > 1 && @p2_parents eq 0)
	{
		$$network_ref{$p2}->add_parent(@p1_parents[0]);
		$$network_ref{$p2}->add_parent(@p1_parents[1]);
		$$network_ref{@p1_parents[0]}->add_child($p2);
		$$network_ref{@p1_parents[1]}->add_child($p2);
		if(check_network($network_ref,3) == 0){return $network_ref;}
	}
	elsif(@p1_parents eq 0 && @p2_parents > 1)
	{
		$$network_ref{$p1}->add_parent(@p2_parents[0]);
		$$network_ref{$p1}->add_parent(@p2_parents[1]);
		$$network_ref{@p2_parents[0]}->add_child($p1);
		$$network_ref{@p2_parents[1]}->add_child($p1);
		if(check_network($network_ref,3) == 0){return $network_ref;}
	}
	elsif(@p1_parents > 1 && @p2_parents > 1)
	{
		my @networks_to_return;
		my $new_network_ref = new_network($network_ref);
		
		## This if statement is here to prevent doing a merge two pairs of parents that share a common individual
		## If this merge is allowed, it would result in a pair of identical parents, which is not possible 
		if(@p1_parents[0] ne @p2_parents[1] && @p1_parents[1] ne @p2_parents[0])
		{
			my @networks = merge_nodes($new_network_ref,@p1_parents[0],@p2_parents[0]);
			foreach (@networks)
			{
				my @new_networks = merge_nodes($_,@p1_parents[1],@p2_parents[1]);
				push(@networks_to_return, @new_networks);
			}
		}
		
		## This if statement is here to prevent doing a merge two pairs of parents that share a common individual
		## If this merge is allowed, it would result in a pair of identical parents, which is not possible 
		if(@p1_parents[0] ne @p2_parents[0] && @p1_parents[1] ne @p2_parents[1])
		{		
			my $new_network_ref = new_network($network_ref);
			my @networks = merge_nodes($new_network_ref,@p1_parents[0],@p2_parents[1]);
			foreach (@networks)
			{
				my @new_networks = merge_nodes($_,@p1_parents[1],@p2_parents[0]);
				push(@networks_to_return, @new_networks);
			}
		}
		
		my @networks;
		foreach(@networks_to_return)
		{
			if(check_network($_,3) == 0){push(@networks,$_)}
		}
		return @networks;
	}
	else
	{
		die "ERROR!!! shouldn't be here\n";
	}
	return -1;
}

sub merge_FS_parents
{
	my $old_network_ref = shift;
	my $network_ref = new_network($old_network_ref);
	my $p1 = shift;
	my $p2 = shift;
	my @p1_parents = $$network_ref{$p1}->parents();
	my @p2_parents = $$network_ref{$p2}->parents();
	
	## Both FS parents must exists		
	if(@p1_parents < 2 || @p2_parents < 2){die "FS must have both parents\n";}
	my @networks_to_return;
	my $new_network_ref = new_network($network_ref);

	## This if statement is here to prevent doing a merge two pairs of parents that share a common individual
	## If this merge is allowed, it would result in a pair of identical parents, which is not possible 
	if(@p1_parents[0] ne @p2_parents[1] && @p1_parents[1] ne @p2_parents[0])
	{
		my @networks = merge_nodes($new_network_ref,@p1_parents[0],@p2_parents[0]);
		foreach (@networks)
		{
			my @new_networks = merge_nodes($_,@p1_parents[1],@p2_parents[1]);
			push(@networks_to_return, @new_networks);
		}
	}
	elsif(@p1_parents[0] ne @p2_parents[0] && @p1_parents[1] ne @p2_parents[1])
	{		
		my $new_network_ref = new_network($network_ref);
		my @networks = merge_nodes($new_network_ref,@p1_parents[0],@p2_parents[1]);
		foreach (@networks)
		{
			my @new_networks = merge_nodes($_,@p1_parents[1],@p2_parents[0]);
			push(@networks_to_return, @new_networks);
		}
	}
	
	my @networks;
	foreach(@networks_to_return)
	{
		if(check_network($_,1) == 0){push(@networks,$_)}
	}
	return @networks;
}


sub merge_avuncular
{
	my $old_network_ref = shift;
	my $network_ref = new_network($old_network_ref);
	my $uncle = shift;
	my $dummy = shift;
	my @dummy_parents = $$network_ref{$dummy}->parents();
	my @parents = $$network_ref{$uncle}->parents();
	
	## If dummy doesn't have parents make its parents the parents of $uncle;
	## else 1 of the parents must be a dummy and the other must be a parent in common with $uncle
	if(@dummy_parents eq 0)
	{
		$$network_ref{$dummy}->add_parent(@parents[0]);
		$$network_ref{$dummy}->add_parent(@parents[1]);
		$$network_ref{@parents[0]}->add_child($dummy);
		$$network_ref{@parents[1]}->add_child($dummy);
		if(check_network($network_ref,2) == 0)
		{
			return $network_ref;
		}
	}
	else
	{
		if($$network_ref{$dummy}->num_real_parents() eq 0)
		{
			## Should not be possible for a dummy to have two dummy parents in phase2
			return -1;
		}
					
		my $dummy_dummy_parent = @dummy_parents[1];
		my $parent_in_common = @dummy_parents[0];
		my $parent_not_in_common;

		## Find the dummy parent
		if(@dummy_parents[0] =~ /Dummy/i)
		{
			$dummy_dummy_parent = @dummy_parents[0];
			$parent_in_common = @dummy_parents[1];
		}
		
		## Find parent in common
		if($parent_in_common eq @parents[0])
		{
			$parent_not_in_common = @parents[1];
		}
		elsif($parent_in_common eq @parents[1])
		{
			$parent_not_in_common = @parents[0];					
		}
		else ## no parent in common
		{
			$parent_not_in_common = "NA";						
		}
		
		## Merge rels dummy grandparent with node_names other parent (the not-in-common parent) if the dummy grandparent doesn't have any parents
		if($parent_not_in_common ne "NA" && $$network_ref{$dummy_dummy_parent}->parents() eq 0)
		{
			my @networks = merge_nodes($network_ref,$dummy_dummy_parent,$parent_not_in_common);
			if(@networks[0] eq -1)
			{
				die "val equals -1 when merging $dummy_dummy_parent and $parent_not_in_common\n";
			}
			my @networks_to_return;
			foreach(@networks)
			{
				if(check_network($_,2) == 0)
				{	
					push(@networks_to_return, $_);
				}
			}
			if(@networks_to_return < 1)
			{
				return -1;
			}
			return @networks_to_return;
		}
	}
	return -1;
}

sub merge_grandparental
{
	my $old_network_ref = shift;
	my $network_ref = new_network($old_network_ref);
	my $GP = shift;
	my $dummy = shift;
	my @dummy_parents = $$network_ref{$dummy}->parents();

	## If dummy doesn't have parents make its parents the parents of $uncle;
	## else 1 of the parents must be a dummy and the other must be a parent in common with $uncle
	if(@dummy_parents == 0)
	{
		my @networks;
		my $new_network_ref = new_network($network_ref);
		my $dummy_parent = add_dummy($new_network_ref);
		$$new_network_ref{$dummy}->add_parent($dummy_parent);
		$$new_network_ref{$dummy_parent}->add_child($dummy);
		
		$$new_network_ref{$dummy}->add_parent($GP);
		$$new_network_ref{$GP}->add_child($dummy);
		push(@networks,$new_network_ref);

		## foreach case where gp has only dummy children with a dummy spouse, add gp and dummy spouse as parents of $dummy
		my @gp_children = $$network_ref{$GP}->children();
		my @gp_dummy_spouses;
		my %invalid_list;
		foreach my $child (@gp_children)
		{
			my @childs_parents = $$network_ref{$child}->parents();
			
			## Get the spouse
			my $gp_spouse;
			if(@childs_parents[0] eq $GP){$gp_spouse = @childs_parents[1]}
			else{$gp_spouse = @childs_parents[0]}

			## If spouse is not dummy, then no need to proceed with it
			if($gp_spouse !~ /Dummy/i){next}

			## If child is dummy, add $gg_spouse to @gg_parent_dummy_spouses
			if($child =~ /Dummy/i)
			{
				push(@gp_dummy_spouses,$gp_spouse);
			}

			## If child is not dummy, add $gp_spouse to %invalid list
			if($child !~ /Dummy/i)
			{
				$invalid_list{$gp_spouse} = 1;
			}
		}
		my @invalids = keys %invalid_list;
		foreach my $gp_spouse (@gp_dummy_spouses)
		{
			my $new_network_ref = new_network($network_ref);
			if(exists $invalid_list{$gp_spouse}){next;}
			## add gp as a parent
			my @dummy_parents = $$new_network_ref{$dummy}->parents();
			my @temp_networks = add_parent($new_network_ref,$GP,$dummy);
			foreach(@temp_networks)
			{
				## Make sure that the pedigree is valid before moving on
				if(check_network($_,2) != 0){next}
				## add gp_spouse as a parent
				my @temp_networks2 = add_parent($_,$gp_spouse,$dummy);
				## add these resulting networks to the possible networks to return
				push(@networks,@temp_networks2);
			}
		}
		my @networks_to_return;
		foreach(@networks)
		{
			if(check_network($_,2) == 0)
			{
				push(@networks_to_return,$_);
			}
		}
		return @networks_to_return;
	}
	else
	{
		if($$network_ref{$dummy}->num_real_parents() eq 0)
		{
			## Should not be possible for a dummy to have two dummy parents in phase2
			return -1;
		}

		## Find the dummy parent					
		my $dummy_dummy_parent = @dummy_parents[1];
		if(@dummy_parents[0] =~ /Dummy/i)
		{
			$dummy_dummy_parent = @dummy_parents[0];
		}
				
		## Merge dummy grandparent with GP if the dummy grandparent doesn't have any parents
		if($$network_ref{$dummy_dummy_parent}->parents() eq 0)
		{
			my @networks = merge_nodes($network_ref,$dummy_dummy_parent,$GP);
			if(@networks[0] eq -1)
			{
				die "val equals -1 when merging $dummy_dummy_parent and $GP\n";
			}
			my @networks_to_return;
			foreach(@networks)
			{
				if(check_network($_,2) == 0)
				{	
					push(@networks_to_return, $_);
				}
			}
			if(@networks_to_return < 1)
			{
				return -1;
			}
			return @networks_to_return;
		}
	}
	return -1;
}

sub merge_grand_avuncular
{
	my $old_network_ref = shift;
	my $network_ref = new_network($old_network_ref);
	my $g_nephew = shift;
	my $g_uncle = shift;
	
	my @networks_to_return;
	my @parents = $$network_ref{$g_nephew}->parents();

	my @g_uncle_parents = $$network_ref{$g_uncle}->parents();
	if(@g_uncle_parents eq 0) ## This should never happen because 
	{
		die "$g_uncle does not have parents when it should : @g_uncle_parents\n";
	}

	my $just_one_parent = 0;
	my $just_one_g_parent = 0;
	my $just_one_parent = expand_on_just_one_parent($network_ref,@parents);
	my $just_one_g_parent = $just_one_parent;
	
	
	foreach my $parent (@parents)
	{
		if($parent !~ /Dummy/i){next;}
		my @grandparents = $$network_ref{$parent}->parents();

		$just_one_g_parent = expand_on_just_one_parent($network_ref,@grandparents);
			
		## If there are no grandparents, add them before continuing
		while (@grandparents < 2)
		{
			my $dummy_parent_name = add_dummy($network_ref);
			my @networks = add_parent($network_ref,$dummy_parent_name,$parent);
			if(@networks > 1){die "IN merge_grand_avuncular; network_ref SHOULD NOT HAVE RETURNED MORE THAN ONE NETWORK\n";}
			$network_ref = shift(@networks);
			@grandparents = $$network_ref{$parent}->parents();
		}

		foreach my $grandparent (@grandparents)
		{
			my @networks;
			my @great_grandparents = $$network_ref{$grandparent}->parents();

			if(@great_grandparents > 0) ## Merge the grand-avuncular parents with great-grandparents
			{
				my $new_network_ref = new_network($network_ref);
				my @network_refs = merge_dummy_parents($new_network_ref,$grandparent,$g_uncle);
				if(@network_refs[0] ne -1)
				{	
					push(@networks,@network_refs);
				}
			}
			elsif(@great_grandparents eq 0) ## Make grand-avuncular parents the parents of the grandparent
			{
				my $new_network_ref = new_network($network_ref);
				foreach my $node (@g_uncle_parents)
				{
					$$new_network_ref{$grandparent}->add_parent($node);
					$$new_network_ref{$node}->add_child($grandparent);
				}
				## Do pedigree check 
				push(@networks,$new_network_ref);
			}
			else
			{
				die "Shouldn't be here\n";
			}
			foreach(@networks)
			{
				if(check_network($_,3) == 0)
				{	
					push(@networks_to_return, $_);
				}
			}
			if($just_one_g_parent == 1){last;}
		}
		if($just_one_parent == 1){last;}
	}
	return @networks_to_return;
}

sub merge_great_grandparent
{
	my $old_network_ref = shift;
	my $network_ref = new_network($old_network_ref);
	my $gg_child = shift;
	my $gg_parent = shift;
	my @networks_to_return;
	my @parents = $$network_ref{$gg_child}->parents();

	## if neither parent has their own parents or other children, then just make one network
	## The same can be said for each parent: (see previous line)
	## *a more general definition would be that if neither parent has a real relative other than through gg_child 
	## then just make one network (implement this later if necessary)
	my $just_one_parent = 0;
	my $just_one_g_parent = 0;
	$just_one_parent = $just_one_g_parent = expand_on_just_one_parent($network_ref,@parents);
	
	foreach my $parent (@parents)
	{
		if($parent !~ /Dummy/i){next;} ## You can only connect through a dummy parent
		my @grandparents = $$network_ref{$parent}->parents();
		
		$just_one_g_parent = expand_on_just_one_parent($network_ref,@grandparents);

		## If there are no grandparents, add them before continuing
		while (@grandparents < 2)
		{
			my $dummy_parent_name = add_dummy($network_ref);
			my @networks = add_parent($network_ref,$dummy_parent_name,$parent);
			if(@networks > 1){die "IN merge_half_avuncular; network_ref SHOULD NOT HAVE RETURNED MORE THAN ONE NETWORK\n";}
			$network_ref = shift(@networks);
			@grandparents = $$network_ref{$parent}->parents();
		}
		foreach my $grandparent (@grandparents)
		{
			## If grandparent has 2 real parents, next
			if($$network_ref{$grandparent}->num_real_parents() eq 2){next;}
			
			my $new_network_ref = new_network($network_ref);
			my @great_grandparents = $$new_network_ref{$grandparent}->parents();			
			
			## If grandparent didn't have any parents, then add dummy parents
			while (@great_grandparents < 2)
			{
				my $dummy_parent_name = add_dummy($new_network_ref);
				my @networks = add_parent($new_network_ref,$dummy_parent_name,$grandparent);
				if(@networks > 1){die "IN merge_great_grandparent; network_ref SHOULD NOT HAVE RETURNED MORE THAN ONE NETWORK\n";}
				$new_network_ref = shift(@networks);
				@great_grandparents = $$new_network_ref{$grandparent}->parents();
			}
			
			## add gg_parent as a parent
			my @networks = add_parent($new_network_ref,$gg_parent,$grandparent);
			
			## foreach case where gg_parent has only dummy children with a dummy spouse, add gg_parent and dummy spouse
			if($$new_network_ref{$grandparent}->num_real_parents() == 0)
			{
				my @gg_parent_children = $$network_ref{$gg_parent}->children();
				my @gg_parent_dummy_spouses;
				my %invalid_list;
				foreach my $child (@gg_parent_children)
				{
					my @childs_parents = $$network_ref{$child}->parents();
					
					## Get the spouse
					my $gg_parent_spouse;
					if(@childs_parents[0] eq $gg_parent){$gg_parent_spouse = @childs_parents[1]}
					else{$gg_parent_spouse = @childs_parents[0]}
					## If spouse is not dummy, then no need to proceed with it
					if($gg_parent_spouse !~ /Dummy/i){next}

					## If child is dummy, add $gg_spouse to @gg_parent_dummy_spouses
					if($child =~ /Dummy/i)
					{
						push(@gg_parent_dummy_spouses,$gg_parent_spouse);
					}

					## If child is not dummy, add $gg_spouse to %invalid list
					if($child !~ /Dummy/i)
					{
						$invalid_list{$gg_parent_spouse} = 1;
					}
				}
				my @invalids = keys %invalid_list;
				foreach my $gg_spouse (@gg_parent_dummy_spouses)
				{
					if(exists $invalid_list{$gg_spouse}){next;}
					my $new_network_ref = new_network($network_ref);
					## add gg_parent as a parent
					my @temp_networks = add_parent($new_network_ref,$gg_parent,$grandparent);
					foreach(@temp_networks)
					{
						## Make sure this pedigree is valid before moving on
						if(check_network($_,3) != 0){next;}	
							## add gg_spouse as a parent
							my @temp_networks2 = add_parent($_,$gg_spouse,$grandparent);
							## add these resulting networks to the possible networks to return
							push(@networks,@temp_networks2);
					}
				}	
			}
			## Do pedigree check
			foreach(@networks)
			{
				if(check_network($_,3) == 0)
				{	
					push(@networks_to_return, $_);
				}
			}
			if($just_one_g_parent == 1){last;}
		}
		if($just_one_parent == 1){last;}
	}
	return @networks_to_return;
}

sub merge_half_avuncular
{
	my $network_ref = shift;
	my $h_uncle = shift;
	my $h_nephew = shift;
	my @networks_to_return;
	my @h_uncle_parents = $$network_ref{$h_uncle}->parents();
	my @h_nephew_parents = $$network_ref{$h_nephew}->parents();
	
	## h_uncle must have at least 1 dummy parent
	if($$network_ref{$h_uncle}->num_real_parents() eq 2){return;}
	
	my $just_one_parent = 0;
	my $just_one_g_parent = 0;
	$just_one_parent = $just_one_g_parent = expand_on_just_one_parent($network_ref,@h_nephew_parents);
	my $just_one_uncles_parent = expand_on_just_one_parent($network_ref,@h_uncle_parents);
	
	foreach my $parent (@h_nephew_parents)
	{
		## Must be a dummy parent
		if($parent != /Dummy/i){next;}
		my @parent_parents = $$network_ref{$parent}->parents();

		$just_one_g_parent = expand_on_just_one_parent($network_ref,@parent_parents);
		
		## If there are no grandparents, add them before continuing
		while (@parent_parents < 2)
		{
			my $dummy_parent_name = add_dummy($network_ref);
			my @networks = add_parent($network_ref,$dummy_parent_name,$parent);
			if(@networks > 1){die "IN merge_half_avuncular; network_ref SHOULD NOT HAVE RETURNED MORE THAN ONE NETWORK\n";}
			$network_ref = shift(@networks);
			@parent_parents = $$network_ref{$parent}->parents();
		}
		
		foreach my $g_parent (@parent_parents)
		{
			## Must be a dummy grandparent
			if($g_parent != /Dummy/i){next;}
			
			foreach my $h_uncle_parent (@h_uncle_parents)
			{
				## Must be a dummy $h_uncle_parent
				if($h_uncle_parent != /Dummy/i){next;}
				
				my @networks = merge_nodes($network_ref,$h_uncle_parent,$g_parent);
				foreach(@networks)
				{
					if(check_network($_,3) == 0)
					{	
						push(@networks_to_return, $_);
					}
				}
				if($just_one_uncles_parent == 1){last;}
			}
			if($just_one_g_parent == 1){last;}
		}
		if($just_one_parent == 1){last;}
	}
	return @networks_to_return;
}

## Add in that the parents can have more children as long as they are shared children
sub expand_on_just_one_parent
{
	my $network_ref = shift;
	my @parents = @_;
	
	if(@parents == 0){return 1;}
	else
	{
		my @p1_ps = $$network_ref{@parents[0]}->parents();
		my @p2_ps = $$network_ref{@parents[1]}->parents();
		my @p1_c = $$network_ref{@parents[0]}->children();
		my @p2_c = $$network_ref{@parents[1]}->children();
		
		## if neither grandparent has their own parents or children with another person, then just make one network
		if(@p1_ps == 0 && @p2_ps == 0 && do_arrays_match(\@p1_c, \@p2_c)){return 1;}
	}
	return 0;

}

sub add_parent
{
	my $old_network_ref = shift;
	my $network_ref = new_network($old_network_ref);
	my $parent = shift;
	my $child = shift;
	my @networks_to_return;

	my @parents = $$network_ref{$child}->parents();
	if(@parents eq 0)
	{
		$$network_ref{$child}->add_parent($parent);
		$$network_ref{$parent}->add_child($child);
		my $dummy = add_dummy($network_ref);
		$$network_ref{$child}->add_parent($dummy);
		$$network_ref{$dummy}->add_child($child);
		push(@networks_to_return, $network_ref);
	}
	elsif(@parents eq 1)
	{
		die "ONLY HAS ONE PARENT\n";	
	}
	elsif(@parents[0] =~ /Dummy/i)
	{
		#print "1. REPLACE @parents[0] with $parent\n";
		my @networks = merge_nodes($network_ref,$parent,@parents[0]);
		push(@networks_to_return, @networks);
		#@parents[0] = $parent;
	}
	elsif(@parents[1] =~ /Dummy/i)
	{
		#print "2. REPLACE @parents[1] with $parent\n";
		my @networks = merge_nodes($network_ref,$parent,@parents[1]);
		push(@networks_to_return, @networks);
		#@parents[1] = $parent;
	}
	else
	{
		die "Can't add $parent to parents of $child; parents full @parents\n";
	}
	return @networks_to_return;
}

## This subroutine will change $network_ref so it only has a first degree relationship between $name and $first_degree_rel
## and it will make a new network_ref with the HAG and beyond relationships
sub separate_PC_and_FS
{
	my $network_ref = shift;
	my $name = shift;
	my $first_degree_rel = shift;
	my $new_network_ref;
	my $old_relationships = $$network_ref{$name}->possible_relationships($first_degree_rel);
	my @new_fs;
	my @new_pc;
	
	foreach(@$old_relationships)
	{
		if($_ eq "PC")
		{
			push(@new_pc,$_);
		}
		elsif($_ eq "FS")
		{
			push(@new_fs,$_);
		}
		else
		{
			die "Should not be any other possible relationship except for PC and FS\n";
		}
	}
	
	$$network_ref{$name}->add_possible_relationships($first_degree_rel,\@new_pc);
	$$network_ref{$first_degree_rel}->add_possible_relationships($name,\@new_pc);

	$new_network_ref = new_network($network_ref);
	$$new_network_ref{$name}->add_possible_relationships($first_degree_rel,\@new_fs);
	$$new_network_ref{$first_degree_rel}->add_possible_relationships($name,\@new_fs);

	## Remove the pc possible relationships from the new network
	$$new_network_ref{$name}->delete_pc($first_degree_rel);
	$$new_network_ref{$first_degree_rel}->delete_pc($name);
	
	## Remove the fs degree possible relationships from the current network
	$$network_ref{$name}->delete_fs($first_degree_rel);
	$$network_ref{$first_degree_rel}->delete_fs($name);

	return $new_network_ref;
}

## This subroutine will change $network_ref so it only has a first degree relationship between $name and $first_degree_rel
## and it will make a new network_ref with the HAG and beyond relationships
sub separate_first_degree_and_HAG
{
	my $network_ref = shift;
	my $name = shift;
	my $first_degree_rel = shift;
	my $new_network_ref;
	my $old_relationships = $$network_ref{$name}->possible_relationships($first_degree_rel);
	my @new_first_degree;
	my @new_rest;
	
	foreach(@$old_relationships)
	{
		if($_ ne "PC" && $_ ne "FS")
		{
			push(@new_rest,$_);
		}
		else
		{
			push(@new_first_degree,$_);
		}
	}
	$$network_ref{$name}->add_possible_relationships($first_degree_rel,\@new_first_degree);
	$$network_ref{$first_degree_rel}->add_possible_relationships($name,\@new_first_degree);

	if(@new_rest > 0)
	{
		$new_network_ref = new_network($network_ref);
		$$new_network_ref{$name}->add_possible_relationships($first_degree_rel,\@new_rest);
		$$new_network_ref{$first_degree_rel}->add_possible_relationships($name,\@new_rest);

		## Remove the 1st degree possible relationships from the new network
		$$new_network_ref{$name}->delete_pc($first_degree_rel);
		$$new_network_ref{$first_degree_rel}->delete_pc($name);
		$$new_network_ref{$name}->delete_fs($first_degree_rel);
		$$new_network_ref{$first_degree_rel}->delete_fs($name);
		
		## Remove the non-1st degree possible relationships from the current network
		$$network_ref{$name}->delete_hag($first_degree_rel);
		$$network_ref{$first_degree_rel}->delete_hag($name);
		$$network_ref{$name}->delete_cgh($first_degree_rel);
		$$network_ref{$first_degree_rel}->delete_cgh($name);
		$$network_ref{$name}->delete_rest($first_degree_rel);
		$$network_ref{$first_degree_rel}->delete_rest($name);
		return $new_network_ref;
	}
	return;
}
#
## This subroutine will change $network_ref so it only has a FS relationship between $name and $fs_rel
## and it will make a new network_ref with the HAG and beyond relationships
sub seperate_HAG_and_CGH
{
	my $network_ref = shift;
	my $name = shift;
	my $fs_rel = shift;
	my $new_network_ref;
	my $old_relationships = $$network_ref{$name}->possible_relationships($fs_rel);
	my @new_fs = ("FS");
	my @new_rest;
	
	foreach(@$old_relationships)
	{
		if($_ ne "FS")
		{
			push(@new_rest,$_);
		}
	}
	$$network_ref{$name}->add_possible_relationships($fs_rel,\@new_fs);
	$$network_ref{$fs_rel}->add_possible_relationships($name,\@new_fs);

	if(@new_rest > 0)
	{
		$new_network_ref = new_network($network_ref);
		$$new_network_ref{$name}->add_possible_relationships($fs_rel,\@new_rest);
		$$new_network_ref{$fs_rel}->add_possible_relationships($name,\@new_rest);

		## remove the FS possible relationship from the new network
		$$new_network_ref{$name}->delete_fs($fs_rel);
		$$new_network_ref{$fs_rel}->delete_fs($name);

		## remove the the non-FS degree relationships from the possible relationships
		$$network_ref{$name}->delete_hag($fs_rel);
		$$network_ref{$fs_rel}->delete_hag($name);
		$$network_ref{$name}->delete_cgh($fs_rel);
		$$network_ref{$fs_rel}->delete_cgh($name);
		$$network_ref{$name}->delete_rest($fs_rel);
		$$network_ref{$fs_rel}->delete_rest($name);
		return $new_network_ref;
	}
	return;
}

## This subroutine will change $network_ref so it only has a HAG relationship between $name and $hag_rel
## and it will make a new network_ref with the CGH and beyond relationships
sub seperate_HAG_and_CGH
{
	my $network_ref = shift;
	my $name = shift;
	my $hag_rel = shift;
	my $new_network_ref;
	my $old_relationships = $$network_ref{$name}->possible_relationships($hag_rel);
	my @new_hag = ("HAG");
	my @new_rest;
	
	foreach(@$old_relationships)
	{
		if($_ ne "HAG")
		{
			push(@new_rest,$_);
		}
	}
	$$network_ref{$name}->add_possible_relationships($hag_rel,\@new_hag);
	$$network_ref{$hag_rel}->add_possible_relationships($name,\@new_hag);

	if(@new_rest > 0)
	{
		$new_network_ref = new_network($network_ref);
		$$new_network_ref{$name}->add_possible_relationships($hag_rel,\@new_rest);
		$$new_network_ref{$hag_rel}->add_possible_relationships($name,\@new_rest);

		## remove the HAG (2nd degree) possible relationship from the new network
		$$new_network_ref{$name}->delete_hag($hag_rel);
		$$new_network_ref{$hag_rel}->delete_hag($name);

		## remove the the non-2nd degree relationships from the possible relationships
		$$network_ref{$name}->delete_cgh($hag_rel);
		$$network_ref{$hag_rel}->delete_cgh($name);
		$$network_ref{$name}->delete_rest($hag_rel);
		$$network_ref{$hag_rel}->delete_rest($name);
		return $new_network_ref;
	}
	return;
}

sub seperate_CGH_and_REST
{
	my $network_ref = shift;
	my $name = shift;
	my $cgh_rel = shift;
	my $new_network_ref;
	my $old_relationships = $$network_ref{$name}->possible_relationships($cgh_rel);
	my @new_cgh = ("CGH");
	my @new_rest;
	
	foreach(@$old_relationships)
	{
		if($_ ne "CGH")
		{
			push(@new_rest,$_);
		}
	}
	$$network_ref{$name}->add_possible_relationships($cgh_rel,\@new_cgh);
	$$network_ref{$cgh_rel}->add_possible_relationships($name,\@new_cgh);

	if(@new_rest > 0)
	{
		$new_network_ref = new_network($network_ref);
		$$new_network_ref{$name}->add_possible_relationships($cgh_rel,\@new_rest);
		$$new_network_ref{$cgh_rel}->add_possible_relationships($name,\@new_rest);

		$$new_network_ref{$name}->delete_cgh($cgh_rel);
		$$new_network_ref{$cgh_rel}->delete_cgh($name);
		$$network_ref{$name}->delete_rest($cgh_rel);
		$$network_ref{$cgh_rel}->delete_rest($name);
		return $new_network_ref;
	}
	return;
}

sub do_arrays_match
{
	my $arr1_ref = shift;
	my $arr2_ref = shift;

	if(@$arr1_ref ne @$arr2_ref)
	{
		return 0;
	}
	foreach (@$arr1_ref)
	{
		my $val = $_;
		if(!grep($_ eq $val,@$arr2_ref))
		{
			return 0;
		}
	}
	foreach (@$arr2_ref)
	{
		my $val = $_;
		if(!grep($_ eq $val,@$arr1_ref))
		{
			return 0;
		}
	}
	return 1;
}

sub new_network
{
	my $network_ref = shift;
	my %new_network;
	
	foreach my $node_name (keys %$network_ref)
	{
		my $new_node = PRIMUS::node_v7::duplicate($$network_ref{$node_name});
		$new_network{$node_name} = $new_node;
	}
	return \%new_network;
}

sub write_cranefoot_files
{
	my $network_ref = shift;
	my $file = shift;
	my $affected_status_ref = shift;
	my $gender_ref = get_genders($network_ref);
	if($gender_ref == 0){return 0;}
	my $output_directory = shift;
	my $network_num = shift;
	my $ages_ref = shift;
	if($network_num eq ""){$network_num = 1}
	if($output_directory eq ""){$output_directory = "../data";}
	
	my %gender_ambigious;
	foreach my $ID (keys %$network_ref)
	{
		my $g = $$network_ref{$ID}->sex();
		#print "ID $ID ($g)\n";
		if($g == 0)
		{
			if($ID !~ /Dummy/i)
			{
				$gender_ambigious{$ID} = 4;
			}
		}
		my @parent_IDs = $$network_ref{$ID}->parents();
		if(@parent_IDs eq 0){next;}
		if(@parent_IDs[0] =~ /Dummy/i && @parent_IDs[1] =~ /Dummy/i)
		{
			$gender_ambigious{@parent_IDs[0]} = 4;
			$gender_ambigious{@parent_IDs[1]} = 4;
		}

	}
	
	open(OUT,">$output_directory/$file.config");
	print OUT "PedigreeFile\t$output_directory/$file.cranefoot.fam\n";
	print OUT "PedigreeName\t$output_directory/$file\n";
	print OUT "SubgraphVariable FID\n";
	print OUT "NameVariable IID\n";
	print OUT "FatherVariable PID\n";
	print OUT "MotherVariable MID\n";
	print OUT "ShapeVariable GENDER\n";
	print OUT "PatternVariable AFFECTED_STATUS\n";
	print OUT "TextVariable IID\n";
	print OUT "VerboseMode off\n";
	my @keys = keys %$ages_ref;
	if(@keys > 0)
	{
		print OUT "AgeVariable AGES\n";
		print OUT "TextVariable AGES\n";
	}
	print OUT "TextVariable TWIN $output_directory/mz_twins\n";
	close(OUT);
	
	## Write the cranefoot fam file
	open(OUT, ">$output_directory/$file.cranefoot.fam");
	if($ages_ref ne "")
	{
		print OUT "FID\tIID\tPID\tMID\tGENDER\tAFFECTED_STATUS\tAGES\n";
	}
	else
	{
		print OUT "FID\tIID\tPID\tMID\tGENDER\tAFFECTED_STATUS\n";
	}

	my @nodes = keys %$network_ref;
	while (my $ID = shift @nodes)
	{
		my $FID;
		my $node_name = $ID;
		my @parent_IDs = $$network_ref{$ID}->parents();
		my @parents; ## without the FID
		foreach(@parent_IDs)
		{
			push (@parents, $_);
		}
		my $g = "0";
		if(exists $$gender_ref{$ID})
		{
			#print "$ID exists\n";
			$g = $$gender_ref{$ID};
			if($g == 1){$g = 7;}
			else{$g = 1;}
		}
		if(exists $gender_ambigious{$ID})
		{
			#print "$ID ambigious\n";
			$g = 0;
		}

		my $a_status = $$affected_status_ref{$ID};
		if($a_status eq $affected_status_value){$a_status = 61;}
		else{$a_status=1;}
		if(@parents eq 0){@parents = (0,0);}
		if(@parents eq 1){push @parents, 0;}
		if(@parents ne 2)
		{
			die "ERROR!!! Incorrect number of parents ". @parents .": @parents\n";
		}
		my $pat;
		my $mat;
		if(exists $$gender_ref{@parent_IDs[0]} && $$gender_ref{@parent_IDs[0]} == $$gender_ref{@parent_IDs[1]} && $$gender_ref{@parent_IDs[0]} ne 0)
		{
			if($verbose > 1)
			{
				print "WARNING2!!! Parents @parents[0] and @parents[1] are the same gender " . $$gender_ref{@parents[1]}." !!!\n" if $verbose > 0;
				print $LOG "WARNING2!!! Parents @parents[0] and @parents[1] are the same gender " . $$gender_ref{@parents[1]}." !!!\n" if $verbose > 0;
			}
		}
		elsif($$gender_ref{@parent_IDs[0]} == 1)
		{
			$pat = @parents[0];
			$mat = @parents[1];
		}
		elsif($$gender_ref{@parent_IDs[0]} == 2)
		{
			$pat = @parents[1];
			$mat = @parents[0];
		}
		elsif($$gender_ref{@parent_IDs[0]} == 0 && $$gender_ref{@parent_IDs[1]} == 1)
		{
			$pat = @parents[1]; 
			$mat = @parents[0]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 0 && $$gender_ref{@parent_IDs[1]} == 2)
		{
			$pat = @parents[0]; 
			$mat = @parents[1]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 0 && $$gender_ref{@parent_IDs[1]} == 0)
		{
			$pat = @parents[1];
			$mat = @parents[0];
		}
		else
		{
			$pat = @parents[1];
			$mat = @parents[0]; 
		}
		

		## test/add the genders to gender_ref
		if(exists $$gender_ref{$pat} && $$gender_ref{$pat} != 0)
		{
			if($$gender_ref{$pat} != 1)
			{
				die "GENDERS FAIL!!! pat: $pat\n";
			}
		}
		elsif($pat ne 0){$$gender_ref{$pat} = 1;}
		if(exists $$gender_ref{$mat})
		{
			if($$gender_ref{$mat} != 2 && $$gender_ref{$mat} != 0)
			{
				die "GENDERS FAIL!!! mat: $mat($$gender_ref{$mat})\n";
			}
		}
		elsif($mat ne 0){$$gender_ref{$mat} = 2;}

		## Set the dummy fill
		if($node_name =~ /Dummy/i){$a_status = 11}
		my $line = "PRIMUS-Network$network_num\t$node_name\t$pat\t$mat\t$g\t$a_status";
		if($ages_ref ne "")
		{
			if(exists $$ages_ref{$node_name}){$line .= "\tAGE=$$ages_ref{$node_name}";}
			else{$line .= "\tAGE=NA"}
		}
		$line =~ s/Dummy/Missing/g;
        $line =~ s/\w+__//g;
		print OUT "$line\n";
	}
	close(OUT);
	system("$bin_dir/$cranefoot_binary $output_directory/$file.config");
	system("rm $output_directory/$file.topology.txt");
}

sub write_fam_file
{
	my $network_ref = shift;
	my $network_name = shift;
	my $affected_status_ref = shift;
	my $gender_ref = get_genders($network_ref);
	if($gender_ref == 0){return 0;}
	my $output_directory = shift;
	my $network_num = shift;
	if($network_num eq ""){$network_num = 1}
	if($output_directory eq ""){$output_directory = "../data";}
	
	open(OUT, ">$output_directory/$network_name.fam") or die "$!\n";
	foreach my $ID (sort keys %$network_ref)
	{
		my $FID;
		my $node_name = $ID;
		my @parent_IDs = $$network_ref{$ID}->parents();
		my @parents; ## without the FID
		foreach(@parent_IDs)
		{
			push (@parents, $_);
		}
		my $g = $$network_ref{$ID}->sex();
		my $a = 0;
		if (exists $$affected_status_ref{$ID})
		{
			$a = $$affected_status_ref{$ID};
		}
		if(@parents eq 0){@parents = (0,0);}
		if(@parents eq 1){push @parents, 0;}
		if(@parents ne 2)
		{
			die "ERROR!!! Incorrect number of parents ". @parents .": @parents\n";
		}
		my $pat;
		my $mat;
		if(exists $$gender_ref{@parent_IDs[0]} && $$gender_ref{@parent_IDs[0]} == $$gender_ref{@parent_IDs[1]} && $$gender_ref{@parent_IDs[0]} ne 0)
		{
			if($verbose > 1)
			{
				print "WARNING3!!! Parents @parents[0] and @parents[1] are the same gender ". $$gender_ref{@parents[1]}." !!!\n" if $verbose > 0;
				print $LOG "WARNING3!!! Parents @parents[0] and @parents[1] are the same gender ". $$gender_ref{@parents[1]}." !!!\n" if $verbose > 0;
			}
		}
		elsif($$gender_ref{@parent_IDs[0]} == 1)
		{
			$pat = @parents[0];
			$mat = @parents[1]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 2)
		{
			$pat = @parents[1];
			$mat = @parents[0]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 0 && $$gender_ref{@parent_IDs[1]} == 2)
		{
			$pat = @parents[0]; 
			$mat = @parents[1]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 0 && $$gender_ref{@parent_IDs[1]} == 1)
		{
			$pat = @parents[1]; 
			$mat = @parents[0]; 
		}
		else
		{
			$pat = @parents[1]; 
			$mat = @parents[0];
		}
		$$gender_ref{$pat} = 1;
		$$gender_ref{$mat} = 2;
		my $line = "Network$network_num\t$node_name\t$pat\t$mat\t$g\t$a";
		$line =~ s/Dummy/Missing/g;
		print OUT "$line\n";
	}
	close(OUT);
}


## tests every node in the network to make sure that every node matches every described relationship
## returns 0 if pass
## value > 0 for fail
sub check_network
{
	my $network_ref = shift;
	if($network_ref eq ""){return 100}
	my $phase = shift;
	foreach my $node_name (sort keys %$network_ref)
	{
		## Check that parents aren't the same
		my @parents = $$network_ref{$node_name}->parents();
		if(@parents[0] ne "" && @parents[0] eq @parents[1])
		{
			print "[check_network] [ERROR] Parents are the same\n";
			return 100;
		}
		my $node = $$network_ref{$node_name};
		my $val = $node->check_pedigree_relationships($network_ref,$phase); 
			## I am running generation check for every node, running it once should be sufficient.
		if($val > 0){return $val;}
		
		my $fail = $$network_ref{$node_name}->are_relationships_missing_in_pedigree($network_ref,$phase-1);
		if($fail >= 1){
			print "[check_network] [ERROR] relationships missing in pedigree (?)\n";
			return 99;
		}
	}	
	my $num_generations = get_num_generations($network_ref);
	if($num_generations > $MAX_GENERATIONS)
	{
		print "[check_network] [ERROR] Network $network_ref has $num_generations generations. Exceeds MAX_GENERATIONS: $MAX_GENERATIONS\n";
		print $LOG "[check_network] [ERROR] Network $network_ref has $num_generations generations. Exceeds MAX_GENERATIONS: $MAX_GENERATIONS\n";
		return 101;
	}
	
	my $sibling_mating = are_sibs_mating($network_ref);
	if($sibling_mating eq 1)
	{
		## Full-sib mating; that is illegal
		print "[check_network] [ERROR] Illegal full-sib mating\n";
		return 102;
	}
	
    ## Age check
    my %age_flags;
    my $age_flags_ref = \%age_flags;
    ## Test ages to rank by age flags
    if(keys %age > 0 && $enforce_age_filtering)
    {
        PRIMUS::get_age_flags::get_age_flags_in_network($network_ref,\%age,$age_flags_ref);
        foreach my $id1 (keys %{ $$age_flags_ref{$network_ref} })
        {
            foreach my $id2 (keys %{ $$age_flags_ref{$network_ref}{$id1} })
            {
                #print "##################  AGE INCOMPATIBLE!!!!!!!!!!!!!!!!!!!!!!!!\n";
				print "[check_network] [ERROR] Age incompatible\n";
                return 105;
            }
        }
    }

	my $final_check = 0;
	$final_check = 1 if $phase > 3;
	my $pass = mito_check($network_ref,$final_check);
	if($pass != 1) 
	{
		print "FAILED MITO CHECK\n" if $verbose > 0;
		return 103;
	}
	my $pass = y_check($network_ref,$final_check);
	if($pass != 1) 
	{	
		print "FAILED Y CHECK\n" if $verbose > 0;
		return 104;
	}
	resolve_sex_with_mito($network_ref) if $phase > 1;
	resolve_sex_with_y($network_ref) if $phase > 1;
	resolve_sex_with_known_parents($network_ref);
	return 0; # passes otherwise 
}

sub resolve_sex_with_known_parents
{
	my $network_ref = shift;
	foreach my $node_name (sort keys %$network_ref)
	{
		my @parents = $$network_ref{$node_name}->parents();
 		next if @parents < 1 || $parents[0] eq "" || $parents[1] eq "";
		my $p1_sex = $$network_ref{$parents[0]}->sex();
		my $p2_sex = $$network_ref{$parents[1]}->sex();

		next if $p1_sex == 0 && $p2_sex == 0;
		if($p1_sex == 0)
		{
			$p1_sex = 1;
			$p1_sex = 2 if $p2_sex == 1;
			$$network_ref{$parents[0]}->sex($p1_sex);
		}
		elsif($p2_sex == 0)
		{
			$p2_sex = 1;
			$p2_sex = 2 if $p1_sex == 1;
			$$network_ref{$parents[1]}->sex($p2_sex);
		}
		elsif($p2_sex == $p1_sex)
		{
			warn "WARNING: $node_name\'s parents $parents[0] ($p1_sex) and $parents[1] ($p2_sex) have same sex\n";
		}
	}
	
}


## For each node, traverse the pedigree identifying direct mito paths between all pairs of people. If there is 1 person with unknown sex, use the predicted mito relationships to assign the sex of that idividual. If there are more than 1 person with unknown sex, then stop exploying that path. If the path contains 1 unkown male (expect the ends) or a male as the most ancetral person, then stop exploying that path.
## Traverse the pedigree for each person ($node):
## 1. Push node into empty array @path and into %visited, and $ancestral_node = "";
## 2. Call traverse($network_ref,\@path,\%visited,$ancestral_node)
sub resolve_sex_with_mito
{
	my $network_ref = shift;
	foreach my $node_name (sort keys %$network_ref)
	{
		next if $node_name =~ /Dummy/i;
		## next if $node_name does not have any mito matching data ### IMPLEMENT THIS
		#print "\n\nRESOLVING WITH $node_name ##########################################\n";
		my %visited;
		my @path;
		my $ancestral_node = "";
		
		push(@path,$node_name);
		$visited{$node_name} = 1;
		traverse_mito($network_ref,\@path,\%visited,$ancestral_node);
		#exit;
	}
}

## Traverse (intial idea of what it would do)
## 1. Check if valid path (anctral node is not male, not male in path except ends, no more than one unknown sex, if only no_match)
## 2. Check if we can resolve unknown sex, and update if possible
## 3. push all first degree relatives into a %visited
## 4. For each non-male parent ($p): push(@{$path,$p); $ancestral_node = $p; traverse($network_ref,$path,\%visited,$ancestral_node)
## 5. For each full-sib ($fs) : push(@{$path,$fs); traverse_down($network_ref,$path,\%visited,$ancestral_node)
## 6. if $ancestral_node == ""; $ancestral_node = $node;
## 7. For each child ($c) : push(@{$path,$c); traverse_down($network_ref,$path,\%visited,$ancestral_node)
sub traverse_mito
{
	my $network_ref = shift;
	my $path_ref = shift;
	my $visited_ref = shift;
	my $ancestral_node = shift;
	
	my $node = $$path_ref[-1];
	my $start = $$path_ref[0];
	$$visited_ref{$node} = 1;
	#print "\nstart: $start\n";
	#print "PATH: @$path_ref\n";
	#print "end: $node\n";

	## Check for valid path
	if(!is_valid_mito_path($network_ref,$path_ref,$ancestral_node))
	{
		#print "invalid1\n";
		return;
	}
	else
	{
		#print "pass1\n\n";
	}

	## Resolved any unkown sex
	my $continue = 1;
	while($continue == 1)
	{
		$continue = resolve_unknown_sex_mito($network_ref,$path_ref,$ancestral_node);
		resolve_sex_with_known_parents($network_ref) if $continue;
		#print "RESOLVED!!!\n" if $continue;
	}
	#print "Done resolving sex\n\n";

	#### Traverse pedigree
	my @parents = $$network_ref{$node}->parents();
	my @children = $$network_ref{$node}->children();
	my @fullsibs = $$network_ref{$node}->get_full_sibs($network_ref);

	#print "$node parents(@parents) sibs(@fullsibs) children(@children)\n";

	## only need to add full-sibs to visited, because they are the only ones that could be checked in a future recursive call
	foreach(@fullsibs){$$visited_ref{$_} = 1}

	## Traverse parents
	foreach my $parent (@parents)
	{
		my $parent_sex = $$network_ref{$parent}->sex() if $parent ne "";
		#print "parent $parent($parent_sex)\n";
		my @new_path = (@$path_ref,$parent);
		traverse_mito($network_ref,\@new_path,$visited_ref,$parent) if $parent_sex != 1;
	}

	## Traverse fullsibs
	foreach my $fullsib (@fullsibs)
	{
		#print "fullsib $fullsib\n";
		my @new_path = (@$path_ref,$fullsib);
		traverse_down_mito($network_ref,\@new_path,$visited_ref,$ancestral_node);
	}

	## Traverse children who are unvisited
	foreach my $child (@children)
	{
		last if ($$network_ref{$node}->sex()) == 1;
		next if exists $$visited_ref{$child};
		#print "child $child\n";
		my @new_path = (@$path_ref,$child);
		$ancestral_node = $node if $ancestral_node eq "";
		traverse_down_mito($network_ref,\@new_path,$visited_ref,$ancestral_node);
	}
}

sub resolve_unknown_sex_mito
{
	my $network_ref = shift;
	my $path_ref = shift;
	my $ancestral_node = shift;
	

	my $start = $$path_ref[0];
	my $end = $$path_ref[-1];
	return 0 if $start =~ /Dummy/i || $end =~ /Dummy/i; ## Can't resolve if I don't have match data, and I don't have mito matches for Dummies
	#print "\nResolve unknown sex @$path_ref?\n";
	#print "Ancestral_node $ancestral_node\n";
	my $match = $$network_ref{$start}->is_mito_match($end);

	#print "$start <-> $end = $match\n";

	if($match == 1 && !$USE_MATCH_MITO)
	{
		return 0; ## don't proceed if it is a match, but it isn't using matches
	}
	elsif($match == 1 && $USE_MATCH_MITO)
	{
		my $change = 0;
		## Set all individuals in the path as female, except the ends, unless one of the ends is ancestral
		if($ancestral_node ne "")
		{
			my $old_sex = 0;
			$old_sex = $$network_ref{$ancestral_node}->sex() if $ancestral_node ne "";
			$change = 1 if $old_sex ne 2;
			warn "WARNING! resolve_unknown_sex_mito is changing the sex of $ancestral_node from $old_sex to 2" if $old_sex == 1;
			$$network_ref{$ancestral_node}->sex(2);
		}
		
		for(my $i = 1; $i < @$path_ref-1; $i++)
		{
			my $node = $$path_ref[$i];
			my $old_sex = $$network_ref{$node}->sex();
			$change = 1 if $old_sex ne 2;
			warn "WARNING! resolve_unknown_sex_mito is changing the sex of $node from $old_sex to 2" if $old_sex == 1;
			$$network_ref{$node}->sex(2);
		}
		#print "here1\n";
		return $change;
	}
	elsif($match == 0 && !$USE_NO_MATCH_MITO)
	{
		return 0; ## don't proceed if it is a no_match, but it isn't using no_matches
	}
	elsif($match == 0 && $USE_NO_MATCH_MITO)
	{
		my $unknown_sex_node = "";
		## If there is only one individual with unknown sex (excluding the ends, unless it is $ancestral_node), then it is male.
		my $num_sex_unknown = 0;
		for(my $i = 1; $i < @$path_ref-1; $i++)
		{
			my $node = $$path_ref[$i];
			my $node_sex = $$network_ref{$node}->sex();
			
			# invalid of there is a male in the path, except for the ends, unless the end is $ancestral_node, which was caught above.
			return 0 if $node_sex == 1; 
			
			# Count the number of individuals in the middle of the path with unknown_sex
			$num_sex_unknown++ if $node_sex == 0;
			$unknown_sex_node = $node if $node_sex == 0;
		}
		#print "unknown sex node: $unknown_sex_node\n";
		#### Check that there aren't too many unknowns
		## If either $end or $start are ancestral and have unknown sex, then it must be counted in $num_sex_unknown
		my $first = @$path_ref[0];
		my $last = @$path_ref[-1];
		my $sex_ancestral_node = 0;
		$sex_ancestral_node = $$network_ref{$ancestral_node}->sex() if $ancestral_node ne "";
		#print "first: $first; ancetral_node $ancestral_node($sex_ancestral_node)\n";
		#print "last: $last; ancetral_node $ancestral_node($sex_ancestral_node)\n";
		#print "unknown sex node: $unknown_sex_node\n";
		if($sex_ancestral_node == 0 && ($ancestral_node eq $first))
		{
			#print "here\n";
			$num_sex_unknown++;
			$unknown_sex_node = $first;
		}
		elsif($sex_ancestral_node == 0 && ($ancestral_node eq $last))
		{
			$num_sex_unknown++;
			$unknown_sex_node = $last;
		}
		#print "unknown sex node: $unknown_sex_node\n";

		return 0 if $num_sex_unknown > 1 || $num_sex_unknown == 0; ## can't resolve unknown sex if there are 0 or >1 unknowns
		my $unknown_sex_node_sex = $$network_ref{$unknown_sex_node}->sex();
		#print "$unknown_sex_node sex = $unknown_sex_node_sex\n";
		$$network_ref{$unknown_sex_node}->sex(1);
		#print "here2\n";
		return 1;
	}
	elsif($match == -1)
	{
		return 0;
	}
	else
	{
		die "SHOULDN'T BE HERE; match = $match\n";
	}
	#print "here3\n";
	return 0;
}

sub is_valid_mito_path
{
	my $network_ref = shift;
	my $path_ref = shift;
	my $ancestral_node = shift;
	
	#print "\nis valid mito path @$path_ref?\n";

	## Check sex of the ancestral node
	my $sex_ancestral_node = 0;
	$sex_ancestral_node = $$network_ref{$ancestral_node}->sex() if $ancestral_node ne "";
	return 0 if $sex_ancestral_node == 1; ## invalid if ancetral node is male
	
	## Count number of individuals with unknown sex and check if there are any known males in the path
	my $num_sex_unknown = 0;
	for(my $i = 1; $i < @$path_ref-1; $i++)
	{
		my $node = $$path_ref[$i];
		my $node_sex = $$network_ref{$node}->sex();
		
		# invalid of there is a male in the path, except for the ends, unless the end is $ancestral_node, which was caught above.
		return 0 if $node_sex == 1; 
		
		# Count the number of individuals in the middle of the path with unknown_sex
		$num_sex_unknown++ if $node_sex == 0;
	}
	return 1 if $USE_MATCH_MITO;

	#### Check that there aren't too many unknowns (ONLY IF NOT USING MATCH)
	my $first = @$path_ref[0];
	my $last = @$path_ref[-1];
	if($sex_ancestral_node == 0 && ($ancestral_node eq $first || $ancestral_node eq $last))
	{
		$num_sex_unknown++ 
	}
	## invalid if more than one unknown sex in path
	return 0 if $num_sex_unknown > 1;
	return 1;
}

## Traverse down
## 1. Check if valid path (anctral node is not male, not male in path except ends, no more than one unknown sex, if only no_match)
## 2. Check if we can resolve unknown sex, and update if possible
## 3. For each child ($c) : push(@{$path,$c); traverse_down($network_ref,$path,\%visited,$ancestral_node)
sub traverse_down_mito
{
	#print "MITO DOWN\n";
	my $network_ref = shift;
	my $path_ref = shift;
	my $visited_ref = shift;
	my $ancestral_node = shift;
	
	my $node = $$path_ref[-1];
	my $start = $$path_ref[0];
	#print "\nstart: $start\n";
	#print "PATH: @$path_ref\n";
	#print "end: $node\n";

	## Check for valid path
	if(!is_valid_mito_path($network_ref,$path_ref,$ancestral_node))
	{
		#print "invalid\n";
		return;
	}
	else
	{
		#print "pass\n\n";
	}

	## Resolved any unkown sex
	my $continue = 1;
	while($continue == 1)
	{
		$continue = resolve_unknown_sex_mito($network_ref,$path_ref,$ancestral_node);
		resolve_sex_with_known_parents($network_ref) if $continue;
		#print "RESOLVED!!!\n" if $continue;
	}
	#print "Done resolving sex\n\n";

	#### Traverse pedigree
	my @children = $$network_ref{$node}->children();

	## Traverse children who are unvisited
	foreach my $child (@children)
	{
		last if ($$network_ref{$node}->sex()) == 1;
		next if exists $$visited_ref{$child};
		#print "child $child\n";
		my @new_path = (@$path_ref,$child);
		$ancestral_node = $node if $ancestral_node eq "";
		traverse_down_mito($network_ref,\@new_path,$visited_ref,$ancestral_node);
	}
}

## Looks to see if there is an intections between expected mito matches in the pedigree and "no matches" from the mito data
## If sex/gender of a sample is not know, then checking does not proceed through that individual
## IF USE_NO_MATCH_MITO:Find the expected matches for each individual and looks if any of those are in the "no matches" list.
## IF USE_MATCH_MITO: Find the expected matches from the pedigree for each individual and check that the merge of that list and the list of "matches" from the mito data is the same size as the expected matches list found from the pedigree.
sub mito_check
{
	my $network_ref = shift;
	my $final_check = shift;
	my $sex_val = 2;

	print "Mito_check\n" if $verbose > 2;
	foreach my $node_name (sort keys %$network_ref)
	{
		my %should_match;
		my %visited;
		my @nodes_to_visit;

		#print "NODE: $node_name\n";

		## Get all the expected matches from the pedigree
		push(@nodes_to_visit,$node_name);
		while(@nodes_to_visit)
		{
			my $curr_node = shift(@nodes_to_visit);
			next if exists $visited{$curr_node};
			
			my $curr_node_sex = $$network_ref{$curr_node}->sex();
			#print "$curr_node ($curr_node_sex)\n";
			

			$visited{$curr_node} = 1;
			
			my @full_sibs = $$network_ref{$curr_node}->get_full_sibs($network_ref);
			#print "$curr_node sibs @full_sibs\n";
			foreach my $sib (@full_sibs)
			{
				my $sib_sex = $$network_ref{$sib}->sex();
				$should_match{$sib} = 1 if $sib ne $node_name;
				push(@nodes_to_visit,$sib);
			}

			## Add the mom for mito, dad for Y, and no one if sex of the individual is not known
			my @parents = $$network_ref{$curr_node}->parents();
			#print "$curr_node parents @parents\n";
			foreach my $parent (@parents)
			{
				my $sex = $$network_ref{$parent}->sex();
				if($sex == $sex_val)
				{
					$should_match{$parent} = 1 if $parent ne $node_name;
					push(@nodes_to_visit,$parent);
				}
			}

			## add appropriate children if individual is the target sex
			if($curr_node_sex == $sex_val)
			{
				my @children = $$network_ref{$curr_node}->children();
				#print "$curr_node children @children\n";
				foreach my $child (@children)
				{
					my $child_sex = $$network_ref{$child}->sex();
					$should_match{$child} = 1 if $child ne $node_name;
					push(@nodes_to_visit,$child);
				}
			}
		}

		## Check if there is any intersection between pedigree matches and genetic matches
		my @no_matches = $$network_ref{$node_name}->get_mito_no_match();
		my @should_matches = keys %should_match;
		#print "should match: @should_matches\n";
		#print "no_matches: @no_matches\n";
		foreach my $node (@no_matches)
		{
			if(exists $should_match{$node})
			{
				#print "\n$node_name matches with $node, but should not\n";
				return 0;
			}
		}


		## ONLY DO WHAT IS BELOW HERE FOR THE FINAL CHECK; IT WON'T PASS UNTIL RECONSTRUCTION IS COMPLETE
		## Don't do the following checks unless you have 100% accuracy when calling a match with mito data
		if($USE_MATCH_MITO == 1 && $final_check == 1)
		{
			## Check that the merge of @matches and should_match is the same size as either of the two sets
			my @matches = $$network_ref{$node_name}->get_mito_match();
			my @should_matches = keys %should_match;

			#print "should mito matches for $node_name: @should_matches\n";
			#print "expected mito matches for $node_name: @matches\n";
			foreach my $node (@matches)
			{
				if(!exists $should_match{$node})
				{
					
					#print "$node_name doesn't matches with $node, but should\n";
					return 0;
				}
			}
		}
	}
	print "done\n" if $verbose > 2;
	return 1;
}

## For each node, traverse the pedigree identifying direct y paths between all pairs of people. If there is 1 person with unknown sex, use the predicted y relationships to assign the sex of that idividual. If there are more than 1 person with unknown sex, then stop exploying that path. 
## Traverse the pedigree for each person ($node):
## 1. Push node into empty array @path and into %visited, and $ancestral_node = "";
## 2. Call traverse($network_ref,\@path,\%visited,$ancestral_node)
sub resolve_sex_with_y
{
	my $network_ref = shift;
	foreach my $node_name (sort keys %$network_ref)
	{
		next if $node_name =~ /Dummy/i || ($$network_ref{$node_name}->sex() == 2); ## Don't bother if dummy or female
		## next if $node_name does not have any y matching data ### IMPLEMENT THIS
		#print "\n\nRESOLVING Y WITH $node_name ##########################################\n";
		my %visited;
		my @path;
		my $ancestral_node = "";
		
		push(@path,$node_name);
		$visited{$node_name} = 1;
		traverse_y($network_ref,\@path,\%visited,$ancestral_node);
	}
}

## Traverse (intial idea of what it would do)
## 1. Check if valid path (anctral node is not male, not male in path except ends, no more than one unknown sex, if only no_match)
## 2. Check if we can resolve unknown sex, and update if possible
## 3. push all first degree relatives into a %visited
## 4. For each non-male parent ($p): push(@{$path,$p); $ancestral_node = $p; traverse($network_ref,$path,\%visited,$ancestral_node)
## 5. For each full-sib ($fs) : push(@{$path,$fs); traverse_down($network_ref,$path,\%visited,$ancestral_node)
## 6. if $ancestral_node == ""; $ancestral_node = $node;
## 7. For each child ($c) : push(@{$path,$c); traverse_down($network_ref,$path,\%visited,$ancestral_node)
sub traverse_y
{
	my $network_ref = shift;
	my $path_ref = shift;
	my $visited_ref = shift;
	my $ancestral_node = shift;
	
	my $node = $$path_ref[-1];
	my $start = $$path_ref[0];
	$$visited_ref{$node} = 1;
	#print "\nstart: $start\n";
	#print "PATH: @$path_ref\n";
	#print "end: $node\n";

	## Check for valid path
	if(!is_valid_y_path($network_ref,$path_ref,$ancestral_node))
	{
		#print "invalid\n";
		return;
	}
	else
	{
		#print "pass\n\n";
	}

	## Resolved any unkown sex
	my $continue = 1;
	while($continue == 1)
	{
		$continue = resolve_unknown_sex_y($network_ref,$path_ref,$ancestral_node);
		resolve_sex_with_known_parents($network_ref) if $continue;
		#print "RESOLVED!!!\n" if $continue;
	}
	#print "Done resolving sex\n\n";

	#### Traverse pedigree
	my @parents = $$network_ref{$node}->parents();
	my @children = $$network_ref{$node}->children();
	my @fullsibs = $$network_ref{$node}->get_full_sibs($network_ref);

	#print "$node parents(@parents) sibs(@fullsibs) children(@children)\n";

	## only need to add full-sibs to visited, because they are the only ones that could be checked in a future recursive call
	foreach(@fullsibs){$$visited_ref{$_} = 1}

	## Traverse parents
	foreach my $parent (@parents)
	{
		my $parent_sex = $$network_ref{$parent}->sex() if $parent ne "";
		#print "parent $parent($parent_sex)\n";
		my @new_path = (@$path_ref,$parent);
		traverse_y($network_ref,\@new_path,$visited_ref,$parent) if $parent_sex != 2;
	}

	## Traverse fullsibs
	foreach my $fullsib (@fullsibs)
	{
		#print "fullsib $fullsib\n";
		my @new_path = (@$path_ref,$fullsib);
		traverse_down_y($network_ref,\@new_path,$visited_ref,$ancestral_node);
	}

	return if ($$network_ref{$node}->sex()) == 2;
	## Traverse children who are unvisited
	foreach my $child (@children)
	{
		next if exists $$visited_ref{$child};
		#print "child $child\n";
		my @new_path = (@$path_ref,$child);
		$ancestral_node = $node if $ancestral_node eq "";
		traverse_down_y($network_ref,\@new_path,$visited_ref,$ancestral_node);
	}
}

sub resolve_unknown_sex_y
{
	my $network_ref = shift;
	my $path_ref = shift;
	my $ancestral_node = shift;
	

	my $start = $$path_ref[0];
	my $end = $$path_ref[-1];
	return 0 if $start =~ /Dummy/i || $end =~ /Dummy/i; ## Can't resolve if I don't have match data, and I don't have y matches for Dummies
	#print "\nResolve unknown sex @$path_ref?\n";
	my $match = $$network_ref{$start}->is_y_match($end);

	#print "$start <-> $end = $match\n";

	if($match == 1 && !$USE_MATCH_Y)
	{
		return 0; ## don't proceed if it is a match, but it isn't using matches
	}
	elsif($match == 1 && $USE_MATCH_Y)
	{
		my $change = 0;
		## Set all individuals in the path as male, including the ends
		for(my $i = 0; $i < @$path_ref; $i++)
		{
			my $node = $$path_ref[$i];
			my $old_sex = $$network_ref{$node}->sex();
			$change = 1 if $old_sex ne 1;
			warn "WARNING! resolve_unknown_sex_y is changing the sex of $node from $old_sex to 1" if $old_sex == 2;
			$$network_ref{$node}->sex(1);
		}
		return $change;
	}
	elsif($match == 0 && !$USE_NO_MATCH_Y)
	{
		return 0; ## don't proceed if it is a no_match, because it isn't using no_matches
	}
	elsif($match == 0 && $USE_NO_MATCH_Y)
	{
		my $unknown_sex_node = "";
		## If there is only one individual with unknown sex (including the ends), then it is female.
		my $num_sex_unknown = 0;
		for(my $i = 0; $i < @$path_ref; $i++)
		{
			my $node = $$path_ref[$i];
			my $node_sex = $$network_ref{$node}->sex();
			
			# invalid of there is a male in the path, except for the ends, unless the end is $ancestral_node, which was caught above.
			return 0 if $node_sex == 1;
			
			# Count the number of individuals in the middle of the path with unknown_sex
			$num_sex_unknown++ if $node_sex == 0;
			$unknown_sex_node = $node if $node_sex == 0;
		}
		#### Check that there aren't too many unknowns 
		return 0 if $num_sex_unknown > 1 || $num_sex_unknown == 0; ## can't resolve unknown sex if there are 0 or >1 unknowns
		$$network_ref{$unknown_sex_node}->sex(2);
		return 1;
	}
	elsif($match == -1)
	{
		return 0;
	}
	else
	{
		die "SHOULDN'T BE HERE; match = $match\n";
	}
	return 1;
}

sub is_valid_y_path
{
	my $network_ref = shift;
	my $path_ref = shift;
	my $ancestral_node = shift;
	
	#print "\nis valid y path @$path_ref?\n";

	## Count number of individuals with unknown sex and check if there are any known females in the path 
	my $num_sex_unknown = 0;
	for(my $i = 0; $i < @$path_ref; $i++) ## From 0 to path size because ends do matter with Y
	{
		my $node = $$path_ref[$i];
		my $node_sex = $$network_ref{$node}->sex();
		
		# invalid of there is a female in the path, including the ancestral and end nodes
		return 0 if $node_sex == 2; 
		
		# Count the number of individuals in the middle of the path with unknown_sex
		$num_sex_unknown++ if $node_sex == 0;
	}
	return 1 if $USE_MATCH_Y;

	#### Check that there aren't too many unknowns (ONLY IF NOT USING MATCH Y)
	## invalid if more than one unknown sex in path
	return 0 if $num_sex_unknown > 1;
	return 1;
}

## Traverse down
## 1. Check if valid path (anctral node is not male, not male in path except ends, no more than one unknown sex, if only no_match)
## 2. Check if we can resolve unknown sex, and update if possible
## 3. For each child ($c) : push(@{$path,$c); traverse_down($network_ref,$path,\%visited,$ancestral_node)
sub traverse_down_y
{
	#print "Y DOWN\n";
	my $network_ref = shift;
	my $path_ref = shift;
	my $visited_ref = shift;
	my $ancestral_node = shift;
	
	my $node = $$path_ref[-1];
	my $start = $$path_ref[0];
	#print "\nstart: $start\n";
	#print "PATH: @$path_ref\n";
	#print "end: $node\n";

	## Check for valid path
	if(!is_valid_y_path($network_ref,$path_ref,$ancestral_node))
	{
		#print "invalid\n";
		return;
	}
	else
	{
		#print "pass\n\n";
	}

	## Resolved any unkown sex
	my $continue = 1;
	while($continue == 1)
	{
		$continue = resolve_unknown_sex_y($network_ref,$path_ref,$ancestral_node);
		resolve_sex_with_known_parents($network_ref) if $continue;
		#print "RESOLVED!!!\n" if $continue;
	}
	#print "Done resolving sex\n\n";

	#### Traverse pedigree
	my @children = $$network_ref{$node}->children();

	return if ($$network_ref{$node}->sex()) == 2;
	## Traverse children who are unvisited
	foreach my $child (@children)
	{
		next if exists $$visited_ref{$child};
		#print "child $child\n";
		my @new_path = (@$path_ref,$child);
		$ancestral_node = $node if $ancestral_node eq "";
		traverse_down_y($network_ref,\@new_path,$visited_ref,$ancestral_node);
	}
}

sub y_check
{
	my $network_ref = shift;
	my $final_check = shift;
	my $sex_val = 1;

	print "Y_check\n" if $verbose > 2;
	foreach my $node_name (sort keys %$network_ref)
	{
		my %should_match;
		my %visited;
		my @nodes_to_visit;

		#print "NODE: $node_name\n";
		my $node_name_sex = $$network_ref{$node_name}->sex();
		next if $node_name_sex == 2;

		## Get all the expected matches from the pedigree
		push(@nodes_to_visit,$node_name);
		while(@nodes_to_visit)
		{
			my $curr_node = shift(@nodes_to_visit);
			next if exists $visited{$curr_node};
			
			my $curr_node_sex = $$network_ref{$curr_node}->sex();
			#print "$curr_node ($curr_node_sex)\n";
			next if $curr_node_sex == 2;


			$visited{$curr_node} = 1;
			
			my @full_sibs = $$network_ref{$curr_node}->get_full_sibs($network_ref);
			#print "$curr_node sibs @full_sibs\n";
			foreach my $sib (@full_sibs)
			{
				my $sib_sex = $$network_ref{$sib}->sex();
				if($sib_sex == 1) # Only add male sibs if we are looking at Y
				{
					$should_match{$sib} = 1 if $sib ne $node_name;
					push(@nodes_to_visit,$sib);
				}
			}

			## Add dad for Y, and no one if sex of the individual is not known
			my @parents = $$network_ref{$curr_node}->parents();
			#print "$curr_node parents @parents\n";
			foreach my $parent (@parents)
			{
				my $sex = $$network_ref{$parent}->sex();
				if($sex == 1)
				{
					$should_match{$parent} = 1 if $parent ne $node_name;
					push(@nodes_to_visit,$parent);
				}
			}

			## add appropriate children if individual is the target sex
			if($curr_node_sex == $sex_val)
			{
				my @children = $$network_ref{$curr_node}->children();
				#print "$curr_node children @children\n";
				foreach my $child (@children)
				{
					my $child_sex = $$network_ref{$child}->sex();
					if($child_sex == 1) # Only add male children if we are looking at Y
					{
						$should_match{$child} = 1 if $child ne $node_name;
						push(@nodes_to_visit,$child);
					}
				}
			}
		}

		## Check if there is any intersection between pedigree matches and genetic matches (correct?)
		## Check if there is any intersection between pedigree should matches and genetic no matches (correct?)
		my @no_matches = $$network_ref{$node_name}->get_y_no_match();
		my @should_matches = keys %should_match;
		#print "should match: @should_matches\n";
		#print "no_matches: @no_matches\n";
		foreach my $node (@no_matches)
		{
			if(exists $should_match{$node})
			{
				#print "\n$node_name matches with $node, but should not\n";
				#exit;
				return 0;
			}
		}

		## ONLY DO WHAT IS BELOW HERE FOR THE FINAL CHECK; IT WON'T PASS UNTIL RECONSTRUCTION IS COMPLETE
		## Don't do the following checks unless you have 100% accuracy when calling a match with y data
		if($USE_MATCH_Y == 1 && $final_check == 1)
		{
			## Check that the merge of @matches and should_match is the same size as either of the two sets
			my @matches = $$network_ref{$node_name}->get_y_match();
			my @should_matches = keys %should_match;

			#print "should y matches for $node_name: @should_matches\n";
			#print "expected y matches for $node_name: @matches\n";
			foreach my $node (@matches)
			{
				if(!exists $should_match{$node})
				{
					
					#print "$node_name Should match $node\n";
					#exit;
					return 0;
				}
			}
		}
	}
	print "done.\n" if $verbose > 2;
	return 1;
}


## Used to detech if there is a sibling mating, which would result in an inbred relationship
sub are_sibs_mating
{
	my $network_ref = shift;
	foreach my $node_name (sort keys %$network_ref)
	{
		## Check if $node_name's has parents in common with his/her spouse's parents 
		my @parents = $$network_ref{$node_name}->parents();
		my @children = $$network_ref{$node_name}->children();
		foreach my $child (@children)
		{
			my @childs_parents = $$network_ref{$child}->parents();
			## If all individuals are real, don't reject the pedigree
			if($child !~ /Dummy/i && $childs_parents[0] !~ /Dummy/i && $childs_parents[1] !~ /Dummy/i && $child !~ /Missing/i && $childs_parents[0] !~ /Missing/i && $childs_parents[1] !~ /Missing/i)
			{
				next;
			}

			foreach(@childs_parents)
			{
				if($_ eq $node_name){next}
				my @childs_grandparents = $$network_ref{$_}->parents();
				
				## Check if both parents match = full sib mating
				if(@parents > 0 && do_arrays_match(\@childs_grandparents,\@parents))
				{
					if($allow_full_sib_dummy_mating eq 0)
					{
						#print "SIB MATING: @childs_grandparents; p: @parents\n";
						return 1;
					}
				}

				## Check if one parent matches = half-sib mating
				if(@parents < 1 || @childs_grandparents < 1 || $allow_half_sib_dummy_mating eq 1){next}
				if(@parents[0] eq @childs_grandparents[0] || @parents[1] eq @childs_grandparents[0] || @parents[0] eq @childs_grandparents[1] || @parents[1] eq @childs_grandparents[1])
				{
					#print "SIB MATING: @childs_grandparents; p: @parents\n";
					return 1;
				}
			}
		}
	}
	return 0;
}

## Checks that individuals are not mating with people from other generations within the family
sub does_network_pass_generation_test
{
	my $network_ref = shift;
	my %nodes;
	foreach(keys %$network_ref){$nodes{$_} = 1;}
	my @keys = keys %nodes;

	while (@keys > 0)
	{
		## Chech that all nodes in the same connected pedigree as $node_name pass generation test
		my $node_name = @keys[0];
		my $val = $$network_ref{$node_name}->pass_generation_check($network_ref);
		if($val eq 0)
		{
			return 0;
		}
		
		## Remove all the nodes that have already been checked and repeat if there are nodes that haven't been checked
		my %names = $$network_ref{$node_name}->get_subpedigree_names($network_ref);
		
		foreach my $name (keys %names)
		{
			delete $nodes{$name};
		}
		@keys = keys %nodes;
	}
	return 1;
}

sub get_genders
{
	my $network_ref = shift;
	my %genders;
	#print "Genders before proceeding:\n";
	foreach my $node (keys %$network_ref)
	{
		$genders{$node} = $$network_ref{$node}->sex();
	#	print "$node $genders{$node}\n";
	}
	#print "\n";
	my $gender_ref = \%genders;
	my @names = keys %$network_ref;
	my %unresolved;	# names of the nodes whose parents are not both resolved
	foreach(@names){$unresolved{$_} = 1 }
	my %resolved; ## Names of the children whose parents are resolved
	
	my $randomly_assign_sex_to_two_unknown_parents = 0;
	my $name_ctr = 0;
	my $continue = 1;

	## Continue cycling through the names in the same order, starting with the first name, until all have been resolved
	my $loop_ctr = 0;
	while($continue)
	{	
		$loop_ctr++;
		if($loop_ctr > 1000000)
		{
			warn "Gender check is in an infinite loop; exiting it\n";
			return 0;
		}
		if($name_ctr > @names-1)
		{
			my @names_resolved = keys %resolved;
			if(@names_resolved == @names)
			{
				## YOU'RE DONE
				return $gender_ref;
			}
			$name_ctr = 0;
			$randomly_assign_sex_to_two_unknown_parents = 1; ## since it is starting over from the begining of the list, allow 
		}
		my $node_name = @names[$name_ctr++]; ## To speed this up use an array of the unresolved samples instead.
		## If it has already been resolved, go to next
		if(exists $resolved{$node_name}){next}

		#my @children = $$network_ref{$node_name}->children();
		#my @parents = $$network_ref{$node_name}->parents();
		my @parent_IDs = $$network_ref{$node_name}->parents();
		if(@parent_IDs eq 0)
		{
			$resolved{$node_name} = 1;
			delete $unresolved{$node_name};
			next;
		}
		my $pat;
		my $mat;
		if(exists $$gender_ref{@parent_IDs[0]} && $$gender_ref{@parent_IDs[0]} eq $$gender_ref{@parent_IDs[1]} && $$gender_ref{@parent_IDs[0]} ne 0)
		{
			if($verbose > 1)
			{
				warn "WARNING1!!! Parents @parent_IDs[0] and @parent_IDs[1] are the same gender ". $$gender_ref{@parent_IDs[1]}." !!!\n";
			}
			return 0;
		}
		elsif($$gender_ref{@parent_IDs[0]} eq 1)
		{
			$pat = @parent_IDs[0];
			$mat = @parent_IDs[1]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} eq 2)
		{
			$pat = @parent_IDs[1];
			$mat = @parent_IDs[0]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 0 && $$gender_ref{@parent_IDs[1]} == 2)
		{
			$pat = @parent_IDs[0]; 
			$mat = @parent_IDs[1]; 
		}
		elsif($$gender_ref{@parent_IDs[0]} == 0 && $$gender_ref{@parent_IDs[1]} == 1)
		{
			$pat = @parent_IDs[1]; 
			$mat = @parent_IDs[0]; 
		}
		else
		{
			#print "parents are both unspecified sex\n";
			if($randomly_assign_sex_to_two_unknown_parents)
			{
				$pat = @parent_IDs[1];
				$mat = @parent_IDs[0];
				$randomly_assign_sex_to_two_unknown_parents = 0;
			}
			else
			{
				next;
			}
		}
		
		if(exists $$gender_ref{$pat})
		{
			if($$gender_ref{$pat} eq 2)
			{
				#print "GENDERS FAIL!!! pat: $pat\n";
				return 0;
			}
			else{$$gender_ref{$pat} = 1;}
		}
		else{$$gender_ref{$pat} = 1;}
		if(exists $$gender_ref{$mat})
		{
			if($$gender_ref{$mat} eq 1)
			{
				#print "GENDERS FAIL!!! mat: $mat\n";
				return 0;
			}
			else{$$gender_ref{$mat} = 2;}
		}
		else{$$gender_ref{$mat} = 2;}
		
		$resolved{$node_name} = 1;
		delete $unresolved{$node_name};
	}
	return $gender_ref;
}

######################################################################################################################
######################################################################################################################
######################################################################################################################


sub build_network_from_fam_file
{
	my $fam_file = shift;
	print "Building network for $fam_file\n" if $verbose > 0;
	print $LOG "Building network for $fam_file\n" if $verbose > 0;
	
	## Read in file
	open(IN,$fam_file);
	my %all_nodes_network;
	my $network_ref = \%all_nodes_network;
	my @network_refs;
	
	## Build pedigree
	while(my $line = <IN>)
	{
		chomp($line);
		my ($FID,$IID,$PID,$MID,$SEX,$PHENOTYPE) = split(/\s+/,$line);
		my $child = "$IID";
		my $mom = "$MID";
		my $dad = "$PID";
		
		if(!exists $$network_ref{$child})
		{
			my $node = new PRIMUS::node_v7($child);
			$$network_ref{$child} = $node;
		}
		if(!exists $$network_ref{$dad} && $PID ne 0)
		{
			my $node = new PRIMUS::node_v7($dad);
			$$network_ref{$dad} = $node;
		}
		if(!exists $$network_ref{$mom} && $MID ne 0)
		{
			my $node = new PRIMUS::node_v7($mom);
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
	
	return @network_refs;	
}


sub write_out_relationships
{
	my $file = shift;
	print "Writing relationships to $file\n" if $verbose > 0;
	print $LOG "Writing relationships to $file\n" if $verbose > 0;
	my @network_refs = @_;
	open(OUT,">$file");
	print OUT "FID1\tIID1\tFID2\tIID2\tRELATIONSHIP\n";
	my %added;
	
	my @PC_likelihoods = (1,0,0,0,0,0);
	my @FS_likelihoods = (0,1,0,0,0,0);
	my @HAG_likelihoods = (0,0,1,0,0,0);
	my @CGH_likelihoods = (0,0,0,1,0,0);
	my @D1C_likelihoods = (0,0,0,2,0,0);
	my @HS1C_likelihoods = (0,0,1,1,0,0);
	my @DHAG_likelihoods = (0,0,2,0,0,0);
	my @PCCGH_likelihoods = (1,0,0,1,0,0);	
	my @DR_likelihoods = (0,0,0,0,1,0);
	my @UN_likelihoods = (0,0,0,0,0,1);
	
	
	
	foreach my $network_ref (@network_refs)
	{
		foreach my $node_name (keys %{$network_ref})
		{
      #my ($FID1, $IID1) = split(/__/,$node_name);
			my $IID1 = $node_name;
			my %rels = $$network_ref{$node_name}->relatives();
			foreach my $rel_name (keys %{$network_ref})
			{
				if($rel_name eq $node_name){next;}
				if(exists $added{"$node_name-$rel_name"}){next;}
				
        #my ($FID2, $IID2) = split(/__/,$rel_name);
				my $IID2 = $rel_name;
				my $relationship = $rels{$rel_name};
				my $rel;
				my @test = @PC_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "PC";}
				my @test = @FS_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "FS";}
				my @test = @HAG_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "HAG";}
				my @test = @CGH_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "CGH";}
				my @test = @D1C_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "D1C";}
				my @test = @HS1C_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "HS1C";}
				my @test = @DR_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "DR";}
				my @test = @UN_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "UN";}
				my @test = @DHAG_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "DHAG";}
				my @test = @PCCGH_likelihoods;
				if(@$relationship[0] eq @test[0] && @$relationship[1] eq @test[1] && @$relationship[2] eq @test[2] && @$relationship[3] eq @test[3] && @$relationship[4] eq @test[4] && @$relationship[5] eq @test[5]){$rel = "PCCGH";}

				if(!exists $added{"$node_name-$rel_name"})
				{
					print OUT "$IID1\t$IID1\t$IID2\t$IID2\t$rel\n";
				}
				$added{"$node_name-$rel_name"} = 1;
				$added{"$rel_name-$node_name"} = 1;
			}
		}
	}
}

return 1;
