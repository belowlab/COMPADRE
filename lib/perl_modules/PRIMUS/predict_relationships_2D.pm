package PRIMUS::predict_relationships_2D;
use strict;
use Data::Dumper;

my $KDE_density_resolution = 1000; 
my $MIN_LIKELIHOOD = 0.1;
my @likelihood_names = qw(PC FS HAG CGH DIS UN MZ);
my $HAG = "";
my $PC;
my $FS;
my $CGH;
my $DR;
my $UN;
my $verbose = 1;
my $curr_bw;
my $curr_kde_type;
#my $ersa_data;

# THIS IS WHAT WE ARE INTERESTED IN REPLACING
sub get_relationship_likelihood_vectors {
	
	my $IBD_file_ref = shift;
	my $IBD_file = $$IBD_file_ref{'FILE'};
	# read in from the reference previously in the program

	print "\nGETTING RELATIONSHIP LIKELIHOOD VECTOR\n\n.genome file: $IBD_file\n";

	$MIN_LIKELIHOOD = shift;
	$verbose = shift;
	my $lib_dir = shift;
	my $output_dir = shift;

	# adding in ersa match file here --ersa_data
	my $match_data = $main::ersa_data_glob;
	print ("ERSA file: $match_data\n\n");

	my $bw = shift;
	my $kde_type = shift;

	my %relationships;
	my %raw_relationship_densities;
	my $original_sum;
	my @fails;
	my $total = 0;
	my $total_possibilities = 0;
	#my $outfile = "$IBD_file\_KDE_likelihood_vectors_prop$MIN_LIKELIHOOD\_bw$bw\_$kde_type.txt\n";
	my $outfile = "$IBD_file\_KDE_likelihood_vectors\n";
	#print "[NEW LOG] OUTFILE : $outfile\n";
	print "[NEW LOG] IBDFILE : $IBD_file\n";
	#my $outfile = "$output_dir/KDE_likelihood_vectors.txt";
	my %possibility_counts;

	if ($MIN_LIKELIHOOD == "") {
		$MIN_LIKELIHOOD = 0.3
	}
	
  	open(PROB_OUT,">$outfile") or die "Can't open likelihood vector output file ($outfile): $!\n";
	open(MZ_OUT,">$output_dir/mz_twins") or die "Can't open mz twin output file ($output_dir/mz_twins): $!\n";
	print MZ_OUT "IID\tTWIN\n";
  	print PROB_OUT "FID1\tIID1\tFID2\tIID2\tMOST_LIKELY_REL\tLIKELIHOOD_VECTOR\tIBD0\tIBD1\tIBD2\tPI_HAT\t-1\tPOSSIBLE_RELS\tLIKELIHOOD_CUTOFF\n";

  	##################### If likelihood vectors were read in as an input file, then bypass looking [at] reference data to save on runtime
  	if (exists $$IBD_file_ref{'likelihood_vectors'} && $$IBD_file_ref{'likelihood_vectors'} == 1) {
		print "LOAD file: $$IBD_file_ref{'FILE'}\n";
		$total_possibilities = load_likelihood_vectors_from_file($$IBD_file_ref{'FILE'},\%raw_relationship_densities,\%relationships);
		close(PROB_OUT);
		close(MZ_OUT);
	  	return (\%relationships,\%raw_relationship_densities, $total_possibilities, @fails);
  	}

  ######################
	if ($bw eq "") {
		#$bw = 9;
		#$type = "clean";
		$bw = 7;
		$kde_type = "noisy";
	}
	if ($kde_type eq "") {
		$kde_type = "noisy";
	}
	if ($bw ne $curr_bw) {
		$HAG ="";
		$curr_bw = $bw;
	}
	if ($kde_type ne $curr_kde_type) {
		$HAG ="";
		$curr_kde_type = $kde_type;
	}

	my $PC_KDE_file = "$lib_dir/KDE_data/PO_KDE_bw17\_$kde_type";
	my $FS_KDE_file = "$lib_dir/KDE_data/FS_KDE_bw2_$kde_type";
	my $HAG_KDE_file = "$lib_dir/KDE_data/2nd_KDE_bw6\_$kde_type";
	my $CGH_KDE_file = "$lib_dir/KDE_data/3rd_KDE_bw4\_$kde_type";
	my $DR_KDE_file = "$lib_dir/KDE_data/UN_KDE_bw6\_$kde_type"; ## This performs better than the 4th degree because the 4th degree has too much overlap with 3rd degree. It would be better to train DR class with a uniform distribution of samples between 4th and 12 degree relatives. Until then, this is the next best thing.
	my $UN_KDE_file = "$lib_dir/KDE_data/UN_KDE_bw1\_$kde_type";
	
	#my $PC_KDE_file = "/nfs/home/grapas2/projects/2011/reconstruct_pedigrees/data/simulations/training_data/PO_KDE_bw17\_$kde_type";
	
	if ($HAG eq "") {
		$HAG = get_data_arrays($HAG_KDE_file);
		$PC = get_data_arrays($PC_KDE_file);
		$FS = get_data_arrays($FS_KDE_file);
		$CGH = get_data_arrays($CGH_KDE_file);
		$DR = get_data_arrays($DR_KDE_file);
		$UN = get_data_arrays($UN_KDE_file);
	}

	if ($MIN_LIKELIHOOD == "") {
		$MIN_LIKELIHOOD = 0.3
	}

	### Load in .genome file 

	open(IN,$IBD_file) or die "Can't open $IBD_file; $!";
	my $header = <IN>;
	while (my $line = <IN>) {
		$line =~ s/^\s+//;
		chomp($line);
		my @temp = split(/\s+/,$line);
		## Multiple by resultion of KDE density values and round to the nearest integer
		my $FID1 = @temp[$$IBD_file_ref{'FID1'}-1]; 
		my $IID1 = @temp[$$IBD_file_ref{'IID1'}-1];
		my $FID2 = @temp[$$IBD_file_ref{'FID2'}-1];
		my $IID2 = @temp[$$IBD_file_ref{'IID2'}-1];
		my $PI_HAT = @temp[$$IBD_file_ref{'PI_HAT'}-1];
		$total++;
		#my $name1 = "$FID1\__$IID1";
		#my $name2 = "$FID2\__$IID2";
		my $name1 = "$IID1";
		my $name2 = "$IID2";
		
		my @vector;

		my $k0 = int(@temp[$$IBD_file_ref{'IBD0'}-1] * $KDE_density_resolution);
		my $k1 = int(@temp[$$IBD_file_ref{'IBD1'}-1] * $KDE_density_resolution);
		my $k2 = int(@temp[$$IBD_file_ref{'IBD2'}-1] * $KDE_density_resolution);
		
		## Make sure all IBD proportions are between 0 and 1, rounding up and down, respectively.
		if ($k2 > 1 * $KDE_density_resolution) {
			$k2 = 1 * $KDE_density_resolution
		};
		if ($k2 < 0) {
			$k2 = 0
		};
		if ($k1 > 1 * $KDE_density_resolution) {
			$k1 = 1 * $KDE_density_resolution
		};
		if ($k1 < 0) {
			$k1 = 0
		};
		if ($k0 > 1 * $KDE_density_resolution) {
			$k0 = 1 * $KDE_density_resolution
		};
		if ($k0 < 0) {
			$k0 = 0
		};

		## This will combine k1 and k2 into k1 so that PRIMUS can do inbred pedigrees that are fully connected by 1st degree relationships, 
		## and at least get the inbred PC relationships correct
		if ($k0 < 0.01 * $KDE_density_resolution && $k1 > .9 * $KDE_density_resolution) {
			$k1 = $k1 + $k2; 
			$k2 = 0;
		}

		my $relationship = @temp[14];
		$relationship =~ s/^AV$/HAG/;
		$relationship =~ s/^GG$/HAG/;
		$relationship =~ s/^HS$/HAG/;
		$relationship =~ s/^GGG$/CGH/;
		$relationship =~ s/^HAV$/CGH/;
		$relationship =~ s/^1C$/CGH/;
		$relationship =~ s/^GAV$/CGH/;
		
		## Hard coded cutoffs to try to catch poor IBD estimates that fall outside the trained KDE regions.
		my @density_vector = ($PC->[$k0][$k1],$FS->[$k0][$k1],$HAG->[$k0][$k1],$CGH->[$k0][$k1],$DR->[$k0][$k1],$UN->[$k0][$k1],0);

		## Check if MZ twins
		if ($k2 > .8 * $KDE_density_resolution) {
			@density_vector = (0,1,0,0,0,0,1);
			@vector = (0,1,0,0,0,0,1);
			my @possibilities = predict_relationship(@vector);
			print MZ_OUT "$name1\tMZ-$name2\n";
			print MZ_OUT "$name2\tMZ-$name1\n";
			$relationships{$name1}{$name2} = \@vector;
		}
		
		my @density_vector = ($PC->[$k0][$k1],$FS->[$k0][$k1],$HAG->[$k0][$k1],$CGH->[$k0][$k1],$DR->[$k0][$k1],$UN->[$k0][$k1]);
		#print "density_vactor: @density_vector\n";
		## Hard coded cutoffs to try to catch poor IBD estimates that fall outside the trained KDE regions.
        if ($k1 > .8 * $KDE_density_resolution) { ## Corrects the 2nd degree splash up near PO
            @density_vector[0] += 1
        }
        
    	## Check of IBD estimates that fall outside the range of the KDE
		if (sum(0,@density_vector) < 0.000001) {
			## MZ test
			if ($k2 > .8 * $KDE_density_resolution) {
				@density_vector[1] = 100;
				@density_vector[6] = 100;
			}
			elsif ($verbose > 0) {
				warn "WARNING!!! Probability vector confidence too low. This is likely due to messy IBD estinates. Applying hard relationship cuttoffs to relationship prediction of $line.\n";
				print "density vector: @density_vector\n";
			
                if ($k1 > .8 * $KDE_density_resolution) {
                    if ($k2 < .01 * $KDE_density_resolution && $k1 < .9 * $KDE_density_resolution) {
                        ## Could be a very inflated but true IBD1 value for 2nd degree relatives 
						@density_vector[2] = 100;
                    }
                    else { # It is possibly noise from a parent/offspring relationship
                        @density_vector[0] = 100;
                    }
                }
                elsif ($k1 < (.1 * $KDE_density_resolution) && ($k0+$k2) > (.9 * $KDE_density_resolution)) {
                    @density_vector[5] = 100;
                }
                elsif ($k1 > (.6 * $KDE_density_resolution) && $k1 < (.8 * $KDE_density_resolution) && ($k0+$k1) > (.9 * $KDE_density_resolution)) {
                    @density_vector[2] = 100;
                }
                else {
                    print "$name1 <-> $name2 = $k0:$k1:$k2 (IBD 0:1:2)\n";
                    print "D_vector: @density_vector\n";
                    die "Unable to predict relationship from IBD estimates for $name1 <-> $name2: $temp[$$IBD_file_ref{'IBD0'}-1]:$temp[$$IBD_file_ref{'IBD1'}-1]:$temp[$$IBD_file_ref{'IBD2'}-1] (IBD 0:1:2); Probability vector confidence too low (@density_vector). Try to get more accurate genome-wide IBD estimates (e.g. adjust for admixture or provide more accurate reference minor allele frequencies when calculating the IBD estimates)\n";
                }
			}
		}
		
		my $original_sum = sum(@density_vector);
		my $original_str = join(',',@density_vector);
		#print ("\nOriginal vector: $original_str\n");
		my @vector = normalize(@density_vector); # divides by sum

		# make a copy for later 
		my @vector_copy = @vector;
		
		###################################
		# density vector updates are done now
		# at this point we can pass the vector AND the set of IDs into the python script that will use ERSA to update it
		## re-update a new density_vector after the ERSA processing has been done 

		
		#print "\nOld vector string ( $name1 $name2 ) : $vector_str\n";

		# Check vector to see if we even need to run ersa
		my $sum01 = $vector[0] + $vector[1];
		print "\n0 1 sum : $sum01\n";

		if ($sum01 < 0.4) {
			
			my $vector_str = join(',',@vector);
			my $python_utility_path = "/data100t1/home/grahame/projects/compadre/primus-ersa-v2/helper3.py";

			my $new_vector = `python3 \"$python_utility_path\" \"$name1\" \"$name2\" \"$vector_str\" \"$match_data\"`;
			chomp($new_vector);

			# re-assign vector to use the new version from the python utility
			@vector = split(',', $new_vector);

			print ("\nNew vector string ( $name1 $name2 ) : $new_vector\n");
		
		}

		else { print ("\nDidn't need to run ERSA\n"); }
		
		###################################


		## If the intial_likelihood_cutoff is dropped low enough, the FS will overlap with HAG (2nd degree), 
		if (@vector[1] > $MIN_LIKELIHOOD && @vector[2] > $MIN_LIKELIHOOD) {
			#if($verbose > -1)
			#{
			#		warn "WARNING!!! Both FS and HAG have sufficiently high likelihoods to be considered. PRIMUS will only reconstruction with the HAG relationship.\n";
			#}
			#@vector[2] += @vector[1];
			#@vector[1] =  0;	
		}
		## If the intial_likelihood_cutoff is dropped low enough, the HAG (2nd degree) will overlap with PC, 
		if (@vector[0] > $MIN_LIKELIHOOD && @vector[2] > $MIN_LIKELIHOOD) {
			if ($verbose > 1) {
				warn "WARNING!!! Both PC and HAG have sufficiently high likelihoods to be considered. PRIMUS will only reconstruction with the PC relationship.\n";
			}
			@vector[0] += @vector[2];
			@vector[2] =  0;
			
		}

		##################################
		# This is the last point at which the vector hashes are updated before being returned.
		# I think that the additional processing starting at line 193 is important and could be affecting why the outside-manufactured estimates aren't computing correctly when it comes to reconstruction

		### Vector added to total relationships data structure

		$relationships{$name1}{$name2} = \@vector;

		# fix density vector with the multiplier from before
		my @density_vector_new = map { $_ * $original_sum } @vector;
		my $new_sum = sum(@density_vector_new);

		my $new_str_normalized = join(',',@vector);
		my $new_str = join(',',@density_vector_new);
		#print ("\nNew vector: $new_str\n");
		#print ("\nNew vector (normalized): $new_str_normalized\n");


		$raw_relationship_densities{$name1}{$name2} = \@density_vector_new;

		##################################

		###########
		# Update 2/25/24
		# Theoretically, the relationship prediction is made AFTER my new vector has been generated  -- it's done in the few lines below, given the vector 

		##Testing stuff
		my @possibilities = predict_relationship(@vector);
		my $rel = get_maximum_relationship(@vector);

		my @possibilities_old = predict_relationship(@vector_copy);
		my $rel_old = get_maximum_relationship(@vector_copy);

		if ($rel ne $rel_old) {
			print "\nNEW relationship prediction ( $name1 $name2 ) : $rel (previously $rel_old)\n";
		}

		my $ibd0 = $k0/$KDE_density_resolution;
		my $ibd1 = $k1/$KDE_density_resolution;
		my $ibd2 = $k2/$KDE_density_resolution;

		### 2/22/24
		### This is where the entire likelihood vector line gets put together and written to outfile 

		print PROB_OUT "$FID1\t$IID1\t$FID2\t$IID2\t$rel\t".join(',',@vector)."\t$ibd0\t$ibd1\t$ibd2\t$PI_HAT\t-1\t".join(',',@possibilities)."\t$MIN_LIKELIHOOD\n";
		
		my $num_possibilities = @possibilities;
		$total_possibilities += $num_possibilities;
		
		if ($verbose > 2) {
			print "line: $line\n";
			print "likelihood vector: $name1 -> $name2 = @vector\n";
			print "densities: @density_vector\n";
		}
		
		## Test if relationship vector matches the expected relationship (if provided)
		if (grep {$_ eq $relationship} @possibilities) {
			if ($verbose > 3) {
				print "PASS.\n";
				print "line: $line\n";
				print "likelihood vector: $name1 -> $name2 = @vector\n";
				print "densities: @density_vector\n";
			}
		}
		else {
			if ($verbose > 3) {
				print "FAIL: $line\n";
				print "$name1 -> $name2 = @vector\n";
				print "densities: @density_vector\n";
				print "@temp[6]:@temp[7]:@temp[8]\n";
				#push(@fails,$line);
				#die "FAIL!!! @k1_vector != $relationship\n";
			}
		}
	}
	close(IN);
	close(MZ_OUT);

	##################################################

	# print "\nPrinting data that's sent to reconstruct_pedigree ...\n\n";
	# print "Relationships : \%relationships\n";
	# print "Relationship densities : \%raw_relationship_densities\n";
	# print "Total possibilities : \%total_possibilities\n\n";


	return (\%relationships,\%raw_relationship_densities, $total_possibilities, @fails);

}

sub load_likelihood_vectors_from_file {

  my $file = shift;
  my $raw_relationships_ref = shift;
  my $relationships_ref = shift;

  my $total_possibilities = 0;

  open(IN,$file) or die "can't open likelihood_vectors_file ($file): $!\n";
  my $header = <IN>; ## remove header

  while (my $line = <IN>) {
    my ($FID1,$IID1,$FID2,$IID2,$REL,$vector,$ibd0,$ibd1,$ibd2,$PI_HAT,@rest) = split(/\s+/,$line);
    my @vector = split(/,/,$vector);
    #my $name1 = "$FID1\__$IID1";
    #my $name2 = "$FID2\__$IID2";
		my $name1 = "$IID1";
		my $name2 = "$IID2";

		## Check if MZ twins
		if (@vector[6] > 0) {
			print MZ_OUT "$name1\tMZ-$name2\n";
			print MZ_OUT "$name2\tMZ-$name1\n";
		}
  
		## If the intial_likelihood_cutoff is dropped low enough, the FS will overlap with HAG (2nd degree); should probably change the reconstruction code to split FS and HAG out from eachother like I do with 2nd and 3rd degree, and 3rd and unrelated 
		if (@vector[1] > $MIN_LIKELIHOOD && @vector[2] > $MIN_LIKELIHOOD) {
			if ($verbose > -1)  {
				warn "WARNING!!! Both FS and HAG have sufficiently high likelihoods to be considered. PRIMUS will only reconstruction with the HAG relationship.\n";
			}
			@vector[2] += @vector[1];
			@vector[1] =  0;
			
		}
		## If the intial_likelihood_cutoff is dropped low enough, the HAG (2nd degree) will overlap with PC, 
		if (@vector[0] > $MIN_LIKELIHOOD && @vector[2] > $MIN_LIKELIHOOD) {
			if ($verbose > -1)  {
				warn "WARNING!!! Both PC and HAG have sufficiently high likelihoods to be considered. PRIMUS will only reconstruction with the PC relationship.\n";
			}
			@vector[0] += @vector[2];
			@vector[2] =  0;
			
		}

		$$raw_relationships_ref{$name1}{$name2} = \@vector;
    	$$relationships_ref{$name1}{$name2} = \@vector;

		##Testing stuff
		my @possibilities = predict_relationship(@vector);
		my $rel = get_maximum_relationship(@vector);

    	#print PROB_OUT "$name1\t$name2\t".join(',',@vector)."\t$rel\t".join(',',@possibilities)."\t$MIN_LIKELIHOOD\n";
		print PROB_OUT "$FID1\t$IID1\t$FID2\t$IID2\t$rel\t".join(',',@vector)."\t$ibd0\t$ibd1\t$ibd2\t$PI_HAT\t-1\t".join(',',@possibilities)."\t$MIN_LIKELIHOOD\n";

    	my $num_possibilities = @possibilities;
    	$total_possibilities += $num_possibilities;
  	}

  	close(IN);
  	return ($total_possibilities);
}

sub sum {

	my @arr = @_;
	my $sum = 0;
	foreach(@arr){$sum+=$_}
	return $sum;
}

sub normalize {

	my @arr = @_;
	my $sum = sum(@arr);
	my @new_arr;
	if ($sum eq 0) {
		return @arr;
	}
	foreach (@arr) {
		push(@new_arr,$_/$sum);
	}
	return @new_arr;
}

sub denormalize {
    my ($target_sum, @arr) = @_;
    
    my @new_arr;
    foreach (@arr) {
        push(@new_arr, $_ * $target_sum);
    }
    return @new_arr;
}

sub predict_relationship {
	
	my @probs = @_;
	my @possibilities;
	
	for (my $i = 0; $i < 6; $i++) {
		if (@probs[$i] > $MIN_LIKELIHOOD) {							
			push(@possibilities, @likelihood_names[$i]);
		}
	}
	return @possibilities;
}

sub get_maximum_relationship {
	
	my @vector = @_;
	my $max_val = 0;
	my $max_pos = -1;
	for (my $i = 0; $i < @vector; $i++) {
		if (@vector[$i] >= $max_val) {
			$max_val = @vector[$i];
			$max_pos = $i;
		}
	}
	if ($max_val < $MIN_LIKELIHOOD) {
		return "FAIL"
	}
	my $max_rel = @likelihood_names[$max_pos];
	return $max_rel;
}

sub get_data_arrays {

	my $file = shift;

	if($verbose > 0) {
		print "Loading KDE look-up-file = $file\n"
	}
	my @data_array;

	for (my $i = 0; $i <= $KDE_density_resolution; $i++) {
		for (my $j = 0; $j <= $KDE_density_resolution; $j++) {
			$data_array[$i][$j] = 0;
		}
	}
	
	open(IN,$file) or die "$file $!";
	
	my $k1 = 0;
	while (my $line = <IN>) {
		my @temp =  split(/\s+/,$line);
		#print "line = " . @temp[0] . " " . @temp[1] . "\n";
		my $k0 = 0;
		foreach (@temp) {
			#print "k1: $k1\t";
			$data_array[$k0][$k1] = $_;
			$k0++;
		}
		$k1++;
		#print "\n\n";
	}
	close(IN);
	return \@data_array;
}

return 1;
