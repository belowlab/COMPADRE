package PRIMUS::compare_fam_files;

use strict;
#use lib "../lib/perl_modules";
use PRIMUS::node_v7;

my %sexes;

sub get_correct_pedigree
{
	#my $missing = shift;
	my $reconstruction = 0;
	
	my $ref_fam_file = shift;
	my $network_dir = shift;
	my $pedigree_name = shift;
	
	my %ref = load_fam_file($ref_fam_file);

	my @correct;
	while(-e "$network_dir/$pedigree_name\_reconstructed_$reconstruction.fam")
	{
		
		my %test = load_fam_file("$network_dir/$pedigree_name\_reconstructed_$reconstruction.fam");
		if(are_fams_same(\%ref,\%test))
		{
			#print " correct = size12-$missing\_reconstructed_$reconstruction.fam\n";
			push(@correct,$reconstruction);
		}
		else
		{
			#print "NOT SAME\n";
		}
		$reconstruction++;
		#last;
	}
	my $num_correct = @correct;
	#print "size12-$missing $num_correct correct: @correct out of $reconstruction\n";
	return ($num_correct,$reconstruction,@correct);
}

sub get_correct_pedigree_simulation
{
	my $missing = shift;
	my $reconstruction = 0;

	my $fam_file_1 = "../data/simulations/size12/size12-$missing.fam";
	my %ref = load_fam_file($fam_file_1);

	my @correct;
	while(-e "../data/simulations/size12/size12-$missing\_reconstructed_$reconstruction.fam")
	{
		
		my %test = load_fam_file("../data/simulations/size12/size12-$missing\_reconstructed_$reconstruction.fam");
		if(are_fams_same(\%ref,\%test))
		{
			#print " correct = size12-$missing\_reconstructed_$reconstruction.fam\n";
			push(@correct,$reconstruction);
		}
		else
		{
			#print "NOT SAME\n";
		}
		$reconstruction++;
		#last;
	}
	my $num_correct = @correct;
	#print "size12-$missing $num_correct correct: @correct out of $reconstruction\n";
	return ($num_correct,$reconstruction,@correct);
}

##################################################################################
## Subs ##########################################################################
##################################################################################

sub are_fams_same
{
	my $fam1 = shift;
	my $fam2 = shift;
	
	my $net_ref1 = make_network_from_fam($fam1);
	my $net_ref2 = make_network_from_fam($fam2);

	my $ctr = 0;
	my %Missing_conversion;
	foreach my $IID (sort {$$fam1{$a} cmp $$fam1{$b} } keys %$fam1)
	{
		## Don't look at dummies directly.
		if($IID =~ /Missing/i){next;}
		#print "\nIID: $IID\n";
		my $val = do_samples_match($fam1,$fam2,$IID,$IID,\%Missing_conversion,$net_ref1,$net_ref2);
		#print "Main loop val: $val\n";
		if($val == 0)
		{
			#print "NOT EQUAL\n";
			return 0;
		}
	}
	#print "EQUAL\n";
	return 1;
}

sub do_samples_match
{
	my $fam1 = shift;
	my $fam2 = shift;
	my $IID1 = shift;
	my $IID2 = shift;
	my $Missing_conversion = shift;
	my $network1 = shift;
	my $network2 = shift;

	## if they are missing samples, check that they haven't already been assigned a conversion
	if($IID1 =~ /Missing/i && $IID2 =~ /Missing/i)
	{
		#print "Missing1: $IID1; Missing2: $IID2\n";
		if(exists $$Missing_conversion{1}{$IID1})
		{
			my $d1 = $$Missing_conversion{1}{$IID1};
			if($d1 ne $IID2)
			{
				#print "(1) Dummies already found and are not equal: $IID1 = $d1 not $IID2\n";
				return 0;
			}
		}
		if(exists $$Missing_conversion{2}{$IID2})
		{
			my $d2 = $$Missing_conversion{2}{$IID2};
			if($d2 ne $IID1)
			{
				#print "(2) Dummies already found and are not equal: $IID2 = $d2 not $IID1\n";
				return 0;
			}
		}
	}

	my $PID1 = $$fam1{$IID1}{"PID"};
	my $MID1 = $$fam1{$IID1}{"MID"};
	my $PID2 = $$fam2{$IID2}{"PID"};
	my $MID2 = $$fam2{$IID2}{"MID"};
	my @children1 = $$network1{$IID1}->children();
	my @children2 = $$network2{$IID2}->children();


	#print "IID1: $IID1 ($PID1 $MID1); IID2: $IID2 ($PID2 $MID2)\n";
	#print "IID1 children: @children1\n";
	#print "IID2 children: @children2\n";
	#print "PID1: $PID1; MID1: $MID1\n";
	#print "PID2: $PID2; MID2: $MID2\n";

	## Check that the children match;
	## This is only necessary to check for the parents of half-siblings
	## This is also not very well tested, and an incomplete test. It only makes sure that the non-dummy children match and that the number of dummy children are the same. It does not compare the dummy children.
	## Also, if you move this to after you check the parents, it will not work because the wrong dummy parent conversion will get set.
	my $val = do_children_match($fam1,$fam2,\@children1,\@children2,$Missing_conversion,$network1,$network2);
	if($val == 0)
	{
		#print "Children don't match: @children1 <=> @children2\n";
		return 0;
	}
	else
	{
		#print "Children match2: @children1 <=> @children2\n";
	}

	my @dummies;
	my @non_dummies;

	if($PID1 =~ /Missing/i){push(@dummies,$PID1);}else{push(@non_dummies,$PID1);}
	if($PID2 =~ /Missing/i){push(@dummies,$PID2);}else{push(@non_dummies,$PID2);}
	if($MID1 =~ /Missing/i){push(@dummies,$MID1);}else{push(@non_dummies,$MID1);}
	if($MID2 =~ /Missing/i){push(@dummies,$MID2);}else{push(@non_dummies,$MID2);}

	## 1. all 4 are dummies => I need to check PID1 to PID2/MID2 and MID1 to PID2/MID2
	## 2. 3 are dummies => fail
	## 3. 2 are dummies => check that the two dummies are the same
	## 4. 1 is Missing => fail
	## 5. 0 are dummies => check that MID1 = MID2 and PID1 = PID2
	if(@dummies == 4)
	{
		my $val1 = do_samples_match($fam1,$fam2,$PID1,$PID2,$Missing_conversion,$network1,$network2);
		my $val2 = do_samples_match($fam1,$fam2,$MID1,$MID2,$Missing_conversion,$network1,$network2);
		
		my $val3 = do_samples_match($fam1,$fam2,$PID1,$MID2,$Missing_conversion,$network1,$network2);
		my $val4 = do_samples_match($fam1,$fam2,$MID1,$PID2,$Missing_conversion,$network1,$network2);

		if($val1 && $val2)
		{
			#print "1+2 $val1 == $val2\n";
			$$Missing_conversion{1}{$PID1} = $PID2;
			$$Missing_conversion{2}{$PID2} = $PID1;
			$$Missing_conversion{1}{$MID1} = $MID2;
			$$Missing_conversion{2}{$MID2} = $MID1;
		}
		elsif($val3 && $val4)
		{
			#print "3+4 $val3 == $val4\n";
			$$Missing_conversion{1}{$PID1} = $MID2;
			$$Missing_conversion{2}{$MID2} = $PID1;
			$$Missing_conversion{1}{$MID1} = $PID2;
			$$Missing_conversion{2}{$PID2} = $MID1;
		}
		else
		{
			#print "Val1/2/3/4: $val1/$val2/$val3/$val4\n";
			return 0;
		}
	}
	elsif(@dummies == 3)
	{
		return 0;
	}
	elsif(@dummies == 2)
	{
		if($PID1 =~ /Missing/i && $PID2 =~ /Missing/i)
		{
			if($MID1 ne $MID2)
			{
				#print "MIDs $MID1 and $MID2 do not match\n";
				return 0;
			}
			if(!do_samples_match($fam1,$fam2,$PID1,$PID2,$Missing_conversion,$network1,$network2))
			{
				return 0;
			}
			$$Missing_conversion{1}{$PID1} = $PID2;
			$$Missing_conversion{2}{$PID2} = $PID1;
		}
		elsif($MID1 =~ /Missing/i && $MID2 =~ /Missing/i)
		{
			if($PID1 ne $PID2)
			{
				#print "PIDs $PID1 and $PID2 do not match\n";
				return 0;
			}
			if(!do_samples_match($fam1,$fam2,$MID1,$MID2,$Missing_conversion,$network1,$network2))
			{
				return 0;
			}
			$$Missing_conversion{1}{$MID1} = $MID2;
			$$Missing_conversion{2}{$MID2} = $MID1;

		}
		elsif($MID1 =~ /Missing/i && $PID2 =~ /Missing/i)
		{
			if($PID1 ne $MID2)
			{
				#print "PID $PID1 and MID $MID2 do not match\n";
				return 0;
			}
			if(!do_samples_match($fam1,$fam2,$MID1,$PID2,$Missing_conversion,$network1,$network2))
			{
				return 0;
			}
			$$Missing_conversion{1}{$MID1} = $PID2;
			$$Missing_conversion{2}{$PID2} = $MID1;

		}
		elsif($PID1 =~ /Missing/i && $MID2 =~ /Missing/i)
		{
			if($MID1 ne $PID2)
			{
				#print "MID $MID1 and PID $PID2 do not match\n";
				return 0;
			}
			if(!do_samples_match($fam1,$fam2,$PID1,$MID2,$Missing_conversion,$network1,$network2))
			{
				return 0;
			}
			$$Missing_conversion{1}{$PID1} = $MID2;
			$$Missing_conversion{2}{$MID2} = $PID1;
		}
		else
		{
			return 0;
		}
	}
	elsif(@dummies == 1)
	{
		return 0;
	}
	else # @dummies == 0
	{
		if(@dummies > 0){die "Dummies should be empty: @dummies\n";}
		if($PID1 eq $PID2 && $MID1 eq $MID2)
		{
			#print "PIDs $PID1 and $PID2 match\n";
			#print "PIDs $MID1 and $MID2 match\n";
			#print "PASS = IID1: $IID1 ($PID1 $MID1); IID2: $IID2 ($PID2 $MID2)\n";
			return 1;
		}
		elsif($PID1 eq $MID2 && $MID1 eq $PID2)
		{
			## Should check that the sex of the samples are unknown before continuing
			if($sexes{$PID1} eq 1 || $sexes{$PID2} eq 1 || $sexes{$MID1} eq 2 || $sexes{$MID2} eq 2)
			{
				#print "Sex checking failed\n";
				#print "IID1: $IID1 ($PID1 $MID1); IID2: $IID2 ($PID2 $MID2) do not match\n";
				return 0;
			}
			
			#print "PIDs $PID1 and $MID2 match\n";
			#print "PIDs $MID1 and $PID2 match\n";
			#print "PASS = IID1: $IID1 ($PID1 $MID1); IID2: $IID2 ($PID2 $MID2)\n";
			return 1;
		}
		else
		{
			#print "IID1: $IID1 ($PID1 $MID1); IID2: $IID2 ($PID2 $MID2) do not match\n";
			return 0;
		}

	}
	

	#print "PASS = IID1: $IID1 ($PID1 $MID1); IID2: $IID2 ($PID2 $MID2)\n";
	return 1;
}

sub load_fam_file
{
	my $ctr = 2000;
	my $file = shift;
	open(IN,$file) or die "cannot open $file; $!\n";
	my %fam;
	my %all_samples;
	while(my $line = <IN>)
	{
		while($line =~ /\s0\s/)
		{
			$line =~ s/\s0\s/\tMissing$ctr\t/;
			$ctr++;
		}
		$line =~ s/\s0\n/\tMissing$ctr\n/;
		chomp($line);
		$ctr++;
		my ($FID,$IID,$PID,$MID,@rest) = split(/\s/,$line);
		#print "$FID,$IID,$PID,$MID\n";
		$fam{$IID}{"FID"} = $FID;
		$fam{$IID}{"PID"} = $PID;
		$fam{$IID}{"MID"} = $MID;
		$fam{$IID}{"sex"} = @rest[0];
		$fam{$IID}{"aff"} = @rest[1];
		$fam{$IID}{"geno"} = @rest[2];
		$all_samples{$IID} = 1;
		$all_samples{$PID} = 1;
		$all_samples{$MID} = 1;
	}
	close(IN);
	
	## If a parent is not and IID, add it
	foreach(keys %all_samples)
	{
		if(!exists $fam{$_})
		{
			$fam{$_}{'FID'} = "?";
			$fam{$_}{'PID'} = "Missing$ctr";
			$ctr++;
			$fam{$_}{'MID'} = "Missing$ctr";
			$ctr++;
		}
	}
	return %fam;
}

sub load_reconstructed_fam_file
{
	my $file = shift;
	open(IN,$file) or die "cannot open $file; $!\n";
	my @fams;
	my %fam;
	my $ctr = 5000;
	while(my $line = <IN>)
	{
		while($line =~ /\s0\s/)
		{
			$line =~ s/\s0\s/\tMissing$ctr\t/;
			$ctr++;
		}
		if($line =~ /\s0\n/)
		{
			$line =~ s/\s0\n/\tMissing$ctr\n/;
			$ctr++;
		}
		chomp($line);
		#print "line: $line\n";
    #my ($network,$name,$PID_long,$MID_long,$sex,@rest) = split(/\s+/,$line);
		my ($network,$name,$PID,$MID,$sex,@rest) = split(/\s+/,$line);
		my $FID = $name;
    my $IID=$name; 
    #if($name !~ /Missing/i && $name !~ /Missing/i)
    #{
      #($FID,$IID) = split(/__/,$name);
      #	$IID = $name

      #}
      #else
      #{
      #$IID = $name;
      #$FID = "?";
      #}
      #my ($PID,$MID,$junk);
      #if($PID_long !~ /Missing/i && $PID_long !~ /Missing/i)
      #{
      #($junk,$PID) = split(/__/,$PID_long);
      #}
      #else
      #{
      #$PID = $PID_long;
      #}
      #if($MID_long !~ /Missing/i && $MID_long !~ /Missing/i)
      #{
      #($junk,$MID) = split(/__/,$MID_long);
      #}
      #else
      #{
      #$MID = $MID_long;
      #}
		#print "$FID,$IID,$PID,$MID\n";
		$fam{$IID}{"FID"} = $FID;
		$fam{$IID}{"PID"} = $PID;
		$fam{$IID}{"MID"} = $MID;
		$sexes{$IID} = $sex;
	}
	close(IN);


	#push(@fams,\%fam);
	#my $network_ref = make_network_from_fam($hash);
	



	return %fam;
}

sub remove_unnecessary_missing
{
	my $hash = shift;

	my $network_ref = make_network_from_fam($hash);
	my $remove = 1;
	while($remove eq 1)
	{
		$remove = 0;
		foreach my $node_name (keys %$network_ref)
		{
			#if(!exists $$hash{$node_name}){next;}
			#my @parents = ($$hash{$node_name}{'PID'},$$hash{$node_name}{'MID'});
			#if(@parents eq 0 || @parents[0] !~ /Missing/ || @parents[1] !~ /Missing/){next;}

			#my @parents0 = ($$hash{@parents[0]}{'PID'},$$hash{@parents[0]}{'MID'});
			#my @parents1 = ($$hash{@parents[1]}{'PID'},$$hash{@parents[1]}{'MID'});
			#if(@parents0 > 0 || @parents1 > 0){next;}
			
			#print "\nnode: $node_name\n";

			if(!exists $$network_ref{$node_name}){next;}
			
			my @children = $$network_ref{$node_name}->children();
			my @parents = $$network_ref{$node_name}->parents();
			#print "parents: @parents\n";
			#print "children: @children\n";
			if($node_name =~ /Missing/i && @children < 1)
			{
				#print "here\n";
				if(@parents > 0)
				{
					$$network_ref{@parents[0]}->remove_child($node_name);
					$$network_ref{@parents[1]}->remove_child($node_name);
				}

				#my @children0 = $$network_ref{@parents[0]}->children();
				#print "children0: @children0\n";
				#my @children1 = $$network_ref{@parents[1]}->children();
				#print "children1: @children1\n";
				
				delete $$network_ref{$node_name};
				delete $$hash{$node_name};
				$remove = 1;
				next;
			}

			#print "parents: @parents\n";
			if(@parents eq 0 || @parents[0] !~ /Missing/ || @parents[1] !~ /Missing/){next;}
			
			my @parents0 = $$network_ref{@parents[0]}->parents();
			my @parents1 = $$network_ref{@parents[1]}->parents();
			#print "parents0: @parents0\n";
			#print "parents1: @parents1\n";
			if(@parents0 > 0 || @parents1 > 0){next;}
			
			my @children0 = $$network_ref{@parents[0]}->children();
			my @children1 = $$network_ref{@parents[1]}->children();
			#print "children0: @children0\n";
			#print "children1: @children1\n";
			if(@children0 > 1 || @children1 > 1){next;}
		
			## Passed all criteria, remove dummy parents
			$remove = 1;
			#print "**REMOVE** parents: @parents\n";
			$$network_ref{$node_name}->delete_parents();
			delete $$network_ref{@parents[0]};
			delete $$network_ref{@parents[1]};
			delete $$hash{@parents[0]};
			delete $$hash{@parents[1]};
		}
	}
}

sub make_network_from_fam
{
	my $hash = shift;
	my %network;
	my $network_ref = \%network;
	
	#print "here\n";

	foreach my $child (keys %$hash)
	{
		#print "child: $child\n";
		my $PID = $$hash{$child}{'PID'};
		my $MID = $$hash{$child}{'MID'};
		my $dad = $$hash{$child}{'PID'};
		my $mom = $$hash{$child}{'MID'};
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
	return $network_ref;
}

sub do_children_match
{
	
	#print "DO CHILDREN MATCH\n";
	my $fam1 = shift;
	my $fam2 = shift;
	my $child1_ref = shift;
	my $child2_ref = shift;
	my $Missing_conversion = shift;
	my $network1 = shift;
	my $network2 = shift;

	if(@$child1_ref ne @$child2_ref)
	{
		return 0;
	}

	my @dummy_samples1;
	my @dummy_samples2;

	foreach my $IID (@$child1_ref)
	{
		if($IID =~ /Missing/i)
		{
			push(@dummy_samples1,$IID);
		}
		elsif(!grep($_ eq $IID,@$child2_ref))
		{
			return 0;
		}
	}
	foreach my $IID (@$child2_ref)
	{
		if($IID =~ /Missing/i)
		{
			push(@dummy_samples2,$IID);
		}
		elsif(!grep($_ eq $IID,@$child1_ref))
		{
			return 0;
		}
	}

	#print "Dummy1: @dummy_samples1; Dummy2: @dummy_samples2\n";
	## Check that they are the same length
	if(@dummy_samples1 != @dummy_samples2)
	{
		return 0;
	}
	
	return 1;

	## For each dummy sample there must be a dummy sample in the other dummy sample set that matches
	foreach my $d1 (@dummy_samples1)
	{
		my $match = 0;
		foreach my $d2 (@dummy_samples2)
		{
			my $val = do_samples_match($fam1,$fam2,$d1,$d2,$Missing_conversion,$network1,$network2);
			if($val == 1)
			{
				$match = 1;
				last;
			}
		}
		if($match == 0)
		{
			return 0;
		}
	}

}


return 1;


