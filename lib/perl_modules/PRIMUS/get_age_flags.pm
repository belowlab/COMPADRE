#! /usr/bin/perl

package PRIMUS::get_age_flags;
use strict;



my $MAX_GENERATION_TIME = 50;
my $MIN_GENERATION_TIME = 13;

my $MAX_FULL_SIBLING_GAP = 25;
my $MAX_HALF_SIBLING_GAP = 25;
my $MAX_AVUNCULAR_GAP = 50;
my $MIN_AVUNCULAR_GAP = -10;
my $MAX_COUSIN_GAP = 30;
my $AGE_FILE;
my $FID_COL;
my $IID_COL;
my $AGE_COL;
my %ages;




#########################
## Subs
########################
sub get_age_flags_in_network
{
	my $network_ref = shift;
	my $ages_ref = shift;
	%ages = %{$ages_ref};
	my $age_flags_ref = shift;
	
	#print "GETTING AGE FLAGS...\n";
	foreach(keys %ages)
	{
        #print "$_ = $ages{$_}\n";
	}
	
	my @ids = keys %$network_ref;
	#print "@ids\n";
	
	for(my $i = 0; $i < @ids - 1; $i++)
	{
		my $IID1 = @ids[$i];
		if($IID1 =~ /Dummy/i){next}
		if(!exists  $ages{$IID1} || $ages{$IID1} !~ /^[\d\.]+$/){next}
        my $age1 = $ages{$IID1};
		for(my $j = $i + 1; $j < @ids; $j++)
		{
			my $IID2 = @ids[$j];
            print "IID2: $IID2\n";
			if($IID2 =~ /Dummy/i){next}
			if(!exists  $ages{$IID2} || $ages{$IID2} !~ /^[\d\.]+$/){next}
			my $age2 = $ages{$IID2};
			my $self = $$network_ref{$IID1};
			my $rel_iid2_to_iid1 = which_relationship_exists_in_pedigree($self,$network_ref,$IID2);

            #print "\n$IID1($age1) -> $IID2($age2) = $rel_iid2_to_iid1\n";

			## Parent offspring tests
			if($rel_iid2_to_iid1 eq "P" && $age2 - $MIN_GENERATION_TIME < $age1)
			{
				#print "FLAG!!!! Parent too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Parent too young";
			}
			if($rel_iid2_to_iid1 eq "P" && $age2 - $MAX_GENERATION_TIME > $age1)
			{
				#print "FLAG!!!! Parent too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Parent too old";
			}
			if($rel_iid2_to_iid1 eq "O" && $age1 - $MIN_GENERATION_TIME < $age2)
			{
				#print "FLAG!!!! Parent too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Parent too young";
			}
			if($rel_iid2_to_iid1 eq "O" && $age1 - $MAX_GENERATION_TIME > $age2)
			{
				#print "FLAG!!!! Parent too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Parent too old";
			}
			
			## Full-sib test
			if($rel_iid2_to_iid1 eq "F" && abs($age1 - $age2) > $MAX_FULL_SIBLING_GAP)
			{
				#print "FLAG!!!! Siblings too far apart\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Siblings too far apart";
			}

			## Half-sib test
			if($rel_iid2_to_iid1 eq "H" && abs($age1 - $age2) > $MAX_HALF_SIBLING_GAP)
			{
				#print "FLAG!!!! Siblings too far apart\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Siblings too far apart";
			}
			
			## Grandparent/child tests
			if($rel_iid2_to_iid1 eq "G" && $age2 - ($MIN_GENERATION_TIME*2) < $age1)
			{
				#print "FLAG!!!! grandParent too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="grandParent too young";
			}
			if($rel_iid2_to_iid1 eq "G" && $age2 - ($MAX_GENERATION_TIME*2) > $age1)
			{
				#print "FLAG!!!! grandParent too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="grandParent too old";
			}
			if($rel_iid2_to_iid1 eq "C" && $age1 - ($MIN_GENERATION_TIME*2) < $age2)
			{
				#print "FLAG!!!! grandParent too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="grandParent too young";
			}
			if($rel_iid2_to_iid1 eq "C" && $age1 - ($MAX_GENERATION_TIME*2) > $age2)
			{
				#print "FLAG!!!! grandParent too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="grandParent too old";
			}
			
			## Avuncular tests
			if($rel_iid2_to_iid1 eq "A" && $age2 - ($MIN_AVUNCULAR_GAP) < $age1)
			{
				#print "FLAG!!!! greatAvuncular too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Avuncular too young";
			}
			if($rel_iid2_to_iid1 eq "A" && $age2 - ($MAX_AVUNCULAR_GAP) > $age1)
			{
				#print "FLAG!!!! greatAvuncular too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Avuncular too old";
			}
			if($rel_iid2_to_iid1 eq "N" && $age1 - ($MIN_AVUNCULAR_GAP) < $age2)
			{
				#print "FLAG!!!! greatAvuncular too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Avuncular too young";
			}
			if($rel_iid2_to_iid1 eq "N" && $age1 - ($MAX_AVUNCULAR_GAP) > $age2)
			{
				#print "FLAG!!!! greatAvuncular too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Avuncular too old";
			}

			## Cousin test
			if($rel_iid2_to_iid1 eq "1C" && abs($age1 - $age2) > $MAX_COUSIN_GAP)
			{
				#print "FLAG!!!! Siblings too far apart\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="Cousins too far apart";
			}

			## Great-grandparent/child tests
			if($rel_iid2_to_iid1 eq "GG" && $age2 - ($MIN_GENERATION_TIME*3) < $age1)
			{
				#print "FLAG!!!! greatgrandParent too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="greatgrandParent too young";
			}
			if($rel_iid2_to_iid1 eq "GG" && $age2 - ($MAX_GENERATION_TIME*3) > $age1)
			{
				#print "FLAG!!!! greatgrandParent too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="greatgrandParent too old";
			}
			if($rel_iid2_to_iid1 eq "GC" && $age1 - ($MIN_GENERATION_TIME*3) < $age2)
			{
				#print "FLAG!!!! greatgrandParent too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="greatgrandParent too young";
			}
			if($rel_iid2_to_iid1 eq "GC" && $age1 - ($MAX_GENERATION_TIME*3) > $age2)
			{
				#print "FLAG!!!! greatgrandParent too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="greatgrandParent too old";
			}
			
			## Great-avuncular tests
			if($rel_iid2_to_iid1 eq "GA" && $age2 - ($MIN_GENERATION_TIME*2) < $age1)
			{
				#print "FLAG!!!! greatAvuncular too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="greatAvuncular too young";
			}
			if($rel_iid2_to_iid1 eq "GA" && $age2 - ($MAX_GENERATION_TIME*2) > $age1)
			{
				#print "FLAG!!!! greatAvuncular too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="greatAvuncular too old";
			}
			if($rel_iid2_to_iid1 eq "GN" && $age1 - ($MIN_GENERATION_TIME*2) < $age2)
			{
				#print "FLAG!!!! greatAvuncular too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="greatAvuncular too young";
			}
			if($rel_iid2_to_iid1 eq "GN" && $age1 - ($MAX_GENERATION_TIME*2) > $age2)
			{
				#print "FLAG!!!! greatAvuncular too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="greatAvuncular too old";
			}

			## Half-avuncular tests
			if($rel_iid2_to_iid1 eq "HA" && $age2 - ($MIN_AVUNCULAR_GAP) < $age1)
			{
				#print "FLAG!!!! greatAvuncular too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="halfAvuncular too young";
			}
			if($rel_iid2_to_iid1 eq "HA" && $age2 - ($MAX_AVUNCULAR_GAP) > $age1)
			{
				#print "FLAG!!!! greatAvuncular too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="halfAvuncular too old";
			}
			if($rel_iid2_to_iid1 eq "HN" && $age1 - ($MIN_AVUNCULAR_GAP) < $age2)
			{
				#print "FLAG!!!! greatAvuncular too young\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="halfAvuncular too young";
			}
			if($rel_iid2_to_iid1 eq "HN" && $age1 - ($MAX_AVUNCULAR_GAP) > $age2)
			{
				#print "FLAG!!!! greatAvuncular too old\n";
				$$age_flags_ref{$network_ref}{$IID1}{$IID2}="halfAvuncular too old";
			}
		}
	}
}

sub load_age_file
{
	print "LOAD_AGE_FILE SHOULD NO LONGER BE USED\n";
	exit; 
	print "Loading ages...\n";

	$AGE_FILE = shift;
	$FID_COL = shift;
	$IID_COL = shift;
	$AGE_COL = shift;

	if($FID_COL eq ""){$FID_COL = 0}
	if($IID_COL eq ""){$IID_COL = 1}
	if($AGE_COL eq ""){$AGE_COL = 2}

	open(IN,$AGE_FILE) or die "Can't open age file $AGE_FILE; $!\n";

	while(my $line = <IN>)
	{
		my @temp = split(/\s+/,$line);
		if(@temp[$AGE_COL] < 0 || @temp[$AGE_COL] eq "NA"){next}
    #$ages{"@temp[$FID_COL]__@temp[$IID_COL]"} = @temp[$AGE_COL];
		$ages{"@temp[$IID_COL]"} = @temp[$AGE_COL];
	}

	return \%ages;
}


sub which_relationship_exists_in_pedigree
{
	my $self = shift;
	my $network_ref = shift;
	my $rel = shift;
	my $phase = 10;
	#if(!exists $$network_ref{$self}){return "UN"}
	if(!exists $$network_ref{$rel}){return "UN"}
	my $self_name = $self->name();


	my @parents = @{$self->{PARENTS} };
	my @children = @{$self->{CHILDREN} };
	
	#print "Self_name: $self_name\n";
	#print "parents: @parents\n";
	#print "children: @children\n";
	
	#print "rel: $rel\n";

	## PC
	if(grep($_ eq $rel, @parents)){return "P";}
	if(grep($_ eq $rel, @children)){return "O";}
	
	if(@parents ne 0)
	{
		## FS
		my @P1_children = $$network_ref{@parents[0]}->children();
		my @P2_children = $$network_ref{@parents[1]}->children();	 
		foreach(@P1_children)
		{
			my $child = $_;
			if(grep($_ eq $child,@P2_children)) #full sibling
			{
				if($child eq $self_name){next;}
				if($child eq $rel){return "F";}
			}
		}
		
		## Check if $rel is a neice/nephew
		my @P1_children = $$network_ref{@parents[0]}->children();
		my @P2_children = $$network_ref{@parents[1]}->children();
		foreach(@P1_children)
		{
			my $child = $_;
			if(grep($_ eq $child,@P2_children)) #full sibling
			{
				if($child eq $self_name){next;}
				my @P1_children = $$network_ref{$child}->children();
				foreach(@P1_children)
				{
					if($_ eq $rel){return "N";}
				}
			}
		}
		
		## check if rel is half-sib
		foreach(@P1_children)
		{
			my $child = $_;
			if(!grep($_ eq $child,@P2_children)) # Half sib
			{
				if($child eq $rel){return "H";}
			}
		}
		foreach(@P2_children)
		{
			my $child = $_;
			#print "child $child\n";
			if(!grep($_ eq $child,@P1_children)) # Half sib
			{
				if($child eq $rel){return "H";}
			}
		}
				
		## check if rel is a grandparent
			my @P1_parents = $$network_ref{@parents[0]}->parents();
		my @P2_parents = $$network_ref{@parents[1]}->parents();
			foreach(@P1_parents)
		{
			if($_ eq $rel){return "G";}
		}
			foreach(@P2_parents)
		{
			if($_ eq $rel){return "G";}
		}
		
		#print "self: $self_name\n";
		#print "parents: @parents\n";
		#print "@parents[0] parents: @P1_parents\n";
		#print "@parents[1] parents: @P2_parents\n";

		## Check if rel is an uncle/aunt
		if(@P1_parents > 0)
		{
			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();		 
				foreach(@G1_1_children)
			{
				my $child = $_;
				#print "child = $child; rel = $rel\n";
				if(grep($_ eq $child,@G1_2_children)) #uncle/aunt
				{
					if($child eq @parents[0]){next;}
					if($child eq $rel){return "A";}
				}
			}
		}
		if(@P2_parents > 0)
		{
			my @G2_1_children = $$network_ref{@P2_parents[0]}->children();
			my @G2_2_children = $$network_ref{@P2_parents[1]}->children();
				
				foreach(@G2_1_children)
			{
				my $child = $_;
				#print "child = $child; rel = $rel\n";
				if(grep($_ eq $child,@G2_2_children)) #Unclde/aunt
				{
					if($child eq @parents[1]){next;}
					if($child eq $rel){return "A";}
				}
			}
		}
	}
	
	## check if rel is grandchild
	foreach(@children)
	{
		my @G_children = $$network_ref{$_}->children();
		foreach(@G_children)
		{
			if($_ eq $rel){return "C";}
		}

	}

	## Phase 3 below

	## Check if rel is great grandchild
	foreach(@children)
	{
		my @G_children = $$network_ref{$_}->children();
		foreach(@G_children)
		{
			my @GG_children = $$network_ref{$_}->children();
			foreach(@GG_children)
			{
			if($_ eq $rel){return "GC";}
			}
		}
	}
	
	if(@parents eq 0){return "UN";} ## I think this might be useless because only dummies don't have parents
	
	my @P1_children = $$network_ref{@parents[0]}->children();
	my @P2_children = $$network_ref{@parents[1]}->children();
	my @P1_parents = $$network_ref{@parents[0]}->parents();
	my @P2_parents = $$network_ref{@parents[1]}->parents();
	my @AUNTS_UNCLES = ();
	
	## check if rel is Half niece/nephew = child of half-sib
	foreach my $h_sib (@P1_children)
	{
		if(!grep($_ eq $h_sib,@P2_children)) # Half sib
		{
			my @h_sib_children = $$network_ref{$h_sib}->children();
			foreach(@h_sib_children)
			{
				if($_ eq $rel){return "HN";}
			}
		}
	}
	foreach my $h_sib(@P2_children)
	{
		if(!grep($_ eq $h_sib,@P1_children)) # Half sib
		{
			my @h_sib_children = $$network_ref{$h_sib}->children();
			foreach(@h_sib_children)
			{
				if($_ eq $rel){return "HN";}
			}
		}
	}
	
	## check if rel is half uncle/aunt = half sibling of parent
	if(@P1_parents > 0)
	{
		my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
		my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
		#print "P1 @P1_parents[0] children: @G1_1_children\n";
		#print "P1 @P1_parents[1] children: @G1_2_children\n";
	 
			foreach(@G1_1_children)
		{
			my $child = $_;
			if(!grep($_ eq $child,@G1_2_children)) # half uncle/aunt
			{
				#print "$child : $self_name\n";
				if($_ eq $rel){return "HA";}
			}
			elsif($child ne @parents[0])
			{
				push(@AUNTS_UNCLES,$child);
			}
		}
		foreach(@G1_2_children)
		{
			my $child = $_;
			if(!grep($_ eq $child,@G1_1_children)) # half uncle/aunt
			{
				#print "$child : $self_name\n";
				if($_ eq $rel){return "HA";}    				
			}
			elsif($child ne @parents[0])
			{
				push(@AUNTS_UNCLES,$child);
			}
		}
	}
	if(@P2_parents > 0)
	{
		my @G2_1_children = $$network_ref{@P2_parents[0]}->children();
		my @G2_2_children = $$network_ref{@P2_parents[1]}->children();
		#print "P2 @P2_parents[0] children: @G2_1_children\n";
		#print "P2 @P2_parents[1] children: @G2_2_children\n";
		
		foreach(@G2_1_children)
		{
			my $child = $_;
			if(!grep($_ eq $child,@G2_2_children)) #half Unclde/aunt
			{
				if($_ eq $rel){return "HA";}    				
			}
			elsif($child ne @parents[1])
			{
				push(@AUNTS_UNCLES,$child);
			}
		}
		foreach(@G2_2_children)
		{
			my $child = $_;
			if(!grep($_ eq $child,@G2_1_children)) #half Unclde/aunt
			{
				if($_ eq $rel){return "HA";}
			}
			elsif($child ne @parents[1])
			{
				push(@AUNTS_UNCLES,$child);
			}
		}
	}
	
	## Check if rel is great grand parent
	my @grandparents;
	push(@grandparents, @P1_parents);
	push(@grandparents, @P2_parents);
	#my @P1_parents = $$network_ref{@parents[0]}->parents();
	#my @P2_parents = $$network_ref{@parents[1]}->parents();
	#print "@parents[0] parents: @P1_parents\n";
	#print "@parents[1] parents: @P2_parents\n";
	
	#print "$self_name grandparents: @grandparents\n";
	foreach my $G_parent (@grandparents)
	{
		my @GG_parents = $$network_ref{$G_parent}->parents();
		foreach(@GG_parents)
		{
			if($_ eq $rel){return "GG";}
		}
	}
		
		
	## Check if rel is first cousin
	my @grandparents0 = $$network_ref{$parents[0]}->parents();
	my @grandparents1 = $$network_ref{$parents[1]}->parents();
	
	my @rel_parents = $$network_ref{$rel}->parents();
	if(@rel_parents ne 0)
	{
		my @rel_grandparents0 = $$network_ref{$rel_parents[0]}->parents();
		my @rel_grandparents1 = $$network_ref{$rel_parents[1]}->parents();

		if(do_arrays_match(\@grandparents0,\@rel_grandparents0) && @grandparents0 > 0){return "1C"}
		elsif(do_arrays_match(\@grandparents0,\@rel_grandparents1) && @grandparents0 > 0){return "1C"}
		elsif(do_arrays_match(\@grandparents1,\@rel_grandparents0) && @grandparents1 > 0){return "1C"}
		elsif(do_arrays_match(\@grandparents1,\@rel_grandparents1) && @grandparents1 > 0){return "1C"}
	} 
		
				
	## Check if rel is great avuncular
	foreach my $G_parent (@grandparents)
	{

			my @GG_parents = $$network_ref{$G_parent}->parents();
			#print "G_parent $G_parent parents: @GG_parents\n";
			if(@GG_parents eq 0){next;}
			my @GG1_children = $$network_ref{@GG_parents[0]}->children();
		my @GG2_children = $$network_ref{@GG_parents[1]}->children();	 
			foreach my $child (@GG1_children)
		{
			if(grep($_ eq $child,@GG2_children)) #full sibling of grandparent
			{
				#print "here\n";
				if($child eq $G_parent){next;}
				if($child eq $rel){return "GA";}
			}
		}
	}
	
	## Check if rel is great nephew/niece
	foreach(@P1_children)
	{
		my $child = $_;
		if(grep($_ eq $child,@P2_children)) #full sibling
		{
			if($child eq $self_name){next;}
			my @FS_children = $$network_ref{$child}->children();
			foreach my $FS_child (@FS_children)
			{
				my @FS_Gchildren = $$network_ref{$FS_child}->children();
				foreach(@FS_Gchildren)
				{
					if($_ eq $rel){return "GN";}
				}
			}
		}
	}	
	return "UN";
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



return 1;
