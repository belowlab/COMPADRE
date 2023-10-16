#! /usr/bin/perl

use lib '../lib/perl_modules';

use strict;
use PRIMUS::node_v7;

my $output_dir = shift;
my $dataset_name = shift;


my $rels_ref = get_possible_relationships($output_dir,"$output_dir/Summary_$dataset_name.txt");
write_table("$output_dir/Summary_$dataset_name\_pairwise_table.txt",$rels_ref);

my %assembled_networks;
my %genders;

sub get_possible_relationships
{
	#print "loading possible relationships...\n";
	
	my $dataset_dir = shift;
	my $summary_file = shift;
	
	my %networks;
	my %id1_is_rels_to_id2;
	
	open(IN,$summary_file) or die "can't open $summary_file; $!\n";
	<IN>;
	while(my $line = <IN>)
	{
		if($line !~ /^\d/){next;}
		my ($num, $name, @temp) = split(/\s+/,$line);
		$networks{$num}{'name'}=$name;
		$networks{$num}{'line'}=$line;
		$networks{$num}{'ids'}=@temp[4];
	}
	close(IN);

	## Check that there are networks in the summary file
	if(!exists $networks{1}){die "No networks generated in $summary_file\n";}

	## For each network
	foreach my $num (sort {$a <=> $b} keys %networks)
	{
		my $net = $networks{$num}{'name'};
		#print "net: $net\n";
		
		my @ids = split(',',$networks{$num}{'ids'});
		#print "ids: @ids\n";
		
		## Check if there are IDs in the network
		if(@ids < 1){warn "$net has no IDs\n";next}

		for(my $i = 0; $i < @ids; $i++)
		{
			for(my $j = 0; $j < @ids; $j++)
			{
				if($i == $j){next;}
				my $IID1 = @ids[$i];
				my $IID2 = @ids[$j];
				my $rel1 = "";
				my $rel2 = "";
				
				my $all_relationships_ref = get_directed_relationships($dataset_dir,$net,$IID1,$rel1,$IID2,$rel2);
				my @rels_of_IID2_to_IID1 = sort {$a cmp $b} keys %$all_relationships_ref;
				#print "$IID2 is @rels_of_IID2_to_IID1 of $IID1\n";
				$id1_is_rels_to_id2{$num}{$IID2}{$IID1} = \@rels_of_IID2_to_IID1;
			}
		}
	}
	foreach(keys %id1_is_rels_to_id2)
	{
		#print "key: $_\n";
	}	
	return \%id1_is_rels_to_id2;
}



sub get_directed_relationships
{
	my $dataset_dir = shift;
	my $net = shift;
	my $IID1 = shift;
	my $rel1 = shift;
	my $IID2 = shift;
	my $rel2 = shift;
	my $results_dir = "$dataset_dir/$net";
	my $num_possible_peds = 1;
	my %all_rels;

	#print "file: $results_dir/$net\_$num_possible_peds.ped\n";
	while(-e "$results_dir/$net\_$num_possible_peds.fam")
	{
		$num_possible_peds++;
	}
	$num_possible_peds--;
	#print "$IID1 <-> $IID2 $net = $num_possible_peds\n";
	
	for(my $i = 1; $i <= $num_possible_peds; $i++)
	{
		my $ped_file = "$results_dir/$net\_$i.fam";
		#print "i: $i\n";

		## Get network_ref
		my $network_ref;
		if(exists $assembled_networks{$net}{$i})
		{
			$network_ref = $assembled_networks{$net}{$i};
		}
		else
		{
			$network_ref = build_network_from_ped_file($ped_file);
			$assembled_networks{$net}{$i} = $network_ref;
		}
		
		## Check relationship
		my $self = $$network_ref{$IID1};
		my $rel = which_relationship_exists_in_pedigree($self,$network_ref,$IID2);
		#print "rel: $rel\n";
		$all_rels{$rel}++;
		if($rel ne $rel2)
		{
			#### Write this next line to file to keep track of all the specific reconstructed pedigrees that don't match a given relationship
			#print "$pop $net $i $IID1 ($rel1) != $IID2 ($rel2)\n";
		}
	}

	return \%all_rels;
}


sub which_relationship_exists_in_pedigree
{
	my $self = shift;
	my $network_ref = shift;
	my $rel = shift;
	my $phase = 10;
	#print "rel: $rel\n";
	my $self_name = $self->name();
	#print "Self_name: $self_name\n";
	#if(!exists $$network_ref{$self}){return "UN"}
	if(!exists $$network_ref{$rel}){return "UN"}


	my @parents = @{$self->{PARENTS} };
	my @children = @{$self->{CHILDREN} };
	
	#print "parents: @parents\n";
	#print "children: @children\n";
	

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


sub build_network_from_ped_file
{
	my $ped_file = shift;
	## Read in file
	open(IN,$ped_file) or die "Cannot open $ped_file; $!\n";
	my %all_nodes_network;
	my $network_ref = \%all_nodes_network;
	my @network_refs;
	
	## Build pedigree
	while(my $line = <IN>)
	{
		chomp($line);
		my ($FID,$IID,$PID,$MID,$SEX,$PHENOTYPE) = split(/\s+/,$line);
		#print "before: $IID:$MID:$PID\n";
		#$IID =~ s/\w+__//;
		#$MID =~ s/\w+__//;
		#$PID =~ s/\w+__//;
		#print "after: $IID:$MID:$PID\n";
		
		my $child = $IID;
		my $mom = $MID;
		my $dad = $PID;
		$genders{$child}=$SEX;
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
	
	## Write out relationships
	foreach my $node_name (keys %$network_ref)
	{
		$$network_ref{$node_name}->make_relative_network_from_pedigree($network_ref);
		my %rels = $$network_ref{$node_name}->relatives();
		foreach(keys %rels)
		{
			#print "$node_name -> $_  = @{$rels{$_} }\n";
		}
	}
	
	return $network_ref;
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
1;

sub write_table
{
	my $relationship_table_file = shift;
	my $rels_ref = shift;
	
	foreach(keys %$rels_ref)
	{
		#print "key: $_\n";
	}
	
	open(TABLE,">$relationship_table_file") or die "can't open $relationship_table_file;$!\n";
	print TABLE "Network\tID1\tID1_sex\trelationship_to_ID2\tID2\tID2_sex\trelationship_to_ID1\n";

	foreach my $net (sort {$a <=> $b} keys %$rels_ref )
	{
		my %id1_is_rels_to_id2 = %{$$rels_ref{$net} };
		#print "net: $net\n";
		my @ids = sort {$a cmp $b} keys %id1_is_rels_to_id2;
		#print "ids: @ids\n";

		for(my $i = 0; $i < @ids - 1; $i++)
		{
			my $id1 = @ids[$i];
			if(!exists $id1_is_rels_to_id2{$id1})
			{
				print "$id1 is missing!!!\n";
				next;
			}
			for(my $j = $i+1 ; $j < @ids; $j++)
			{
				my $id2 = @ids[$j];
				if(!exists $id1_is_rels_to_id2{$id1}{$id2})
				{
					print "$id2 is missing!!!\n";
					next;
				}
				my $gender1 = $genders{$id1};
				my $gender2 = $genders{$id2};
				if($gender1 eq 1){$gender1 = "M"}
				elsif($gender1 eq 2){$gender1 = "F"}
				else{$gender1 = "0"}
				
				if($gender2 eq 1){$gender2 = "M"}
				elsif($gender2 eq 2){$gender2 = "F"}
				else{$gender2 = "0"}


				my @one2two_rels = @{$id1_is_rels_to_id2{$id1}{$id2} };
				my @two2one_rels = @{$id1_is_rels_to_id2{$id2}{$id1} };
				if(@one2two_rels < 1){@one2two_rels[0] = "?"}
				if(@two2one_rels < 1){@two2one_rels[0] = "?"}

				if(@one2two_rels[0] eq "UN" || @one2two_rels[0] eq "?"){next}

				#print "$net\t$id1\t$gender1\t".join(",",@one2two_rels)."\t$id2\t$gender2\t".join(",",@two2one_rels)."\n";
				print TABLE "$net\t$id1\t$gender1\t".join(",",@one2two_rels)."\t$id2\t$gender2\t".join(",",@two2one_rels)."\n";
			}
		}
	}

	print TABLE "\n\n\n";
print TABLE "
# RELATIONSHIP KEY:
# 1st degree:
# P = parent
# O = offspring
# F = full-sib
#
# 2nd degree:
# A = uncle/aunt
# N = niece/nephew
# H = half-sib
# G = grandparent
# C = grandchild
#
# 3rd degree:
# GA = great-uncle/aunt
# GN = great-niece/nephew
# HA = half-uncle/aunt
# HN = half-niece/nephew
# 1C = first cousin
# GG = great-grandparent
# GC = great=grandchild\n";
	close(TABLE);
}

