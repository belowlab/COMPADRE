package PRIMUS::node_v7;
use strict;

## v2 is changing the format the relationships are stored from Hash of values ('pc','fs',ect) to a hash of array references to a vector of likelihoods
## This requires changing the functions that make new nodes, store relatives, add relatives, and return relatives
## This will require additional functions to store CGH relationships
## Also requires changing the checking functions to allow for the same individual to be in more than one relationship category.
## I still need to figure out how to build all the networks with the likelihood vectors and not just limit myself to one of the relationships
## 		maybe I can split the likelihood vector on the fly and add another network with the more distant relationship
## 		if I have three relationship (e.g. Hag, CGH, DR) then I would first split hag and then put CGH and DR in another network, when I get to them
##		I would then split CGH and DR
## v7 is very minor changes and is cleaned up for the first public release 

my $verbose;
my $PRINT_FAIL = 0;
my $MIN_PROBABILITY = 0.1;
my $MAX_GEN_GAP = 0;	
my @likelihood_names = qw(PC FS HAG CGH DR UN MZ);

##################################################
## the object constructor (simplistic version)  ##
##################################################
sub new {
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self  = {};
	$self->{NAME} = shift;
	$self->{SEX} = shift;
	$self->{PARENTS} = [];
	$self->{CHILDREN}  = [];
	$self->{RELATIVES} = {};
	$self->{POSSIBLE_RELATIONSHIPS} = {};
	$self->{MITO} = {};
	$self->{NO_MITO_MATCH} = [];
	$self->{MITO_MATCH} = [];
	$self->{Y} = {};
	$self->{NO_Y_MATCH} = [];
	$self->{Y_MATCH} = [];
	$self->{PC} = {};
	$self->{FS} = {};
	$self->{UNRES_first_degree} = {}; #unresolved_first_degree
	$self->{HAG} = {};
	$self->{UNRES_HAG} = {}; #unresolved_HAG
	$self->{CGH} = {};
	$self->{UNRES_CGH} = {}; #unresolved_CGH
	bless($self); # but see below
	bless ($self, $class);
	return $self;
}

sub duplicate {
	my $old_node = shift;
	my $new_node = {};
	$new_node->{NAME} = $old_node->{NAME};
	$new_node->{SEX}  = $old_node->{SEX};
	@{ $new_node->{PARENTS} } = @{$old_node->{PARENTS} };
	@{ $new_node->{CHILDREN} } = @{$old_node->{CHILDREN} };
	%{ $new_node->{RELATIVES} } = %{$old_node->{RELATIVES} };
	%{ $new_node->{POSSIBLE_RELATIONSHIPS} } = %{$old_node->{POSSIBLE_RELATIONSHIPS} };
	%{ $new_node->{MITO} } = %{$old_node->{MITO} };
	@{ $new_node->{NO_MITO_MATCH} } = @{$old_node->{NO_MITO_MATCH} };
	@{ $new_node->{MITO_MATCH} } = @{$old_node->{MITO_MATCH} };
	@{ $new_node->{NO_Y_MATCH} } = @{$old_node->{NO_Y_MATCH} };
	%{ $new_node->{Y} } = %{$old_node->{Y} };
	%{ $new_node->{PC} } = %{ $old_node->{PC} };
	%{ $new_node->{FS} } = %{ $old_node->{FS} }; 
	%{ $new_node->{UNRES_first_degree} } = %{ $old_node->{UNRES_first_degree} };
	%{ $new_node->{HAG} } = %{ $old_node->{HAG} };
	%{ $new_node->{UNRES_HAG} } = %{ $old_node->{UNRES_HAG} };
	%{ $new_node->{CGH} } = %{ $old_node->{CGH} };
	%{ $new_node->{UNRES_CGH} } = %{ $old_node->{UNRES_CGH} };
	bless($new_node);		
	return $new_node;
}

##############################################
## methods to access per-object data        ##
##                                          ##
## With args, they set the value.  Without  ##
## any, they only retrieve it/them.         ##
##############################################
sub set_min_likelihood
{
	$MIN_PROBABILITY = shift;
}

sub set_verbose
{
	$verbose = shift;
}

sub name {
	my $self = shift;
	if (@_) { $self->{NAME} = shift }
	return $self->{NAME};
}

sub sex {
	my $self = shift;
	if (@_) { $self->{SEX} = shift }
	return 0 if $self->{SEX} eq "";
	return $self->{SEX};
}

sub relatives {
	my $self = shift;
	if (@_) 
	{
		%{ $self->{RELATIVES} } = @_;
		foreach my $rel (keys %{ $self->{RELATIVES} })
		{
			my $relationship_likelihood_vector = $self->get_relationship($rel);
			my @possible_relationships = get_possible_relationships(@$relationship_likelihood_vector);
			add_possible_relationsips($rel,\@possible_relationships);
		}
	}
	return %{$self->{RELATIVES} };
}

sub possible_relationships {
		my $self = shift;
		my $rel = shift;
		return ${$self->{POSSIBLE_RELATIONSHIPS} }{$rel};
}

sub all_possible_relationships {
		my $self = shift;
		return %{$self->{POSSIBLE_RELATIONSHIPS} };
}

sub add_possible_relationships {
	my $self = shift;
	my $rel = shift;
	my $possible_relationship_ref = shift;
	${ $self->{POSSIBLE_RELATIONSHIPS} }{$rel} = $possible_relationship_ref;
}

sub get_possible_relationships{
	my $self = shift;
	my @likelihoods = @_;
	my @possibilities;
	
	for(my $i = 0; $i < 6; $i++)
	{
		if(@likelihoods[$i] > $MIN_PROBABILITY)
		{							
			push(@possibilities, @likelihood_names[$i]);
		}
	}
	return @possibilities;
}

sub add_relative {
	my $self = shift;
	my $self_name = $self->{NAME};
	my $rel = shift;
	my @likelihood_vector = @_;
	
	if(exists ${ $self->{RELATIVES} }{$rel}) # I am no longer sure when this is used or why I coded it this way.
	{
		my $old_likelihood_vector = ${ $self->{RELATIVES} }{$rel};
		#die "Can't add $rel (@likelihood_vector) as relative of $self_name because it already exists @$old_likelihood_vector\n";
		#print "$self_name -> $rel = [@likelihood_vector] + @$old_likelihood_vector\n";
		for(my $i = 0; $i < 6; $i++)
		{
			@likelihood_vector[$i] += @$old_likelihood_vector[$i];
		} 
	}
	${ $self->{RELATIVES} }{$rel} = \@likelihood_vector;		
	
	my @possible_relationships = $self->get_possible_relationships(@likelihood_vector);
	$self->add_possible_relationships($rel,\@possible_relationships);
	
}

sub add_y_match {
	my $self = shift;
	my $rel = shift;
	my $match = shift;
	${ $self->{Y} }{$rel} = $match;
	push(@{ $self->{NO_Y_MATCH} }, $rel) if $match eq "0";
	push(@{ $self->{Y_MATCH} }, $rel) if $match eq "1";
}

sub add_mito_match {
	my $self = shift;
	my $rel = shift;
	my $match = shift;
	${ $self->{MITO} }{$rel} = $match;
	#print "adding $self->{NAME} <-> $rel = $match\n";
	push(@{ $self->{NO_MITO_MATCH} }, $rel) if $match eq "0";
	push(@{ $self->{MITO_MATCH} }, $rel) if $match eq "1";
	#print "get no match: @{ $self->{NO_MITO_MATCH} }\n";
}

sub is_y_match{
	my $self = shift;
	my $rel = shift;
	if(!exists ${ $self->{Y} }{$rel})
	{
		return -1;
	}
	return ${ $self->{Y} }{$rel};
}

sub is_mito_match{
	my $self = shift;
	my $rel = shift;
	if(!exists ${ $self->{MITO} }{$rel})
	{
		return -1;
	}
	return ${ $self->{MITO} }{$rel};
}

sub get_y_match {
	my $self = shift;
	return () if $self->{Y_MATCH} eq "";
	return @{$self->{Y_MATCH} };
}

sub get_mito_match {
	my $self = shift;
	return () if $self->{MITO_MATCH} eq "";
	return @{$self->{MITO_MATCH} };
}

sub get_y_no_match {
	my $self = shift;
	#print "get no match: @{ $self->{NO_MITO_MATCH} }\n";
	return () if $self->{NO_Y_MATCH} eq "";
	#return () if !exists $self->{NO_MITO_MATCH};
	return @{$self->{NO_Y_MATCH} };
}

sub get_mito_no_match {
	my $self = shift;
	#print "get no match: @{ $self->{NO_MITO_MATCH} }\n";
	return () if $self->{NO_MITO_MATCH} eq "";
	#return () if !exists $self->{NO_MITO_MATCH};
	return @{$self->{NO_MITO_MATCH} };
}

sub parents {
	my $self = shift;
	if (@_) { @{ $self->{PARENTS} } = @_ }
	return @{$self->{PARENTS} };
}

sub delete_children {
	my $self = shift;
	@{ $self->{CHILDREN} } = ();
}

sub delete_parents {
	my $self = shift;
	@{ $self->{PARENTS} } = ();
}

sub children {
	my $self = shift;
	#print "\@_: @_\n";
	if (@_) { @{ $self->{CHILDREN} } = @_ }
	return @{ $self->{CHILDREN} };
}

sub get_full_sibs {
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	my @parents = $self->parents();
	if(@parents < 1){return}
	my @P1_children = $$network_ref{@parents[0]}->children();
	my @P2_children = $$network_ref{@parents[1]}->children();

	my @full_sibs;
	foreach(@P1_children)
	{
		my $child = $_;
		if(grep($_ eq $child,@P2_children)) #full sibling
		{
			if($child eq $self_name){next;}
			push(@full_sibs,$child);
		}
	}
	return @full_sibs;
}

sub pc {
	my $self = shift;
	if (@_) { %{ $self->{PC} } = @_ }
	return %{ $self->{PC} };
}

sub delete_pc {
	my $self = shift;
	my $rel = shift;
	delete ${ $self->{PC} }{$rel}; 
}

sub fs {
	my $self = shift;
	if (@_) { %{ $self->{FS} } = @_ }
	return %{ $self->{FS} };
}

sub delete_fs {
	my $self = shift;
	my $rel = shift;
	delete ${ $self->{FS} }{$rel}; 
}

sub hag {
	my $self = shift;
	if (@_) { %{ $self->{HAG} } = @_ }
	return %{ $self->{HAG} };
}

sub delete_hag {
	my $self = shift;
	my $rel = shift;
	delete ${ $self->{HAG} }{$rel}; 
}

sub unres_hag {
	my $self = shift;
	if (@_) { %{ $self->{UNRES_HAG} } = @_ }
	return %{ $self->{UNRES_HAG} };
}

sub unres_first_degree {
	my $self = shift;
	if (@_) { %{ $self->{UNRES_first_degree} } = @_ }
	return %{ $self->{UNRES_first_degree} };
}

sub cgh {
	my $self = shift;
	if (@_) { %{ $self->{CGH} } = @_ }
	return %{ $self->{CGH} };
}

sub delete_cgh {
	my $self = shift;
	my $rel = shift;
	delete ${ $self->{CGH} }{$rel};
}

sub delete_rest {
	my $self = shift;
	my $rel = shift;
	delete ${ $self->{DR} }{$rel}; 
	delete ${ $self->{UN} }{$rel}; 
}

sub unres_cgh {
	my $self = shift;
	if (@_) { %{ $self->{UNRES_CGH} } = @_ }
	return %{ $self->{UNRES_CGH} };
}

sub add_parent {
	my $self = shift;
	my $parent = shift;
	
	if(grep ($_ eq $parent, @{ $self->{PARENTS} } ))
	{
		return 1;
	}		
	
	if(@{ $self->{PARENTS} } >= 2)
	{
		die "Can't add $parent to parents of $self->{NAME}; parents full @{$self->{PARENTS} }; use the merge_nodes function instead\n";
		
		## Doing what is below messes things up when the dummy node has more than just this one child
		my @parents = $self->parents();
		if(@parents[0] =~ /Dummy/i)
		{
			#print "REPLACE @parents[0] with $parent\n";
			@parents[0] = $parent;
		}
		elsif(@parents[1] =~ /Dummy/i)
		{
			#print "REPLACE @parents[1] with $parent\n";
			@parents[1] = $parent;
		}
		else
		{
			die "Can't add $parent to parents of $self->{NAME}; parents full @{$self->{PARENTS} }\n";
		}
		$self->parents(@parents);
	}
	else
	{
		push (@{ $self->{PARENTS} }, $parent);
	}
	return 1;
}

sub replace_parent {
	my $self = shift;
	my $parent = shift;
	my $parent_to_replace = shift;
	
	if ($parent_to_replace eq "" || !grep ($_ eq $parent_to_replace, @{ $self->{PARENTS} } ))
	{
		die "ERROR!!! Can't replace parent $parent_to_replace with $parent in @{ $self->{PARENTS} }; does not exist\n";
	}
	if($parent_to_replace !~ /Dummy/i){die "ERROR!!! Replacing non Dummy parent $parent_to_replace with $parent\n";}
	
	#print "REPLACE $parent_to_replace with $parent\n";
	my @parents = $self->parents();
	if(@parents[0] eq $parent_to_replace)
	{
		@parents[0] = $parent;
	}
	elsif(@parents[1] eq $parent_to_replace)
	{
		@parents[1] = $parent;
	}
	else
	{
		die "ERROR!  CHECK IT OUT\n";
	}
	$self->parents(@parents);
	return 1;
}

sub add_child {
	my $self = shift;
	my $child = shift;

	if(!grep($_ eq $child,@{ $self->{CHILDREN} } ) )
	{
		push (@{ $self->{CHILDREN} }, $child);
		if($self->name() eq "Dummy1" && $child eq "Dummy5")
		{
			#print "ADDING\n";
			return 1000;
		}
		return 1;
	}
	return 0;
}

sub replace_child {
	my $self = shift;
	my $child = shift;
	my $child_to_replace = shift;
	
	if (!grep ($_ eq $child_to_replace, @{ $self->{CHILDREN} } ))
	{
		die "ERROR!!! Can't replace child $child_to_replace with $child in". @{ $self->{CHILDREN} }."; does not exist\n";
	}

	#print "REPLACE child $child_to_replace with $child\n";
	my @children = $self->children();
	for(my $i = 0; $i < @children; $i++)
	{
		if(@children[$i] =~ $child_to_replace)
		{
			@children[$i] = $child;
			last;
		}
	}
	$self->children(@children);
	return 1;
}

sub remove_child {
	my $self = shift;
	my $self_name = $self->name();
	#print "self_name: $self_name\n";
	my $child = shift;
		
	if (!grep ($_ eq $child, @{ $self->{CHILDREN} } ))
	{
		print "WARNING!!! Can't remove child $child in @{ $self->{CHILDREN} }; does not exist\n";
		return 0;
	}
	
	#print "REMOVING child $child\n";
	my @children = $self->children();
	my @new_children;
	for(my $i = 0; $i < @children; $i++)
	{
		if(@children[$i] ne $child)
		{
			#print "push @children[$i]\n";
			push(@new_children,@children[$i]);
		}
	}
	#print "new children: @new_children\n";
		@{ $self->{CHILDREN} } = @new_children;
	#$self->children(@new_children);
	return 1;
}

sub get_relationship {
	my $self = shift;
	my $relative = shift;
	my %rels = $self->relatives();
	return $rels{$relative};
}


## add first_degree if it doesn't already exists
sub update_unres_first_degree {
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	my %PCs = % {$self->{PC} };
	my %FSs = % {$self->{FS} };
	my %unres_first_degrees;

	my %subpedigree = ();

	my @PC_names = keys %PCs;
	my @FS_names = keys %FSs;
	#print "node: $self_name; pc_rels: @PC_names; fs_names: @FS_names\n";

	foreach(keys %PCs)
	{
		#print "checking $_\n";
		my $resolved_PC = $self->does_relationship_exist_in_pedigree($network_ref,$_,"PC",4);
		#print "resolved: $resolved_PC\n";
		if($resolved_PC == 0)
		{
			$unres_first_degrees{$_} = $PCs{$_};
		}
	}
	foreach(keys %FSs)
	{
		#print "checking $_\n";
		my $resolved_FS = $self->does_relationship_exist_in_pedigree($network_ref,$_,"FS",4);
		#print "resolved: $resolved_FS\n";
		if($resolved_FS == 0)
		{
			$unres_first_degrees{$_} = $FSs{$_};
		}
	}

	my @unres_first_degree_names = keys %unres_first_degrees;
	if(@unres_first_degree_names < 1)
	{
		%{ $self->{UNRES_first_degree} } = ();
	}

	$self->unres_first_degree(%unres_first_degrees);
}

## return first degree if it doesn't already exists
sub get_unres_first_degree {
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	my %unres_first_degrees = %{$self->{UNRES_first_degree} };
	my @unres_first_degree_names = keys %unres_first_degrees;
	foreach(keys %unres_first_degrees)
	{
		my $resolved = $self->does_relationship_exist_in_pedigree($network_ref,$_,"FS",4);
		if($resolved eq 0)
		{
			delete ($self->{UNRES_first_degree}{$_} );
			return $_;
		}
		else
		{
			delete ($self->{UNRES_first_degree}{$_} );
		}
	}
	return "";
}

## add hag if it doesn't already exists
sub update_unres_hag {
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	my %HAGs = % {$self->{HAG} };
	my %unres_HAGs;

	my %subpedigree = ();

	foreach(keys %HAGs)
	{
		my $resolved = $self->does_relationship_exist_in_pedigree($network_ref,$_,"HAG",4);
		## This messes things up if half siblings need to be merged on a common dummy parent if the two dummy parents are already children of the real grandparents.
		#if($resolved eq 0)
		#{
		#	$resolved = $self->does_relationship_exist_in_pedigree($network_ref,$_,"CGH",4);
		#}
		if($resolved == 0)
		{
			$unres_HAGs{$_} = $HAGs{$_};
		}
	}
	my @unres_hag_names = keys %unres_HAGs;
	#print "$self_name: " . @unres_hag_names ."\n";
	if(@unres_hag_names < 1)
	{
		%{ $self->{UNRES_hag} } = ();
	}
	$self->unres_hag(%unres_HAGs);
}


## return hag if it doesn't already exists
## and if a CGH doesn't already exist in the pedigree. 
sub get_unres_hag {
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	my %unres_HAGs = %{$self->{UNRES_HAG} };
	my @unres_hag_names = keys %unres_HAGs;
	#if(@unres_hag_names eq 0)
	#{
	#	$self->update_unres_hag($network_ref);
	#}
	#%unres_HAGs = %{$self->{UNRES_HAG} };
	#print "unres_hags $self_name: @unres_hag_names\n";
	foreach(keys %unres_HAGs)
{
	my $resolved = $self->does_relationship_exist_in_pedigree($network_ref,$_,"HAG",4);
	## This messes things up if half siblings need to be merged on a common dummy parent if the two dummy parents are already children of the real grandparents.
	#if($resolved eq 0)
	#{
	#	$resolved = $self->does_relationship_exist_in_pedigree($network_ref,$_,"CGH",4);
	#}
	#print "\n\nunres_hags $self_name: $_ : $resolved\n";
	if($resolved eq 0)
	{
		delete ($self->{UNRES_HAG}{$_} );
		return $_;
	}
	else
	{
		delete ($self->{UNRES_HAG}{$_} );
	}
}
return "";
}

## add CGH if it doesn't already exists
## and if a hag doesn't already exist in the pedigree
sub update_unres_cgh {
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	#print "self_name: $self_name\n";
	my %CGHs = % {$self->{CGH} };
	my %unres_CGHs;

	my %subpedigree = ();

	foreach(keys %CGHs)
	{
		my $resolved = $self->does_relationship_exist_in_pedigree($network_ref,$_,"CGH",4);
		if($resolved == 0)
		{
			$resolved = $self->does_relationship_exist_in_pedigree($network_ref,$_,"HAG",4);
		}
		#print "\n\nunres_cghs_update $self_name: $_ : $resolved\n";
		if($resolved == 0)
		{
			$unres_CGHs{$_} = $CGHs{$_};
		}
	}
	my @unres_cgh_names = keys %unres_CGHs;
	if(@unres_cgh_names < 1)
	{
		%{ $self->{UNRES_cgh} } = ();
	}
	$self->unres_cgh(%unres_CGHs);
}

## return CGH if it doesn't already exists
## and if a hag doesn't already exist in the pedigree
sub get_unres_cgh {
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	my %unres_CGHs = %{$self->{UNRES_CGH} };
	my @unres_cgh_names = keys %unres_CGHs;
	#print "unres_cghs $self_name: @unres_cgh_names\n";
		if(@unres_cgh_names == 0)
	{
		#$self->update_unres_cgh($network_ref);
	}
		#%unres_CGHs = %{$self->{UNRES_CGH} };
	#my @unres_cgh_names = keys %unres_CGHs;
	#print "unres_cghs $self_name: @unres_cgh_names\n";

		foreach(keys %unres_CGHs)
	{
		my $resolved = $self->does_relationship_exist_in_pedigree($network_ref,$_,"CGH",4);
		if($resolved eq 0)
		{
			$resolved = $self->does_relationship_exist_in_pedigree($network_ref,$_,"HAG",4);
		}
		#print "\n\nunres_cghs $self_name: $_ : $resolved\n";
		if($resolved eq 0)
		{
			delete ($self->{UNRES_CGH}{$_} );
			return $_;
		}
		else
		{
			delete ($self->{UNRES_CGH}{$_} );
		}
	}
	return "";
}
	
sub num_real_parents
{
	my $self = shift;
	my $num = 0;
	my @parents = @{$self->{PARENTS} };
	#print "parents: @parents\n";
	foreach( @parents )
	{
		#print "########parent: $_\n";
		if($_ !~ /Dummy/i){$num++;}
	}
	return $num;
}

## Is relatioship supported by the current pedigree, returns 0 or 1
sub does_relationship_exist_in_pedigree
{
#if(0 == $self->does_relationship_exist_in_pedigree($network_ref, $_,"", $phase_test) )
	my $self = shift;
	my $network_ref = shift;
	my $rel = shift;
	my $relationship = shift;
	my $phase = shift;
	my $self_name = $self->name();
	
	## If the call did not specify a relationship, then you need to test all relationships in the likelihood vector
	## It will start
	if ($relationship eq "")
	{
		my $relationship_likelihood_vector = $self->get_relationship($rel);
		my $relationships_ref = $self->possible_relationships($rel);
		foreach(@$relationships_ref)
		{ 
			if($_ eq "PC" && $phase < 1){return 99;}
			if($_ eq "FS" && $phase < 1){return 99;}
			if($_ eq "HAG" && $phase < 2){return 199;}
			if($_ eq "CGH" && $phase < 3){return 299;}
			if($_ eq "DR" && $phase < 4){return 399;}
			if($_ eq "UN" && $phase < 4){return 499;}
			my $return_val = $self->does_relationship_exist_in_pedigree($network_ref, $rel, $_, $phase);
			if($return_val ne 0)
			{
				return $return_val;
			}
		}
		return 0;
	}

	my @parents = @{$self->{PARENTS} };
	my @children = @{$self->{CHILDREN} };
		
	if($relationship eq "PC")
	{
		if(grep($_ eq $rel, @parents)){return 1;}
		if(grep($_ eq $rel, @children)){return 1;}
	}
	elsif($relationship eq "FS")
	{
		if(@parents < 1){return 0;}
		my @P1_children = $$network_ref{@parents[0]}->children();
		my @P2_children = $$network_ref{@parents[1]}->children();	 
		foreach(@P1_children)
		{
			my $child = $_;
			if(grep($_ eq $child,@P2_children)) #full sibling
			{
				if($child eq $self_name){next;}
				if($child eq $rel){return 1;}
			}
		}
	}
	elsif($relationship eq "HAG")
	{
		if($phase < 2)
		{
			return 599;
		}
		## Check if $rel is a neice/nephew
		if(@parents ne 0)
		{
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
						if($_ eq $rel){return 2;}
					}
				}
			}
			
			## check if rel is half-sib
			foreach(@P1_children)
			{
				my $child = $_;
				if(!grep($_ eq $child,@P2_children)) # Half sib
				{
					if($child eq $rel){return 3;}
				}
			}
			foreach(@P2_children)
			{
				my $child = $_;
				#print "child $child\n";
				if(!grep($_ eq $child,@P1_children)) # Half sib
				{
					if($child eq $rel){return 4;}
				}
			}
						
			## check if rel is a grandparent
			my @P1_parents = $$network_ref{@parents[0]}->parents();
			my @P2_parents = $$network_ref{@parents[1]}->parents();
			foreach(@P1_parents)
			{
				if($_ eq $rel){return 5;}
			}
			foreach(@P2_parents)
			{
				if($_ eq $rel){return 6;}
			}
			
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
						if($child eq $rel){return 7;}
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
						if($child eq $rel){return 8;}
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
				if($_ eq $rel){return 9;}
			}

		}
	}
	elsif($relationship eq "CGH")
	{
		if($phase < 3)
		{
			return 699;
		}
		## Check if rel is great grandchild
		foreach(@children)
		{
			my @G_children = $$network_ref{$_}->children();
			foreach(@G_children)
			{
				my @GG_children = $$network_ref{$_}->children();
				foreach(@GG_children)
				{
				if($_ eq $rel){return 10;}
				}
			}
		}
		
		if(@parents eq 0){return 0;} ## I think this might be useless because only dummies don't have parents
			
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
					if($_ eq $rel){return 11;}
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
					if($_ eq $rel){return 12;}
				}
			}
		}
		
		## check if rel is half uncle/aunt = half sibling of parent
		if(@P1_parents > 0)
		{
			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
	 
			foreach(@G1_1_children)
			{
				my $child = $_;
				if(!grep($_ eq $child,@G1_2_children)) # half uncle/aunt
				{
					#print "$child : $self_name\n";
					if($_ eq $rel){return 13;}
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
					if($_ eq $rel){return 14;}					
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
			
			foreach(@G2_1_children)
			{
				my $child = $_;
				if(!grep($_ eq $child,@G2_2_children)) #half Unclde/aunt
				{
					if($_ eq $rel){return 15;}					
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
					if($_ eq $rel){return 16;}
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
		
		#print "$self_name grandparents: @grandparents\n";
		foreach my $G_parent (@grandparents)
		{
			my @GG_parents = $$network_ref{$G_parent}->parents();
			foreach(@GG_parents)
			{
				if($_ eq $rel){return 17;}
			}
		}
		
		## Check if rel is first cousin
		my @grandparents0 = $$network_ref{$parents[0]}->parents();
		my @grandparents1 = $$network_ref{$parents[1]}->parents();
		my @rel_parents = $$network_ref{$rel}->parents();
		my @rel_grandparents0 = $$network_ref{$rel_parents[0]}->parents();
		my @rel_grandparents1 = $$network_ref{$rel_parents[1]}->parents();
		if(do_arrays_match(\@grandparents0,\@rel_grandparents0) && @grandparents0 > 0){return 17.5}
		elsif(do_arrays_match(\@grandparents0,\@rel_grandparents1) && @grandparents0 > 0){return 17.5}
		elsif(do_arrays_match(\@grandparents1,\@rel_grandparents0) && @grandparents1 > 0){return 17.5}
		elsif(do_arrays_match(\@grandparents1,\@rel_grandparents1) && @grandparents1 > 0){return 17.5}
					
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
					if($child eq $rel){return 18;}
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
						if($_ eq $rel){return 19;}
					}
				}
			}
		}	
	}
	elsif($relationship eq "DR")
	{
		return 799;
		## PRIMUS doesn't check that DR are indeed DR because those relations are difficult to predict at this point
		if($phase < 4)
		{
			return 99;
		}
	
		## Check for more distant relationships
		if($self->is_relative_of($rel))
		{
			return 20;
		}
	}
	elsif($relationship eq "UN")
	{
		return 899;
		if($phase < 4)
		{
			return 99;
		}
		## Check for the unrelateds
			my %all_relatives = $self->get_all_relatives($network_ref);
		if(!exists $all_relatives{$rel})
		{
			return 21;
		}
	}
	
	return 0;
}

## Once two individuals are in the pedigree, they should not create any relationships that are not described.
## This subroutine tests that all relationship described by the subpedigree match the relationship described by the IBD file	
sub check_pedigree_relationships
{	
	my $self = shift;
	my $network_ref = shift;
	my $phase = shift;
	my $self_name = $self->name();
	
	if($self->pass_generation_check($network_ref) eq 0)
	{
		#print "$self_name is related to itself in a different generation\n";
		return 100;
	}
	if($self_name =~/Dummy/i){next;}
		my @parents = @{$self->{PARENTS} };
	my @children = @{$self->{CHILDREN} };
		#my %rels = % {$self->{RELATIVES} };
		my %rels = % {$self->{POSSIBLE_RELATIONSHIPS} };
		my $num_fails = 0;
		
	## If $self is ancestor or descendant of itself, then fail
		
		
	## Check if rel is child of self
	foreach my $child (@children)
	{
	if($child =~ /Dummy/i){next;}
		if(!grep($_ eq "PC", @{$rels{$child} }))
	{
		if($PRINT_FAIL){print "FAIL!!! $child not child of $self_name\n";}
		$num_fails++;
	}
		#if($rels{$_} ne "PC"){print "FAIL!!! $_ not child of $self_name\n";$num_fails++;}
	}		
	
	if ($phase >= 2)
	{
		## Check if rel is grandchild of self
		#print "$self_name children: @children\n";
		foreach my $child (@children)
		{
			#print "$self_name child: $child\n";
			my @G_children = $$network_ref{$child}->children();
			foreach my $grandchild (@G_children) # These are neices/nephews
			{
				if($grandchild =~ /Dummy/i){next;}
			
				#print "grandchild = $grandchild; grandpa = $self_name\n";
				if(!grep($_ eq "HAG", @{$rels{$grandchild} }))
				{
					if($PRINT_FAIL){print "FAIL!!! $grandchild not grandchild of $self_name\n";}
					$num_fails++;
				}
			}
		}
	}

	if ($phase >= 3)
	{
		## Check if rel is great grandchild self
		foreach my $c (@children)
		{
			my @G_children = $$network_ref{$c}->children();
			foreach my $gc (@G_children) # These are neices/nephews
			{
				my @GG_children = $$network_ref{$gc}->children();
				foreach my $ggc(@GG_children) # These are neices/nephews
				{
					if($ggc =~ /Dummy/i){next;}
					if(!grep($_ eq "CGH", @{$rels{$ggc} }))
					{
						if($PRINT_FAIL){print "FAIL1!!! $ggc not great-grandchild of $self_name\n";}
						$num_fails++;
					}
				}
			}
		}
	}
		
	if(@parents ne 0)
	{
		## Check parents
		foreach my $p (@parents)
		{
			if($p =~/Dummy/i){next;}
			if(!grep($_ eq "PC", @{$rels{$p} }))
			{
				if($PRINT_FAIL){print "FAIL!!! $p not PC with $self_name. it is " . $rels{$p} ."\n";}
				$num_fails++;
			}
		}
		
		## Check full siblings
		my @P1_children = $$network_ref{@parents[0]}->children();
		my @P2_children = $$network_ref{@parents[1]}->children();
		#print "@parents[0] children: @P1_children\n";
		#print "@parents[1] children: @P2_children\n";
		my @FS;
		my @HS;
		my @UNCLES_AUNTS;
	 
		foreach(@P1_children)
		{
			my $child = $_;
			if(grep($_ eq $child,@P2_children)) #full sibling
			{
				if($child eq $self_name){next;}
				push(@FS,$child);
				if($child =~/Dummy/i){next;}
				if(!grep($_ eq "FS", @{$rels{$child} }))
				{
					if($PRINT_FAIL){print "FAIL!!! $child not FS with $self_name\n";}
					$num_fails++;
				}
			}
		}
		
		## Don't do the phase 2 checks if we are not in Phase 2
		if ($phase < 2){
			return $num_fails;
		}

		
		## Check half siblings
		foreach my $child (@P1_children)
		{
			if(!grep($_ eq $child, @P2_children)) # Half sib
			{
				push(@HS,$child);
				if($child =~ /Dummy/i){next;}
				if(!grep($_ eq "HAG", @{$rels{$child} }))
				{
					if($PRINT_FAIL){print "FAIL!!! $child not HS with $self_name\n";}
					$num_fails++;
				}
			}
		}
		foreach my $child (@P2_children)
		{
			if(!grep($_ eq $child, @P1_children)) # Half sib
			{
				push(@HS,$child);
				if($child =~ /Dummy/i){next;}
				if(!grep($_ eq "HAG", @{$rels{$child} }))
				{
					if($PRINT_FAIL){print "FAIL!!! $child not HS with $self_name\n";}
					$num_fails++;
				}
			}
		}

		## Check Nieces/nephews
		foreach my $fs (@FS)
		{
			my @FS_children = $$network_ref{$fs}->children();
			foreach my $child (@FS_children) # These are neices/nephews
			{
				if($child =~ /Dummy/i){next;}
					if(!grep($_ eq "HAG", @{$rels{$child} }))
				{
					if($PRINT_FAIL){print "FAIL!!! $child not neices/nephew with $self_name\n";}
					$num_fails++;
				}
			}
		}

		## Check grandparents
		my @P1_parents = $$network_ref{@parents[0]}->parents();
		my @P2_parents = $$network_ref{@parents[1]}->parents();
		#print "@parents[0] parents: @P1_parents\n";
		#print "@parents[1] parents: @P2_parents\n";
		foreach my $gg (@P1_parents)
		{
			if($gg =~/Dummy/i){next;}
			if(!grep($_ eq "HAG", @{$rels{$gg} }))
			{
				if($PRINT_FAIL){print "FAIL!!! $gg not GG with $self_name\n";}
				$num_fails++;
			}
		}
		foreach my $gg (@P2_parents)
		{
			if($gg =~/Dummy/i){next;}
			if(!grep($_ eq "HAG", @{$rels{$gg} }))
			{	
				if($PRINT_FAIL){print "FAIL!!! $gg not GG with $self_name\n";}
				$num_fails++;
			}
		}
		
		
		## Check if $rel is uncle/aunt of $self
		if(@P1_parents > 0)
		{
			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
	 
			foreach(@G1_1_children)
			{
				my $child = $_;
				if(grep($_ eq $child,@G1_2_children)) #uncle/aunt
				{
					#print "$child : $self_name\n";
				if($child eq @parents[0]){next;}
					push(@UNCLES_AUNTS,$child);
					if($child =~ /Dummy/i){next;}
					if(!grep($_ eq "HAG", @{$rels{$child} }))
					{
						if($PRINT_FAIL){print "FAIL!!! $self_name not neice or nephew1 of $child\n";}
						$num_fails++;
					}
				}
			}
		}
		if(@P2_parents > 0)
		{
			#print "$self_name parent @parents[1] parents @P2_parents\n";
			my @G2_1_children = $$network_ref{@P2_parents[0]}->children();
			my @G2_2_children = $$network_ref{@P2_parents[1]}->children();
			#print "P2 @P2_parents[0] children: @G2_1_children\n";
			#print "P2 @P2_parents[1] children: @G2_2_children\n";
	 
			foreach(@G2_1_children)
			{
				my $child = $_;
				if(grep($_ eq $child,@G2_2_children)) #Unclde/aunt
				{
				if($child eq @parents[1]){next;}
					push(@UNCLES_AUNTS,$child);					
					if($child =~/Dummy/i){next;}
					if(!grep($_ eq "HAG", @{$rels{$child} }))
					{
						if($PRINT_FAIL){print "FAIL!!! $self_name not neice or nephew2 of $child\n";}
						$num_fails++;
					}
				}
			}
		}
		
		
		## Don't do the phase 3 checks if we are not in Phase 3
		if ($phase < 3){
			return $num_fails;
		}
		
		## Check half uncle/aunts
		if(@P1_parents > 0)
		{
			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
			#print "P1 @P1_parents[0] children: @G1_1_children\n";
			#print "P1 @P1_parents[1] children: @G1_2_children\n";
	 
			foreach(@G1_1_children)
			{
				my $child = $_;
				if(!grep($_ eq $child,@G1_2_children)) #uncle/aunt
				{
					#print "$child : $self_name\n";
					if($child =~ /Dummy/i || $child eq @parents[0]){next;}
					if(!grep($_ eq "CGH", @{$rels{$child} }))
					{
						if($PRINT_FAIL){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";}
						$num_fails++;
					}
				}
			}
			foreach(@G1_2_children)
			{
				my $child = $_;
				if(!grep($_ eq $child,@G1_1_children)) #uncle/aunt
				{
					#print "$child : $self_name\n";
					if($child =~ /Dummy/i || $child eq @parents[0]){next;}
					if(!grep($_ eq "CGH", @{$rels{$child} }))
					{
						if($PRINT_FAIL){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";}
						$num_fails++;
					}
				}
			}
		}
		if(@P2_parents > 0)
		{
			my @G2_1_children = $$network_ref{@P2_parents[0]}->children();
			my @G2_2_children = $$network_ref{@P2_parents[1]}->children();
			#print "P2 @P2_parents[0] children: @G2_1_children\n";
			#print "P2 @P2_parents[1] children: @G2_2_children\n";
		
			foreach my $child (@G2_1_children)
			{
				if(!grep($_ eq $child,@G2_2_children)) #Unclde/aunt
				{
					if($child =~/Dummy/i || $child eq @parents[1]){next;}
					if(!grep($_ eq "CGH", @{$rels{$child} }))
					{
						if($PRINT_FAIL){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";}
						$num_fails++;
					}
				}
			}
			foreach my $child (@G2_2_children)
			{
				if(!grep($_ eq $child,@G2_1_children)) #Unclde/aunt
				{
					if($child =~/Dummy/i || $child eq @parents[1]){next;}
					if(!grep($_ eq "CGH", @{$rels{$child} }))
					{
						if($PRINT_FAIL){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";}
						$num_fails++;
					}
				}
			}
		}
			
		## Check if $rel is half neice/nephew of $self
		foreach my $hs (@HS)
		{
			my @HS_children = $$network_ref{$hs}->children();
			foreach my $hs_child (@HS_children) # These are neices/nephews
			{
				if($hs_child =~ /Dummy/i){next;}
				if(!grep($_ eq "CGH", @{$rels{$hs_child} }))
				{
					if($PRINT_FAIL){print "FAIL!!! $hs_child not half neices/nephew with $self_name\n";}
					$num_fails++;
				}
			}
		}
			
			
		## Check if $rel is GREAT-grandparent of $self
		my @grandparents;
		push(@grandparents, @P1_parents);
		push(@grandparents, @P2_parents);
		#my @P1_parents = $$network_ref{@parents[0]}->parents();
		#my @P2_parents = $$network_ref{@parents[1]}->parents();
		#print "@parents[0] parents: @P1_parents\n";
		#print "@parents[1] parents: @P2_parents\n";
		
		#print "$self_name grandparents: @grandparents\n";
		foreach my $gp (@grandparents)
		{
			my @GG_parents = $$network_ref{$gp}->parents();
			foreach my $gg_parent (@GG_parents)
				{
					if($gg_parent =~/Dummy/i){next;}
					if(!grep($_ eq "CGH", @{$rels{$gg_parent} }))
					{
						if($PRINT_FAIL){print "FAIL!!! $gg_parent not GG-parent of $self_name\n";}
						$num_fails++;
					}
				}
		}
		
		## CHECK if COUSINS
		foreach my $aa (@UNCLES_AUNTS)
		{
			my @cousins = $$network_ref{$aa}->children();
			foreach my $c1 (@cousins)
			{
				if($c1 =~ /Dummy/i){next;}
				if(($c1 =~ /NA21379/ && $self_name =~ /NA21517/) || ($c1 =~ /NA21517/ && $self_name =~ /NA21379/) )
				{
					print "$self_name -> $c1 = @{$rels{$c1} }\n";
				}
				if(!grep($_ eq "CGH", @{$rels{$c1} }))
				{
					if($PRINT_FAIL){print "FAIL!!! $c1 not cousin with $self_name\n";}
					if(($c1 =~ /NA21379/ && $self_name =~ /NA21517/) || ($c1 =~ /NA21517/ && $self_name =~ /NA21379/) )
					{
						print "FAIL!!! $c1 not cousin with $self_name\n";
					}
					$num_fails++;
				}
			}
		}
		
		## CHECK IF REL IS GREAT UNCLE/AUNT of self
		my @G_parents;
		push(@G_parents, @P1_parents);
		push(@G_parents, @P2_parents);
		
		foreach my $G_parent (@G_parents)
		{
			my @GG_parents = $$network_ref{$G_parent}->parents();
			if(@GG_parents eq 0){next;}
			
			my @GG1_children = $$network_ref{@GG_parents[0]}->children();
			my @GG2_children = $$network_ref{@GG_parents[1]}->children();
			#print "GG1 @GG_parents[0] children: @GG1_children\n";
			#print "GG2 @GG_parents[1] children: @GG2_children\n";
	 
			foreach my $child (@GG1_children)
			{
				if(grep($_ eq $child,@GG2_children)) #great uncle/aunt
				{
					#print "$child : $self_name\n";
					if($child =~ /Dummy/i || $child eq $G_parent){next;}
					if(!grep($_ eq "CGH", @{$rels{$child} }))
					{
						if($PRINT_FAIL){print "FAIL!!! $self_name not great neice or nephew of $child\n";}
						$num_fails++;
					}
				}
			}
		}
				
		## Check if rel is great neice or nephew of self
		foreach my $fs (@FS)
		{
			my @FS_children = $$network_ref{$fs}->children();
			foreach my $fs_child (@FS_children) # These are neices/nephews
			{
				my @FS_Gchildren = $$network_ref{$fs_child}->children();
				foreach my $fs_gchild (@FS_Gchildren) # These are neices/nephews
				{
					if($fs_gchild =~ /Dummy/i){next;}
					if(!grep($_ eq "CGH", @{$rels{$fs_gchild} }))
					{
						if($PRINT_FAIL){print "FAIL!!! $fs_gchild not great neices/nephew with $self_name\n";}
						$num_fails++;
					}
				}
			}
		}
	
		########################################################
		## Don't do the phase 4 checks if we are not in Phase 4
		if ($phase < 4){
			return $num_fails;
		}
		
		## Check for more distant relationships
		my %all_relatives = $self->get_all_relatives($network_ref);
		my %current_rels = $self->relatives();
		
		foreach my $temp_rel (keys %all_relatives)
		{
			if(!exists $current_rels{$temp_rel} && $temp_rel ne $self_name)
			{
				if($temp_rel =~ /Dummy/i){next;}
				#if($rels{$_} ne "DR"){print "FAIL!!! $_ not DR with $self_name\n";$num_fails++;}
			}
		}
		
		## Check for the unrelateds
		my %all_pedigree = $self->get_subpedigree_names($network_ref);
		
		foreach my $temp_rel (keys %all_pedigree)
		{
			if(!exists $all_relatives{$temp_rel})
			{
				if($temp_rel =~ /Dummy/i){next;}
			#if(exists $rels{$_} && $rels{$_} ne "UN"){print "FAIL!!! $_ not unrelated to $self_name\n";$num_fails++;}
			}
		}
	}
	return $num_fails;
}

## once two individuals are part of the same subpedigree, they should match all described relationships of a lower relationship not being connecged in this phase; 
## this subroutine tests that the subpedigree contains all relationship that the IBD file requires
sub are_relationships_missing_in_pedigree ### yes = 1; no = 0;
{
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	my $phase_test = shift;
	#print "phase $phase_test\n";
	my %names = $self->get_subpedigree_names($network_ref);
	my %relatives = % {$self->{RELATIVES} };
	foreach(keys %relatives)
	{
		## Don't worry about the relationship if both individuals aren't in the same subpedigree yet
		## THIS DOES NOT WORK, BECAUSE IT WILL LET PEDIGREES THROUGH WITH UNRESOLVED RELATIONSHIPS
		## need to make this check that all relationships between sub pedigrees can be DR or UN
		#if(!exists $names{$_}){next;}
		
		if(0 == $self->does_relationship_exist_in_pedigree($network_ref, $_,"", $phase_test) )
		{
			#print "$self_name + $_ -> [" . $relatives{$_} . "] is missing\n";
			return 1;
		}
	}
	return 0;
}

## This method could potentially be very expensive (we will see)	
## I can store this data as a hash or array in each node and update every time I combine two 
## subpedigree by just merging the two hashes/arrays
sub get_subpedigree_names
{
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	my @parents = @{$self->{PARENTS} };
	my @children = @{$self->{CHILDREN} };
	
	my @names_to_explore;
	push(@names_to_explore,@parents);
	push(@names_to_explore,@children);
	
	my %subpedigree_names = ();
	$subpedigree_names{$self_name} = 1;
	while(my $name = shift(@names_to_explore) )
	{
		#print "NAME $name\n";
		if(exists $subpedigree_names{$name}){next;}
		$subpedigree_names{$name} = 1;
		my @temp_names = $$network_ref{$name}->parents();
		push(@names_to_explore,@temp_names);
		my @temp_names = $$network_ref{$name}->children();
		push(@names_to_explore,@temp_names);
	}
	return %subpedigree_names;
}


## Go through all relatives and identify their relationship
sub make_relative_network_from_pedigree
{	
	my $self = shift;
	my $network_ref = shift;
	my $phase = 10;
	my $self_name = $self->name();
	if($self_name =~/Dummy/i){next;}
	my @parents = @{$self->{PARENTS} };
	my @children = @{$self->{CHILDREN} };
	my %rels = % {$self->{RELATIVES} };
	my $num_fails = 0;
	
	my @PC_likelihoods = (1,0,0,0,0,0);
	my @FS_likelihoods = (0,1,0,0,0,0);
	my @HAG_likelihoods = (0,0,1,0,0,0);
	my @CGH_likelihoods = (0,0,0,1,0,0);
	my @D1C_likelihoods = (0,0,0,2,0,0);
	my @HS1C_likelihoods = (0,0,1,1,0,0);
	my @DR_likelihoods = (0,0,0,0,1,0);
	my @UN_likelihoods = (0,0,0,0,0,1);
	
	## Check if rel is child of self
	foreach(@children)
	{
		if($_ =~ /Dummy/i){next;}
		#if($rels{$_} ne "PC"){print "FAIL!!! $_ not child with $self_name\n";$num_fails++;}
		$$network_ref{$self_name}->add_relative($_,@PC_likelihoods);
	}		
	
	if ($phase >= 2)
	{
		## Check if rel is grandchild of self
		foreach(@children)
		{
			my @G_children = $$network_ref{$_}->children();
			foreach(@G_children) # These are neices/nephews
			{
				if($_ =~ /Dummy/i){next;}
				#if($rels{$_} ne "HAG"){print "FAIL!!! $_ not grandchild of $self_name\n";$num_fails++;}
				$$network_ref{$self_name}->add_relative($_,@HAG_likelihoods);
			}
		}
	}

	if ($phase >= 3)
	{
		## Check if rel is great grandchild self
		foreach(@children)
		{
			my @G_children = $$network_ref{$_}->children();
			foreach(@G_children) # These are neices/nephews
			{
				my @GG_children = $$network_ref{$_}->children();
				foreach(@GG_children) # These are neices/nephews
				{
					if($_ =~ /Dummy/i){next;}
				#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not great-grandchild of $self_name\n";$num_fails++;}
				$$network_ref{$self_name}->add_relative($_,@CGH_likelihoods);
				}
			}
		}
	}

	#print "parents: @parents\n";
	if(@parents ne 0)
	{
		## Check parents
		foreach(@parents)
		{
			if($_ =~/Dummy/i){next;}
			#if($rels{$_} ne "PC"){print "FAIL!!! $_ not PC with $self_name. it is " . $rels{$_} ."\n";$num_fails++;}
			$$network_ref{$self_name}->add_relative($_,@PC_likelihoods);
		}
		
		## Check full siblings
		my @P1_children = $$network_ref{@parents[0]}->children();
		my @P2_children = $$network_ref{@parents[1]}->children();
		#print "@parents[0] children: @P1_children\n";
		#print "@parents[1] children: @P2_children\n";
		my @FS;
		my @HS;
		my @UNCLES_AUNTS;
	 
		foreach(@P1_children)
		{
			my $child = $_;
			if(grep($_ eq $child,@P2_children)) #full sibling
			{
				if($child eq $self_name){next;}
				push(@FS,$child);
				if($child =~/Dummy/i){next;}
				$$network_ref{$self_name}->add_relative($child,@FS_likelihoods);
			}
		}
		
		## Don't do the phase 2 checks if we are not in Phase 2
		if ($phase < 2){
			return $num_fails;
		}

		
		## Check half siblings
		foreach(@P1_children)
		{
			my $child = $_;
			if(!grep($_ eq $child, @P2_children)) # Half sib
			{
				push(@HS,$child);
				if($child =~ /Dummy/i){next;}
				$$network_ref{$self_name}->add_relative($child,@HAG_likelihoods);
			}
		}
		foreach(@P2_children)
		{
			my $child = $_;
			if(!grep($_ eq $child, @P1_children)) # Half sib
			{
			push(@HS,$child);
				if($child =~ /Dummy/i){next;}
			$$network_ref{$self_name}->add_relative($child,@HAG_likelihoods);
			}
		}

		## Check Nieces/nephews
		foreach(@FS)
		{
			my @FS_children = $$network_ref{$_}->children();
			foreach(@FS_children) # These are neices/nephews
			{
				if($_ =~ /Dummy/i){next;}
				#if($rels{$_} ne "HAG"){print "FAIL!!! $_ not neices/nephew with $self_name\n";$num_fails++;}
				$$network_ref{$self_name}->add_relative($_,@HAG_likelihoods);
			}
		}


		## Check grandparents
		my @P1_parents = $$network_ref{@parents[0]}->parents();
		my @P2_parents = $$network_ref{@parents[1]}->parents();
		#print "$self_name parents:\n";
		#print "@parents[0] parents: @P1_parents\n";
		#print "@parents[1] parents: @P2_parents\n";
		foreach(@P1_parents)
		{
			if($_ =~/Dummy/i){next;}
			#if($rels{$_} ne "HAG"){print "FAIL!!! $_ not GG with $self_name\n";$num_fails++;}
			$$network_ref{$self_name}->add_relative($_,@HAG_likelihoods);
		}
		foreach(@P2_parents)
		{
			if($_ =~/Dummy/i){next;}
			#if($rels{$_} ne "HAG"){print "FAIL!!! $_ not GG with $self_name\n";$num_fails++;}
			$$network_ref{$self_name}->add_relative($_,@HAG_likelihoods);
		}
		
		
		## Check if $rel is uncle/aunt of $self
		if(@P1_parents > 0)
		{
			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
			#print "P1 @P1_parents[0] children: @G1_1_children\n";
			#print "P1 @P1_parents[1] children: @G1_2_children\n";
	 
			foreach(@G1_1_children)
			{
				my $child = $_;
				if(grep($_ eq $child,@G1_2_children)) #uncle/aunt
				{
					#print "$child : $self_name\n";
					if($child eq @parents[0]){next;}
				push(@UNCLES_AUNTS,$child);
					if($child =~ /Dummy/i){next;}
					$$network_ref{$self_name}->add_relative($child,@HAG_likelihoods);
				}
			}
		}
		if(@P2_parents > 0)
		{
			#print "$self_name parent @parents[1] parents @P2_parents\n";
			my @G2_1_children = $$network_ref{@P2_parents[0]}->children();
			my @G2_2_children = $$network_ref{@P2_parents[1]}->children();
			#print "P2 @P2_parents[0] children: @G2_1_children\n";
			#print "P2 @P2_parents[1] children: @G2_2_children\n";
	 
			foreach(@G2_1_children)
			{
				my $child = $_;
				if(grep($_ eq $child,@G2_2_children)) #Unclde/aunt
				{
					if($child eq @parents[1]){next;}
					push(@UNCLES_AUNTS,$child);
					if($child =~/Dummy/i){next;}
					$$network_ref{$self_name}->add_relative($child,@HAG_likelihoods);
				}
			}
		}
		
		
		## Don't do the phase 3 checks if we are not in Phase 3
		if ($phase < 3){
			return $num_fails;
		}
		
		## Check half uncle/aunts
		if(@P1_parents > 0)
		{
			my @G1_1_children = $$network_ref{@P1_parents[0]}->children();
			my @G1_2_children = $$network_ref{@P1_parents[1]}->children();
			#print "P1 @P1_parents[0] children: @G1_1_children\n";
			#print "P1 @P1_parents[1] children: @G1_2_children\n";
	 
			foreach(@G1_1_children)
			{
				my $child = $_;
				if(!grep($_ eq $child,@G1_2_children)) #uncle/aunt
				{
					#print "$child : $self_name\n";
					if($child =~ /Dummy/i || $child eq @parents[0]){next;}
					#if($rels{$child} ne "CGH"){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";$num_fails++;}
					$$network_ref{$self_name}->add_relative($child,@CGH_likelihoods);
				}
			}
			foreach(@G1_2_children)
			{
				my $child = $_;
				if(!grep($_ eq $child,@G1_1_children)) #uncle/aunt
				{
					#print "$child : $self_name\n";
					if($child =~ /Dummy/i || $child eq @parents[0]){next;}
					#if($rels{$child} ne "CGH"){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";$num_fails++;}
					$$network_ref{$self_name}->add_relative($child,@CGH_likelihoods);
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
				if(!grep($_ eq $child,@G2_2_children)) #Unclde/aunt
				{
					if($child =~/Dummy/i || $_ eq @parents[1]){next;}
					#if($rels{$child} ne "CGH"){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";$num_fails++;}
					$$network_ref{$self_name}->add_relative($child,@CGH_likelihoods);
				}
			}
			foreach(@G2_2_children)
			{
				my $child = $_;
				if(!grep($_ eq $child,@G2_1_children)) #Unclde/aunt
				{
					if($child =~/Dummy/i || $_ eq @parents[1]){next;}
					#if($rels{$child} ne "CGH"){print "FAIL!!! $child is not half uncle/aunt of $self_name\n";$num_fails++;}
					$$network_ref{$self_name}->add_relative($child,@CGH_likelihoods);
				}
			}
		}
			
		## Check if $rel is half neice/nephew of $self
		foreach(@HS)
		{
			my @HS_children = $$network_ref{$_}->children();
			foreach(@HS_children) # These are neices/nephews
			{
				if($_ =~ /Dummy/i){next;}
				#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not half neices/nephew with $self_name\n";$num_fails++;}
				$$network_ref{$self_name}->add_relative($_,@CGH_likelihoods);
			}
		}
			
			
		## Check if $rel is GREAT-grandparent of $self
		my @grandparents;
		push(@grandparents, @P1_parents);
		push(@grandparents, @P2_parents);
		#my @P1_parents = $$network_ref{@parents[0]}->parents();
		#my @P2_parents = $$network_ref{@parents[1]}->parents();
		#print "@parents[0] parents: @P1_parents\n";
		#print "@parents[1] parents: @P2_parents\n";
		
		#print "$self_name grandparents: @grandparents\n";
		foreach(@grandparents)
		{
			my @GG_parents = $$network_ref{$_}->parents();
			foreach(@GG_parents)
			{
				if($_ =~/Dummy/i){next;}
				#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not GG-parent of $self_name\n";$num_fails++;}
				$$network_ref{$self_name}->add_relative($_,@CGH_likelihoods);
			}
		}
		
		## CHECK if COUSINS
		if($self_name eq "E" && @UNCLES_AUNTS[0] eq "F")
		{
			#print "$self_name\'s uncles_aunts: @UNCLES_AUNTS#######################################\n";
		}
		foreach(@UNCLES_AUNTS)
		{
			my @cousins = $$network_ref{$_}->children();
			foreach(@cousins)
			{
				if($_ =~ /Dummy/i){next;}
					#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not cousin with $self_name\n";$num_fails++;}
					$$network_ref{$self_name}->add_relative($_,@CGH_likelihoods);
			}
		}
		
		## CHECK IF REL IS GREAT UNCLE/AUNT of self
		my @G_parents;
		push(@G_parents, @P1_parents);
		push(@G_parents, @P2_parents);
		
		foreach my $G_parent (@G_parents)
		{
			my @GG_parents = $$network_ref{$G_parent}->parents();
			if(@GG_parents eq 0){next;}
			
			my @GG1_children = $$network_ref{@GG_parents[0]}->children();
			my @GG2_children = $$network_ref{@GG_parents[1]}->children();
			#print "GG1 @GG_parents[0] children: @GG1_children\n";
			#print "GG2 @GG_parents[1] children: @GG2_children\n";
	 
			foreach my $child (@GG1_children)
			{
				if(grep($_ eq $child,@GG2_children)) #great uncle/aunt
				{
					#print "$child : $self_name\n";
					if($child =~ /Dummy/i || $child eq $G_parent){next;}
					#if($rels{$child} ne "CGH"){print "FAIL!!! $self_name not great neice or nephew of $child\n";$num_fails++;}
					$$network_ref{$self_name}->add_relative($child,@CGH_likelihoods);
				}
			}
		}
				
		## Check if rel is great neice or nephew of self
		foreach(@FS)
		{
			my @FS_children = $$network_ref{$_}->children();
			foreach(@FS_children) # These are neices/nephews
			{
				my @FS_Gchildren = $$network_ref{$_}->children();
				foreach(@FS_Gchildren) # These are neices/nephews
				{
					if($_ =~ /Dummy/i){next;}
					#if($rels{$_} ne "CGH"){print "FAIL!!! $_ not great neices/nephew with $self_name\n";$num_fails++;}
					$$network_ref{$self_name}->add_relative($_,@CGH_likelihoods);
				}
			}
		}	
	} ## Oct 2, 2012 moved this bracket here from after "get all unrelateds" section below => needs to be tested


	## Check for more distant relationships
	my %all_relatives = $self->get_all_relatives($network_ref);
	my %current_rels = $self->relatives();

	foreach(keys %all_relatives)
	{
		if(!exists $current_rels{$_} && $_ ne $self_name)
		{
			if($_ =~ /Dummy/i){next;}
			$$network_ref{$self_name}->add_relative($_,@DR_likelihoods);
		}
	}

	## get all unrelateds
	my %all_pedigree = $self->get_subpedigree_names($network_ref);

	foreach(keys %all_pedigree)
	{
		if(!exists $all_relatives{$_})
		{
			if($_ =~ /Dummy/i){next;}
			$$network_ref{$self_name}->add_relative($_,@UN_likelihoods);
		}
	}
	return $num_fails;
}

sub is_related_to
{
	my $self = shift;
	my $network_ref = shift;
	my $rel_name = shift;
	my %relatives = $self->get_all_relatives($network_ref);
	my $are_related = 0;
		
	if(exists $relatives{$rel_name})
	{
		$are_related = 1;
	}
	
	return $are_related;
}

sub get_all_relatives
{
	my $self = shift;
	my $network_ref = shift;
	my $self_name = $self->name();
	my @parents_to_explore = @{$self->{PARENTS} };
	my @children_to_explore = @{$self->{CHILDREN} };
	my %relatives; 
	$relatives{$self_name} = 1;
	#print "\nFinding rels $self_name\n";
	while(my $name = shift(@parents_to_explore) )
	{
		if(exists $relatives{$name})
		{
			next;
		}
		$relatives{$name} = 1;
		my @parents = $$network_ref{$name}->parents();
		push(@parents_to_explore,@parents);
		
		my @children = $$network_ref{$name}->children();
		foreach(@children)
		{
			if($_ ne "$name"){push(@children_to_explore,@children);}
		}
	}
	while(my $name = shift(@children_to_explore) )
	{
		if(exists $relatives{$name})
		{
			next;
		}
		$relatives{$name} = 1;
		my @children = $$network_ref{$name}->children();
		push(@children_to_explore,@children);
	}		
	
	return %relatives;
}

## Check if anyone is mating with someone else outside the permitted generation gap
sub pass_generation_check
{
	#print "generation_check\n";
	my $self = shift;
	my $network_ref = shift;
	
	my $self_name = $self->name();
	my @nodes_to_explore;
	my %relatives;
	$relatives{$self_name} = 0;
	push(@nodes_to_explore,$self_name);
	
	#print "\nFinding rels $self_name\n";
	while(my $name = shift(@nodes_to_explore) )
	{
		#print "$name\n";
		if(!exists $$network_ref{$name})
		{
			return 0;
		}
		my @parents = $$network_ref{$name}->parents();
		foreach my $parent (@parents)
		{
			my $expected_generation = $relatives{$name} + 1;
			if(exists $relatives{$parent})
			{
				if($relatives{$parent} > $MAX_GEN_GAP + $expected_generation || $relatives{$parent} < $expected_generation - $MAX_GEN_GAP)
				{
					#print "$parent is generation ". $relatives{$parent} . " and " . $expected_generation . "\n";
					return 0;
				}
				next;
			}
			push(@nodes_to_explore,$parent);
			$relatives{$parent} = $relatives{$name} + 1;
		}
		
		my @children = $$network_ref{$name}->children();
		foreach my $child(@children)
		{
			my $expected_generation = $relatives{$name} - 1;
			if(exists $relatives{$child})
			{
				if($relatives{$child} > $MAX_GEN_GAP + $expected_generation || $relatives{$child} < $expected_generation - $MAX_GEN_GAP)
				{
					#print "$child is generation ". $relatives{$child} . " and " . $expected_generation . "\n";
					return 0;
				}
				next;
			}
			push(@nodes_to_explore,$child);
			$relatives{$child} = $relatives{$name} - 1;
		}
		#print "nodes to explore: @nodes_to_explore\n";
	}

	## Make sure that no one is mating with someone outside the permitted generation gap
	foreach my $name (keys %relatives)
	{
		my ($parent1,$parent2) = $$network_ref{$name}->parents();
		if($relatives{$parent1} > $MAX_GEN_GAP + $relatives{$parent2} || $relatives{$parent2} > $MAX_GEN_GAP + $relatives{$parent1})
		{
			print "$parent1 ($relatives{$parent1}) mating with $parent2 ($relatives{$parent2})\n";
			return 0;
		}
	}
	return 1;
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
