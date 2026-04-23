#!/opt/homebrew/Caskroom/mambaforge/base/envs/compadre_env_test/bin/perl

#####################################
# This file contains tests for the clique algorithms used in PRIMUS.
#
# These tests check things like the BronKerbosch algorithm (add more later) 
# to make sure that networks are being identified correctly.
#
# Tests can be run using the prove utility from the root of the COMPADRE 
# directory
#####################################

use strict;
use warnings;
use Test::More tests => 31;
use Test::Deep;
use lib 'lib/perl_modules';
use lib 't/lib';

# Import IMUS to access its functions and globals
use PRIMUS::IMUS;
use Types::IMUS_types;

#####################################
# Test Helpers
#####################################

# All setup functions return a tuple of ($config, $state, $network_ref) for use with refactored functions.
# This list of individuals will start as the initial candidate pool for the 
# BronKerbosch algorithm.
sub setup_k3_network {
    # Complete graph: 3 nodes all connected
    # Edges: 1-2, 1-3, 2-3 (all > threshold)
    
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {
            'ID01;ID02' => 0.25,
            'ID01;ID03' => 0.28,
            'ID02;ID03' => 0.26,
        }
    );
    
    return ($config, $state, { ID01 => 1, ID02 => 1, ID03 => 1 });
}

sub setup_k4_network {
    # Complete graph: 4 nodes all connected
    
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {
            'ID01;ID02' => 0.25,  'ID01;ID03' => 0.28,  'ID01;ID04' => 0.30,
            'ID02;ID03' => 0.26,  'ID02;ID04' => 0.27,
            'ID03;ID04' => 0.29,
        }
    );
    
    return ($config, $state, { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1 });
}

sub setup_disconnected_network {
    # No edges: 4 isolated nodes
    
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {}
    );  # Empty = no relationships
    
    return ($config, $state, { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1 });
}

sub setup_bipartite_network {
    # Bipartite graph: 2 groups (5 nodes each) - for finding independent set
    # Group A: nodes 1-5 (unrelated to each other)
    # Group B: nodes 6-10 (unrelated to each other)
    # Between-group: related (score > threshold) - so they don't connect in complement graph
    
    my %edges = (
        # Between groups only: related (score > threshold)
        # Within-group pairs left undefined (treated as unrelated <= threshold)
        'ID01;ID06' => 0.25,  'ID01;ID07' => 0.28,  'ID01;ID08' => 0.30,  'ID01;ID09' => 0.26,  'ID01;ID10' => 0.27,
        'ID02;ID06' => 0.26,  'ID02;ID07' => 0.27,  'ID02;ID08' => 0.25,  'ID02;ID09' => 0.29,  'ID02;ID10' => 0.28,
        'ID03;ID06' => 0.30,  'ID03;ID07' => 0.26,  'ID03;ID08' => 0.27,  'ID03;ID09' => 0.25,  'ID03;ID10' => 0.29,
        'ID04;ID06' => 0.28,  'ID04;ID07' => 0.29,  'ID04;ID08' => 0.26,  'ID04;ID09' => 0.27,  'ID04;ID10' => 0.25,
        'ID05;ID06' => 0.27,  'ID05;ID07' => 0.25,  'ID05;ID08' => 0.29,  'ID05;ID09' => 0.28,  'ID05;ID10' => 0.26,
    );
    
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => \%edges
    );
    
    return ($config, $state, { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1, ID05 => 1, ID06 => 1, ID07 => 1, ID08 => 1, ID09 => 1, ID10 => 1 });
}

sub setup_empty_network {

    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {}
    );  # No relationships, empty hash

    return ($config, $state, { ID01 =>1, ID02 => 1, ID03 => 1, ID04 => 1});
}


#####################################
# Test BronKerbosch Algorithm
#####################################

# Test 1 Make sure the BronKerbosch code is identifying 
# 1 clique for Graph K3 (3 nodes, all connected). We check to 
# see how many cliques are identified and then what ifs and in 
# the clique
{
    my ($config, $state, $network_ref) = setup_k3_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;
    
    # Call the actual BronKerbosh from IMUS
    PRIMUS::IMUS::BronKerbosh($config, $state, \@maximal_cliques, \%R, \%P, \%X, \$num_visited);
    
    # All scores > threshold means complement graph has no edges
    # So each node is its own clique
    is(scalar(@maximal_cliques), 3, "K3: Three maximal cliques found (one per node)");
}

# Test 2: Complete Graph K4 (4 nodes, all connected)
{
    my ($config, $state, $network_ref) = setup_k4_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;
    
    PRIMUS::IMUS::BronKerbosh($config, $state, \@maximal_cliques, \%R, \%P, \%X, \$num_visited);

    
    # All scores > threshold means complement graph has no edges
    is(scalar(@maximal_cliques), 4, "K4: Four maximal cliques found (one per node)");
}

# Test 3: Disconnected Graph (no edges, should find multiple cliques)
{
    my ($config, $state, $network_ref) = setup_disconnected_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;
    
    
    PRIMUS::IMUS::BronKerbosh($config, $state, \@maximal_cliques, \%R, \%P, \%X, \$num_visited);
    
 
    
    # NOTE: With completely unrelated nodes (no edges), the algorithm 
    # finds all nodes in one clique due to how inverse neighbors work.
    # This is edge case behavior - the algorithm is designed for highly 
    # connected genetic networks, not completely disconnected graphs.
    is(scalar(@maximal_cliques), 1, "Disconnected graph: All unrelated nodes found in one clique");
}

# Test 4,5,6: checking how the algorithm performs for a bipartite 
# graph. In this case our graph has 2 sets of 5 nodes where there # are only connectinos between sets and not within sets. We 
# expect the function to return 2 sets. Nodes 1-5 have no 
# connections to each other and nodes 6-10 have no connections to 
# each other. Therefore we should get 2 sets where 1 contains 
#nodes 1-5 and the other contains nodes 6-10.
{
    my ($config, $state, $network_ref) = setup_bipartite_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;


    PRIMUS::IMUS::BronKerbosh($config, $state, \@maximal_cliques, \%R, \%P, \%X, \$num_visited);

    # Just verify we're finding the expected 2 cliques from the bipartite structure
    is (scalar(@maximal_cliques), 2, "Bipartite graph: Two maximal cliques found");
    
    # Verify correct nodes are in each clique (the unrelated sets)
    my @clique1_nodes = sort keys %{ $maximal_cliques[0] };
    my @clique2_nodes = sort keys %{ $maximal_cliques[1] };
    
    my @group_a_expected = qw(ID01 ID02 ID03 ID04 ID05);
    my @group_b_expected = qw(ID06 ID07 ID08 ID09 ID10);
    
    # Check if clique 1 matches group A or group B
    my $clique1_is_group_a = (@clique1_nodes == @group_a_expected && 
                              join(",", @clique1_nodes) eq join(",", @group_a_expected));
    
    if ($clique1_is_group_a) {
        cmp_deeply(\@clique1_nodes, \@group_a_expected, "Bipartite graph: Clique 1 contains unrelated set (ID01-ID05)");
        cmp_deeply(\@clique2_nodes, \@group_b_expected, "Bipartite graph: Clique 2 contains unrelated set (ID06-ID10)");
    } else {
        cmp_deeply(\@clique1_nodes, \@group_b_expected, "Bipartite graph: Clique 1 contains unrelated set (ID06-ID10)");
        cmp_deeply(\@clique2_nodes, \@group_a_expected, "Bipartite graph: Clique 2 contains unrelated set (ID01-ID05)");
    }
}

# Test 7,8: We need to check the case where no individuals are related to each other. This means the %id_id_score hash is empty. The complement graph should have all individuals connected to each other, so we should get one big clique.
{

    my ($config, $state, $network_ref) = setup_empty_network();

    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;

    PRIMUS::IMUS::BronKerbosh($config, $state, \@maximal_cliques, \%R, \%P, \%X, \$num_visited);

    is(scalar(@maximal_cliques), 1, "Empty graph test: One maximal clique identified when all individuals are unrelated (complement graph fully connected)");

    is(scalar(keys %{ $maximal_cliques[0] }), 4, "Empty graph test: Maximal clique contains all individuals when all are unrelated");
}

# Test 9: Hypothesis test - verify how undefined hash values behave in comparison
# This test validates that undefined hash lookups correctly treated as <= THRESHOLD
# when checking for inverse neighbors (complement graph edges)
{
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {
            'ID01;ID02' => 0.25,
        }
    );
    
    my @maximal_cliques;
    my %R = ();
    my %P = ( ID01 => 1, ID02 => 1, ID03 => 1 );
    my %X = ();
    my $num_visited = 0;
    
    PRIMUS::IMUS::BronKerbosh($config, $state, \@maximal_cliques, \%R, \%P, \%X, \$num_visited);
    
    # In complement graph:
    # - Edge 1-2 does NOT exist (score > threshold)
    # - All other edges DO exist (score <= threshold, including undefined)
    # So valid cliques: {1,3}, {2,3}, or {1,2,3} if they're isolated in complement
    # We expect node 3 to connect with both 1 and 2 in complement (one or both cliques)
    
    my $has_node_3 = 0;
    for my $clique_ref (@maximal_cliques) {
        $has_node_3++ if exists $clique_ref->{ID03};
    }
    
    is($has_node_3 > 0, 1, "Undefined scores: Node ID03 found in clique (treated as unrelated to ID01 and ID02)");
}

############################
# Test select_pivot function
############################

{
    # Test 10: Check to make sure that there are no neighbors if all individuals are related to each other. In this case the complement graph has no edges, so every node should have 0 neighbors in the complement graph. We can check that the pivot selection correctly identifies that there are no neighbors.
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {
            'ID01;ID02' => 0.25,
            'ID01;ID03' => 0.28,
            'ID02;ID03' => 0.26,
        }
    );

    my %P = ( ID01 => 1, ID02 => 1, ID03 => 1 );
    my %X = ();

    my ($pivot, %neighbors) = PRIMUS::IMUS::select_pivot($config, $state, \%P, \%X);

    my $no_neighbors_found = (scalar(keys %neighbors) == 0);    

    is($no_neighbors_found, 1, "select_pivot test: No neighbors found when all individuals are related (complete graph)");
}

{
    # Test 11: Test to make sure that if all individuals are unrelated to each other then 
    # returned group of neighbors includes all other individuals in the candidate set.
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {}
    );

    my %P = ( ID01 => 1, ID02 => 1, ID03 => 1 );
    my %X = ();

    my ($pivot, %neighbors) = PRIMUS::IMUS::select_pivot($config, $state, \%P, \%X);

    my $two_neighbors_found = (scalar(keys %neighbors) == 2);

    is ($two_neighbors_found, 1, "select_pivot test: Two neighbors found when all individuals are unrelated (no edges)");
}

{
    # Test 12: Test to make sure that the correct pivot is being returned.
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {
            'ID02;ID03' => 0.25
        }
    );

    my %P = ( ID01 => 1, ID02 => 1, ID03 => 1 );
    my %X = ();

    my ($pivot, %neighbors) = PRIMUS::IMUS::select_pivot($config, $state, \%P, \%X);

    my $correct_pivot = ($pivot eq 'ID01');
    my $two_neighbors_found = (scalar(keys %neighbors) == 2);

    is($correct_pivot && $two_neighbors_found, 1, "select_pivot test: Correct pivot (ID01) with 2 neighbors when ID01 is unrelated to both ID02 and ID03");
}

{
    # Test 13: Test to make sure that if only 1 individual is in the candidate set then the 
    # correct pivot is returned and the correct number of neighbors are returned
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {}
    );

    my %P = ( ID01 => 1);
    my %X = ();
    my ($pivot, %neighbors) = PRIMUS::IMUS::select_pivot($config, $state, \%P, \%X);

    my $correct_pivot = ($pivot eq 'ID01');
    my $no_neighbors_found = (scalar(keys %neighbors) == 0);

    is($correct_pivot && $no_neighbors_found, 1, "select_pivot test: Single node pivot with no neighbors when only one individual in candidate set");
}

# Test 14-15: select_pivot with tied candidates (multiple nodes with same max neighbor count)
# This tests that select_pivot correctly identifies a max-degree node when there are ties
{
    my $config = PRIMUS::IMUS::Config->new();
    
    my $state = PRIMUS::IMUS::State->new(
        id_id_scores => {
            'ID01;ID04' => 0.25,  # ID01 related to ID04
            'ID02;ID04' => 0.25,  # ID02 related to ID04
            'ID03;ID04' => 0.25,  # ID03 related to ID04
            # ID01;ID02, ID01;ID03, ID02;ID03 undefined = unrelated
        }
    );
    
    my %P = ( ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1 );
    my %X = ();
    
    my ($pivot, %neighbors) = PRIMUS::IMUS::select_pivot($config, $state, \%P, \%X);
    
    # Verify pivot is one of the tied max-degree nodes
    my $is_max_degree_pivot = ($pivot eq 'ID01' || $pivot eq 'ID02' || $pivot eq 'ID03');
    my $two_neighbors_found = (scalar(keys %neighbors) == 2);

    is($is_max_degree_pivot && $two_neighbors_found, 1, "select_pivot: Pivot is one of max-degree nodes with 2 neighbors (tied at 2 neighbors)");
    
    
    # Verify neighbors are from the correct set {ID01, ID02, ID03}
    my $valid_neighbors = 1;
    for my $neighbor (keys %neighbors) {
        if (!($neighbor =~ /^ID0[1-3]$/ && $neighbor ne $pivot)) {
            $valid_neighbors = 0;
            last;
        }
    }
    is($valid_neighbors, 1, "select_pivot: All neighbors are from the tied set (excluding self)");
}

############################
# Test King_method function
############################

# Test 16-18: King_method with bipartite network - verify greedy correctness. There are 2 equal sized independent sets generated by the setup function. We check the following things: 1. Only one of the two groups are returned, 2. the one set consist of 5 individuals. The set consist of either individuals ID1-ID5 or ID6-ID10.
{
    my ($config, $state, $network_ref) = setup_bipartite_network();
    
    my @maximal_cliques = ();
    
    # Call King_method with bipartite network
    PRIMUS::IMUS::King_method($config, $state, \@maximal_cliques, $network_ref);
    
    # King_method returns one independent set
    is(scalar(@maximal_cliques), 1, "King_method: Returns a single independent set");
    
    # In bipartite structure, independent set should contain one complete group (5 nodes)
    my $selected_set_ref = $maximal_cliques[0];
    my $set_size = scalar(keys %$selected_set_ref);
    is($set_size, 5, "King_method bipartite: Selected set contains 5 nodes (one complete group)");
    
    # Verify selected set is either Group A (ID01-ID05) or Group B (ID06-ID10)
    my @selected_ids = sort keys %$selected_set_ref;
    my $is_group_a = join(",", @selected_ids) eq "ID01,ID02,ID03,ID04,ID05";
    my $is_group_b = join(",", @selected_ids) eq "ID06,ID07,ID08,ID09,ID10";
    
    ok($is_group_a || $is_group_b, "King_method bipartite: Selected set is either Group A or Group B");
}

# Test 19: King_method independence property - verify no edges within selected set. 
{
    my ($config, $state, $network_ref) = setup_bipartite_network();
    
    my @maximal_cliques = ();
    
    # Call King_method with bipartite network
    PRIMUS::IMUS::King_method($config, $state, \@maximal_cliques, $network_ref);
    
    my $selected_set_ref = $maximal_cliques[0];
    my @selected_ids = keys %$selected_set_ref;
    
    # Verify independence: no pair in selected set should have score > threshold
    my $is_independent = 1;
    for (my $i = 0; $i < @selected_ids; $i++) {
        for (my $j = $i + 1; $j < @selected_ids; $j++) {
            my $id1 = $selected_ids[$i];
            my $id2 = $selected_ids[$j];
            
            # Create canonical key (lower ID first)
            my $key = ($id1 lt $id2) ? "$id1;$id2" : "$id2;$id1";
            
            # Check if this pair has a score in id_id_scores
            if (exists $state->{id_id_scores}->{$key}) {
                my $score = $state->{id_id_scores}->{$key};
                # If score > threshold, independence is violated
                if ($score > $config->{threshold}) {
                    $is_independent = 0;
                    last;
                }
            }
            # If score undefined or <= threshold, independence holds for this pair
        }
        last unless $is_independent;
    }
    
    ok($is_independent, "King_method: No two nodes in selected set are related (independence property satisfied)");
}

############################
# Test get_maximum_clique function
############################

# Test 20: get_maximum_clique with single clique (edge case)
{
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => ['trait_file_1'],
        trait_files => {'trait_file_1' => 'size'},
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Create 1 clique with 2 individuals
    my %clique1 = (ID01 => 1, ID02 => 1);
    my @cliques = (\%clique1);
    
    # Create trait data: 1 trait with values for both individuals
    my %trait_data = (ID01 => 5, ID02 => 3);
    my @trait_refs = (\%trait_data);
    
    # Call get_maximum_clique
    my $result = PRIMUS::IMUS::get_maximum_clique($config, $state, @cliques, \@trait_refs);
    
    is($result, 0, "get_maximum_clique single clique: Returns index 0");
}

# Test 21: get_maximum_clique with three cliques and high preference
{
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => ['trait_file_1'],
        trait_files => {'trait_file_1' => 'high_qtrait'},
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Create 3 cliques with different trait averages
    my %clique1 = (ID01 => 1, ID02 => 1);  # avg = 5
    my %clique2 = (ID03 => 1, ID04 => 1);  # avg = 7.5 (sum=15)
    my %clique3 = (ID05 => 1, ID06 => 1);  # avg = 15
    my @cliques = (\%clique1, \%clique2, \%clique3);
    
    # Trait data for each clique
    my %trait_data = (
        ID01 => 5, ID02 => 5,        # Clique 1: sum=10, avg=5
        ID03 => 7.5, ID04 => 7.5,    # Clique 2: sum=15, avg=7.5
        ID05 => 15, ID06 => 15,      # Clique 3: sum=30, avg=15
    );
    my @trait_refs = (\%trait_data);
    
    # Call get_maximum_clique
    my $result = PRIMUS::IMUS::get_maximum_clique($config, $state, @cliques, \@trait_refs);
    
    is($result, 2, "get_maximum_clique high preference: Selects clique with highest average (index 2)");
}

# Test 22: get_maximum_clique with three cliques and low preference
{
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => ['trait_file_1'],
        trait_files => {'trait_file_1' => 'low_qtrait'},
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Create 3 cliques (same as Test 21)
    my %clique1 = (ID01 => 1, ID02 => 1);  # avg = 5
    my %clique2 = (ID03 => 1, ID04 => 1);  # avg = 7.5 (sum=15)
    my %clique3 = (ID05 => 1, ID06 => 1);  # avg = 15
    my @cliques = (\%clique1, \%clique2, \%clique3);
    
    # Trait data for each clique
    my %trait_data = (
        ID01 => 5, ID02 => 5,        # Clique 1: sum=10, avg=5
        ID03 => 7.5, ID04 => 7.5,    # Clique 2: sum=15, avg=7.5
        ID05 => 15, ID06 => 15,      # Clique 3: sum=30, avg=15
    );
    my @trait_refs = (\%trait_data);
    
    # Call get_maximum_clique
    my $result = PRIMUS::IMUS::get_maximum_clique($config, $state, @cliques, \@trait_refs);
    
    is($result, 0, "get_maximum_clique low preference: Selects clique with lowest average (index 0)");
}

# Test 23: get_maximum_clique with multiple traits and priority ordering
{
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => ['trait_file_1', 'trait_file_2'],
        trait_files => {
            'trait_file_1' => 'high_qtrait',
            'trait_file_2' => 'low_qtrait',
        },
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Create 3 cliques
    my %clique1 = (ID01 => 1, ID02 => 1);
    my %clique2 = (ID03 => 1, ID04 => 1);
    my %clique3 = (ID05 => 1, ID06 => 1);
    my @cliques = (\%clique1, \%clique2, \%clique3);
    
    # Trait data: 2 traits
    my %trait_data_1 = (
        ID01 => 5, ID02 => 5,        # Clique 1: high avg = 5
        ID03 => 8, ID04 => 8,        # Clique 2: high avg = 8 (WINNER for first trait)
        ID05 => 3, ID06 => 3,        # Clique 3: high avg = 3
    );
    my %trait_data_2 = (
        ID01 => 10, ID02 => 10,      # Clique 1: low avg = 10
        ID03 => 8, ID04 => 8,        # Clique 2: low avg = 8
        ID05 => 9, ID06 => 9,        # Clique 3: low avg = 9
    );
    my @trait_refs = (\%trait_data_1, \%trait_data_2);
    
    # Call get_maximum_clique
    my $result = PRIMUS::IMUS::get_maximum_clique($config, $state, @cliques, \@trait_refs);
    
    is($result, 1, "get_maximum_clique multiple traits: First trait priority selects index 1 (highest high value)");
}

# Test get_highest_degree_node function
############################

# Test 24-25: get_highest_degree_node with single maximum degree
{
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => ['trait_file_1'],
        trait_files => {'trait_file_1' => 'size'},
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Set up trait data
    my %trait_data = (ID01 => 10, ID02 => 5, ID03 => 8);
    $state->{trait_refs} = [\%trait_data];
    
    # Create degrees: ID01 has degree 3 (highest), others have 1-2
    my %degrees = (ID01 => 3, ID02 => 2, ID03 => 1);
    
    # Call get_highest_degree_node
    my ($node, $degree) = PRIMUS::IMUS::get_highest_degree_node($config, $state, \%degrees);
    
    is($node, 'ID01', "get_highest_degree_node single max: Returns node with highest degree (ID01 with degree 3)");
    
    is($degree, 3, "get_highest_degree_node single max: Returns correct degree value (3)");
}

# Test 26: get_highest_degree_node with tied degrees using trait tiebreaker
{
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => ['trait_file_1'],
        trait_files => {'trait_file_1' => 'high_qtrait'},  # Higher trait wins
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Set up trait data: ID01 and ID02 are tied in degree, but ID02 has better traits
    my %trait_data = (ID01 => 5, ID02 => 10, ID03 => 3);
    $state->{trait_refs} = [\%trait_data];
    
    # Create degrees: ID01 and ID02 both have degree 2 (tie)
    my %degrees = (ID01 => 2, ID02 => 2, ID03 => 1);
    
    # Call get_highest_degree_node
    my ($node, $degree) = PRIMUS::IMUS::get_highest_degree_node($config, $state, \%degrees);
    
    is($node, 'ID02', "get_highest_degree_node tied degrees: Tiebreaker selects node with better traits (ID02 with trait 10 > 5)");
}

# Test 27: get_highest_degree_node with low trait preference (lower value should win)
{
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => ['trait_file_1'],
        trait_files => {'trait_file_1' => 'low_qtrait'},  # Lower trait wins (prefer to remove)
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Set up trait data: ID01 and ID03 are tied in degree
    # With low_qtrait, ID03 with lower value (3) should be preferred for removal
    my %trait_data = (ID01 => 8, ID02 => 4, ID03 => 3);
    $state->{trait_refs} = [\%trait_data];
    
    # Create degrees: ID01 and ID03 both have degree 2 (tie)
    my %degrees = (ID01 => 2, ID02 => 1, ID03 => 2);
    
    # Call get_highest_degree_node
    my ($node, $degree) = PRIMUS::IMUS::get_highest_degree_node($config, $state, \%degrees);
    
    is($node, 'ID03', "get_highest_degree_node low preference tiebreaker: Selects node with lower trait value for removal (ID03 with trait 3)");
}

# Test get_connected_components function
############################

# Test 28: Single connected component - all 5 nodes connected
{
    # Build a fully connected network (all nodes connected to each other)
    my %neighbors = (
        ID01 => { ID02 => 1, ID03 => 1, ID04 => 1, ID05 => 1 },
        ID02 => { ID01 => 1, ID03 => 1, ID04 => 1, ID05 => 1 },
        ID03 => { ID01 => 1, ID02 => 1, ID04 => 1, ID05 => 1 },
        ID04 => { ID01 => 1, ID02 => 1, ID03 => 1, ID05 => 1 },
        ID05 => { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1 },
    );
    
    my $network_ref = { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1, ID05 => 1 };
    
    # Call get_connected_components
    my @components = PRIMUS::IMUS::get_connected_components($network_ref, \%neighbors);
    
    is(scalar(@components), 1, "Connected component single: Returns 1 component for fully connected network (5 nodes)");
}

# Test 29: Multiple components - 3 isolated nodes + 2 connected nodes
{
    # Build network with 3 isolated nodes and 2 that are connected to each other
    my %neighbors = (
        ID01 => { ID02 => 1 },        # ID01 connected to ID02
        ID02 => { ID01 => 1 },        # ID02 connected to ID01
        ID03 => {},                   # ID03 isolated
        ID04 => {},                   # ID04 isolated
        ID05 => {},                   # ID05 isolated
    );
    
    my $network_ref = { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1, ID05 => 1 };
    
    # Call get_connected_components
    my @components = PRIMUS::IMUS::get_connected_components($network_ref, \%neighbors);
    
    is(scalar(@components), 4, "Connected components multiple: Returns 4 components (1 pair + 3 isolated)");
}

# Test 30: All isolated nodes - 5 nodes with no edges
{
    # Build network where all nodes are isolated (no connections)
    my %neighbors = (
        ID01 => {},
        ID02 => {},
        ID03 => {},
        ID04 => {},
        ID05 => {},
    );
    
    my $network_ref = { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1, ID05 => 1 };
    
    # Call get_connected_components
    my @components = PRIMUS::IMUS::get_connected_components($network_ref, \%neighbors);
    
    is(scalar(@components), 5, "Connected components all isolated: Returns 5 components (one per node)");
}

# Test 31: Bipartite structure - 2 groups with connections only within groups
{
    # Build bipartite network: Group A (ID01-ID03) connected within, Group B (ID04-ID05) connected within
    # No connections between groups
    my %neighbors = (
        ID01 => { ID02 => 1, ID03 => 1 },     # Group A: ID01 connects to ID02, ID03
        ID02 => { ID01 => 1, ID03 => 1 },     # Group A: ID02 connects to ID01, ID03
        ID03 => { ID01 => 1, ID02 => 1 },     # Group A: ID03 connects to ID01, ID02
        ID04 => { ID05 => 1 },                # Group B: ID04 connects to ID05
        ID05 => { ID04 => 1 },                # Group B: ID05 connects to ID04
    );
    
    my $network_ref = { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1, ID05 => 1 };
    
    # Call get_connected_components
    my @components = PRIMUS::IMUS::get_connected_components($network_ref, \%neighbors);
    
    is(scalar(@components), 2, "Connected components bipartite: Returns 2 components (one for each group)");
}