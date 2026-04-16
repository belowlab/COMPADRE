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
use Test::More tests => 15;
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