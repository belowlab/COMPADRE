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
use Test::More tests => 12;
use Test::Deep;
use lib 'lib/perl_modules';
use lib 't/lib';

# Import IMUS to access its functions and globals
use PRIMUS::IMUS;

#####################################
# Test Helpers
#####################################

sub reset_imus_globals {
    # Reset all IMUS package globals to clean state for a test
    no strict 'refs';
    
    %PRIMUS::IMUS::id_id_scores = ();
    $PRIMUS::IMUS::THRESHOLD = 0.1;
    @PRIMUS::IMUS::trait_refs = ();
    @PRIMUS::IMUS::trait_order = ();
    %PRIMUS::IMUS::trait_files = ();
    %PRIMUS::IMUS::networks = ();
    %PRIMUS::IMUS::id_network = ();
}

# All setup functions return a hash that list all individuals in the graph space. 
# This list of individuals will start as the initial candidate pool for the 
# BronKerbosch algorithm.
sub setup_k3_network {
    # Complete graph: 3 nodes all connected
    # Edges: 1-2, 1-3, 2-3 (all > threshold)
    
    reset_imus_globals();
    
    no strict 'refs';
    %PRIMUS::IMUS::id_id_scores = (
        'ID01;ID02' => 0.25,
        'ID01;ID03' => 0.28,
        'ID02;ID03' => 0.26,
    );
    
    return { ID01 => 1, ID02 => 1, ID03 => 1 };
}

sub setup_k4_network {
    # Complete graph: 4 nodes all connected
    
    reset_imus_globals();
    
    no strict 'refs';
    %PRIMUS::IMUS::id_id_scores = (
        'ID01;ID02' => 0.25,  'ID01;ID03' => 0.28,  'ID01;ID04' => 0.30,
        'ID02;ID03' => 0.26,  'ID02;ID04' => 0.27,
        'ID03;ID04' => 0.29,
    );
    
    return { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1 };
}

sub setup_disconnected_network {
    # No edges: 4 isolated nodes
    
    reset_imus_globals();
    
    no strict 'refs';
    %PRIMUS::IMUS::id_id_scores = ();  # Empty = no relationships
    
    return { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1 };
}

sub setup_bipartite_network {
    # Bipartite graph: 2 groups (5 nodes each) - for finding independent set
    # Group A: nodes 1-5 (unrelated to each other)
    # Group B: nodes 6-10 (unrelated to each other)
    # Between-group: related (score > threshold) - so they don't connect in complement graph
    
    reset_imus_globals();
    
    no strict 'refs';
    
    my %edges = (
        # Between groups only: related (score > threshold)
        # Within-group pairs left undefined (treated as unrelated <= threshold)
        'ID01;ID06' => 0.25,  'ID01;ID07' => 0.28,  'ID01;ID08' => 0.30,  'ID01;ID09' => 0.26,  'ID01;ID10' => 0.27,
        'ID02;ID06' => 0.26,  'ID02;ID07' => 0.27,  'ID02;ID08' => 0.25,  'ID02;ID09' => 0.29,  'ID02;ID10' => 0.28,
        'ID03;ID06' => 0.30,  'ID03;ID07' => 0.26,  'ID03;ID08' => 0.27,  'ID03;ID09' => 0.25,  'ID03;ID10' => 0.29,
        'ID04;ID06' => 0.28,  'ID04;ID07' => 0.29,  'ID04;ID08' => 0.26,  'ID04;ID09' => 0.27,  'ID04;ID10' => 0.25,
        'ID05;ID06' => 0.27,  'ID05;ID07' => 0.25,  'ID05;ID08' => 0.29,  'ID05;ID09' => 0.28,  'ID05;ID10' => 0.26,
    );
    
    %PRIMUS::IMUS::id_id_scores = %edges;
    
    return { ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1, ID05 => 1, ID06 => 1, ID07 => 1, ID08 => 1, ID09 => 1, ID10 => 1 };
}

sub setup_empty_network {

    reset_imus_globals();
    
    no strict 'refs';

    %PRIMUS::IMUS::id_id_scores = ();  # No relationships, empty hash

    return { ID01 =>1, ID02 => 1, ID03 => 1, ID04 => 1};
}


#####################################
# Test BronKerbosch Algorithm
#####################################

# Test 1 and 2: Make sure the BronKerbosch code is identifying 
# 1 clique for Graph K3 (3 nodes, all connected). We check to 
# see how many cliques are identified and then what ifs and in 
# the clique
{
    my $network_ref = setup_k3_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;
    
    # Call the actual BronKerbosh from IMUS
    PRIMUS::IMUS::BronKerbosh(\@maximal_cliques, \%R, \%P, \%X, \$num_visited, \%PRIMUS::IMUS::id_id_scores);
    
    # All scores > threshold means complement graph has no edges
    # So each node is its own clique
    is(scalar(@maximal_cliques), 3, "K3: Three maximal cliques found (one per node)");
}

# Test 3: Complete Graph K4 (4 nodes, all connected)
{
    my $network_ref = setup_k4_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;
    
    PRIMUS::IMUS::BronKerbosh(\@maximal_cliques, \%R, \%P, \%X, \$num_visited, \%PRIMUS::IMUS::id_id_scores);

    
    # All scores > threshold means complement graph has no edges
    is(scalar(@maximal_cliques), 4, "K4: Four maximal cliques found (one per node)");
}

# Test 4: Disconnected Graph (no edges, should find multiple cliques)
{
    my $network_ref = setup_disconnected_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;
    
    
    PRIMUS::IMUS::BronKerbosh(\@maximal_cliques, \%R, \%P, \%X, \$num_visited, \%PRIMUS::IMUS::id_id_scores);
    
 
    
    # NOTE: With completely unrelated nodes (no edges), the algorithm 
    # finds all nodes in one clique due to how inverse neighbors work.
    # This is edge case behavior - the algorithm is designed for highly 
    # connected genetic networks, not completely disconnected graphs.
    is(scalar(@maximal_cliques), 1, "Disconnected graph: All unrelated nodes found in one clique");
}

# Test 5: checking how the algorithm performs for a bipartite 
# graph. In this case our graph has 2 sets of 5 nodes where there # are only connectinos between sets and not within sets. We 
# expect the function to return 2 sets. Nodes 1-5 have no 
# connections to each other and nodes 6-10 have no connections to 
# each other. Therefore we should get 2 sets where 1 contains 
#nodes 1-5 and the other contains nodes 6-10.
{
    my $network_ref = setup_bipartite_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;


    PRIMUS::IMUS::BronKerbosh(\@maximal_cliques, \%R, \%P, \%X, \$num_visited, \%PRIMUS::IMUS::id_id_scores);

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

# We need to check the case where no individuals are related to each other. This means the %id_id_score hash is empty. The complement graph should have all individuals connected to each other, so we should get one big clique.
{

    my $network_ref = setup_empty_network();

    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;

    PRIMUS::IMUS::BronKerbosh(\@maximal_cliques, \%R, \%P, \%X, \$num_visited, \%PRIMUS::IMUS::id_id_scores);

    is(scalar(@maximal_cliques), 1, "Empty graph test: One maximal clique identified when all individuals are unrelated (complement graph fully connected)");

    is(scalar(keys %{ $maximal_cliques[0] }), 4, "Empty graph test: Maximal clique contains all individuals when all are unrelated");
}

# Test 6: Hypothesis test - verify how undefined hash values behave in comparison
# This test validates that undefined hash lookups correctly treated as <= THRESHOLD
# when checking for inverse neighbors (complement graph edges)
{
    reset_imus_globals();
    
    no strict 'refs';
    
    # Setup: only define one edge (1;2 > threshold)
    # All other pairs undefined, should be treated as <= threshold
    %PRIMUS::IMUS::id_id_scores = (
        'ID01;ID02' => 0.25,
    );
    
    my @maximal_cliques;
    my %R = ();
    my %P = ( ID01 => 1, ID02 => 1, ID03 => 1 );
    my %X = ();
    my $num_visited = 0;
    
    PRIMUS::IMUS::BronKerbosh(\@maximal_cliques, \%R, \%P, \%X, \$num_visited, \%PRIMUS::IMUS::id_id_scores);
    
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

# Test 8-10: select_pivot with tied candidates (multiple nodes with same max neighbor count)
# This tests that select_pivot correctly identifies a max-degree node when there are ties
{
    reset_imus_globals();
    
    no strict 'refs';
    
    # Setup: 4 nodes where ID01, ID02, ID03 each have 2 unrelated neighbors
    # (they form a triangle in the complement graph), and ID04 is only related to all three
    # Candidates P:
    #   ID01 unrelated to: ID02, ID03 (2 neighbors in P)
    #   ID02 unrelated to: ID01, ID03 (2 neighbors in P)
    #   ID03 unrelated to: ID01, ID02 (2 neighbors in P)
    #   ID04 unrelated to: none (0 neighbors in P; related to all others)
    # Expected: select_pivot should return one of ID01/ID02/ID03 with 2 neighbors
    
    %PRIMUS::IMUS::id_id_scores = (
        'ID01;ID04' => 0.25,  # ID01 related to ID04
        'ID02;ID04' => 0.25,  # ID02 related to ID04
        'ID03;ID04' => 0.25,  # ID03 related to ID04
        # ID01;ID02, ID01;ID03, ID02;ID03 undefined = unrelated
    );
    
    my %P = ( ID01 => 1, ID02 => 1, ID03 => 1, ID04 => 1 );
    my %X = ();
    
    my ($pivot, %neighbors) = PRIMUS::IMUS::select_pivot(\%P, \%X, \%PRIMUS::IMUS::id_id_scores);
    
    # Verify pivot is one of the tied max-degree nodes
    my $is_max_degree_pivot = ($pivot eq 'ID01' || $pivot eq 'ID02' || $pivot eq 'ID03');
    is($is_max_degree_pivot, 1, "select_pivot: Pivot is one of max-degree nodes (tied at 2 neighbors)");
    
    # Verify neighbor count is exactly 2 (the maximum)
    is(scalar(keys %neighbors), 2, "select_pivot: Pivot has 2 neighbors (maximum degree)");
    
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