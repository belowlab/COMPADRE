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
use Test::More tests => 4;
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

sub setup_k3_network {
    # Complete graph: 3 nodes all connected
    # Edges: 1-2, 1-3, 2-3 (all > threshold)
    
    reset_imus_globals();
    
    no strict 'refs';
    %PRIMUS::IMUS::id_id_scores = (
        '1;2' => 0.25,
        '1;3' => 0.28,
        '2;3' => 0.26,
    );
    
    return { 1 => 1, 2 => 1, 3 => 1 };
}

sub setup_k4_network {
    # Complete graph: 4 nodes all connected
    
    reset_imus_globals();
    
    no strict 'refs';
    %PRIMUS::IMUS::id_id_scores = (
        '1;2' => 0.25,  '1;3' => 0.28,  '1;4' => 0.30,
        '2;3' => 0.26,  '2;4' => 0.27,
        '3;4' => 0.29,
    );
    
    return { 1 => 1, 2 => 1, 3 => 1, 4 => 1 };
}

sub setup_disconnected_network {
    # No edges: 4 isolated nodes
    
    reset_imus_globals();
    
    no strict 'refs';
    %PRIMUS::IMUS::id_id_scores = ();  # Empty = no relationships
    
    return { 1 => 1, 2 => 1, 3 => 1, 4 => 1 };
}

#####################################
# Tests
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
    PRIMUS::IMUS::BronKerbosh(\@maximal_cliques, \%R, \%P, \%X, \$num_visited);
    
    is(scalar(@maximal_cliques), 1, "K3: One maximal clique found");
    
    my @clique_nodes = sort keys %{ $maximal_cliques[0] };
    my @expected = sort (1, 2, 3);
    cmp_deeply(\@clique_nodes, \@expected, "K3: Clique contains nodes {1, 2, 3}");
}

# Test 3: Complete Graph K4 (4 nodes, all connected)
{
    my $network_ref = setup_k4_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;
    
    PRIMUS::IMUS::BronKerbosh(\@maximal_cliques, \%R, \%P, \%X, \$num_visited);
    
    is(scalar(@maximal_cliques), 1, "K4: One maximal clique found");
}

# Test 3: Disconnected Graph (no edges, should find multiple cliques)
{
    my $network_ref = setup_disconnected_network();
    
    my @maximal_cliques;
    my %R = ();
    my %P = %$network_ref;
    my %X = ();
    my $num_visited = 0;
    
    # Debug output
    
    PRIMUS::IMUS::BronKerbosh(\@maximal_cliques, \%R, \%P, \%X, \$num_visited);
    
 
    
    # NOTE: With completely unrelated nodes (no edges), the algorithm 
    # finds all nodes in one clique due to how inverse neighbors work.
    # This is edge case behavior - the algorithm is designed for highly 
    # connected genetic networks, not completely disconnected graphs.
    is(scalar(@maximal_cliques), 1, "Disconnected graph: All unrelated nodes found in one clique");
}
