#!/opt/homebrew/Caskroom/mambaforge/base/envs/compadre_env_test/bin/perl

#####################################
# This file contains tests for the data loading used in IMUS.
#
# These tests check things like load_data and load_samples
#
#
# Tests can be run using the prove utility from the root of the COMPADRE 
# directory
#####################################

use strict;
use warnings;
use Test::More tests => 1;
use Test::Deep;
use Path::Tiny; # This library we allow us to make temp files that we can use for testing
use lib 'lib/perl_modules';
use lib 't/lib';

use PRIMUS::IMUS;
use Types::IMUS_types;

# Test 1: Test that the load data function is properly reading in the data and storing it in the state object
{
    # Lets make a temporary file that gets automatically cleaned up
    my $temp = Path::Tiny->tempfile(SUFFIX => '.genome');

    # Example space separated string where listing 
    # relationships between indiviudals ID1, ID2, and ID3.
    my $test_file_str = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n" . "FAM1 ID01 FAM2 ID02 OT 0 1.0 0.0 0.0 0.5 -1 0.75 0.3 2.7\n" . "FAM1 ID01 FAM3 ID03 OT 0 0.25 0.5 0.25 0.8 -1 0.7 0.8 1.2\n" . "FAM2 ID02 FAM3 ID03 OT 0 0.4 0.6 0.0 0.7 -1 0.6 1.0 NA\n";

    $temp->spew($test_file_str); # write to teh temp file
    # Create the config object
    my $config = PRIMUS::IMUS::Config->new(
        id1_column => 1,
        id2_column => 3,
        relatedness_column => 9,
        fid1_column => 0,
        fid2_column => 2,
    );

    # The state object will be what we have to check to make 
    # sure the function is working properly
    my $state = PRIMUS::IMUS::State->new();

    # We can use the $temp file in the load_data function
    PRIMUS::IMUS::load_data($config, $state, $temp);

    is(scalar keys %{$state->{id_id_scores}}, 3, "There should be 3 pairs of individuals in the id_id_scores hash");
    
}