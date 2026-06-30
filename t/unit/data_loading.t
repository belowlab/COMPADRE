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
use Test::More tests => 23;
use Test::Deep;
use Path::Tiny; # This library we allow us to make temp files that we can use for testing
use lib 'lib/perl_modules';
use lib 't/lib';

use PRIMUS::IMUS;
use Types::IMUS_types;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($OFF);

# Test 1-5: Test that the load data function is properly reading in the data and storing it in the state object when a valid .genome file is provided
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
    
    # Test that all three expected keys are in id_id_all_info
    # For the following 2 test we have the following pairs from the genome file:
    # ID01;ID02, ID01;ID03, ID02;ID03 but they are now being mapped to integers so 
    # we expect the keys to be 1;2 1;3 2;3. Since the file is read sequentially, 
    # then we can expect the ids to be mapped to the same integer
    my @expected_keys = qw(1;2 1;3 2;3);
    my $all_keys_in_all_info = 1;
    foreach my $key (@expected_keys) {
        if (!exists $state->{id_id_all_info}->{$key}) {
            $all_keys_in_all_info = 0;
            last;
        }
    }
    ok($all_keys_in_all_info, "All expected keys (1;2, 1;3, 2;3) are present in id_id_all_info");
    
    # Test that all three expected keys are in id_id_scores
    my $all_keys_in_scores = 1;
    foreach my $key (@expected_keys) {
        if (!exists $state->{id_id_scores}->{$key}) {
            $all_keys_in_scores = 0;
            last;
        }
    }
    ok($all_keys_in_scores, "All expected keys (1;2, 1;3, 2;3) are present in id_id_scores");
    
    # Test that the header line is correctly stored in the state object. (We expect it to have a newline character)
    my $expected_header = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n";
    is($state->{outfile_header}, $expected_header, "Header line correctly parsed and stored in state object");
    
    # We need to test here and make sure that all three individuals were loaded into 
    # their own network singleton by this function
    is(scalar keys %{$state->{networks}}, 3, "There should be 3 singleton networks for the 3 unique individuals in the .genome file");
}

# Test 6: Test invalid and malformed .genome files. Empty/whitespace-only line causes 
# die
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.genome');
    
    # File with header, valid line, then empty line
    my $test_file_str = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n" .
                        "FAM1 ID01 FAM2 ID02 OT 0 1.0 0.0 0.0 0.5 -1 0.75 0.3 2.7\n" .
                        "\n";  # empty line
    
    $temp->spew($test_file_str);
    
    my $config = PRIMUS::IMUS::Config->new(
        id1_column => 1,
        id2_column => 3,
        relatedness_column => 9,
        fid1_column => 0,
        fid2_column => 2,
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Catch the die
    my $error;
    eval {
        PRIMUS::IMUS::load_data($config, $state, $temp);
        1;
    } or $error = $@;
    
    like($error, qr/Detected an empty line in the genome file/, "Dies when empty line is encountered");
}

# Test 7: Insufficient columns causes die. This happens because we check and make sure 
# that the largest required index is a valid value in the line
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.genome');
    
    # File with header, then line missing PI_HAT column (only 8 fields instead of 14)
    my $test_file_str = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n" .
                        "FAM1 ID01 FAM2 ID02 OT 0 1.0 0.0\n";  # only 8 fields
    
    $temp->spew($test_file_str);
    
    my $config = PRIMUS::IMUS::Config->new(
        id1_column => 1,
        id2_column => 3,
        relatedness_column => 9,
        fid1_column => 0,
        fid2_column => 2,
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Catch the die
    my $error;
    eval {
        PRIMUS::IMUS::load_data($config, $state, $temp);
        1;
    } or $error = $@;
    
    like($error, qr/insufficient columns/, "Dies when line has insufficient columns");
}

# Test 8-10: Non-numeric PI_HAT value results in warning and continues
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.genome');
    
    # File with header and line where PI_HAT is invalid (not numeric, not 'nan')
    my $test_file_str = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n" .
                        "FAM1 ID01 FAM2 ID02 OT 0 1.0 0.0 0.0 XYZ -1 0.75 0.3 2.7\n";
    
    $temp->spew($test_file_str);
    
    my $config = PRIMUS::IMUS::Config->new(
        id1_column => 1,
        id2_column => 3,
        relatedness_column => 9,
        fid1_column => 0,
        fid2_column => 2,
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Should NOT die, but should continue processing
    my $error;
    eval {
        PRIMUS::IMUS::load_data($config, $state, $temp);
        1;
    } or $error = $@;
    

    ok(!$error, "Non-numeric PI_HAT does not cause die");
    
    # Verify the line was stored in id_id_all_info (all pairs stored regardless of validity)
    ok(!exists $state->{id_id_all_info}->{"ID01;ID02"}, "Line with invalid PI_HAT is not stored in id_id_all_info");
    
    # Verify the line is NOT stored in id_id_scores (invalid PI_HAT means failed threshold check)
    ok(!exists $state->{id_id_scores}->{"ID01;ID02"}, "Line with invalid PI_HAT is not stored in id_id_scores");
}

# Test 11: PI_HAT value of 'nan' causes die
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.genome');
    
    # File with header and line where PI_HAT is 'nan'
    my $test_file_str = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n" .
                        "FAM1 ID01 FAM2 ID02 OT 0 1.0 0.0 0.0 nan -1 0.75 0.3 2.7\n";
    
    $temp->spew($test_file_str);
    
    my $config = PRIMUS::IMUS::Config->new(
        id1_column => 1,
        id2_column => 3,
        relatedness_column => 9,
        fid1_column => 0,
        fid2_column => 2,
    );
    
    my $state = PRIMUS::IMUS::State->new();
    
    # Should die when PI_HAT is 'nan'
    my $error;
    eval {
        PRIMUS::IMUS::load_data($config, $state, $temp);
        1;
    } or $error = $@;
    
    like($error, qr/Terminating program since PI_HAT of nan/, "Dies when PI_HAT value is 'nan'");
}

{
    # Test when the user provides a higher threshold to make sure the pair doesn't end up in the id_id_scores hash
    my $temp = Path::Tiny->tempfile(SUFFIX => '.genome');

    # Example space separated string where listing 
    # relationships between indiviudals ID1, ID2, and ID3.The PI_HAT value for ID2;ID3 
    # is very low so these individuals should not be included in the id_id_scores hash
    my $test_file_str = "FID1 IID1 FID2 IID2 RT EZ Z0 Z1 Z2 PI_HAT PHE DST PPC RATIO\n" . "FAM1 ID01 FAM2 ID02 OT 0 1.0 0.0 0.0 0.5 -1 0.75 0.3 2.7\n" . "FAM1 ID01 FAM3 ID03 OT 0 0.25 0.5 0.25 0.8 -1 0.7 0.8 1.2\n" . "FAM2 ID02 FAM3 ID03 OT 0 0.4 0.6 0.0 0.03 -1 0.6 1.0 NA\n";

    $temp->spew($test_file_str); # write to teh temp fi

    my $config = PRIMUS::IMUS::Config->new(
        id1_column => 1,
        id2_column => 3,
        relatedness_column => 9,
        fid1_column => 0,
        fid2_column => 2,
        threshold => 0.2
    );
    
    my $state = PRIMUS::IMUS::State->new();

    PRIMUS::IMUS::load_data($config, $state, $temp);

    my $only_2_pairs_in_scores = scalar keys %{$state->{id_id_scores}} == 2;
    my $id02_id03_not_in_scores = !exists $state->{id_id_scores}->{"ID02;ID03"};

    ok($only_2_pairs_in_scores && $id02_id03_not_in_scores, "Only pairs with PI_HAT above the threshold are included in id_id_scores");
}

#################################
## Following test check the 
## load_samples function
#################################
# Test 7: load_samples with .fam format (FID IID ... → extract IID)
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.fam');
    
    # .fam format: FID IID PID MID SEX PHENO ...
    my $test_file_str = "FAM1 ID01 0 0 1 -9\n" .
                        "FAM2 ID02 0 0 2 -9\n" .
                        "FAM3 ID03 0 0 1 -9\n";
    
    $temp->spew($test_file_str);
    
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new();
    
    # Load the samples
    PRIMUS::IMUS::load_samples($config, $state, $temp);
    
    # Verify that 3 individuals were initialized with unique networks
    is(scalar keys %{$state->{id_network}}, 3, ".fam format: 3 samples loaded with unique networks");
}

# Test 8: load_samples with single column format (IID only)
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.txt');
    
    # Single column format: just IID
    my $test_file_str = "ID04\n" .
                        "ID05\n" .
                        "ID06\n";
    
    $temp->spew($test_file_str);
    
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new();
    
    # Load the samples
    PRIMUS::IMUS::load_samples($config, $state, $temp);
    
    # Verify that 3 individuals were initialized
    is(scalar keys %{$state->{id_network}}, 3, "Single column format: 3 samples loaded");
}

# Test 9: load_samples skips header lines
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.fam');
    
    # File with header and data
    my $test_file_str = "FID IID PID MID SEX PHENO\n" .
                        "FAM1 ID07 0 0 1 -9\n" .
                        "FAM2 ID08 0 0 2 -9\n";
    
    $temp->spew($test_file_str);
    
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new();
    
    # Load the samples
    PRIMUS::IMUS::load_samples($config, $state, $temp);
    
    # Verify that only 2 data rows were loaded (header skipped)
    is(scalar keys %{$state->{id_network}}, 2, "Header line skipped correctly");
}

# Test 10: load_samples adds new ids not in original networks
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.fam');
    
    # Because the program reads through line by line, they will 
    # be encountered sequentially. 
    my $test_file_str = "FAM1 ID09 0 0 1 -9\n" .
                        "FAM2 ID10 0 0 2 -9\n" .
                        "FAM3 ID11 0 0 1 -9\n";
    
    $temp->spew($test_file_str);
    
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new();
    
    # Pre-populate with one ID
    $state->{id_network}{"ID09"} = 0;
    push @{ $state->{networks}{0} }, "ID09";
    $state->{network_ctr} = 1;
    
    # Load samples (some already in network, some new)
    PRIMUS::IMUS::load_samples($config, $state, $temp);
    
    # Verify that all 3 samples are now in networks
    is(scalar keys %{$state->{id_network}}, 3, "All 3 samples have networks after load_samples");
    
    # Verify all 3 individuals are actually in the networks array
    my $id09_in_networks = grep { $_ eq "ID09" } @{ $state->{networks}{0} };
    my $id10_in_networks = grep { $_ eq "ID10" } @{ $state->{networks}{1} };
    my $id11_in_networks = grep { $_ eq "ID11" } @{ $state->{networks}{2} };
    
    ok($id09_in_networks && $id10_in_networks && $id11_in_networks, "All 3 individuals are in their respective networks arrays");
}

###################################
## Thest functions check the 
## load_trait_data subroutine
###################################

# Test 18-19: load_trait_data binary trait conversion (1 and 2 → 0 and 1)
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.txt');
    
    # Binary trait file: header + 3 individuals with trait values of 1 and 2
    # Values of 1 should convert to 0, and 2 should convert to 1
    my $test_file_str = "FAM1 ID01 1\n" .
                        "FAM2 ID02 2\n" .
                        "FAM3 ID03 1\n";
    
    $temp->spew($test_file_str);
    my $temp_str = $temp->stringify;
    # Create config with pre-configured column indices
    # (column detection happens in set_values2, not load_trait_data)
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => [$temp_str],
        trait_files => {$temp_str => 'btrait'},
        trait_fid_columns => { $temp_str => 0},
        trait_id_columns => { $temp_str => 1},
        trait_data_columns => { $temp_str => 2},
        exclude_value => -9,
    );
    
    # Create state with id_network populated 
    # (load_trait_data needs this to know which IDs to process)
    my $state = PRIMUS::IMUS::State->new();
    $state->{id_network}{"ID01"} = 0;
    $state->{id_network}{"ID02"} = 1;
    $state->{id_network}{"ID03"} = 2;
    
    # Call load_trait_data
    PRIMUS::IMUS::load_trait_data($config, $state);
    
    # Verify trait_refs array has one hash reference
    is(scalar @{$state->{trait_refs}}, 1, "One trait hash stored in trait_refs");
    
    # Get the trait hash from trait_refs
    my $trait_hash = $state->{trait_refs}->[0];
    
    # Verify binary conversion: 1 → 0, 2 → 1
    my $binary_conversion_correct = 
        $trait_hash->{ID01} == 0 && 
        $trait_hash->{ID02} == 1 && 
        $trait_hash->{ID03} == 0;
    ok($binary_conversion_correct, "Binary trait values correctly converted (1→0, 2→1)");
}

# Test 20: load_trait_data exclude_value conversion to NA
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.txt');
    
    # Trait file with some individuals having -9 (exclude_value) and others with valid values
    my $test_file_str = "FAM1 ID01 100\n" .
                        "FAM2 ID02 -9\n" .
                        "FAM3 ID03 50\n" .
                        "FAM4 ID04 -9\n";
    
    $temp->spew($test_file_str);
    my $temp_str = $temp->stringify;
    
    # Create config with exclude_value set to -9 (default)
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => [$temp_str],
        trait_files => {$temp_str => 'qtrait'},
        trait_fid_columns => { $temp_str => 0},
        trait_id_columns => { $temp_str => 1},
        trait_data_columns => { $temp_str => 2},
        exclude_value => -9,
    );
    
    # Create state with id_network populated
    my $state = PRIMUS::IMUS::State->new();
    $state->{id_network}{"ID01"} = 0;
    $state->{id_network}{"ID02"} = 1;
    $state->{id_network}{"ID03"} = 2;
    $state->{id_network}{"ID04"} = 3;
    
    # Call load_trait_data
    PRIMUS::IMUS::load_trait_data($config, $state);
    
    # Get the trait hash from trait_refs
    my $trait_hash = $state->{trait_refs}->[0];
    
    # Verify exclude_value conversion to NA and normal values preserved
    my $exclude_values_correct = 
        $trait_hash->{ID01} == 100 &&
        $trait_hash->{ID02} eq "NA" &&
        $trait_hash->{ID03} == 50 &&
        $trait_hash->{ID04} eq "NA";
    ok($exclude_values_correct, "Quantitative trait values matching with excluded IDs (-9) converted to NA string");
}

# Test 21: load_trait_data binary trait with excluded individuals
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.txt');
    
    # Binary trait file with a mix of valid values (1, 2) and excluded values (-9)
    my $test_file_str = "FAM1 ID01 2\n" .
                        "FAM2 ID02 -9\n" .
                        "FAM3 ID03 1\n" .
                        "FAM4 ID04 -9\n" .
                        "FAM5 ID05 2\n";
    
    $temp->spew($test_file_str);
    my $temp_str = $temp->stringify;
    
    # Create config with exclude_value set to -9 and binary trait type
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => [$temp_str],
        trait_files => {$temp_str => 'btrait'},
        trait_fid_columns => { $temp_str => 0},
        trait_id_columns => { $temp_str => 1},
        trait_data_columns => { $temp_str => 2},
        exclude_value => -9,
    );
    
    # Create state with id_network populated
    my $state = PRIMUS::IMUS::State->new();
    $state->{id_network}{"ID01"} = 0;
    $state->{id_network}{"ID02"} = 1;
    $state->{id_network}{"ID03"} = 2;
    $state->{id_network}{"ID04"} = 3;
    $state->{id_network}{"ID05"} = 4;
    
    # Call load_trait_data
    PRIMUS::IMUS::load_trait_data($config, $state);
    
    # Get the trait hash from trait_refs
    my $trait_hash = $state->{trait_refs}->[0];
    
    # Verify binary conversion AND exclude_value handling:
    # 2→1, 1→0 for valid values; -9→"NA" for excluded values
    my $binary_and_exclude_correct = 
        $trait_hash->{ID01} == 1 &&
        $trait_hash->{ID02} eq "NA" &&
        $trait_hash->{ID03} == 0 &&
        $trait_hash->{ID04} eq "NA" &&
        $trait_hash->{ID05} == 1;
    ok($binary_and_exclude_correct, "Binary trait values converted (2→1, 1→0) with excluded individuals as NA");
}

# Test 22-23: load_trait_data dies when binary trait contains invalid value (not 1 or 2)
{
    my $temp = Path::Tiny->tempfile(SUFFIX => '.txt');
    
    # Binary trait file with an invalid value (3 is not 1 or 2)
    my $test_file_str = "FID IID TRAIT\n" .
                        "FAM1 ID01 1\n" .
                        "FAM2 ID02 3\n" .
                        "FAM3 ID03 2\n";
    
    $temp->spew($test_file_str);
    my $temp_str = $temp->stringify;
    
    # Create config with binary trait type
    my $config = PRIMUS::IMUS::Config->new(
        trait_order => [$temp_str],
        trait_files => {$temp_str => 'btrait'},
        trait_fid_columns => { $temp_str => 0},
        trait_id_columns => { $temp_str => 1},
        trait_data_columns => { $temp_str => 2},
        exclude_value => -9,
    );
    
    # Create state with id_network populated
    my $state = PRIMUS::IMUS::State->new();
    $state->{id_network}{"ID01"} = 0;
    $state->{id_network}{"ID02"} = 1;
    $state->{id_network}{"ID03"} = 2;
    
    # Verify that load_trait_data dies when encountering invalid binary value
    my $died = 0;
    my $error_msg = "";
    eval {
        PRIMUS::IMUS::load_trait_data($config, $state);
    };
    if ($@) {
        $died = 1;
        $error_msg = $@;
    }
    
    ok($died, "load_trait_data dies when binary trait contains value other than 1 or 2");
    like($error_msg, qr/Binary trait file.*contains a binary value other than 1 or 2/, "Error message mentions binary value requirement");
}