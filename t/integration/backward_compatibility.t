#!/usr/bin/env perl
use Test::More;
use IPC::Cmd qw(can_run);
use lib 't/lib';
use File::Temp qw(tempdir);
use CompadreTestHelpers;
use CompareFiles;
use File::Spec;

# We have to define a port that can be used by compadre
my $port = CompadreTestHelpers::get_free_port();


# Making a temporary output directory
my $temp_output_dir = tempdir(CLEANUP => 0);
# Next three lines define where the input files and comparison files are located    
my $fixtures = File::Spec->catfile('t', 'fixtures');
my $truth_set_dir = File::Spec->catfile($fixtures, 'expected_outputs', 'backwards_compat');
my $inputs_dir = File::Spec->catfile($fixtures, 'input');

# This represents the total number of test being run. We may have to update 
# this at some point
my $num_test = 1;

SKIP: { 
    skip "the plink binary was not found in the users PATH. PLINK is required for the testing suite. Please install plink and then rerun the tests", $num_test unless can_run('plink');

    # We are going to first make sure that we can run through all of the compadre test suite
    my $result = run_compadre(
        File::Spec->catfile($inputs_dir, "input"),
        File::Spec->catfile($inputs_dir, 'segments.txt'),
        $temp_output_dir,
        $port
    );
    #TEST1: we want to make sure that the program runs 
    # successfully without crashing. We can do this just by 
    # checking the response code 
    ok($result->{success}, "compadre ran successfully") or diag("compadre failed with error: $result->{stderr}");

    #TEST2: now we want to make sure that the correct output 
    # files were created. Because COMPADRE makes the file 
    # names then we can check and make sure that these exact 
    # files are made.
    my ($files_found, $err_msg) = verify_output_files_made($temp_output_dir,
        "input_cleaned.genome_maximum_independent_set_PLINK",
        "input_cleaned.genome_maximum_independent_set_KING",
        "input_cleaned.genome_maximum_independent_set_PRIMUS",
        "input_cleaned.genome_maximum_independent_set",
        "input_cleaned.genome_networks",
        "input_cleaned.genome_network1",
        "input_cleaned.genome_network1/input_cleaned.genome_network1_1.pdf",
        "input_cleaned.genome_network1/Summary_input_cleaned.genome_network1.txt",
        "Summary_input_cleaned.genome.txt",
    );
    ok($files_found, "all expected output files were created") or diag($err_msg);
    

    #TEST3: check and make sure that the ids in the unrelated_set file are as should be expected.
    my ($success_code, $err_message) = compare_independent_set_files(
        File::Spec->catfile($truth_set_dir, "input_cleaned.genome_maximum_independent_set"),
        File::Spec->catfile($temp_output_dir, "input_cleaned.genome_maximum_independent_set")
    );

    ok($success_code, "checking IMUS unrelated set ids") or diag($err_message);

    #TEST4: check and make sure that the ids in the unrelated_set (PLINK version) file are as should be expected.
    my ($success_code_plink, $err_message_plink) = compare_independent_set_files(
        File::Spec->catfile($truth_set_dir, "input_cleaned.genome_maximum_independent_set_PLINK"),
        File::Spec->catfile($temp_output_dir, "input_cleaned.genome_maximum_independent_set_PLINK")
    );
    ok($success_code_plink, "checking PLINK unrelated set ids") or diag($err_message_plink);

    # TEST5: check and make sure that the ids in the unrelated set (KING version) file are as should be expected.
    my ($success_code_king, $err_message_king) = compare_independent_set_files(
        File::Spec->catfile($truth_set_dir, "input_cleaned.genome_maximum_independent_set_KING"),
        File::Spec->catfile($temp_output_dir, "input_cleaned.genome_maximum_independent_set_KING")
    );
    ok($success_code_king, "checking KING unrelated set ids") or diag($err_message_king);

    # TEST6: check and make sure that hte ids in the unrelated set (PRIMUS version) file are as should be expected.
    my ($success_code_primus, $err_message_primus) = compare_independent_set_files(
        File::Spec->catfile($truth_set_dir, "input_cleaned.genome_maximum_independent_set_PRIMUS"),
        File::Spec->catfile($temp_output_dir, "input_cleaned.genome_maximum_independent_set_PRIMUS")
    );
    ok($success_code_primus, "checking PRIMUS unrelated set ids") or diag($err_message_primus);

    # TEST7: check and make sure that the network files are 
    # correct

    cleanup_test_output($temp_output_dir);
}

done_testing();