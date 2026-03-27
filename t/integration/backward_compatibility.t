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
my $num_test = 9;

# Paths to the reference data required from compadre_data.zip
my $onekg_ref    = File::Spec->catfile('lib', '1KG', '1KG_reference');
my $hapmap3_dir  = File::Spec->catfile('lib', 'hapmap3');
my $kde_data_dir = File::Spec->catfile('lib', 'KDE_data');

my $data_url = 'https://github.com/belowlab/COMPADRE/releases/download/pre-release-0.2.0/compadre_data.zip';
my $data_msg = "Required reference data is missing. Download $data_url and extract it into the lib/ directory.";

SKIP: { 
    skip "the plink binary was not found in the users PATH. PLINK is required for the testing suite. Please install plink and then rerun the tests", $num_test unless can_run('plink');
    skip "1KG reference files not found at $onekg_ref (.bed/.bim/.fam). $data_msg", $num_test
        unless (-e "$onekg_ref.bed" && -e "$onekg_ref.bim" && -e "$onekg_ref.fam");
    skip "hapmap3 directory not found at $hapmap3_dir. $data_msg", $num_test
        unless -d $hapmap3_dir;
    skip "KDE_data directory not found at $kde_data_dir. $data_msg", $num_test
        unless -d $kde_data_dir;

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
    
    # Debug output for understanding what happened
    diag("COMPADRE stdout:\n$result->{stdout}") if $result->{stdout};
    diag("COMPADRE stderr:\n$result->{stderr}") if $result->{stderr};
    diag("COMPADRE exit code: $result->{exit_code}");
    diag("Output directory: $temp_output_dir");
    diag("Files in output directory: " . join(", ", glob("$temp_output_dir/*")));

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
    my ($success_code_genome, $err_message_genome) = compare_genome_file_content(
        File::Spec->catfile($truth_set_dir, "input_cleaned.genome_networks"),
        File::Spec->catfile($temp_output_dir, "input_cleaned.genome_networks")
    );
    ok($success_code_genome, "checking genome network file content") or diag($err_message_genome);

    # TEST8: we now need to make sure the original network file created by the 
    # plink --genome command is the same
    my ($success_code_genome_plink, $err_message_genome_plink) = compare_genome_file_content(
        File::Spec->catfile($truth_set_dir, "prePRIMUS", "input_cleaned.genome"),
        File::Spec->catfile($temp_output_dir, "input_prePRIMUS", "input_cleaned.genome")
    );
    ok($success_code_genome_plink, "checking original genome file from PLINK") or diag($err_message_genome_plink);

    # TEST9: we now need to check the network summary file nad 
    # make sure that each line in this file matches. The 
    # ordering of the file here does matter so we don't need to 
    # perform an sorting
    my ($success_code_summary, $err_message_summary) = compare_network_summary_files(
        File::Spec->catfile($truth_set_dir, "Summary_input_cleaned.genome_network1.txt"),
        File::Spec->catfile($temp_output_dir, "input_cleaned.genome_network1", "Summary_input_cleaned.genome_network1.txt")
    );

    ok($success_code_summary, "checking network summary file content") or diag($err_message_summary);

    cleanup_test_output($temp_output_dir);
}

done_testing();