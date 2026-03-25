#!/usr/bin/env perl
use Test::More;
use IPC::Cmd qw(can_run);
use lib 't/lib';
use File::Temp qw(tempdir);
use CompadreTestHelpers;
use File::Spec;

# We have to define a port that can be used by compadre
my $port = CompadreTestHelpers::get_free_port();


# Making a temporary output directory
my $temp_output_dir = tempdir(CLEANUP => 0);
# Next three lines define where the input files and comparison files are located    
my $fixtures = File::Spec->catfile('t', 'fixtures');
my $truth_set_dir = File::Spec->catfile($fixtures, 'truth_sets');
my $inputs_dir = File::Spec->catfile($fixtures, 'input');

# This represents the total number of test being run. We may have to update 
# this at some point
my $num_test = 1;

SKIP: { 
    skip "the plink binary was not found in the users PATH. PLINK is required for the testing suite. Please install plink and then rerun the tests", $num_test unless can_run('plink');

    my $result = run_compadre(
        File::Spec->catfile($inputs_dir, "input"),
        File::Spec->catfile($inputs_dir, 'segments.txt'),
        $temp_output_dir,
        $port
    );

    ok($result->{success}, "compadre ran successfully") or diag("compadre failed with error: $result->{stderr}");
    cleanup_test_output($temp_output_dir);
}

done_testing();