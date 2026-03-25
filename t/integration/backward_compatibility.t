use Test::More;
use lib 't/lib';
use CompadreTestHelpers;
use File::Spec

# We have to define a port that can be used by compadre
my $port = get_free_port();
# Making a temporary output directory
my $temp_output_dir = tempdir(CLEANUP => 0);
# Next three lines define where the input files and comparison files are located    
my $fixtures = File::Spec->catfile('t', 'fixtures');
my $truth_set_dir = File::Spec->catfile($fixtures, 'truth_sets');
my $input_file = File::Spec->catfile($fixtures, 'input');