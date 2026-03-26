package CompareFiles;

use strict;
use warnings;
use Exporter 'import';
use File::Spec;

our @EXPORT = qw(verify_output_files_made compare_independent_set_files);

=head1 Functions

=cut

=head2 extract_individual_ids($file_path) 

Open provided file and extract the ids into an array 
and then sort the values

Parameters:
    - $file_path: the path to the file that contains the ids. The function assumes that the first column of the file contains the ids. We assume that this file has a header line with FID and IID as the first two columns.

Returns: a sorted array of the ids extracted from the file
=cut

sub extract_individual_ids {
    my ($file_path) = @_;
    open my $fh, '<', $file_path or die "Could not open file '$file_path': $!";
    my @ids;
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^(IID|FID)/i; # Skip header lines that start with IID or FID (case-insensitive)
        my @fields = split /\t/, $line;
        push @ids, $fields[0]; # Assuming the first column contains the individual ID
    }
    close $fh;
    return sort @ids;
}

=head2 compare_independent_set_files($true_file, $test_file)

Compare the ids in the truth set to the ids in the output produced by the test. We assume that ids are in the first column of the file and that the file is tab delimited

Parameters:
    - $true_file: the path to the file that contains the true ids. We assume that this file has a header line with FID and IID as the first two columns.
    - $test_file: the path to the file that contains the test ids. We assume that this file has a header line with FID and IID as the first two columns.
    
Returns: a list where the first value is either 1 or 0 (indicating success or failure), and the second value is an error message if any ids do not match.

=cut

sub compare_independent_set_files {
    my ($true_file, $test_file) = @_;
    # Lets extract the ids that we need to compare from 
    # each file
    my @true_ids = extract_individual_ids($true_file);
    my @test_ids = extract_individual_ids($test_file);

    my $success_code = 1;
    my $err_message = "";

    if (scalar(@true_ids) != scalar(@test_ids)) {
        $success_code = 0;
        $err_message = "The number of ids in the true file ($true_file) is different from the number of ids in the test file ($test_file). True file has " . scalar(@true_ids) . " ids, while test file has " . scalar(@test_ids) . " ids.";
    } else {
        for (my $i = 0; $i < scalar(@true_ids); $i++) {
            if ($true_ids[$i] ne $test_ids[$i]) {
                $success_code = 0;
                $err_message = "The ids in the true file ($true_file) and the test file ($test_file) do not match. The first mismatch is at index $i: true id is '$true_ids[$i]', while test id is '$test_ids[$i]'.";
                last;
            }
        }
    }

    return ($success_code, $err_message);
}

=head2 verify_output_files_made($output_dir, @files)

Verify that all expected output files exist in the output directory.

Parameters:
  - $output_dir: directory path where output files should be
  - @files: list of expected file/directory names (relative to $output_dir)

Returns: a list where the first value is either 1 or 0 (indicating success or failure), and the second value is an error message if any files are missing.

Dies with detailed error message listing missing files if any are not found.

=cut

sub verify_output_files_made {
    my ($output_dir, @files) = @_;
    
    my @missing;
    my $return_code = 1;
    my $err_message = "";
    
    foreach my $file (@files) {
        my $full_path = File::Spec->catfile($output_dir, $file);
        
        unless (-e $full_path) {
            push @missing, $file;
        }
    }
    
    if (@missing) {
        my $missing_list = join("\n  - ", @missing);
        $err_message = "Expected output files not found in $output_dir:\n  - $missing_list";
        $return_code = 0;
    }
    
    return ($return_code, $err_message);
}

1;