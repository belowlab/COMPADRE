package CompareFiles;

use strict;
use warnings;
use Exporter 'import';
use File::Spec;

our @EXPORT = qw(verify_output_files_made extract_individual_ids);

sub extract_individual_ids {
    my ($file_path) = @_;
    open my $fh, '<', $file_path or die "Could not open file '$file_path': $!";
    my @ids;
    while (my $line = <$fh>) {
        chomp $line;
        my @fields = split /\t/, $line;
        push @ids, $fields[0]; # Assuming the first column contains the individual ID
    }
    close $fh;
    return sort @ids;
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