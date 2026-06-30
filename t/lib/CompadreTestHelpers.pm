package CompadreTestHelpers;

use strict;
use warnings;
use Exporter 'import';
use File::Temp qw(tempdir);
use File::Spec;
use IPC::Run qw(run);
use IO::Socket::INET;

our @EXPORT = qw(get_free_port run_compadre cleanup_test_output verify_output_exists);

=head1 NAME

CompadreTestHelpers - Utility functions for COMPADRE test suite

=head1 SYNOPSIS

  use lib 't/lib';
  use CompadreTestHelpers;

  my $port = get_free_port();
  my $output_dir = tempdir(CLEANUP => 0);
  
  my $result = run_compadre(
    './examples/input',
    './examples/segments.txt',
    $output_dir,
    $port
  );

  if ($result->{success}) {
    verify_output_exists($output_dir, qw(
      input_cleaned.genome_maximum_independent_set_PLINK
      input_cleaned.genome_networks
    ));
  }

  cleanup_test_output($output_dir);

=head1 FUNCTIONS

=cut

=head2 get_free_port()

Find and return an available TCP port for socket communication.

Starts at port 6000 and increments until an available port is found.
Useful for avoiding conflicts when tests run in parallel.

Returns: integer port number

=cut

sub get_free_port {
    my $socket;
    my $port = 6000;
    my $max_attempts = 100;
    
    for (my $i = 0; $i < $max_attempts; $i++) {
        $socket = IO::Socket::INET->new(
            LocalHost => 'localhost',
            LocalPort => $port,
            Proto     => 'tcp',
            Listen    => 1,
            Reuse     => 0,
            Timeout   => 5,
        );
        
        last if $socket;
    
        $port++;
    }
    
    close($socket) if $socket;
    
    # If we tried up to port number 6099 then at the final iteration we 
    # will get to 6100 right as we exit the for loop, so we check if we 
    # exceeded the max attempts and if so we die with an error message
    if ($port >= 6000 + $max_attempts) {
        die "Could not find free port after $max_attempts attempts";
    }
    return $port;
}

=head2 run_compadre($input_file, $segment_file, $output_dir, $port)

Execute COMPADRE process with given parameters.

Parameters:
  - $input_file: Path to input file (without extension, e.g., './examples/input')
  - $segment_file: Path to segments file (e.g., './examples/segments.txt')
  - $output_dir: Directory for output files
  - $port: Port number for socket communication

Returns: hashref with keys:
  - success: boolean (1 if exit code 0)
  - exit_code: integer exit code
  - stdout: captured standard output
  - stderr: captured standard error
  - run_time: seconds elapsed

Dies on fatal errors (command not found, etc.)

=cut

sub run_compadre {
    my ($input_file, $segment_file, $output_dir) = @_;
    
    die "Input file not provided" unless defined $input_file;
    die "Segment file not provided" unless defined $segment_file;
    die "Output directory not provided" unless defined $output_dir;
    
    my @command = (
        'perl',
        'bin/run_COMPADRE.pl',
        '--file', $input_file,
        '--segment_data', $segment_file,
        '--output', $output_dir,
        '--genome',
        '--verbose', '1',
    );
    
    my ($stdout, $stderr);
    my $start_time = time();
    
    eval {
        run \@command, \undef, \$stdout, \$stderr or die "Command failed: $?";
    };
    
    my $end_time = time();
    my $exception = $@;
    my $exit_code = $? >> 8;
    
    return {
        success => !$exception && $exit_code == 0,
        exit_code => $exit_code,
        stdout => $stdout || '',
        stderr => $stderr || '',
        run_time => $end_time - $start_time,
        exception => $exception,
    };
}

=head2 cleanup_test_output($dir)

Recursively remove a test output directory and all its contents.

Parameters:
  - $dir: directory path to remove

Uses File::Temp's cleanup behavior, safely handles nested subdirectories.
Warns but doesn't die if directory doesn't exist.

=cut

sub cleanup_test_output {
    my ($dir) = @_;
    
    return unless defined $dir;
    return unless -d $dir;
    
    # Use system rm for robustness with complex directory structures
    system('rm', '-rf', $dir);
    
    if ($? != 0) {
        warn "Failed to remove test output directory: $dir";
    }
}


=head2 verify_output_exists($dir, @files)

Verify that all expected output files exist in the output directory.

Parameters:
  - $dir: output directory path
  - @files: list of expected file/directory names (relative to $dir)

Returns: 1 if all files exist

Dies with detailed error message if any files are missing.

=cut

sub verify_output_exists {
    my ($dir, @files) = @_;
    
    die "Output directory not provided" unless defined $dir;
    die "No files specified to verify" unless @files;
    
    my @missing;
    
    foreach my $file (@files) {
        my $full_path = File::Spec->catfile($dir, $file);
        
        unless (-e $full_path) {
            push @missing, $file;
        }
    }
    
    if (@missing) {
        my $missing_list = join("\n  - ", @missing);
        die "Expected output files not found in $dir:\n  - $missing_list";
    }
    
    return 1;
}

1;

=head1 AUTHOR

COMPADRE Integration Test Suite

=cut
