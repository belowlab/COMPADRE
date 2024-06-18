#!/usr/bin/perl
use strict;
use warnings;
use File::Find;

# Directory to search for .ps files
my $directory = '/data100t1/home/grahame/projects/compadre/unified-simulations/analysis/primus-results/oct23/sim106-3/uniform3_size20_sim106-3.genome_network1';

# Subroutine to process each file found
sub process_file {
    if (/\.(ps)$/i) {
        my $file = $File::Find::name;
        (my $output_file = $file) =~ s/\.ps$/.pdf/i;
        my $command = "ps2pdf -dFirstPage=2 $file $output_file";
        system($command) == 0 or warn "Failed to execute: $command\n";
    }
}

# Find all files in the directory
find(\&process_file, $directory);
