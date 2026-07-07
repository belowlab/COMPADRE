#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use Time::HiRes qw(gettimeofday tv_interval);
use Scalar::Util qw(looks_like_number);

# Add lib/perl_modules to path programmatically
use lib File::Spec->catdir('lib', 'perl_modules');
use PRIMUS::IMUS;
use Log::Log4perl;

# Initialize Log4perl so it doesn't complain about missing configuration
Log::Log4perl->init(\"
    log4perl.rootLogger = ERROR, screen
    log4perl.appender.screen = Log::Log4perl::Appender::Screen
    log4perl.appender.screen.stderr = 1
    log4perl.appender.screen.layout = Log::Log4perl::Layout::SimpleLayout
");

# Get cross-platform RSS memory footprint
sub get_process_rss_kb {
    # Get the memory on a linux machine
    if ( -d "/proc" && -f "/proc/self/status" ) {
        if ( open( my $fh, '<', "/proc/self/status" ) ) {
            while ( my $line = <$fh> ) {
                if ( $line =~ /^VmRSS:\s+(\d+)/ ) {
                    close($fh);
                    return int($1);
                }
            }
            close($fh);
        }
    }
    # Fallback to macOS posix shell ps to get memory
    my $pid = $$;
    my $rss = `ps -o rss= -p $pid`;
    if (defined $rss && $rss =~ /(\d+)/) {
        return int($1);
    }
    return 0;
}

my $genome_file = $ARGV[0];
if (!defined $genome_file || !-f $genome_file) {
    die "Usage: perl -Ilib/perl_modules benchmarks/profile_perl_memory.pl <path_to_genome_file>\n";
}

# Configuration representing standard PLINK .genome format
my $config = {
    threshold          => 0.09375, # Standard PRIMUS relatedness threshold
    relatedness_column => 9,       # PI_HAT column (0-indexed)
    id1_column         => 1,       # IID1
    id2_column         => 3,       # IID2
    fid1_column        => 0,       # FID1
    fid2_column        => 2,       # FID2
};

my $state = {
    id_id_scores   => {},
    id_id_all_info => {},
    id_network     => {},
    networks       => {},
    network_ctr    => 0,
};

print "========================================\n";
print "Starting load_data memory benchmark\n";
print "Target File: $genome_file\n";
print "========================================\n";

# Measure baseline memory
my $baseline_mem = get_process_rss_kb();
printf("Baseline Process Memory (RSS): %d KB\n", $baseline_mem);

# Measure execution time
my $t0 = [gettimeofday];

# Load the data
PRIMUS::IMUS::load_data($config, $state, $genome_file);

my $elapsed = tv_interval($t0);
my $post_load_mem = get_process_rss_kb();
my $delta_mem = $post_load_mem - $baseline_mem;

# Count loaded entries
my $total_loaded = keys %{$state->{id_id_all_info}};
my $related_loaded = keys %{$state->{id_id_scores}};

printf("\nData Loading Complete in %.4f seconds.\n", $elapsed);
printf("Post-Load Process Memory (RSS): %d KB\n", $post_load_mem);
printf("Memory Growth (Delta RSS):       %d KB (%.2f MB)\n", $delta_mem, $delta_mem / 1024);
printf("Total Unique Pairs Stored:       %d\n", $total_loaded);
printf("Related Pairs Stored (>= 0.09):  %d\n", $related_loaded);

if ($total_loaded > 0) {
    my $bytes_per_pair = ($delta_mem * 1024) / $total_loaded;
    printf("Approx Memory per Pair:          %.2f bytes\n", $bytes_per_pair);
} else {
    print "No pairs loaded.\n";
}
print "========================================\n";
