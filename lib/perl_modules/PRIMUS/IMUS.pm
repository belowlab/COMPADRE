#! /bin/env perl

#####################################
## Written by Jeff Staples
## PhD Student
## Genome Sciences, UW
## May 8, 2012
## contact: grapas2@uw.edu
#####################################

#	*******************************************************************************
#	*                                                                             *
#	*  Program: PRIMUS			                                      *
#	*  by Jeffrey Staples, Deborah A. Nickerson, and Jennifer E. Below            *
#	*  University of Washington                                                   *
#	*                                                                             *
#	*******************************************************************************

package PRIMUS::IMUS;
use strict;
use Getopt::Long qw(GetOptionsFromArray);
use PRIMUS::predict_relationships_2D;
use Types::IMUS_types;
use Log::Log4perl;

my $useage =
"\n\nUSAGE: $0\t  -input [IBD_file]  -output_dir [output_dir]  -threshold [num; default = 0.1]  -[high|low|mean|tails]/[b|q]trait [trait_file]

\n\nType '$0 -help' for explanation of input\n\n";

my $do_IMUS = 1;
my $do_PR   = 1;
my $verbose;
my $lib_dir;
my $LOG = Log::Log4perl->get_logger(__PACKAGE__);

my $samples_file;
my %arg;
my $relatedness_file;
my $relatedness_file_name;
my $ped_file         = "none";
my $missingness_file = "none";
my $output_dir;
my $ersa_data;
my $network_ctr;        # the next network # to add to %networks hash
my $LOWEST_MAX_NETWORK_SIZE = 60;
our $THRESHOLD = .1;    # default is .1
my $MIN_LIKELIHOOD = .1;
my $EXCLUDE_VALUE  = 0;
my $OUTFILE_HEADER;
my $IBD_file_ref;
my $RELATEDNESS_COLUMN      = -1;
my $ID1_COLUMN              = -1;
my $ID2_COLUMN              = -1;
my $FID1_COLUMN             = -1;
my $FID2_COLUMN             = -1;
my %TRAIT_FID_COLUMNS       = -1;
my %TRAIT_ID_COLUMNS        = -1;
my %TRAIT_DATA_COLUMNS      = -1;
my $PRINT_ALTERNATE_RESULTS = 1;

### Here 'networks' could also be considered related groups or families;
## Hash of arrays: Key = network; Value = array (list) of IIDs in that network
my %networks;
## Regular hash; Key = IID1; Value = network that the individual
my %id_network;
my %id_id_all_info;
my @trait_refs;
my @trait_order;
my %trait_files; ## Key = file name; Value = trait type (quantitative or binary)
my %child_parents;

# @purpose Main entry point for IMUS (Independent set of Maximally Unrelated Samples)
# @param @_ - Hash of configuration parameters (passed to set_values2)
# @return (list) - (output_file_path, list_of_unrelated_ids)
# @side_effects Loads data from files, modifies all global variables, writes output files
# @calls load_samples, load_data, load_trait_data, colapse_networks, breakup_large_networks, write_out_independent_set, compare_alternative_methods

sub run_IMUS {

    $LOG->info("Running IMUS");
    my $config = PRIMUS::IMUS::Config->new();
    my $state = PRIMUS::IMUS::State->new();
    #parseCommandLine(@_);
    reset_values();
    set_values2($config, @_);

    $LOG->info("Relatedness_file: $config->{relatedness_file}");
    $LOG->info("Threshold: $config->{threshold}");
    $LOG->info("Selection criteria are based on the following:");
    foreach (@trait_order) {
        $LOG->info("\t$_ ($trait_files{$_})");
    }

    $LOG->debug("IDENTIFYING FAMILY NETWORKS IN DATA");
    $LOG->debug("Writing network files to $config->{output_dir}/");

    $LOG->debug("Loading data...");
    load_samples($config->{samples_file}) if defined $config->{samples_file} && $config->{samples_file} ne "";

    # He we create an empty hash for the id_id_scores
    my %id_id_scores = ();
    load_data( $config->{relatedness_file}, \$state->{id_id_scores});
    #TODO: Need to update this function to not use global values
    load_trait_data();

    $LOG->debug("done.");

    $LOG->debug("colapsing networks...");
    colapse_networks( \%id_id_scores );

    $LOG->debug("done.");

    foreach my $key ( keys %networks ) {

        #print "key: $key; @{ $networks{$key} }\n";
    }

    write_out_networks( $relatedness_file_name, \%id_id_scores );

    if ( !$do_IMUS ) { return 1; }

    $LOG->info("IDENTIFYING A MAXIMUM UNRELATED SET");

    $LOG->debug("Checking for large networks...");
    breakup_large_networks( \%id_id_scores );

    $LOG->debug("done.");

    my $num_networks = keys %networks;

    $LOG->info("# of family networks: $num_networks");

    $LOG->info("Writing out unrelated set");
    my %PRIMUS_unrelated_set =
      write_out_independent_set( $relatedness_file_name, \%id_id_scores );

    $LOG->info("done.");

    $LOG->info("Testing alternative methods...");
    compare_alternative_methods( $relatedness_file_name,
        \%PRIMUS_unrelated_set, %networks, \%id_id_scores );

    $LOG->info("done.");

    $LOG->info("unrelated_file: $relatedness_file_name\_maximum_independent_set")
    $LOG->info("unrelated_set size: " . ( keys %PRIMUS_unrelated_set ));
    $LOG->info("unrelated_file: $relatedness_file_name\_maximum_independent_set");
    $LOG->info("unrelated_set size: " . ( keys %PRIMUS_unrelated_set ));

    return ( "$output_dir/$relatedness_file_name\_maximum_independent_set",
        ( keys %PRIMUS_unrelated_set ) );
}

#############################################################################################
#### SUBROUTINES ############################################################################
#############################################################################################
# @purpose Reset all global variables to initial state between test runs or invocations
# @param (none)
# @return (void)
# @side_effects Clears all package global variables including %id_id_scores, %networks, etc.

sub reset_values {
    $do_IMUS = 1;
    $do_PR   = 1;
    $verbose = 1;
    $lib_dir = "";

    %arg                     = ();
    $relatedness_file        = "";
    $relatedness_file_name   = "";
    $ped_file                = "none";
    $missingness_file        = "none";
    $output_dir              = "";
    $network_ctr             = 0;  # the next network # to add to %networks hash
    $LOWEST_MAX_NETWORK_SIZE = 60;
    $THRESHOLD               = .1; # default is .1
    $MIN_LIKELIHOOD          = .1;
    $EXCLUDE_VALUE           = 0;
    $OUTFILE_HEADER          = "";
    $IBD_file_ref            = "";
    $RELATEDNESS_COLUMN      = -1;
    $ID1_COLUMN              = -1;
    $ID2_COLUMN              = -1;
    $FID1_COLUMN             = -1;
    $FID2_COLUMN             = -1;
    %TRAIT_FID_COLUMNS       = -1;
    %TRAIT_ID_COLUMNS        = -1;
    %TRAIT_DATA_COLUMNS      = -1;
    $PRINT_ALTERNATE_RESULTS = 1;

    ### Here 'networks' could also be considered related groups or families;
    ## Hash of arrays: Key = network; Value = array (list) of IIDs in that network
    %networks = ();
    ## Regular hash; Key = IID1; Value = network that the individual belongs to
    %id_network = ();

    %id_id_all_info = ();
    @trait_refs     = ();
    @trait_order    = ();
    %trait_files =
      ();    ## Key = file name; Value = trait type (quantitative or binary)
    %child_parents = ();

}

# @purpose Print help/usage information to console
# @param (none)
# @return (void) - prints to STDOUT
# @side_effects Clears screen and displays command-line usage

sub help {
    system("clear");
    print "
Command-line options:
USAGE: $0\t  -input [IBD_file]  -output_dir [output_dir]  -threshold [num; default = 0.1]  -[high|low|mean|tails]/[b|q]trait [trait_file]

Synopsys:
-input [input_file_name]
-threshold [value]; (default = 0.1)
-output_dir [path_to_existing_directory]; (default = [input_file_name]_results/)
-trait size
-[high|low|mean|tails|'user_specified_value']/[q|p]trait [name]
-ped [name]
-missing_value [value]

Quick start help:

perl $0 -input ../example_data/example.genome

Should run the example data file and output the results in ../example_data/example.genome_results/


For more information on options and input/output files, please read the documentation, which is available at 
http://sourceforge.net/projects/primus-beta/

For questions, feature requests, and bug reports contact:
Jeff Staples - grapas2\@uw.edu
Piper Below - below\@uw.edu

\n\n";
}

# @param $name1 - First identifier
# @param $name2 - Second identifier
# @return (scalar) - Sorted key string "SMALLER;LARGER" using lexicographic order
# @notes Uses string comparison (lt) for consistent ordering; ensures single storage of pairs

sub sort_key {
    my ( $name1, $name2 ) = @_;
    return ( $name1 lt $name2 ) ? "$name1\;$name2" : "$name2\;$name1";
}

# @purpose Parse command-line arguments from @ARGV
# @param (none, uses @ARGV directly)
# @return (void)
# @side_effects Sets %arg hash with parsed options
# @deprecated - use set_values2 instead

# sub parseCommandLine {
#     my $trait_ctr = 1;
#     for ( my $i = 0 ; $i <= $#ARGV ; $i++ ) {
#         if ( $ARGV[$i] =~ /^-/ ) {
#             if ( $ARGV[$i] =~ /trait/ ) {
#                 if ( $ARGV[ $i + 1 ] ne "NA" && $ARGV[ $i + 1 ] !~ /^-/ ) {
#                     ## Check if user is trying to weight on a qtrait first; if so, default to weighting on size, then qtrait
#                     if ( $ARGV[$i] =~ /qtrait/ && @trait_order eq 0 ) {
#                         push( @trait_order, "size" );
#                     }
#                     push( @trait_order, $ARGV[ $i + 1 ] );
#                     $trait_files{ $ARGV[ $i + 1 ] } =
#                       $ARGV[$i];    ## Key is file and value is trait type
#                 }
#             }
#             $arg{ $ARGV[$i] } = $ARGV[ $i + 1 ];
#         }
#     }

#     ## Append the size trait to the end of traits if it is not already in it.
#     my @arr = %trait_files;
#     if ( !grep( /size/i, @trait_order ) ) {
#         push( @trait_order, "size" );
#     }
#     $trait_files{"size"} =
#       "-size";    ## Key is supposed to be file and value is trait type

#     if ( exists $arg{-help} ) {
#         help();
#         exit;
#     }
#     die("\n\nInput file required $useage") if ( !( $arg{-input} ) );

# }

# @purpose Parse and validate configuration options passed as array reference
# @param @_ - Array reference of command-line style arguments
# @return (void)
# @side_effects Sets numerous global variables: $THRESHOLD, $output_dir, %trait_files, etc.
# @calls GetOptionsFromArray, die on validation errors
# @uses Getopt::Long module

sub set_values2 {

    my ($config, @args) = @_;

    GetOptionsFromArray(
        \@args,

        # Diagnostic options
        'verbose=i' => sub { $config->{verbose} = $_[1]},
        'help|?'    => sub { help() },

        # Settings
        "rel_threshold=f"         => sub { $config->{threshold} = $_[1] },
        "int_likelihood_cutoff=f" => sub { $config->{min_likelihood} = $_[1] },
        "do_IMUS=i"               => sub { $config->{do_IMUS} = $_[1] },
        "do_PR=i"                 => sub { $config->{do_PR} = $_[1] },
        "missing_val=f"           => sub { $config->{exclude_value} = $_[1] },
        "output_dir=s"            => sub { $config->{output_dir} = $_[1] },
        "samples_file=s"          => sub { $config->{samples_file} = $_[1] },

        #"ersa_data=s" => \$ersa_data,
        "trait_order=s"     => sub { $config->{trait_order} = $_[1] },
        "lib=s"             => sub { $config->{lib_dir} = $_[1] },
        "traits=s"          => sub { $config->{trait_files} = $_[1] },
        "log_file_handle=s" => sub { $config->{log_file_handle} = $_[1] },
        # Need to adjust this section to wrok with the configuration object
        "ibd_estimates=s"   => sub {
            my %ibds = %{ @_[1] };
            $IBD_file_ref     = \%ibds;
            $relatedness_file = $ibds{'FILE'};
            if ( $relatedness_file =~ /\/([^\/]+)$/
              ) ## if there is a path to the file, ignore all the path, and select the name of the file
            {
                $relatedness_file_name = $1;
            }
            else {
                $relatedness_file_name = $relatedness_file;
            }
            $FID1_COLUMN        = $ibds{'FID1'} - 1;
            $ID1_COLUMN         = $ibds{'IID1'} - 1;
            $FID2_COLUMN        = $ibds{'FID2'} - 1;
            $ID2_COLUMN         = $ibds{'IID2'} - 1;
            $RELATEDNESS_COLUMN = $ibds{'PI_HAT'} - 1;

        },

    ) or die "Failed to parse options for IMUS\n";

    if ( !-d $config->{output_dir} ) {
        mkdir("$config->{output_dir}") or die "Can't make $config->{output_dir}; $!\n";
    }

    # Initialize trait_order and trait_files if not yet set
    if ( !defined $config->{trait_order} || ref($config->{trait_order}) ne 'ARRAY' ) {
        $config->{trait_order} = [];
    }
    if ( !defined $config->{trait_files} || ref($config->{trait_files}) ne 'HASH' ) {
        $config->{trait_files} = {};
    }

    ## If traits were empty, then add size
    if ( !exists $config->{trait_files}{'size'} ) {
        $config->{trait_files}{'size'} = 'size';
        push( @{ $config->{trait_order} }, 'size' );
    }
    
    # Initialize trait column tracking in config if not present
    if ( !defined $config->{trait_fid_columns} ) {
        $config->{trait_fid_columns} = {};
    }
    if ( !defined $config->{trait_id_columns} ) {
        $config->{trait_id_columns} = {};
    }
    if ( !defined $config->{trait_data_columns} ) {
        $config->{trait_data_columns} = {};
    }

    foreach my $file ( keys %{ $config->{trait_files} } ) {

        #print "Trait file: $file\n";
        ## Check that trait files exist
        if ( $file eq "NA" || $file =~ /^-/ ) {
            delete $config->{trait_files}{$file};
            next;
        }
        elsif ( $file =~ /size/i ) {
            next;
        }
        elsif ( !-e $file ) {
            die "Trait file $file does not exist\n";
        }

        my $trait_type = $config->{trait_files}{$file};

        open( IN, $file ) or die "Trait file $file failed to open; $!\n";
        my $header = <IN>;
        chomp($header);
        close(IN);
        $header =~ s/^\s+//;

        ## If there are only two columns, assume columns "ID  Trait_value"
        my @temp = split( /\s+/, $header );
        if ( @temp == 3 ) {
            $config->{trait_fid_columns}{$file}  = 0;
            $config->{trait_id_columns}{$file}   = 1;
            $config->{trait_data_columns}{$file} = 2;
        }
        elsif ( $file =~ /.ped$/i ) {
            $config->{trait_fid_columns}{$file}  = 0;
            $config->{trait_id_columns}{$file}   = 1;
            $config->{trait_data_columns}{$file} = 5;
        }
        else    ## get columns
        {
            $LOG->info("For $file:");
            $config->{trait_fid_columns}{$file} =
              get_correct_column( "Enter the column # containing the FIDs:  ",
                $header );
            $config->{trait_id_columns}{$file} =
              get_correct_column( "Enter the column # containing the IDs:  ",
                $header );
            $config->{trait_data_columns}{$file} = get_correct_column(
                "Enter the column # containing the $trait_type data:  ",
                $header );
        }
    }
}

# sub set_values {
#     $relatedness_file = $arg{-input};
#     if ( $relatedness_file =~ /\/([^\/]+)$/
#       ) ## if there is a path to the file, ignore all the path, and select the name of the file
#     {
#         $relatedness_file_name = $1;
#     }
#     else {
#         $relatedness_file_name = $relatedness_file;
#     }
#     if ( !-e $relatedness_file ) {
#         die "input_file $relatedness_file does not exist\n";
#     }

#     if ( exists $arg{-BKcutoff} ) {
#         print "ERROR! -BKcutoff invalid with versions 2.7 and above\n";
#         exit;

#         #$LOWEST_MAX_NETWORK_SIZE = $arg{-BKcutoff};
#     }

#     if ( exists $arg{-threshold} ) {
#         $THRESHOLD = $arg{-threshold};
#     }

#     if ( exists $arg{-missing_value} ) {
#         $EXCLUDE_VALUE = $arg{-missing_value};
#     }

#     if ( exists $arg{-output_dir} ) {
#         $output_dir = $arg{-output_dir};
#     }
#     else {
#         $output_dir = "$relatedness_file\_results";
#     }
#     if ( !-d $output_dir ) {
#         mkdir("$output_dir") or die "Can't make $output_dir; $!\n";
#     }

#     ## Determine the type of intput file or column containing relatedness
#     open( IN, $relatedness_file );
#     my $header = <IN>;
#     chomp($header);
#     $header =~ s/^\s+//;
#     close IN;
#     my @h_elements = split( /\s+/, $header );

#     for ( my $col = 0 ; $col < @h_elements ; $col++ ) {
#         if ( @h_elements[$col] =~ /PI_HAT/i || @h_elements[$col] =~ /Kinship/i )
#         {
#             $RELATEDNESS_COLUMN = $col;
#         }
#         if ( @h_elements[$col] =~ /IID1/i
#             || ( @h_elements[$col] =~ /ID1/i && @h_elements[$col] =~ /FID1/i ) )
#         {
#             $ID1_COLUMN = $col;
#         }
#         if ( @h_elements[$col] =~ /IID2/i
#             || ( @h_elements[$col] =~ /ID2/i && @h_elements[$col] =~ /FID2/i ) )
#         {
#             $ID2_COLUMN = $col;
#         }
#         if ( @h_elements[$col] =~ /FID1/i || @h_elements[$col] =~ /FID$/i ) {
#             $FID1_COLUMN = $col;
#         }
#         if ( @h_elements[$col] =~ /FID2/i || @h_elements[$col] =~ /FID$/i ) {
#             $FID2_COLUMN = $col;
#         }
#     }
#     if ( $RELATEDNESS_COLUMN eq -1 ) {
#         print "Unable to identify the relatedness column in input file\n";
#         $RELATEDNESS_COLUMN = get_correct_column(
#             "Enter the column # containing the relatedness data:  ", $header );
#     }
#     if ( $ID1_COLUMN eq -1 ) {
#         print "Unable to identify the ID1 column in input file\n";
#         $ID1_COLUMN =
#           get_correct_column( "Enter the column # containing ID1:  ", $header );
#     }
#     if ( $ID2_COLUMN eq -1 ) {
#         print "Unable to identify the ID2 column in input file\n";
#         $ID2_COLUMN =
#           get_correct_column( "Enter the column # containing ID2:  ", $header );
#     }
#     if ( $FID1_COLUMN eq -1 ) {
#         print "Unable to identify the FID1 column in input file\n";
#         $FID1_COLUMN =
#           get_correct_column( "Enter the column # containing FID1:  ",
#             $header );
#     }
#     if ( $FID2_COLUMN eq -1 ) {
#         print "Unable to identify the FID2 column in input file\n";
#         $FID2_COLUMN =
#           get_correct_column( "Enter the column # containing FID2:  ",
#             $header );
#     }

#     foreach my $file ( keys %trait_files ) {
#         print "Trait file: $file\n";
#         ## Check that trait files exist
#         if ( $file eq "NA" || $file =~ /^-/ ) {
#             delete $trait_files{$file};
#             next;
#         }
#         elsif ( $file =~ /size/i ) {
#             next;
#         }
#         elsif ( !-e $file ) {
#             die "Trait file $file does not exist\n";
#         }

#         my $trait_type = $trait_files{$file};

#         open( IN, $file ) or die "Trait file $file failed to open; $!\n";
#         my $header = <IN>;
#         chomp($header);
#         close(IN);
#         $header =~ s/^\s+//;

#         ## If there are only two columns, assume columns "ID  Trait_value"
#         my @temp = split( /\s+/, $header );
#         if ( @temp == 3 ) {
#             $TRAIT_FID_COLUMNS{$file}  = 0;
#             $TRAIT_ID_COLUMNS{$file}   = 1;
#             $TRAIT_DATA_COLUMNS{$file} = 2;
#         }
#         elsif ( $file =~ /.ped$/i ) {
#             $TRAIT_FID_COLUMNS{$file}  = 0;
#             $TRAIT_ID_COLUMNS{$file}   = 1;
#             $TRAIT_DATA_COLUMNS{$file} = 5;
#         }
#         else    ## get columns
#         {
#             print "For $file:\n";
#             $TRAIT_FID_COLUMNS{$file} =
#               get_correct_column( "Enter the column # containing the FIDs:  ",
#                 $header );
#             $TRAIT_ID_COLUMNS{$file} =
#               get_correct_column( "Enter the column # containing the IDs:  ",
#                 $header );
#             $TRAIT_DATA_COLUMNS{$file} = get_correct_column(
#                 "Enter the column # containing the $trait_type data:  ",
#                 $header );
#         }
#     }
# }

## Method used for setting the data column of interest
# @purpose Prompt user to select correct column from file header
# @param $question (scalar) - Question to display to user
# @param $header (scalar) - Header line from file with column names
# @return (scalar) - 0-based column index selected by user
# @side_effects Reads from STDIN, displays menu of options

sub get_correct_column {
    my $question = shift;
    my $header   = shift;
    chomp($header);
    my $column = -1;

    my @h_elements = split( /\s+/, $header );

    while ( $column eq -1 ) {
        for ( my $col = 1 ; $col <= @h_elements ; $col++ ) {
            $LOG->info("$col: @h_elements[$col-1]");
        }
        $LOG->info("$question");
        ## READ IN RESPONSE, and check if valid. set it as the column
        my $answer = <STDIN>;
        chomp $answer;

        if ( $answer >= 1 && $answer <= @h_elements ) {
            $column = $answer;
        }
        else {
            $LOG->warn("Invalid response\n\n");
        }
    }

    return $column -
      1;   # needs to be zero based, not one based like printed out to the user.
}

# @purpose Load sample IDs from file into %id_network and %networks
# @param $file (scalar) - Path to sample file
# @return (void)
# @side_effects Populates %id_network (IID => network_ctr) and %networks arrays
# @notes Supports .fam format (FID IID ...) or single column format

sub load_samples {
    my $file = shift;
    $LOG->info("Loading all samples from $file...");

    if ( !open( IN, $file ) ) {
        my $msg = "ERROR!!! Samples file $file cannot be read; $!\n";
        $LOG->warn($msg);
        die $msg;
    }
    while ( my $line = <IN> ) {

        # remove leading whitespace and newline characters
        $line =~ s/^\s+//;
        chomp($line);

        # if line is empty or we think it is a header, skip it
        next if $line eq "" || $line =~ /^FID/;    # Skip empty and header

        my @temp = split( /\s+/, $line );
        my $iid;
        if ( @temp >= 2 ) {
            $iid = $temp[1];    # Standard .fam format: FID IID ...
        }
        else {
            $iid = $temp[0];    # Single column format
        }

        if ( !exists $id_network{$iid} ) {
            $id_network{$iid} = $network_ctr;
            push @{ $networks{$network_ctr} }, "$iid";
            $network_ctr++;
        }
    }
    close(IN);
}

# @purpose Load relatedness/IBD data from file into %id_id_scores and %id_id_all_info
# @param $file (scalar) - Path to relatedness file
# @param $id_id_scores_ref (hashref) - Reference to %id_id_scores hash to populate
# @return (void)
# @side_effects Populates hash referenced by $id_id_scores_ref (key => PI_HAT value) and %id_id_all_info
# @side_effects Populates %networks and %id_network with IID pairs
# @uses global $THRESHOLD, $RELATEDNESS_COLUMN, $ID1_COLUMN, $ID2_COLUMN, $FID1_COLUMN, $FID2_COLUMN
# @notes Only stores pairs with score > THRESHOLD in hash

sub load_data {
    my ( $file, $id_id_scores_ref ) = (@_);
    open( IN, $file )
      or die "ERROR!!! Relatedness input file $file cannot be read in; $!\n";
    $OUTFILE_HEADER = <IN>;    ## skip header
    while ( my $line = <IN> ) {
        $line =~ s/^\s+//;
        chomp($line);
        if ( $line =~ /^FID/ ) { $OUTFILE_HEADER = $line; next; } ## skip header
        my @temp = split( /\s+/, $line );
        my $FID1 = @temp[$FID1_COLUMN];
        my $FID2 = @temp[$FID2_COLUMN];
        my $IID1 = @temp[$ID1_COLUMN];
        my $IID2 = @temp[$ID2_COLUMN];

        #if($IID1 eq "."){$IID1 = $FID1;}
        #if($IID2 eq "."){$IID2 = $FID2;}

        #my $name1 = "$FID1**$IID1";
        #my $name2 = "$FID2**$IID2";
        my $name1 = "$IID1";
        my $name2 = "$IID2";

        # Lets sort the names so that we dont have to store both
        # copies of the key
        my $key = sort_key( $name1, $name2 );

        my $PI_HAT = @temp[$RELATEDNESS_COLUMN];
        if ( $PI_HAT > $THRESHOLD ) {
            $id_id_scores_ref->{$key} = $PI_HAT;
        }

        $id_id_all_info{$key} = $line;

        if ( !exists $id_network{$name1} ) {
            $id_network{$name1} = $network_ctr;
            push @{ $networks{$network_ctr} }, "$name1";
            $network_ctr++;
        }
        if ( !exists $id_network{$name2} ) {
            $id_network{$name2} = $network_ctr;
            push @{ $networks{$network_ctr} }, "$name2";
            $network_ctr++;
        }
    }
    close(IN);
}

# @purpose Load trait data from files into @trait_refs array of hashes
# @param (none)
# @return (void)
# @side_effects Populates @trait_refs with trait hash references for each trait in @trait_order
# @uses global @trait_order, %trait_files, @trait_refs, %trait_refs
# @notes Special handling for 'size' trait (always added); supports binary and quantitative traits

sub load_trait_data {
    foreach my $file (@trait_order) {

        #print "Loading trait data: $file\n";
        my $trait_type = $trait_files{$file};
        my %trait_hash;

        ## input the trait size values of 1 for each ID from the input file
        if ( $file =~ /size/i ) {
            foreach my $ID ( keys %id_network ) {
                $trait_hash{$ID} = 1;
            }
            push( @trait_refs, \%trait_hash )
              ; ## The hash references are loaded onto this array in the order of selection
            next;
        }
        open( IN, $file ) or die "Trait file $file failed to open; $!\n";
        while ( my $line = <IN> ) {
            $line =~ s/^\s+//;
            chomp($line);
            my @temp = split( /\s+/, $line );

            my $FID_COLUMN = $TRAIT_FID_COLUMNS{$file};
            my $IID_COLUMN = $TRAIT_ID_COLUMNS{$file};
            my $FID        = @temp[$FID_COLUMN];
            my $IID        = @temp[$IID_COLUMN];

            #my $ID = "$FID**$IID";
            my $ID          = "$IID";
            my $DATA_COLUMN = $TRAIT_DATA_COLUMNS{$file};
            my $trait_val   = @temp[$DATA_COLUMN];

            if ( $trait_val == $EXCLUDE_VALUE ) {
                $trait_val = "NA";
            }
            elsif ( $trait_type =~ /btrait/i ) {
                if ( $trait_val == 2 ) {
                    $trait_val = 1;
                }
                elsif ( $trait_val == 1 ) {
                    $trait_val = 0;
                }
                else {
                    die
"Binary trait file $file contains a binary value other than 1 or 2 in the data column $DATA_COLUMN:\n\t$line\n";
                }
            }
            $trait_hash{$ID} = $trait_val;
        }
        close(IN);

        if ( $trait_type =~ /(mean|tail|\d+)/ ) {
            if ( $trait_type =~ /btrait/i ) {
                die
"INVALID TRAIT OPTION: $trait_type for file $file.\nOnly options 'high' and 'low' are valid with binary traits\n\n";
            }

            my $selection_val = fold_trait_data( $trait_type, \%trait_hash );
            $LOG->info("Selection value for $file: " . $selection_val);
        }
        push( @trait_refs, \%trait_hash )
          ; ## The hash references are loaded onto this array in the order of selection
    }
}

## Method used for doing weighting of tails, mean, or user specified value. It will adjust all the trait values to be the distance from the mean or the user specified value. This is NOT the most effective way to do the tail weighting of skewed traits where the mean is not actually halfwaye between the highest and lowest value. A better way would be to find the highest and lowest value, and for each value change it to the min(max-val,val-min). This will set it to the distance from the max or min whichever it is closer to.

# @purpose Adjust trait data by distance from fold value (mean, tail weighting, or user value)
# @param $trait_type (scalar) - Type of folding: 'mean', 'tail', or numeric value
# @param $hash_ref (hashref) - Hash of ID => trait_value pairs
# @return (scalar) - The fold value used (mean or user-specified)
# @side_effects Modifies values in $hash_ref to be absolute distance from fold_value
# @notes Used for weighting selection criteria; improves quantitative trait handling

sub fold_trait_data {
    my $trait_type = shift;
    my $hash_ref   = shift;

    my $fold_value = 0;
    if ( $trait_type =~ /(mean|tail)/ ) {
        ## set fold value to the mean of all trait values
        my $ctr = 0;
        my $sum = 0;
        foreach ( keys %$hash_ref ) {
            $ctr++;
            $sum = $sum + $$hash_ref{$_};
        }
        $fold_value = $sum / $ctr;
    }
    if ( $trait_type =~ /\d+/ ) {
        ## up out the fold_value from the trait type
        $trait_type =~ /-?(-?\.?\d+)/;
        $fold_value = $1;
    }

    ## subtract $fold_value from each trait value and then take the absolute value, resulting in possitive distance from the fold value
    foreach ( keys %$hash_ref ) {
        my $old_val = $$hash_ref{$_};
        my $new_val = abs( $old_val - $fold_value );
        $$hash_ref{$_} = $new_val;
    }
    return $fold_value;
}

## Network processing subroutines
# @purpose Merge networks containing related individuals using predict_relationships_2D
# @param $id_id_scores (hashref) - Reference to %id_id_scores hash (pair_key => PI_HAT)
# @return (void)
# @side_effects Merges %networks entries where individuals are related; updates %id_network
# @uses global $THRESHOLD, %id_network, %networks, $do_PR, $IBD_file_ref, $MIN_LIKELIHOOD, $lib_dir, $output_dir
# @calls PRIMUS::predict_relationships_2D::get_relationship_likelihood_vectors
# @notes Core network construction algorithm; combines individuals into families

sub colapse_networks {
    my ($id_id_scores) = @_;

    # print all unique entries
    #foreach my $network (sort {$a <=> $b} keys %networks)
    #{
    #	my @temp = @{ $networks{$network} };
    #}
    my $relationships_ref;
    if ($do_PR) {
        eval {
            ($relationships_ref) =
              PRIMUS::predict_relationships_2D::get_relationship_likelihood_vectors(
                $IBD_file_ref, $MIN_LIKELIHOOD, 0, $lib_dir, $output_dir );
        };
        if ($@) {
            $LOG->warn("ERROR running predict_relationships_2d.pm: $@");
            die;
        }
    }

    foreach my $id_pair ( keys %$id_id_scores ) {
        my ( $id1, $id2 ) = split( /\;/, $id_pair );
        my $score       = $id_id_scores->{$id_pair};
        my $id1_network = $id_network{$id1};
        my $id2_network = $id_network{$id2};

        #if unrelated, do nothing
        if ( !$do_PR ) {
            if ( $score <= $THRESHOLD ) {
                next;
            }
        }
        else {
            #$id1 =~ s/\*\*/__/;
            #$id2 =~ s/\*\*/__/;
            foreach ( keys %$relationships_ref ) {

                #print "key: $_\n";
            }
            if ( !exists $$relationships_ref{$id1}{$id2} ) {
                my $temp = $id1;
                $id1 = $id2;
                $id2 = $temp;
            }

            #print "id1: $id1; id2: $id2\n";
            my @vector = @{ $$relationships_ref{$id1}{$id2} };
            my @possibilities =
              PRIMUS::predict_relationships_2D::predict_relationship(@vector);
            ## Check that there is at least one.
            if ( @possibilities < 1 ) {

            }
            ## If that one is UN then treat as unrelated
            if ( @possibilities == 1 ) {
                if ( $possibilities[0] eq "UN" ) {
                    next;
                }
            }
        }

        #if already in the same network, do nothing
        if ( $id1_network eq $id2_network ) {
            next;
        }

        # combine $id1_network and $id2_network
        my @network1  = @{ $networks{$id1_network} };
        my @network2  = @{ $networks{$id2_network} };
        my @new_array = ( @network1, @network2 );
        $networks{$id1_network} = [@new_array];
        delete $networks{$id2_network};

        foreach my $id (@network2) {
            my $old_network = $id_network{$id};
            if ( $old_network ne $id2_network ) {
                $LOG->warn(
"ERROR!!! $id claims to be in $old_network; actually in $id2_network");
                exit;
            }
            $id_network{$id} = $id1_network;
        }
    }
}

# @purpose Write formatted network data to output files
# @param $file (scalar) - Base filename for output files
# @return (void)
# @side_effects Creates output files: ${file}_networks, ${file}_unrelated_samples.txt, ${file}_network*.genome, ${file}_network*.dot
# @uses global %networks, %id_id_all_info, %id_id_scores, $THRESHOLD, $output_dir
# @calls write_out_dot_file for networks > 4 nodes

sub write_out_networks {
    my ( $file, $id_id_scores ) = @_;
    my $num_networks = keys %networks;

    open( NETWORKS_OUT, ">$output_dir/$file\_networks" )
      or die "Can't write to $output_dir/$file\_networks; $!\n";
    open( UNREL_OUT, ">$output_dir/$file\_unrelated_samples.txt" )
      or die "Can't write to $output_dir/$file\_unrelateds; $!\n";
    print NETWORKS_OUT "Network\t$OUTFILE_HEADER";
    print UNREL_OUT "FID\tIID\n";
    my $network_ctr = 1;
    foreach my $network ( sort { $a <=> $b } keys %networks ) {

        my @temp = @{ $networks{$network} };
        if ( @temp < 2 ) {
            my $name = @temp[0];
            $name =~ s/\*\*/\t/;
            print UNREL_OUT "$name\n";
            next;
        }
        if ( @temp > 4 ) {
            write_out_dot_file( $file, $network, $network_ctr, $id_id_scores );
        }
        ## This will write out the .genome file for each network (if a .genome file was read in)
        open( GENOMES_OUT, ">$output_dir/$file\_network$network_ctr.genome" )
          or die "Can't write to $output_dir/$file\_networks; $!\n";
        print GENOMES_OUT "$OUTFILE_HEADER";

        for my $i ( 0 .. @temp - 1 ) {
            for ( my $j = 0 ; $j < $i ; $j++ ) {
                my $id1  = $temp[$i];
                my $id2  = $temp[$j];
                my $key  = sort_key( $id1, $id2 );
                my $info = $id_id_all_info{$key};
                if ( $info ne "" ) {
                    print GENOMES_OUT "$info\n";
                }
                my $score = $id_id_scores->{$key};
                if ( $score > $THRESHOLD ) {
                    print NETWORKS_OUT "$network_ctr\t$info\n";
                }
            }
        }
        $network_ctr++;
    }
    close(NETWORKS_OUT);
    close(UNREL_OUT);

}

# @purpose Calculate connectivity of a network (ratio of edges to possible edges)
# @param $network_ref (hashref) - Network nodes (nodes => value)
# @param $id_id_scores (hashref) - Reference to %id_id_scores hash (pair_key => PI_HAT)
# @return (scalar) - Connectivity value 0-1 (0=isolated, 1=fully connected)
# @side_effects None (read-only operation on parameters)
# @uses global $THRESHOLD
# @notes Used to determine algorithm parameters for large networks

sub get_connectedness {
    my ( $network_ref, $id_id_scores ) = (@_);
    my $num_connections;
    my %visited;
    my @ids             = keys %$network_ref;
    my $num_ids         = keys %$network_ref;
    my $max_connections = 0;
    for ( my $i = 0 ; $i < $num_ids ; $i++ ) {
        for ( my $j = $i + 1 ; $j < $num_ids ; $j++ ) {
            my $key        = @ids[$i];
            my $key2       = @ids[$j];
            my $sorted_key = sort_key( $key, $key2 );
            $max_connections++;
            my $PI_HAT = $id_id_scores->{$sorted_key};
            if ( $PI_HAT > $THRESHOLD ) {
                $num_connections++;
            }
        }
    }
    if ( $max_connections eq 0 ) {
        return 1;
    }
    my $connectedness = $num_connections / $max_connections;
    return $connectedness;
}

sub breakup_large_networks {

    my ($id_id_scores) = (@_);

    foreach my $network ( sort { $a <=> $b } keys %networks ) {
        my @temp = @{ $networks{$network} };
        my %P    = map { $_ => 1 } @temp;

        if ( keys %P > $LOWEST_MAX_NETWORK_SIZE ) {
            my $connectedness = get_connectedness( \%P, $id_id_scores );
            my $size          = keys %P;
            my $MAX_NETWORK_SIZE =
              get_max_network_size( $connectedness, $size );
            if ( $size > $MAX_NETWORK_SIZE ) {
                $LOG->warn("WARNING!!! " . @temp
                  . " NODES WITH CONNECTIVITY OF $connectedness WILL TAKE TOO LONG TO RUN; USING NEXT BEST SOLUTION.");
                breakup_large_network( \%P, $MAX_NETWORK_SIZE, $id_id_scores );
                @{ $networks{$network} } = keys %P;
            }
        }
    }

}

sub breakup_large_network {

    my ( $P_ref, $MAX_NETWORK_SIZE, $id_id_scores ) = (@_);
    my %degrees;
    my %neighbors;

    load_degrees_and_neighbors( $P_ref, \%degrees, \%neighbors, $id_id_scores );

    while ( keys %$P_ref > $MAX_NETWORK_SIZE ) {
        my ( $node_to_remove, $degree ) = get_highest_degree_node( \%degrees );
        if ( $degree > 0 ) {
            delete $$P_ref{$node_to_remove};
            delete $degrees{$node_to_remove};
            reduce_neighbors( $node_to_remove, \%neighbors, \%degrees );
            delete $neighbors{$node_to_remove};

            # Get connected components of new graph
            my @component_refs =
              get_connected_components( $P_ref, \%neighbors );
            %{$P_ref} = %{ $component_refs[0] };
            shift(@component_refs);

# Add smaller connected components to end of $networks{$network_ctr}
#foreach my $hash_ref (sort {keys %{ $component_refs[$a] } <=> keys %{ $component_refs[$b]} } @component_refs)
            foreach my $hash_ref (@component_refs) {
                my @temp = keys %$hash_ref;
                foreach (@temp) {
                    delete $degrees{$_};
                    delete $neighbors{$_};
                }

                ## if still too big, call this routine recursively
                if ( keys %$hash_ref > $LOWEST_MAX_NETWORK_SIZE ) {
                    my $connectedness =
                      get_connectedness( $hash_ref, $id_id_scores );
                    my $LOCAL_MAX_NETWORK_SIZE =
                      get_max_network_size($connectedness);
                    $LOG->info("Connectedness: $connectedness");
                    $LOG->info("LOCAL_MAX_NETWORK_SIZE: $LOCAL_MAX_NETWORK_SIZE");
                    if ( keys %$hash_ref > $LOCAL_MAX_NETWORK_SIZE ) {
                        ## DECEND RECURSIVELY
                        breakup_large_network( $hash_ref,
                            $LOCAL_MAX_NETWORK_SIZE, $id_id_scores );
                    }
                }
                $network_ctr++;
                @{ $networks{$network_ctr} } = keys %$hash_ref;
            }

# proceed with larger network %P until it is smaller than its $MAX_NETWORK_SIZE for its connectivity
        }
        else {
            # It should never get here
            die "ERROR!!! PRUNING NODE WITHOUT RELATIVES!!!\n";
        }
        my $connectedness = get_connectedness( $P_ref, $id_id_scores );
        $MAX_NETWORK_SIZE = get_max_network_size($connectedness);

    }
}

# @purpose Find all connected components in a graph
# @param $network_ref (hashref) - Nodes in graph (nodes => value)
# @param $neighbors_ref (hashref) - Adjacency info (node => hashref of neighbors)
# @return (array) - List of hashrefs representing connected component node sets
# @side_effects None (read-only operation)
# @algorithm BFS from each unvisited node; O(V+E) complexity

sub get_connected_components {
    my $network_ref   = shift;
    my $neighbors_ref = shift;
    my @component_refs;
    my %remaining_ids = %$network_ref;

    ## This loop will cycle through all connected components;
    while ( keys %remaining_ids > 0 ) {
        ## For this node, identify the connected component, and remove all nodes in the connected component from %remaining_ids
        my @temp = keys %remaining_ids;
        my $node = @temp[0];
        my %connected_component;
        my %waiting_to_visit;
        $waiting_to_visit{$node} = 1;

        while ( keys %waiting_to_visit > 0 ) {
            foreach my $node ( keys %waiting_to_visit ) {
                $connected_component{$node} = 1;
                delete $waiting_to_visit{$node};
                delete $remaining_ids{$node};
                my @neighbors = keys %{ $neighbors_ref->{$node} };
                foreach (@neighbors) {
                    if ( !exists $connected_component{$_} ) {
                        $waiting_to_visit{$_} = 1;
                    }
                }
            }
        }
        push( @component_refs, \%connected_component );
    }
    return @component_refs;
}

sub write_out_independent_set {
    my ( $file, $id_id_scores ) = (@_);
    my %unrelated_set;
    open( UNIQUE_OUT, ">$output_dir/$file\_maximum_independent_set_PRIMUS" );
    print UNIQUE_OUT "FID\tIID\n";
    ## get most unrelateds in each network with Bron-Kerbosch Algorithm
    foreach my $network ( sort { $a <=> $b } keys %networks ) {
        my @temp = @{ $networks{$network} };

        my %P = map { $_ => 1 } @temp;
        my %R;
        my %X;
        my @maximal_cliques;
        my %degrees;
        my %neighbors;

        my $nodes_visited = 0;
        if ( @temp > 45 ) {
            $LOG->info("Running BronKerbosh for network $network (size = " . @temp . ")");
        }
        BronKerbosh( \@maximal_cliques, \%R, \%P, \%X, \$nodes_visited,
            $id_id_scores );

        my $maximum_clique = get_maximum_clique(@maximal_cliques, \@trait_refs);
        my @maximum_ids    = keys %{ $maximal_cliques[$maximum_clique] };
        write_out_maximum_clique_ids(@maximum_ids);
        foreach (@maximum_ids) { $unrelated_set{$_} = 1; }
    }
    close(UNIQUE_OUT);
    return %unrelated_set;
}

# @purpose Calculate degree (neighbor count) for each node and store neighbor sets
# @param $network_ref (hashref) - Network nodes (nodes => value)
# @param $degree_ref (hashref) - Output: nodes => degree count
# @param $neighbors_ref (hashref) - Output: nodes => neighbor hashref
# @param $id_id_scores (hashref) - Reference to %id_id_scores hash (pair_key => PI_HAT)
# @return (void)
# @side_effects Modifies $degree_ref and $neighbors_ref hashrefs
# @uses get_actual_neighbors to find related nodes (original graph)
# @notes Used for network analysis and degree-based pruning

sub load_degrees_and_neighbors {

    my ( $network_ref, $degree_ref, $neighbors_ref, $id_id_scores ) = (@_);

    foreach my $node ( keys %$network_ref ) {
        my %neighbors =
          get_actual_neighbors( $node, $network_ref, $id_id_scores );

        $$neighbors_ref{$node} =
          \%neighbors;    #get_actual_neighbors($node,$network_ref);
        my @temp   = keys %{ $neighbors_ref->{$node} };    # %neighbors;
        my $degree = @temp;
        $$degree_ref{$node} = $degree;
    }
}

# @purpose Decrement degree count for neighbors when node is removed
# @param $node (scalar) - Node being removed
# @param $neighbors_ref (hashref) - Adjacency info (node => neighbors)
# @param $degree_ref (hashref) - Degree counts (node => count)
# @return (void)
# @side_effects Decrements $degree_ref for all neighbors of $node; removes edges
# @notes Used during network pruning to maintain accurate degree counts

sub reduce_neighbors {
    my $node          = shift;
    my $neighbors_ref = shift;
    my $degree_ref    = shift;
    foreach my $neighbor ( keys %{ $$neighbors_ref{$node} } ) {
        $$degree_ref{$neighbor}--;
        delete $$neighbors_ref{$neighbor}{$node};
    }
}

# @purpose Select node with highest degree (most related to others)
# @param $degrees_ref (hashref) - Node => degree_count
# @return (list) - ($max_node, $max_degree)
# @side_effects None (read-only)
# @algorithm Breaks ties using weighted_comparison of trait values
# @uses get highest-degree node; breaks ties by trait values to prefer nodes to keep

sub get_highest_degree_node {
    my $degrees_ref = shift;
    my $max_degree  = -1;
    my $max_node;
    my @max_trait_values = qw(NA);

    foreach my $node ( keys %$degrees_ref ) {
        my $degree = $$degrees_ref{$node};

        if ( $degree > $max_degree ) {
            $max_degree = $degree;
            $max_node   = $node;
        }
        elsif ( $degree == $max_degree ) {
            my @temp_values;
            for ( my $i = 0 ; $i < @trait_order ; $i++ ) {
                @temp_values[$i] = $trait_refs[$i]{$node};
            }

            my $max = weighted_comparison( \@max_trait_values, \@temp_values )
              ;    ## return 1 for max, and 2 for temp
            ## We want the lower of the two to remove; therefore if max has higher trait values, switch max to temp
            if ( $max == 1 || @max_trait_values[0] eq "NA" ) {
                $max_degree       = $degree;
                $max_node         = $node;
                @max_trait_values = @temp_values;
            }
        }
    }
    return ( $max_node, $max_degree );
}

# @purpose Select the best clique/result based on trait scores; supports both selecting from multiple cliques and comparing method results
# @param @maximal_cliques - Variable number of clique hashrefs (array of hashes) OR independent sets (hashes) when comparing methods
# @param $trait_refs (array reference) - Reference to @trait_refs array, where each element is a hashref mapping ID => trait_value
# @return (scalar) - Index of best clique/method (0=PLINK, 1=KING, 2=PRIMUS for comparisons; array index for clique selection)
# @side_effects None (read-only)
# @algorithm Calculates trait score sum for each clique/set; uses weighted_comparison to select best
# @uses $trait_refs_ref, @trait_order, %trait_files to evaluate quality
# @notes Handles both: (1) selecting best from multiple cliques within a network, (2) comparing independent sets from different methods

sub get_maximum_clique {
    my ($trait_refs) = pop;
    my (@maximal_cliques) = @_;
    my @max_values      = qw(NA);
    my $maximum_clique;

    for ( my $i = 0 ; $i < @maximal_cliques ; $i++ ) {
        my @temp_values;
        my $ctr = 0;
        foreach my $id ( keys %{ $maximal_cliques[$i] } ) {
            $ctr++;
            for ( my $i = 0 ; $i < @trait_order ; $i++ ) {
                my $trait_type = $trait_files{ $trait_order[$i] };
                my $trait_val  = $trait_refs->[$i]{$id};
                if ( exists $trait_refs->[$i]{$id} ) {
                    @temp_values[$i] += $trait_val;
                }
            }
        }

        ## Average  non-binary trait values
        for ( my $i = 0 ; $i < @trait_order ; $i++ ) {
            my $trait_type = $trait_files{ $trait_order[$i] };
            if ( $trait_type =~ /qtrait/i && $ctr ne 0 ) {
                @temp_values[$i] = @temp_values[$i] / $ctr;
            }
        }

        ## Compare
        my $max = weighted_comparison( \@max_values, \@temp_values );    ## returns 1 for max  and 2 for temp

        ## Update max if necessary
        if ( $max eq 2 ) {
            $maximum_clique = $i;
            @max_values     = @temp_values;
        }
    }
    return $maximum_clique;
}

# @purpose Compare two sets of trait values and return index of better one
# @param $a_ref (arrayref) - First set of trait values
# @param $b_ref (arrayref) - Second set of trait values
# @return (scalar) - 1 if $a better, 2 if $b better
# @side_effects None (read-only)
# @algorithm Lexicographically compares trait values in @trait_order; applies high/low/size preferences
# @uses @trait_order, %trait_files for trait type information
# @globals %trait_files, @trait_order

sub weighted_comparison {
    my $a_ref = shift;
    my $b_ref = shift;

    if ( $$a_ref[0] eq "NA" ) {
        return 2;
    }

    ## This will loop through the traits in order and returns the proper number of array once a difference is found.
    ##If no differences, the first array # is returned.
    for ( my $i = 0 ; $i < @$a_ref ; $i++ ) {
        my $trait_type = $trait_files{ $trait_order[$i] };
        if (   $trait_type =~ /high/i
            || $trait_type =~ /ped/i
            || $trait_type =~ /size/i
            || $trait_type =~ /tail/i )
        {
            if ( $$a_ref[$i] > $$b_ref[$i] ) {
                return 1;
            }
            if ( $$a_ref[$i] < $$b_ref[$i] ) {
                return 2;
            }
        }
        elsif ($trait_type =~ /low/i
            || $trait_type =~ /mean/i
            || $trait_type =~ /\d+/i )
        {
            if ( $$a_ref[$i] < $$b_ref[$i] ) {
                return 1;
            }
            if ( $$a_ref[$i] > $$b_ref[$i] ) {
                return 2;
            }
        }
    }
    return 1;
}

# @purpose Output individual IDs from clique to file in PLINK format
# @param @clique - Array of individual IDs
# @return (void)
# @side_effects Prints to UNIQUE_OUT filehandle
# @notes Format: FID IID (splits ** delimited IDs)

sub write_out_maximum_clique_ids {
    my @clique = @_;

    foreach my $ID (@clique) {
        my ( $FID, $IID ) = split( /\*\*/, $ID );
        print UNIQUE_OUT "$FID\t$IID\n";
    }

}

### BronKerbosh Algorithm subroutines
# @param $maximal_cliques_ref (arrayref) - Will be populated with maximal cliques found
# @param $R_ref (hashref) - Current clique being built (nodes => 1)
# @param $P_ref (hashref) - Candidate set of nodes that could extend R (nodes => 1)
# @param $X_ref (hashref) - Already processed nodes to avoid duplicate cliques (nodes => 1)
# @param $num_visited_ref (scalarref) - Reference to counter tracking recursive calls
# @return (scalar) - Number of nodes visited in this recursion level
# @uses %id_id_scores (global - REFACTOR: should be parameter)
# @uses $THRESHOLD (global - for determining neighbors)
# @algorithm Bron-Kerbosch with pivot optimization; finds maximal cliques in complement graph

sub BronKerbosh {

    # my $maximal_cliques_ref = shift;
    # my $R_ref               = shift;
    # my $P_ref               = shift;
    # my $X_ref               = shift;
    # my $num_visited_ref     = shift;

    my ( $maximal_cliques_ref, $R_ref, $P_ref, $X_ref, $num_visited_ref,
        $id_id_scores )
      = (@_);

    $$num_visited_ref++;
    my $nodes_visited = 1;

    if ( keys %$P_ref == 0 and %$X_ref == 0 ) {
        push @{$maximal_cliques_ref}, {%$R_ref};
        return;
    }

    ## Select pivot node from P union X that maximizes the cardinality of P intersection N(u);
    my ( $u, %u_neighbors ) = select_pivot( $P_ref, $X_ref, $id_id_scores );

    foreach my $v ( keys %$P_ref ) {
        ## Use the pivot node: continue only if $v is not a neighbor of the pivot $u
        if ( exists $u_neighbors{$v} ) { next; }

        my %temp_R = %$R_ref;
        $temp_R{$v} = $$P_ref{$v};    # Load temp_R = R union v

        my %temp_P;
        get_inverse_neighbors( $v, $P_ref, \%temp_P, $id_id_scores )
          ;                           # Load temp_P = P intersection N(v)
        my @temp_arr = keys %temp_P;

        my %temp_X;
        get_inverse_neighbors( $v, $X_ref, \%temp_X, $id_id_scores )
          ;                           # Load temp_X = X intersect N(v)

        BronKerbosh( $maximal_cliques_ref, \%temp_R, \%temp_P, \%temp_X,
            $num_visited_ref, $id_id_scores );
        $$X_ref{$v} = $$P_ref{$v};
        delete( $$P_ref{$v} );
    }
    return $nodes_visited;
}

# @param $P_ref (hashref) - Candidate set of nodes (nodes => 1)
# @param $X_ref (hashref) - Already processed nodes (nodes => 1)
# @param $id_id_scores (hashref) - Reference to %id_id_scores hash (pair_key => PI_HAT)
# @return (list) - ($pivot_node, %pivot_neighbors)
# @return_detail $pivot_node - Node from P with most neighbors in P
# @return_detail %pivot_neighbors - Inverse neighbors of pivot (score <= THRESHOLD)
# @uses get_inverse_neighbors()
# @notes Pivot selection maximizes |P ∩ N(u)| to reduce search space

## Select pivot node from P union X that maximizes the cardinality of P intersection N(u);
sub select_pivot {

    my ( $P_ref, $X_ref, $id_id_scores ) = (@_);

    my $u;
    my $max_size = -1;
    my %u_neighbors;

    foreach my $key ( keys %$P_ref ) {
        my %neighbors;
        get_inverse_neighbors( $key, $P_ref, \%neighbors, $id_id_scores );
        if ( $max_size < keys %neighbors ) {
            $max_size    = keys %neighbors;
            $u           = $key;
            %u_neighbors = %neighbors;
        }
    }
    return ( $u, %u_neighbors );
}

# @purpose Find all related neighbors of a node (where score > THRESHOLD)
# @param $v (scalar) - Node ID to find related neighbors for
# @param $hash_ref (hashref) - Candidate nodes to check against (nodes => value)
# @param $id_id_scores (hashref) - Reference to %id_id_scores hash (pair_key => PI_HAT)
# @return (hash) - Hash of related neighbors (nodes => value)
# @side_effects None (read-only operation on parameters)
# @uses global $THRESHOLD
# @notes Related = kinship coefficient > THRESHOLD (typically > 0.1); opposite of get_inverse_neighbors

sub get_actual_neighbors {
    my ( $v, $hash_ref, $id_id_scores ) = (@_);
    my %neighbors;

    foreach my $n_v ( keys %$hash_ref ) {
        if ( $n_v eq $v ) {
            next;
        }
        my $key   = sort_key( $v, $n_v );
        my $score = $id_id_scores->{$key};
        if ( $score > $THRESHOLD ) {
            $neighbors{$n_v} = $$hash_ref{$n_v};
        }
    }

    return %neighbors;
}

# @purpose Find all unrelated neighbors of a node (where score <= THRESHOLD)
# @param $v (scalar) - Node ID to find unrelated neighbors for
# @param $hash_ref (hashref) - Candidate nodes to check against (nodes => value)
# @param $neighbors (scalar or hashref) - Output: pass undef for auto-vivification or hashref to populate
# @param $id_id_scores (hashref) - Reference to %id_id_scores hash (pair_key => PI_HAT)
# @return (hashref) - Reference to neighbors hash (auto-vivified if needed), modified in place
# @side_effects Modifies $neighbors hash (or creates new one via auto-vivification) with unrelated nodes
# @uses global $THRESHOLD
# @notes Unrelated = kinship coefficient <= THRESHOLD (typically <= 0.1); KEY function for Bron-Kerbosch algorithm

sub get_inverse_neighbors {

    my ( $v, $hash_ref, $neighbors, $id_id_scores ) = (@_);

    foreach my $n_v ( keys %$hash_ref ) {

        if ( $n_v eq $v ) {
            next;
        }
        my $key   = sort_key( $v, $n_v );
        my $score = $id_id_scores->{$key};
        if ( $score <= $THRESHOLD ) {
            $$neighbors{$n_v} = $$hash_ref{$n_v};
        }
    }
    return $neighbors;
}

### Write the .dot file that can be read into a graph visualization program like Graphviz
# @purpose Write network graph in Graphviz DOT format for visualization
# @param $file (scalar) - Base filename
# @param $network (scalar) - Network ID to output
# @param $network_ctr (scalar) - Counter for output file naming
# @return (void)
# @side_effects Creates ${file}_network${network_ctr}.dot file
# @uses %networks, %id_id_all_info, %id_id_scores for edge information
# @notes Black edges = related pairs (score > THRESHOLD); nodes for isolated individuals

sub write_out_dot_file {
	my ($file, $network, $network_ctr, $id_id_scores) = (@_);

    open( GRAPH_OUT, ">$output_dir/$file\_network$network_ctr.dot" );
    print GRAPH_OUT "graph network$network_ctr {\n";
    print GRAPH_OUT "\tnode [shape=circle];\n\n";

    my @temp = @{ $networks{$network} };

    ## Write out all connections
    for my $i ( 0 .. @temp - 1 ) {
        for ( my $j = 0 ; $j < $i ; $j++ ) {
            my $id1   = $temp[$i];
            my $id2   = $temp[$j];
            my $key   = sort_key( $id1, $id2 );
            my $info  = $id_id_all_info{$key};
            my $score = $id_id_scores->{$key};

            ## Change the ** delimiter between FID and IID to and _
            #$id1 =~ s/\*\*/_/;
            #$id2 =~ s/\*\*/_/;

            if ( $score > $THRESHOLD ) {
                my $color = "black";
                print GRAPH_OUT "\t\"$id1\" -- \"$id2\" [color=$color];\n";
            }
            elsif ( $network eq "all" ) {
                print GRAPH_OUT "\tnode \"$id1\";\n";
                print GRAPH_OUT "\tnode \"$id2\";\n";
            }
        }
    }
    print GRAPH_OUT "}";
    close(GRAPH_OUT);
}

### Subroutine is necessary if processing more than one unrelated file in a single run. Call this inbetween runs.
# @purpose Determine maximum network size based on connectivity
# @param $connectedness (scalar) - Network connectivity value 0-1
# @return (scalar) - Maximum recommended node count for that connectivity level
# @side_effects None (read-only)
# @notes Used to determine BronKerbosh feasibility; higher connectivity = larger networks feasible

sub get_max_network_size {
    my $connectedness = shift;
    if ( $connectedness <= .15 ) { return 60; }
    if ( $connectedness <= .2 )  { return 65; }
    if ( $connectedness <= .25 ) { return 70; }
    if ( $connectedness <= .35 ) { return 100; }
    if ( $connectedness <= .45 ) { return 130; }
    if ( $connectedness <= .55 ) { return 170; }
    if ( $connectedness <= .65 ) { return 230; }
    if ( $connectedness <= .75 ) { return 330; }
    if ( $connectedness <= .85 ) { return 500; }
    if ( $connectedness <= .95 ) { return 500; }
    return 500;
}

############ ALTERNATIVE METHOD SUBROUTINES #################################
# @purpose Run KING method and alternative PLINK method, compare results with PRIMUS
# @param $file (scalar) - Base filename
# @param $PRIMUS_unrelated_set (hashref) - PRIMUS result
# @param %networks_alt (hash) - Networks for alternative methods
# @return (void)
# @side_effects Creates output files for KING, PLINK methods; copies best result to final output
# @calls write_out_independent_set_KING, breakup_large_networks_PLINK, write_out_independent_set_PLINK
# @uses get_maximum_clique to determine best method

sub compare_alternative_methods {
    ### Run alternative methods

    my $file = shift;
    my $PRIMUS_unrelated_set = shift;
    my $id_id_scores = pop;
    my %networks_alt = @_;

    my %KING_network =
      write_out_independent_set_KING( $file, %networks_alt, $id_id_scores );

    breakup_large_networks_PLINK( \%networks_alt, $id_id_scores );

    my %PLINK_network = write_out_independent_set_PLINK( $file, %networks_alt );

    ### Compare KING, PLINK, and PRIMUS
    my $val = get_maximum_clique( \%PLINK_network, \%KING_network,
        $PRIMUS_unrelated_set, \@trait_ref );

    if ( $val eq 0 ) {
        ## PLINK performed the best
        system(
"cp $output_dir/$file\_maximum_independent_set_PLINK $output_dir/$file\_maximum_independent_set"
        );
    }
    elsif ( $val eq 1 ) {
        ## KING performed the best
        system(
"cp $output_dir/$file\_maximum_independent_set_KING $output_dir/$file\_maximum_independent_set"
        );
    }
    elsif ( $val eq 2 ) {
        ## PRIMUS performed the best
        system(
"cp $output_dir/$file\_maximum_independent_set_PRIMUS $output_dir/$file\_maximum_independent_set"
        );
    }
    else {
        die "ERROR!!!!!\n";
    }

    if ( $PRINT_ALTERNATE_RESULTS eq 0 ) {
        system("rm $output_dir/$file\_maximum_independent_set_PLINK");
        system("rm $output_dir/$file\_maximum_independent_set_KING");
        system("rm $output_dir/$file\_maximum_independent_set_PRIMUS");
    }
}

# @purpose Find maximum independent set using KING method (greedy algorithm)
# @param $file (scalar) - Base filename for output
# @param %networks_king (hash) - Networks to process
# @return (hash) - IID => 1 for selected unrelated individuals
# @side_effects Creates ${file}_maximum_independent_set_KING
# @calls King_method for each network
# @uses get_maximum_clique for clique selection

sub write_out_independent_set_KING {
    my $file         = shift;
    my $id_id_scores = pop;
    my %networks_king = @_;
    my %unrelated_set;
    open( UNIQUE_OUT, ">$output_dir/$file\_maximum_independent_set_KING" );
    print UNIQUE_OUT "FID\tIID\ta_status\n";

    ## get most unrelateds in each network
    foreach my $network ( sort { $a <=> $b } keys %networks_king ) {
        my @temp = @{ $networks_king{$network} };
        my %P    = map { $_ => 1 } @temp;
        my @maximal_cliques;

        King_method( \@maximal_cliques, \%P, $id_id_scores );

        my $maximum_clique = get_maximum_clique(@maximal_cliques, \@trait_ref);
        my @maximum_ids    = keys %{ $maximal_cliques[$maximum_clique] };
        write_out_maximum_clique_ids(@maximum_ids);
        foreach (@maximum_ids) { $unrelated_set{$_} = 1; }
    }
    close(UNIQUE_OUT);
    return %unrelated_set;
}

# @purpose KING method: greedy independent set selection based on inverse neighbors
# @param $maximal_cliques_ref (arrayref) - Output array (single clique)
# @param $network_ref (hashref) - Network nodes (nodes => 1)
# @return (void)
# @side_effects Pushes result clique into $maximal_cliques_ref
# @algorithm Greedy: iteratively select highest degree unrelated node
# @uses load_inverse_degrees_and_neighbors, load_degrees_and_neighbors, get_highest_degree_node
# @notes Alternative to BronKerbosh; faster but may produce smaller independent set

sub King_method {

    my ( $maximal_cliques_ref, $network_ref, $id_id_scores ) = (@_);

    my %actual_degrees;
    my %actual_neighbors;
    my %inverse_degrees;
    my %inverse_neighbors;
    my %i_set = ();

    load_inverse_degrees_and_neighbors( $network_ref, \%inverse_degrees,
        \%inverse_neighbors, $id_id_scores );
    load_degrees_and_neighbors( $network_ref, \%actual_degrees,
        \%actual_neighbors, $id_id_scores );

    my @sorted_ids = (
        sort { $inverse_degrees{$b} cmp $inverse_degrees{$a} }
          keys %inverse_degrees
    );

    while ( keys %inverse_degrees ) {
        my ( $id, $degree ) = get_highest_degree_node( \%inverse_degrees );
        delete $inverse_degrees{$id};
        my $add = 1;
        foreach my $test_id ( keys %i_set ) {
            if ( exists $actual_neighbors{$id}{$test_id} ) {
                $add = 0;
                last;
            }
        }
        if ( $add eq 1 ) {
            $i_set{$id} = 1;
        }
    }
    push( @$maximal_cliques_ref, \%i_set );
}

# @purpose Calculate inverse degree (unrelated neighbor count) for each node
# @param $network_ref (hashref) - Network nodes (nodes => value)
# @param $degree_ref (hashref) - Output: nodes => inverse_degree_count
# @param $neighbors_ref (hashref) - Output: nodes => unrelated_neighbors_hashref
# @return (void)
# @side_effects Modifies $degree_ref and $neighbors_ref hashrefs
# @uses get_inverse_neighbors to find unrelated nodes (complement graph)
# @notes Used by KING method; counts nodes with score <= THRESHOLD

sub load_inverse_degrees_and_neighbors {

    my ( $network_ref, $degree_ref, $neighbors_ref, $id_id_scores ) = (@_);

    foreach my $node ( keys %$network_ref ) {

        my $neighbors =
          get_inverse_neighbors( $node, $network_ref, undef, $id_id_scores );

        $$neighbors_ref{$node} = $neighbors;
        my @temp   = keys %{ $neighbors_ref->{$node} };    # %neighbors;
        my $degree = @temp;
        $$degree_ref{$node} = $degree;
    }
}

# @purpose Find maximum independent set using PLINK method
# @param $file (scalar) - Base filename for output
# @param %networks (hash) - Networks to process
# @return (hash) - IID => 1 for selected individuals (already singleton networks)
# @side_effects Creates ${file}_maximum_independent_set_PLINK
# @notes PLINK method outputs one individual per network (very conservative/small independent set)

sub write_out_independent_set_PLINK {
    my $file     = shift;
    my %networks = @_;
    my %unrelated_set;
    open( UNIQUE_OUT, ">$output_dir/$file\_maximum_independent_set_PLINK" );
    print UNIQUE_OUT "FID\tIID\n";

    ## get most unrelateds in each network
    foreach my $network ( sort { $a <=> $b } keys %networks ) {
        my @node = @{ $networks{$network} };

        if ( @node > 1 ) {
            die "ERROR!!! " . @node . " NODES. Should only be 1\n";
        }
        write_out_maximum_clique_ids(@node);
        foreach (@node) { $unrelated_set{$_} = 1; }
    }
    close(UNIQUE_OUT);
    return %unrelated_set;
}

# @purpose Break up large networks for PLINK method (one individual per network)
# @param $networks_ref (hashref) - Networks to modify
# @return (void)
# @side_effects Modifies $networks_ref by splitting networks
# @uses breakup_large_network_PLINK to split networks
# @algorithm Removes highest-degree nodes until max 1 individual per network

sub breakup_large_networks_PLINK {
    my ( $networks_ref, $id_id_scores ) = (@_);
    foreach my $network ( sort { $a <=> $b } keys %$networks_ref ) {

        my @temp = @{ $$networks_ref{$network} };

        #print "\n\nNEXT NETWORK $network: @temp\n";

        my %P = map { $_ => 1 } @temp;

        if ( keys %P > 0 ) {
            breakup_large_network_PLINK( $networks_ref, \%P, $id_id_scores );
            @{ $$networks_ref{$network} } = keys %P;
        }
    }
}

# @purpose Recursively split network down to single individuals (PLINK method)
# @param $networks_ref (hashref) - Main networks hash to add split components to
# @param $P_ref (hashref) - Current network nodes (nodes => 1)
# @return (void)
# @side_effects Modifies $P_ref, modifies $networks_ref by adding components
# @uses get_node_to_remove, reduce_neighbors, get_connected_components recursively
# @globals $network_ctr, %networks
# @algorithm Removes highest-degree/lowest-trait node until one node per network

sub breakup_large_network_PLINK {

    my ( $networks_ref, $P_ref, $id_id_scores ) = (@_);
    my %degrees;
    my %neighbors;

    load_degrees_and_neighbors( $P_ref, \%degrees, \%neighbors, $id_id_scores );

    while ( keys %$P_ref > 1 ) {
        my @P_nodes   = keys %$P_ref;
        my @degrees   = keys %degrees;
        my @neighbors = keys %neighbors;

        my ( $node_to_remove, $degree ) =
          get_node_to_remove( \%degrees, \%neighbors );
        if ( $degree > 0 ) {
            delete $$P_ref{$node_to_remove};
            reduce_neighbors( $node_to_remove, \%neighbors, \%degrees );
            delete $degrees{$node_to_remove};
            delete $neighbors{$node_to_remove};

            # Get connected components of new graph
            my @component_refs =
              get_connected_components( $P_ref, \%neighbors );
            %{$P_ref} = %{ $component_refs[0] };
            shift(@component_refs); #remove the first element from the $hash_ref

            foreach my $hash_ref (@component_refs) {
                my @temp = keys %$hash_ref;
                foreach (@temp) {
                    delete $degrees{$_};
                    delete $neighbors{$_};
                }

                ## if still too big, call this routine recursively
                if ( keys %$hash_ref > 1 ) {
                    breakup_large_network_PLINK( $networks_ref, $hash_ref,
                        $id_id_scores );
                }
                $network_ctr++;
                @{ $$networks_ref{$network_ctr} } = keys %$hash_ref;
            }

            # proceed with network %P until it is smaller than $MAX_NETWORK_SIZE
        }
        else {
            # I should never get here
            die "ERROR!!! PRUNING NODE WITHOUT RELATIVES!!!\n";
        }

        #exit;
    }

    # Replace $networks{$network} with keys %P
}

# @purpose Select node to remove from network (PLINK method)
# @param $degrees_ref (hashref) - Node => degree count (inverse neighbors)
# @param $neighbors_ref (hashref) - Adjacency information
# @return (list) - ($node_to_remove, $degree)
# @side_effects None (read-only)
# @algorithm Selects node with lowest trait values from highest degree's neighbor
# @uses get_highest_degree_node, weighted_comparison
# @notes PLINK uses inverse degree; breaks ties by trait preference

sub get_node_to_remove {
    my $degrees_ref   = shift;
    my $neighbors_ref = shift;
    my @nodes         = sort keys %$degrees_ref;
    my $node          = @nodes[0];
    my $degree        = $$degrees_ref{$node};
    if ( $degree == 0 ) {
        die "ERROR REMOVING A NODE WITHOUT NEIGHBORS $node\n";
    }
    my @neighbors       = sort keys %{ $neighbors_ref->{$node} };  # %neighbors;
    my $neighbor        = @neighbors[0];
    my $neighbor_degree = $$degrees_ref{$neighbor};
    my @node_values;
    my @neighbor_values;

    for ( my $i = 0 ; $i < @trait_order ; $i++ ) {
        @node_values[$i] = $trait_refs[$i]{$node};
    }
    for ( my $i = 0 ; $i < @trait_order ; $i++ ) {
        @neighbor_values[$i] = $trait_refs[$i]{$neighbor};
    }

    my $max = weighted_comparison( \@node_values, \@neighbor_values )
      ;    ## return 1 for node, and 2 for neighbor
    if ( $max == 2 ) {
        return ( $node, $degree );
    }
    return ( $neighbor, $neighbor_degree );
}

1;
