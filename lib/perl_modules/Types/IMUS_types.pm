#####################################
## IMUS Configuration and State Classes
## Data class containers for IMUS.pm
#####################################

package PRIMUS::IMUS::Config;
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

# @purpose Data class for IMUS configuration parameters
# @param %args - Named parameters for configuration
# @return (object) - Blessed hash reference

sub new {
    my ($class, %args) = @_;

    my $threshold_val = $args{threshold} // 0.1; # Default threshold value

    if (!looks_like_number($threshold_val) || $threshold_val < 0 || $threshold_val > 1) {
        die "ERROR: Threshold value '$threshold_val' is not a valid number. Values should be between 0 and 1.\n";
    }
    return bless {
        # Threshold settings
        threshold            => ($args{threshold} // 0.1) + 0, # We are forcing the threshold value to be a number
        min_likelihood       => ($args{min_likelihood} // 0.1) + 0,
        exclude_value        => $args{exclude_value} // 0,
        lowest_max_network_size => $args{lowest_max_network_size} // 60,
        
        # I/O and output
        output_dir            => $args{output_dir},
        outputfile_header     => $args{outputfile_header} // "",
        relatedness_file      => $args{relatedness_file},
        relatedness_file_name => $args{relatedness_file_name} // "relatedness_output.txt",
        samples_file          => $args{samples_file},
        ersa_data             => $args{ersa_data} // "none",
        ped_file              => $args{ped_file} // "none",
        missingness_file      => $args{missingness_file} // "none",
        ibd_file_ref          => $args{ibd_file_ref} // "",
        
        # Trait configuration
        trait_order          => $args{trait_order} // [],
        trait_files          => $args{trait_files} // {},
        
        # Column indices
        relatedness_column   => $args{relatedness_column} // -1,
        id1_column           => $args{id1_column} // -1,
        id2_column           => $args{id2_column} // -1,
        fid1_column          => $args{fid1_column} // -1,
        fid2_column          => $args{fid2_column} // -1,
        trait_fid_columns    => $args{trait_fid_columns} // {},
        trait_id_columns     => $args{trait_id_columns} // {},
        trait_data_columns   => $args{trait_data_columns} // {},
        
        # Genome file column index mapping
        # Maps field names to their 0-based column indices in standard PLINK genome files
        genome_file_columns  => $args{genome_file_columns} // {
            FID1   => 0,
            IID1   => 1,
            FID2   => 2,
            IID2   => 3,
            RT     => 4,
            EZ     => 5,
            Z0     => 6,
            Z1     => 7,
            Z2     => 8,
            PI_HAT => 9,
            PHE    => 10,
            DST    => 11,
            PPC    => 12,
            RATIO  => 13,
        },
        
        # Processing flags
        do_IMUS              => $args{do_IMUS} // 1,
        do_PR                => $args{do_PR} // 1,
        print_alternate      => $args{print_alternate} // 1,
        verbose              => $args{verbose} // 0,
        
        # Other settings
        lib_dir              => $args{lib_dir},
        log_file_handle      => $args{log_file_handle} // "",
    }, $class;
}

1;

package PRIMUS::IMUS::State;
use strict;
use warnings;

# @purpose Data class for IMUS runtime state
# @param %args - Named parameters for state initialization
# @return (object) - Blessed hash reference

sub new {
    my ($class, %args) = @_;
    return bless {
        # Network data structures
        networks             => $args{networks} // {},
        id_network           => $args{id_network} // {},
        id_id_scores         => $args{id_id_scores} // {},
        id_id_all_info       => $args{id_id_all_info} // {},
        
        # Trait data
        trait_refs           => $args{trait_refs} // [],
        
        # Relationship tracking
        child_parents        => $args{child_parents} // {},
        
        # Counter for network numbering
        network_ctr          => $args{network_ctr} // 0,

        # Header for the output file
        outfile_header       => $args{outfile_header} // "",
        # IID mapping
        iid_map               => $args{id_map} // {},
        #FID mapping
        fid_map               => $args{fid_map} // {},
        #Keep track of the ids read in. We are going to treat this like a set
        ids_loaded            => $args{ids_loaded} // {},
    }, $class;
}

1;