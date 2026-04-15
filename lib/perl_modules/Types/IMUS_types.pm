#####################################
## IMUS Configuration and State Classes
## Data class containers for IMUS.pm
#####################################

package PRIMUS::IMUS::Config;
use strict;
use warnings;

# @purpose Data class for IMUS configuration parameters
# @param %args - Named parameters for configuration
# @return (object) - Blessed hash reference

sub new {
    my ($class, %args) = @_;
    return bless {
        # Threshold settings
        threshold            => $args{threshold} // 0.1,
        min_likelihood       => $args{min_likelihood} // 0.1,
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
    }, $class;
}

1;