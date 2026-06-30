package Types::State;

use Moo;
use strict;
use warnings;

# ==============================================================================
# IMUS State
# ==============================================================================

# Key = network ID; Value = array of IIDs in that network
has 'networks' => (
    is      => 'rw',
    default => sub { {} },
);

# Key = IID; Value = network ID the individual belongs to
has 'id_network' => (
    is      => 'rw',
    default => sub { {} },
);

# Key = "IID1;IID2"; Value = PI_HAT score
has 'id_id_scores' => (
    is      => 'rw',
    default => sub { {} },
);

# Key = "IID1;IID2"; Value = original raw line from input genome file
has 'id_id_all_info' => (
    is      => 'rw',
    default => sub { {} },
);

# List of hashes containing loaded trait data
has 'trait_refs' => (
    is      => 'rw',
    default => sub { [] },
);

# Order of traits for prioritizing maximum independent set selection
has 'trait_order' => (
    is      => 'rw',
    default => sub { [] },
);

# Key = trait filename; Value = trait type (quantitative/binary)
has 'trait_files' => (
    is      => 'rw',
    default => sub { {} },
);

# Key = child IID; Value = array of parent IIDs
has 'child_parents' => (
    is      => 'rw',
    default => sub { {} },
);

# Counter for assigning new network IDs
has 'network_ctr' => (
    is      => 'rw',
    default => 0,
);

# Header string parsed from input relatedness file
has 'outfile_header' => (
    is      => 'rw',
    default => "",
);

# Key = numeric mapped ID; Value = original IID string
has 'iid_map' => (
    is      => 'rw',
    default => sub { {} },
);

# Key = numeric mapped ID; Value = original FID string
has 'fid_map' => (
    is      => 'rw',
    default => sub { {} },
);

# Key = original IID string; Value = numeric mapped ID
has 'ids_loaded' => (
    is      => 'rw',
    default => sub { {} },
);

# ==============================================================================
# Pedigree Reconstruction (PR) State
# ==============================================================================

# Key = IID; Value = sex (1 = male, 2 = female, etc.)
has 'gender' => (
    is      => 'rw',
    default => sub { {} },
);

# Key = IID; Value = age
has 'age' => (
    is      => 'rw',
    default => sub { {} },
);

# Key = IID; Value = affection status
has 'affected_status' => (
    is      => 'rw',
    default => sub { {} },
);

# Counter for assigning IDs to dummy parents/ancestors
has 'dummy_ctr' => (
    is      => 'rw',
    default => 1,
);

1;
