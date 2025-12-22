package CLI;

use strict;
use GetOpt::Long::Descriptive;
use FindBin;
use Pod::Usage;

sub new_cli {
    my ($cli) = @_;

    my ($opt, $usage) = describe_options(
      'run_COMPADRE.pl [options] -p file | -i FILE=file [options] | -h',
      [],# This blank list allow for spaces
      ["For usage and documentation:"],
      ["help|h", "Brief help message", { shortcircuit => 1 }],
      [], 
      ["Required options (one of either of the following):"], # This line will serve as a section heading
      ["plink_ibd|p", "Specify path to a .genome IBD estimates files produced by PLINK"],
      ["input|i=s", "Specify path to an IBD estimates file and additional column information"],
      ["(or provide both of the following):"],
      ["file=s", "Path to PLINK formatted data without the file extensions. Behaves the same as in PLINK (this flag requires the '--genome' flag)"],
      ["genome", "Read in --file and calculate IBD estimates using PLINK"],
      []
      ["COMPADRE options (new):"],
      ["segment_data=s", "Shared pairwise IBD segment data in a format reaable by ERSA (see full documentation for more details)."],
      ["port_number", "Port number used for additinoal Python computation", { default => 6000 }],
      ["run_padre", "Run PADRE after standard pedigree reconstruction is complete"],
      [],
      ["Gereral options:"],
      ["rel_threshold|t=f", "Set the minimum level of relatedness for two people to be considered related", {default => 0.1}],
      ["degree_rel_cutoff=i", "Set the maximum degree of relatedness for two people to be considered related", {default => 3}],
      ["output_dir|o=s", "Specify path to the output directory for all results"], # TODO: add default value
      ["verbose|v=i", "Set verbosity level (0=none; 1=default; 2=more; 3=max)", {default=>1}],
      [],
      ["prePRIMUS IBD estimation options:"],
      ["--plink_ex=s", "Path to the plink executable file (searches environment variables by default)"],
      ["ref_pops", "Comma separated list of 1000 Genomes populations used for reference allele freqs (overrides default method)"],
      ["no_automatic", "Turun off automatic selection of the HapMap3 populations for reference allele freqs (On by default)."], #TODO: Also add default here. Uncertain if this is a bool or string
      ["remove_AIMs", "Automatically remove ancetsry informative markers (off by default)."], #TODO: Add default for same reason as above
      ["internal_ref", "Use the dataset provided in --file to get reference allele frequencies"],
      ["alt_ref_stem", "Path to PLINK formatted data (no file extensions) used for allele frequencies"],
      ["keep_inter_files", "Keep intermediate files used to create teh IBD estimates with prePRIMUS"],
      ["min_pihat_threshold", "Set a minimum pi-hat threshold that will be used in the plink --genome calculation"],
      ["max_memory", "Specify amoutn of memory to be used in PLINK prePRIMUS command (in MB)"],
      [],
      ["Identification of maximum unrelated set (IMUS) options:"],
      ["--no_IMUS", "Don't identify a maximum unrelated set (by default IMUS is run)"],
      ["missing_val", "Set value that denotes missing dat in IBD file"],
      ["size|s", "Specify to weight on set size (Default unless a bianry trait is specified first)"],
      ["high_btrait", "File with FID, IID, and binary trait to weight for the higher value"],
      ["low_btrait", "File with FID, IID, and binary trait to weight for the lower value"],
      ["high_qtrait", "File with FID, IID, and quantitative trait to weight for higher values"],
      ["low_qtrait", "File with FID, IID, and quantitative trait to weight for lower values"],
      ["mean_qtrait", "File with FID, IID, and quantitative trait to weight towards the mean value"],
      ["tails_qtrait", "File with FID, IID, andquantitative trait to weight against the middle values"],
      [],
      ["Pedigree reconstruction (PR) options:"],
      ["no_PR", "Don't reconstruct pedigress (runs pedigree reconstruction by default)"],
      ["max_gens", "Max number of generations to be sampled in reconstructed pedigrees. (default = no limit)"],
      ["max_gen_gap=i", "Max number of generations between two people that have a child (default = 0)"],
      ["age_file=s", "Specify path to the file containing the age of each sample."],
      ["ages=s", "Like --age_file but requires FILE=[file], optional specification of file columns."],
      ["--sex_file=s", "Specify path to the file containing the sex of each sample"],
      ["sexes=s", "Like --sex-file but requires FILE=[file], optional specification of file columns"],
      ["mito_matches=s", "Path to mito matching status for each pair of samples. (requires File=[file])"],
      ["y_matches", "Path to y matching status for each pair of samples (requires FILE=[file])"],
      ["MT_error_rate", "Proportion of the MT sequence that must not match to be call a non-match"],
      ["Y_error_rate", "Proportion of the Y sequence that must not match to be called a non-match"],
      ["affection_file=s", "Speccify path to the file containing the affection status of each sample"],
      ["affections", "Like --afection_file. Need File=[file], optional specification of file columns"],
      ["int_likelihood_cutoff=f", "Initial minimum likelihood for a relationship to reconstruction (default=0.1)", {default => 0.1}],
      [],
      ["ERSA options:"],
      ["min_cm=f", "minimum segment size to consider (default=2.5)", {default => 2.5}],
      ["max_cm=f", "maximum segment size to consider for estimating the exponential disribution of segment sizes in the population. (default=10)", {default=> 10.0}],
      ["max_meioses=i", "maximum number of meioses to consider. (default=40)",  {default => 40}],
      ["rec_per_meioses=f", "expected number of recombination events per meioses (default=35.2548101)", {default => 35.2548101}],
      ["ascertained_chromosome=i", "chromsome number of ascertained disease locus"],
      ["ascertained_position=i", "chromosomal position of ascertained disease locus"],
      ["control_files", "GERMLINE or Beagle fibd output files(s) for population control"],
      ["control_sample_size", "Sample size of control population. Used with --control_files"],
      ["--exp_mean=f", "Mean of exp distribution of shared segment size in population (default=3.197036753)", {default => 3.19703}],
      ["pois_mean=f", "Mean of the Poisson distribution of the number of segments shared between a pair of individuals in the population. (default=13.73)", {default => 13.73}],
      ["pair_file=s", "Restrict pairwise comparisons to the ID pairs specified in this file."],
      ["single_pair=s", "Restrict pairwise comparisons to the pairs specified in this flag (id1:id2)"],
      ["number_of_ancestors", "Restrict relationships to [1] one parent (hal-sibs/cousins), [2] two parents (full-sibs/cousins), or [0] (parent-offspring/grandparent-granchild)."],
      ["number_of_chromosomes=i", "Number of chromosomes (default=22)", {default=22}],
      ["parent_offspring_option=s", "Option to evaluate potential parent-offspring and sibling relationships based on total proprtion of the genome that is shared ibd1 (default=true)", {default=>"true"}],
      ["parent_offspring_zscore=f", "Z-score for rejecting a sibling relationship in favor of a parent-offspring relationship (default=2.33, alpha=0.01)", {default=> 2.33}],
      ["adjust_pop_dist", "Option to adjust the population distribution of shared segments downward for segments that could no be detected due to recent ancestry. (default=false)", {default="false"}],
      ["confidence_level=f", "Confidence level for confidence interval around the estimated degree of relationship. (default=0.95)", {default=>0.95}],
      ["mask_common_shared_regions=s", "excludes chromosomal regions that are commonly shared from evaluation. Used only when the control_files or mask_region_file parameter is specified (default=false)", {default=>"false"}],
      ["mask_region_cross_length", "length in base pairs that a shared segment must extend past a masked segment in order to avoid truncation. Used only when mask_common_shared_regions parameter is specified (default=1000000)", {default => 1000000}],
      ["mask_region_file=s", "File containing chromosomal regions to exclude from evaluation."],
      ["mask_region_threshold", "Threhsold for the ratio of observed vs. expected segment sharing in controls before a region will be masked."],
      ["mask_region_sim_count", "number of simulations performed of the null distribution of shared sgemnet locations in controls; results written to output_file.sim"],
      ["recombination_files", "file containing genetic distances for all chromsomes. This parameter must be specified with Beagle fibd input files."],
      ["beagle_markers_files", "Beagle marker files (one file required for each chromosome, wildcards required, ex: chr*beagle.marker). Each filename must begin with the chromosome numae followed by a period. This parameter must be specified with Beagle fibd input files."]
      ["show_docs=s", "Show the complete documentation file", {default="false"}]
    )

    if ($opt->help) {
        print $usage->text;
        exit 0;
    }

    # If the user passes the flags to view the documentation then we need to invoke the pod2usage function for teh docs.pod
    if ($opt->show_docs == "true") {
      my $docs_path = $FindBin::Bin/docs.pod
      if (! -e $docs_path) {
        die "The documentation file was not found"
      } 
      pod2usage(
        -input   => $docs_path,  # <--- The key parameter
        -exitval => 0,
        -verbose => 2
    );
    }
}

sub validate_args {
  return 1;
}
