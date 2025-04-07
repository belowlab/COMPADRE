#! /usr/bin/perl

# Script to run PADRE

# Author          : Jeffrey Staples
# Created On      : Jan 22 2016
# Last Modified By: Jeffrey Staples
# Last Modified On: Jan 22 2016
# Update Count    : 0
# Status          : Unknown, Use with caution!

################ Common stuff ################
use strict;

my @commandline_options = @ARGV;

use Getopt::Long 2.13;
use PADRE;
use File::Path qw(make_path);

my $package = "PADRE";
my $lib_dir; 
my $bin_dir; 

my $pod2usage = sub {
        # Load Pod::Usage only if needed.
        require Pod::Usage;
        Pod::Usage->import;
        &pod2usage;
    };


################ Command line parameters ################

## NOTE: COLUMNS NUMBERS ARE 1 BASED IN THIS FILE; IT WILL BE ADJUSTED TO 0 BASED AS IT GETS PASSED TO OTHER FILES

## REQUIRED PARAMETERS
my $dataset_name;

## Command line options.
my $verbose = 1; 
my $study_name = "";
my $output_dir = "";
my $log_file;
my $LOG;

## PADRE variables
my $run_PADRE = 0;
my $ersa_model_output;
my $ersa_results;
my $project_summary_file;
my $PADRE_multiple_test_correct = 0;
my $relatedness_threshold = .09375; ## Values is halfway between the expected mean pi_hat for 3rd degree and 4th degree
my $degree_rel_cutoff = 3;

## Development options (not shown with -help).
my $debug = 0;			# debugging
my $test = 0;			# test mode.
my $cluster = 0;		# run on cluster

## Process command line options.
apply_options();
open($LOG,">$log_file") if($LOG eq "");

print $LOG "Commandline options used: @commandline_options\n";

################ Print all files and settings ################

print_files_and_settings() if $verbose > 0;

#################### RUN PROGRAMS ###########################
	my $results = PADRE::run_PADRE_project_summary($project_summary_file,$ersa_model_output,$ersa_results,$degree_rel_cutoff,$output_dir,$PADRE_multiple_test_correct);
	print "PADRE results: $results\n";

exit 0;

##################################################################################
########################## Subroutines ###########################################
##################################################################################
sub print_files_and_settings
{
	print "\nFILES AND COLUMNS\n";
	print $LOG "\nFILES AND COLUMNS\n";
	

	print "\nSETTINGS\n";
	print "Verbose: $verbose\n";
	print "Relatedness threshold: $relatedness_threshold\n";
}

sub apply_options 
{
    my $help = 0;		
    my $ident = 0;		
    my $man = 0;		
    
    # Process options.
    if ( @ARGV > 0 ) 
    {
	GetOptions(
		# Diagnostic options
		'verbose|v=i' => \$verbose,
		'test' => \$test,
		'help|h|?'  => \$help,
		'man'	  => \$man,
		'debug'	  => \$debug,
		'cluster' => \$cluster,

		# Settings
		"degree_rel_cutoff|d=i" => \$degree_rel_cutoff,
    "rel_threshold|threshold|t=f" => \$relatedness_threshold,
		"output_dir|o=s" => \$output_dir,

		# Files and directories
		"bin=s" => \$bin_dir,
		"lib=s" => \$lib_dir,

		## PADRE options
		"ersa_model_output=s" => \$ersa_model_output,
		"ersa_results=s" => \$ersa_results,
		"project_summary_file=s" => \$project_summary_file,
		"PADRE_multiple_test_correct" => \$PADRE_multiple_test_correct,		
		)
		or $pod2usage->(2);
    	}
    	if ( $man or $help ) {
		$pod2usage->(1) if $help;
		$pod2usage->(VERBOSE => 2) if $man;
    	}
	
	## PADRE requires ERSA likelihood file and results file as well as a PRIMUS summary file.
	if($ersa_model_output ne "" || $ersa_results ne "")
	{
		if($ersa_model_output ne "" && $ersa_results ne "")
		{
			if($project_summary_file ne "")
			{
				$run_PADRE = 1;
			}
      else
			{
				die "INVALID INPUTS: PADRE requires a PRIMUS project summary file; either provide one on the commandline or generate one by running PRIMUS.\n";
			}
		}
		else
		{
			die "INVALID INPUTS: PADRE requires but the model_output_fule and the results files output by ERSA.\n";
		}
	}

	## Set output directory if not passed in
  #if($output_dir eq "" && exists $ibd_estimates{'FILE'}){$output_dir = "$ibd_estimates{'FILE'}\_PRIMUS"}
	

	## If the relatedness threshold or degree_rel_cutoff changed, then adjust them both here
	if($degree_rel_cutoff != 3)
	{
		## If degree_rel_cutoff is unchanged, adjust it to account for the change in the $relatedness_threshold
		if($degree_rel_cutoff == 3)
		{
			$degree_rel_cutoff = 3 if $relatedness_threshold <= 0.1;
			$degree_rel_cutoff = 2 if $relatedness_threshold > 0.1 && $relatedness_threshold <= 0.2;
			$degree_rel_cutoff = 1 if $relatedness_threshold > 0.2;
		}
		elsif($degree_rel_cutoff == 2)
		{
			$relatedness_threshold = 0.1875;
		}
		elsif($degree_rel_cutoff == 1)
		{
			$relatedness_threshold = 0.375;
		}
		else
		{
			die "INVALID value for -d|--degree_rel_cutoff; Valid options are are integers from 1 to 3, representing 1st through 3rd degree relatives\n";  
		}
	}

	$log_file = "$output_dir/PRIMUS_output.log";
}


__END__

################ Documentation ################

=head1 NAME

run_PADRE.pl - Run PADRE on the results from PADRE and ERSA

=head1 SYNOPSIS

B<run_PADRE.pl> [I<options>] --ersa_model_output F<file> --ersa_results F<file> --project_summary F<file> --degree_rel_cutoff <1|2|3> | -h 


=head1 DESCRIPTION

B<run_PADRE.pl> takes in the two output files from ERSA, the project level summary file from PRIMUS in the relative positition to the other PRIMUS results for that project, and the degree of relatedness used to reconstruct the PRIMUS pedigrees. 


=head1 OPTIONS

 For usage and documentation:
   -h, --help		Brief help message
   --man        Full documentation
 
 General options:
   -t, --rel_threshold	Set the minimum level of relatedness for two people to be considered related (default=0.1)
   --degree_rel_cutoff	Set the maximum degree of relatedness for two people to be considered related (default=3; i.e. 3rd degree relatives)
   -o, --output_dir	Specify path to the output directory for all results(default=[PATH_TO_IBD_FILE]_PRIMUS/)
   -v, --verbose	Set verbosity level (0=none; 1=default; 2=more; 3=max)

 PADRE options:
   --ersa_model_output	Path to the model_output_file generated by ERSA
   --ersa_results	Path to the standard ERSA results file
   --project_summary	Path to the PRIMUS generated project level summary file (usualy in *_PRIMUS/ dir)
   --degree_rel_cutoff	Specify the minimum degree of relatedness used to generated pedigrees (default = 3)
   --PADRE_multiple_test_correct  Apply Bonferroni correction for the p-value threshold to determine statistical signficance of a relationship (default = no correction)


=head1 DOCUMENTATION

=head3 For usage and documentation:

=over 8

=item B<-h, --help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head3 PADRE options:

=over 8
=item B<OVERVIEW>

PADRE is an algorithm uses the pairwise relationship estimate likelihoods and results generated by ERSA to connect the pedigrees reconstructed by PRIMUS. This tool will generate two output files. One will contain the network numbers, pedigree numbers, and names of founder in each pedigree that is most likely to be connected to the other pedigree at the degree of relationship specified in the final column. There will be one line per pair of networks. Individuals who are unrelated at the relatedness threshold specified are treated as their own network. The other results file provides the same information but for every related pair of genotyped individuals in all the networks.

=item B<--ersa_model_output> F<path/file>

ERSA will generate this file when the B<--model_output_file> F<path/file> option specified on the commandline when running ERSA. For each pair of individuals in the dataset, it will contain the likelihood that they are related as 1st through 40th degree relatives sharing 0, 1, or two parents in common at each degree. PRIMUS uses these likelihhods and the likelihood that they are unrelated (obtained from the file specified with the B<--ersa_results> option) to find the most likely way each family network is connected.

=item B<--ersa_results> F<path/file>

ERSA's main output file (specified with the --output_file option) needs to be provided. PADRE uses the likelihood that a pair of individuals are unrelated in its algorithm. Provide the path to the file here.

=item B<--project_summary> F<path/Summary_*.genome.txt>

The project level summary file is in the main output directory which we call the Project level directory (default is *_PRIMUS). In the main ouput directory, there are two summary files. One is a summary containing the pairwise relationship table (*_pairwuse_table.txt), and this is NOT the file you want to provide with this option. You want to provide the other summary file.

Note: you can rerun PADRE in the same run as the pedigree reconstructions. In which case, you do not need to include this option.

=item B<--degree_rel_cutoff [1|2|3]>

If you do not reconsruct the pedigrees in the same run as PADRE, then you need to specify what degree of relatedness you used as your cutoff during reconstruction. Currently, there is no way for PADRE to know if you used something other than the default relationship cutoff (default is 3), so you must specify if you used a different degree. 

Alternative, if you reconstructed using the --threshold or -t option, then you can use the same option and value you used during reconstructionto correctly run PADRE. 


=item B<--PADRE_multiple_test_correct [1|2|3]>

Apply Bonferroni correction for the p-value threshold to determine statistical signficance of a relationship (default = no correction)


=head1 DEFINITIONS

=over 8

=item F<DIR>

A full path or a relative path to a directory.

=item F<FILE>

A full path or a relative path to an input files.

=item F<file_stem>

A full path or a relative path to an input file, but without including the file extensions (e.g. the file stem for "foo_bar.ped" is "foo_bar").

=item I<NUM>

A floating point number.

=back

=head1 OUTPUT FILES

=over 4

=item B<PADRE estimated pairwise relationships (PADRE_sample_relatedness.txt)>

PADRE estimated pairwise relationships based on the most likely pedigree connections. 

=item B<PADRE estimated network connections (PADRE_network_connections.txt)>

PADRE estimated most likely connectes between networks. 

=item B<PADRE estimated network connections image (PADRE_network_connections.dot)>

Graph image of the PADRE estimated most likely connectes between networks. 

=back

=head1 QUICK START

There is one example dataset in the F<example_data/> directory. 

Example. Run PADRE on the ERSA and PRIMUS results

=over 8

../bin/run_PADRE.pl --ersa_model_output test_ERSA_models --ersa_results test_ERSA_results --project_summary test.genome_PRIMUS/Summary_test.genome.txt --degree_rel 1
B<./run_PADRE.pl> B<--ersa_model_output> F<../example_data/test_ERSA_models> B<--ersa_results> F<../example_data/test_ERSA_results> B<--project_summary> F<../example_data/test.genome_PRIMUS/Summary_test.genome.txt> B<--degree_rel> 1

This command will run PADRE on the example dataset and write the results to F<../example_data/test.genome_PRIMUS/>. 
 
=back


=cut






