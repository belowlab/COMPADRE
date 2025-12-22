#! /usr/bin/perl

# Script to run PRIMUS

# Author          : Jeffrey Staples
# Created On      : May 30 2014
# Last Modified By: Jeffrey Staples
# Last Modified On: May 30 2014
# Update Count    : 0
# Status          : Unknown, Use with caution!

################ Common stuff ################
use strict;
use File::Find;
use File::Basename;
use File::Spec;
#use warnings;
use IO::Socket::INET; # added for socket compadre helper connectionuse IPC::Open2;
use IPC::Open2; # also for new compadre helper

my @commandline_options = @ARGV;

#use lib '../lib/perl_modules';
use Getopt::Long 2.13;
use PRIMUS::IMUS qw(run_IMUS);
use PRIMUS::reconstruct_pedigree_v7;
use PRIMUS::predict_relationships_2D;
use PRIMUS::prePRIMUS_pipeline_v7;
use PRIMUS::PRIMUS_plus_ERSA;
use File::Path qw(make_path);

use Cwd qw(abs_path);
my $script_path = abs_path($0);
my $project_root = dirname(dirname($script_path));

my @onekg_pops = ("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI");

my $package = "PRIMUS";
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
my %ibd_estimates = ("FID1",1,"IID1",2,"FID2",3,"IID2",4,"IBD0",7,"IBD1",8,"IBD2",9,"PI_HAT",10); ## Defaults to PLINK .genome file format
my %mito_matches = ("FID1",1,"IID1",2,"FID2",3,"IID2",4,"MATCH",5,"MATCH_VAL",1);
my %y_matches = ("FID1",1,"IID1",2,"FID2",3,"IID2",4,"MATCH",5,"MATCH_VAL",1);
my $dataset_name;

## Command line options
my $verbose = 1; 
my $study_name = "";
my $output_dir = "";
my $ersa_data = "";
my $port_number = 6000;
my $run_padre = 0;  # New flag for PADRE
my $reference_pop = "";
my $log_file;
my $LOG;

# Store PID globally for signal handlers
our $compadre_pid;
our $cleanup_called = 0; 
our $port_number_glob;

sub cleanup_compadre {
    my $sig = shift;
    
    # Prevent multiple cleanup calls
    return if $cleanup_called;
    $cleanup_called = 1;
    
    print "\nReceived signal $sig, cleaning up COMPADRE helper...\n" if $verbose > 0;
    
    if (defined $compadre_pid) {
        # Skip graceful shutdown for forced termination - just kill it
        if ($sig eq 'INT' || $sig eq 'TERM' || $sig eq 'QUIT') {
            if (kill(0, $compadre_pid)) {
                print "Force terminating COMPADRE helper (PID: $compadre_pid)\n" if $verbose > 0;
                kill('KILL', $compadre_pid);  # Use KILL immediately for forced shutdown
                sleep(1);
            }
        } else {
            # For normal END cleanup, try graceful first
            eval {
                my $shutdown_ack = send_to_compadre_helper("close", $port_number);
                print "COMPADRE graceful shutdown: $shutdown_ack\n" if $verbose > 0;
            };
            
            sleep(1);
            
            if (kill(0, $compadre_pid)) {
                print "Force terminating COMPADRE helper (PID: $compadre_pid)\n" if $verbose > 0;
                kill('TERM', $compadre_pid);
                sleep(1);
                kill('KILL', $compadre_pid) if kill(0, $compadre_pid);
            }
        }
    }
    
    exit(0) unless $sig eq 'END';
}

# Set up signal handlers
$SIG{INT}  = \&cleanup_compadre;  # CTRL+C
$SIG{TERM} = \&cleanup_compadre;  # Termination
$SIG{QUIT} = \&cleanup_compadre;  # CTRL+\

END {
    # This runs when the script exits for any reason
    cleanup_compadre('END') if defined $compadre_pid;
}

#######################################################################################

## PRIMUS variables
my %sexes = ("FID",1,"IID",2,"SEX",3,"MALE",1,"FEMALE",2);
my %affections = ("FID",1,"IID",2,"AFFECTION",3,"AFFECTION_VALUE",2);
my %ages = ("FID",1,"IID",2,"AGE",3);
my %traits;
my @trait_order;
my $relatedness_threshold = .09375; ## Values is halfway between the expected mean pi_hat for 3rd degree and 4th degree
my $degree_rel_cutoff = 3;
my $max_memory = 0;
my $initial_likelihood_cutoff = 0.3; ## Arbitrary cutoff for the relationship likelihood vectors
my $max_generations = "none"; ## The maximum number of generations in a pedigree allowed during reconstruction
my $max_generation_gap = 0; ## What is maximum different in generations that 2 people can mate
my $missing_data_value = 0;
my $get_max_unrelated_set = 1; ## Defaults to determining the maximum unrelated set
my $reconstruct_pedigrees = 1; ## Defaults to reconstructing pedigrees
my $generate_likelihood_vectors_only = 0; ## Defaults to reconstructing pedigrees

## prePRIMUS variables
my $run_prePRIMUS = 0;
my $data_stem = "";
my $plink_path = "plink";
my $no_automatic_IBD = 0;
my $remove_AIMs = 0;
my $keep_AIMs = 0;
my $no_mito = 1;
my $no_y = 1;
my $use_mito_match = 0;
my $use_y_match = 0;
my $MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH = 0.01;
my $Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH = 0.01;
my $internal_ref = 0;
my $alt_ref_stem = "";
my $keep_prePRIMUS_intermediate_files = 1;
my $no_PCA_plot = 0;
my @ref_pops;
my $rerun = 0;
my $min_pihat_threshold = 0;

## PADRE variables
my $run_PRIMUS_plus_ERSA = 0;
my $ersa_model_output;
my $ersa_results;
my $project_summary_file;
my $PADRE_multiple_test_correct = 0;

## Development options (not shown with -help)
my $debug = 0;			# debugging
my $test = 0;			# test mode.
my $cluster = 0;		# run on cluster
my $public_html_dir = "";

## ERSA options (new)
my %ersa_settings = (
	'min_cm' => 2.5,
	'max_cm' => 10,
	'max_meioses' => 40,
	'rec_per_meioses' => 35.2548101,
	'ascertained_chromosome' => 'no_ascertainment',
	'ascertained_position' => -1,
	'control_sample_size' => 'None',
	'exp_mean' => 3.197036753,
	'pois_mean' => 13.73,
	'number_of_ancestors' => 'None',
	'number_of_chromosomes' => 22,
	'parent_offspring_option' => 'true',
	'parent_offspring_zscore' => 2.33,
	'adjust_pop_dist' => 'false',
	'confidence_level' => 0.95,
	'mask_common_shared_regions' => 'false',
	'mask_region_cross_length' => 1000000,
	'mask_region_threshold' => 4.0,
	'mask_region_sim_count' => 0,
	'control_files' => 'None',
	'mask_region_file' => 'None',
	'recombination_files' => 'None',
	'beagle_markers_files' => 'None',
);


## Process command line options
apply_options();

if(-d $output_dir) # The output directory exist we are going to rename
{
	system("mv $output_dir $output_dir\_OLD");	
} 
# Then we need to create the directory to output to
make_path($output_dir) if !-d $output_dir;

open($LOG,">$log_file") if($LOG eq "");

print $LOG "Commandline options used: @commandline_options\n";

################ Print all files and settings ################

print_files_and_settings() if $verbose > 0;

if($max_generations eq 'none'){$max_generations = 1000}

# if(-d $output_dir)
# {
# 	system("mv $output_dir $output_dir\_OLD");
# }


#################### RUN PROGRAMS ###########################
if($generate_likelihood_vectors_only)
{
  generate_likelihood_vectors_only();
  print "done.\n";
  exit;
}

if($run_prePRIMUS){
	run_prePRIMUS();
}

## Run IMUS to get family networks and maximum unrelated set (unless turned off)
if($reconstruct_pedigrees || $get_max_unrelated_set)
{
	my @IMUS_commands = ("--do_IMUS",$get_max_unrelated_set,"--do_PR",$reconstruct_pedigrees,"--ibd_estimates",\%ibd_estimates,"--verbose",$verbose,"--trait_order",\@trait_order,"--traits",\%traits,"--output_dir",$output_dir,"--rel_threshold",$relatedness_threshold,"--lib",$lib_dir,"--int_likelihood_cutoff",$initial_likelihood_cutoff,"--log_file_handle",$LOG);
	print "IMUS_commands: @IMUS_commands\n" if $verbose > 1;
	print $LOG "IMUS_commands: @IMUS_commands\n" if $verbose > 1;
	if(!PRIMUS::IMUS::run_IMUS(@IMUS_commands)){die "IMUS FAILED TO COMPLETE\n\n"}
}

## Run pedigree reconstruction
if($reconstruct_pedigrees){

	check_names();
	run_PR();
}

## Run PRIMUS+ERSA
if($run_PRIMUS_plus_ERSA)
{
	my $results = PRIMUS::PRIMUS_plus_ERSA::run_PRIMUS_plus_ERSA_project_summary($project_summary_file,$ersa_model_output,$ersa_results,$degree_rel_cutoff,$output_dir,$PADRE_multiple_test_correct);
	print "PADRE results: $results\n";
}


exit 0;

##################################################################################
########################## Subroutines ###########################################
##################################################################################

sub generate_likelihood_vectors_only
{
  print "Generating likelihood vector file...\n";
  if(!-d $output_dir)
  {
    make_path($output_dir);
  }
  PRIMUS::predict_relationships_2D::get_relationship_likelihood_vectors(\%ibd_estimates,$initial_likelihood_cutoff ,$verbose,$lib_dir,$output_dir);
}

sub run_prePRIMUS
{
	## Test that all these paths exists
	my %paths = ("LIB",$lib_dir,"PLINK",$plink_path);
	PRIMUS::prePRIMUS_pipeline_v7::set_paths(\%paths);

	# Create the prePRIMUS directory with specific permissions
	my $preprimus_dir = "$output_dir/$study_name\_prePRIMUS";
	make_path($preprimus_dir, { mode => 0755 }) if !-d $preprimus_dir;
	
	my @IBD_commands = ("--verbose",$verbose,"--study_name",$study_name,"--output_dir",$preprimus_dir,"--lib",$lib_dir,"--file",$data_stem,"--rerun",$rerun,"--ref_pops_ref",\@ref_pops,"--remove_AIMs",$remove_AIMs,"--keep_AIMs",$keep_AIMs,"--internal_ref",$internal_ref,"--alt_ref",$alt_ref_stem,"--no_PCA_plot",$no_PCA_plot,"--keep_intermediate_files",$keep_prePRIMUS_intermediate_files,"--no_automatic_IBD",$no_automatic_IBD,"--rel_threshold",$relatedness_threshold,"--log_file_handle",$LOG,"--MT_error_rate",$MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH,"--Y_error_rate",$Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH,"--no_mito",$no_mito,"--no_y",$no_y, "--min_pihat_threshold", $min_pihat_threshold, "--max_memory",$max_memory);

	## Run the PLINK IBD pipeline
	my ($temp_genome_file,$temp_sex_file,$temp_mt_match_file,$temp_y_match_file) = PRIMUS::prePRIMUS_pipeline_v7::run_prePRIMUS_main(@IBD_commands);
	$ibd_estimates{'FILE'} = $temp_genome_file;
	if(!exists $sexes{'FILE'} && -e $temp_sex_file)
	{
		$sexes{'FILE'} = $temp_sex_file;
		$sexes{'SEX'} = 4;
	}
	if(!exists $mito_matches{'FILE'} && -e $temp_mt_match_file)
	{
		$mito_matches{'FILE'}=$temp_mt_match_file;
	}
	if(!exists $y_matches{'FILE'} && -e $temp_y_match_file)
	{
		$y_matches{'FILE'}=$temp_y_match_file;
	}
	$dataset_name = PRIMUS::prePRIMUS_pipeline_v7::get_file_name_from_stem($temp_genome_file);
	print "Dataset name: $dataset_name\n" if $verbose > 0;
	print $LOG "Dataset name: $dataset_name\n" if $verbose > 0;
	print "PrePRIMUS done.\n" if $verbose > 0;
	print $LOG "PrePRIMUS done.\n" if $verbose > 0;
}

sub run_PR {

	## Load networks
	my %networks;
	open(IN,"$output_dir/$dataset_name\_networks") or die "can't open $output_dir/$dataset_name\_networks; $!\n\n";
	<IN>; ## get rid of header
	while(my $line = <IN>)
	{
		my ($network,@rest) = split(/\s+/,$line);
		$networks{$network} = 1;
	}
	close(IN);
	
	## Load additional attributes files into simple hash references
	my $sexes_ref; 
	my $affections_ref; 
	my $ages_ref; 
	if(exists $sexes{'FILE'})
	{
		if(!-e $sexes{'FILE'})
		{
			warn "$sexes{'FILE'} does not exists; continuing without sex info.\n";
			delete $sexes{'FILE'};
		}
		else
		{
			$sexes_ref = get_sample_info($sexes{'FILE'},$sexes{'FID'},$sexes{'IID'},$sexes{'SEX'});
			foreach my $id (keys %$sexes_ref)
			{
				if($$sexes_ref{$id} eq $sexes{'MALE'}){$$sexes_ref{$id} = 1}
				elsif($$sexes_ref{$id} eq $sexes{'FEMALE'}){$$sexes_ref{$id} = 2}
				else{$$sexes_ref{$id} = 0}
				if($verbose > 1)
				{
					print "$id sex = $$sexes_ref{$id}\n";
					print $LOG "$id sex = $$sexes_ref{$id}\n";
				}
			}
		}
	}
	if(exists $affections{'FILE'}){$affections_ref = get_sample_info($affections{'FILE'},$affections{'FID'},$affections{'IID'},$affections{'AFFECTION'})}
	if(exists $ages{'FILE'}){$ages_ref = get_sample_info($ages{'FILE'},$ages{'FID'},$ages{'IID'},$ages{'AGE'})}
	
	## Set the project summary file in case PRIMUS+ERSA is also called
	$project_summary_file = "$output_dir/Summary_$dataset_name.txt";

	## Submit each network to be reconstructed
	foreach my $net (sort {$a <=> $b} keys %networks)
	{
		my $network_name = "$dataset_name\_network$net";
		my $network_dir = "$output_dir/$network_name";
		if(!-d $network_dir)
		{
			make_path($network_dir);
		}
		system("mv $output_dir/$network_name.genome $network_dir");
		
		my %network_ibd_estimates = %ibd_estimates;
		$network_ibd_estimates{'FILE'} = "$network_dir/$network_name.genome";


		if($cluster)
		{	
			#my $new_network_dir = "../$network_dir";
			my $new_network_dir = "$network_dir";
			my @PR_commands = ("--network",$network_name,"--verbose",$verbose,"--output_dir",$new_network_dir,"--int_likelihood_cutoff",$initial_likelihood_cutoff,"--network_num",$net,"--affection_status_value",$affections{'AFFECTION_VALUE'},"--sex_file","$network_dir/sex.temp","--age_file","$network_dir/age.temp","--affection_file","$network_dir/affection.temp","--max_gen",$max_generations, "--lib",$lib_dir, "--bin",$bin_dir,"--no_mito",$no_mito,"--no_y",$no_y,"--use_mito_match",$use_mito_match,"--use_y_match",$use_y_match,"--degree_rel_cutoff",$degree_rel_cutoff);
			foreach my $key (keys %network_ibd_estimates)
			{
				push(@PR_commands, "--$key");
				#$network_ibd_estimates{$key} = "../$network_ibd_estimates{$key}" if $key =~ /FILE/;
				push(@PR_commands, "$network_ibd_estimates{$key}");
				#print "--$key = $network_ibd_estimates{$key}\n";
			}
			foreach my $key (keys %mito_matches) ## MIGHT NOT WORK YET
			{
				push(@PR_commands, "--MITO_$key");
				#$mito_matches{$key} = "../$mito_matches{$key}" if $key =~ /FILE/;
				push(@PR_commands, "$mito_matches{$key}");
				#print "--$key = $network_ibd_estimates{$key}\n";
			}
			foreach my $key (keys %y_matches) ## MIGHT NOT WORK YET
			{
				push(@PR_commands, "--Y_$key");
				#$y_matches{$key} = "../$y_matches{$key}" if $key =~ /FILE/;
				push(@PR_commands, "$y_matches{$key}");
				#print "--$key = $network_ibd_estimates{$key}\n";
			}
			
			open(SEX,">$network_dir/sex.temp");
			foreach my $key (keys %$sexes_ref){print SEX "$key\t$$sexes_ref{$key}\n";}
			close(SEX);
			open(AGE,">$network_dir/age.temp");
			foreach my $key (keys %$ages_ref){print AGE "$key\t$$ages_ref{$key}\n";}
			close(AGE);
			open(AFFECTION,">$network_dir/affection.temp");
			foreach my $key (keys %$affections_ref){print AFFECTION "$key\t$$affections_ref{$key}\n";}
			close(AFFECTION);
			#print "pr: @PR_commands\n";
			my $file = "$network_dir/run_dataset_cluster_$network_name.sh";
			print "writing $file\n" if $verbose > 0;
			print $LOG "writing $file\n" if $verbose > 0;
			open(OUT,">$file") or die "Can't open $file;$!\n\n";
			print OUT "\#\$ -S /bin/bash\n";
			print OUT "\#\$ -l h_vmem=8G -l mem_requested=8G\n";
			#print OUT "\#\$ -l mem_requested=12G\n";
			print OUT "\#\$ -e $bin_dir/$new_network_dir/run_dataset_cluster_$network_name.e\n";
			print OUT "\#\$ -q nick-grad.q\n";
			print OUT "\#\$ -o $bin_dir/$new_network_dir/run_dataset_cluster_$network_name.o\n";
			print OUT "module load modules modules-init modules-gs\n";
			print OUT "module load modules\n";
			print OUT "module load modules-init\n";
			print OUT "#module load perl/5.14.2\n\n";
			print OUT "cd $bin_dir\n";
			print OUT "pwd\n";
			print OUT "GSITSHOST=`/bin/hostname`\n";
			print OUT "GSITSPWD=`/bin/pwd`\n";
			print OUT "GSITSDATE=`/bin/date`\n"; 
			print OUT "echo \"**** JOB STARTED ON \$GSITSHOST AT \$GSITSDATE\"\n";
			print OUT "echo \"**** JOB RUNNING IN \$GSITSPWD\"\n\n";
			print OUT "PERL5LIB=\"\"\n"; ## REMOVE THIS 
			print OUT "./run_reconstruction.pl @PR_commands\n";
			#print OUT "rm $network_dir/sex.temp\n";
			#print OUT "rm $network_dir/age.temp\n";
			#print OUT "rm $network_dir/affection.temp\n";
			print OUT "./make_dataset_summary.pl $output_dir $dataset_name\n";
			print OUT "./make_dataset_pairwise_summary.pl $output_dir $dataset_name";
			close(OUT);
			system("qsub -js 12 $file");
		}
		else
		{
			my @PR_commands = ("--network",$network_name,"--ibd_estimates",\%network_ibd_estimates,"--y_hash",\%y_matches,"--mito_hash",\%mito_matches,"--verbose",$verbose,"--output_dir",$network_dir,"--int_likelihood_cutoff",$initial_likelihood_cutoff,"--max_gen",$max_generations,"--sex_ref",$sexes_ref,"--age_ref",$ages_ref,"--affection_ref",$affections_ref,"--network_num",$net,"--affection_status_value",$affections{'AFFECTION_VALUE'}, "--lib",$lib_dir,"--bin",$bin_dir,"--log_file_handle",$LOG,"--no_mito",$no_mito,"--no_y",$no_y,"--use_mito_match",$use_mito_match,"--use_y_match",$use_y_match,"--degree_rel_cutoff",$degree_rel_cutoff);
			eval
			{
				PRIMUS::reconstruct_pedigree_v7::reconstruct_pedigree(@PR_commands);
			};
			if($@)
			{
				warn "Reconstruction failed for $network_name: $@\n";
			};
		}
	}
	
	if($verbose > 0){print "Writing dataset summary file $output_dir/Summary_$dataset_name.txt\n"}
	if($verbose > 0){print $LOG "Writing dataset summary file $output_dir/Summary_$dataset_name.txt\n"}
	system("$bin_dir/make_dataset_summary.pl $output_dir $dataset_name");
	system("./make_dataset_pairwise_summary.pl $output_dir $dataset_name");

	print "done.\n" if $verbose > 0;
	print $LOG "done.\n" if $verbose > 0;

	# convert ps to pdf files 
	find(\&process_file, $output_dir);

	#print "converted .ps to .pdf.\n" if $verbose > 0;
	#print $LOG "converted .ps to .pdf.\n" if $verbose > 0;

	###############################################################################################

	if ($run_padre) {

		print "\nInitializing PADRE\n" if $verbose > 0;
		print $LOG "Initializing PADRE\n" if $verbose > 0;

		print "Requesting ERSA output path from COMPADRE helper... (port $port_number_glob)\n" if $verbose > 1;
		print $LOG "Requesting ERSA output path from COMPADRE helper... (port $port_number_glob)\n" if $verbose > 1;

		my $ersa_output_path_prefix = PRIMUS::predict_relationships_2D::send_to_compadre_helper("padre", $port_number_glob);

		print "Received ERSA output path prefix: $ersa_output_path_prefix\n" if $verbose > 1;
		print $LOG "Received ERSA output path prefix: $ersa_output_path_prefix\n" if $verbose > 1;

		# Get the exact ones for regular and model files
		my $ersa_model_output_path = "$ersa_output_path_prefix.model";
		my $ersa_output_path = "$ersa_output_path_prefix.out";

		my $padre_command = "perl $project_root/padre/bin/run_PADRE.pl --ersa_model_output $ersa_model_output_path --ersa_results $ersa_output_path --project_summary $output_dir/Summary_$dataset_name.txt --degree_rel_cutoff $degree_rel_cutoff --output_dir $output_dir/padre_results";

		system($padre_command) == 0
			or warn "Failed to run PADRE: $padre_command\n";

		print "PADRE complete.\n\n" if $verbose > 0;
		print $LOG "PADRE complete.\n\n" if $verbose > 0;

	}

	my $shutdown_ack = PRIMUS::predict_relationships_2D::send_to_compadre_helper("close", $port_number_glob);
	print "COMPADRE socket shutdown message: $shutdown_ack\n" if $verbose > 0;	

	## Write out pairwise Summary file based on the results in the Summary file and possible pedigrees
	#my $rels_ref = PRIMUS::get_pairwise_summary::get_possible_relationships($output_dir,"$output_dir/Summary_$dataset_name.txt");
	#PRIMUS::get_pairwise_summary::write_table("$output_dir/Summary_$dataset_name\_pairwise_table.txt",$rels_ref);
}

# Convert .ps from cranefoot to .pdf 
sub process_file 
{
    if (/\.(ps)$/i) {
        my $file = $File::Find::name;
        (my $output_file = $file) =~ s/\.ps$/.pdf/i;
        my $command = "ps2pdf -dFirstPage=2 $file $output_file";
        system($command) == 0 or warn "Failed to execute: $command\n";
    }
}

sub print_files_and_settings {

	print "\nFILES AND COLUMNS\n";
	print $LOG "\nFILES AND COLUMNS\n";
	
	print "LOG FILE: $log_file\n";

	print "Data stem: $data_stem\n" if $data_stem ne "";
	print $LOG "Data stem: $data_stem\n" if $data_stem ne "";
	print "IBD file: $ibd_estimates{'FILE'} (FID1=".($ibd_estimates{'FID1'})."; IID1=".($ibd_estimates{'IID1'})."; FID2=".($ibd_estimates{'FID2'})."; IID2=".($ibd_estimates{'IID2'})."; IBD0=".($ibd_estimates{'IBD0'})."; IBD1=".($ibd_estimates{'IBD1'})."; IBD2=".($ibd_estimates{'IBD2'})."; PI_HAT/RELATEDNESS=".($ibd_estimates{'PI_HAT'}).")\n" if exists $ibd_estimates{'FILE'};
	print $LOG "IBD file: $ibd_estimates{'FILE'} (FID1=".($ibd_estimates{'FID1'})."; IID1=".($ibd_estimates{'IID1'})."; FID2=".($ibd_estimates{'FID2'})."; IID2=".($ibd_estimates{'IID2'})."; IBD0=".($ibd_estimates{'IBD0'})."; IBD1=".($ibd_estimates{'IBD1'})."; IBD2=".($ibd_estimates{'IBD2'})."; PI_HAT/RELATEDNESS=".($ibd_estimates{'PI_HAT'}).")\n" if exists $ibd_estimates{'FILE'};
	
	print "Dataset results dir: $output_dir\n";
	print $LOG "Dataset results dir: $output_dir\n";

	if(!exists $ages{'FILE'}){print "Age file: none\n";print $LOG "Age file: none\n";}
	elsif(!-e $ages{'FILE'}){print "Age file: $ages{'FILE'} does not exists\n"; $pod2usage->(2);print $LOG "Age file: $ages{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print "Age file: $ages{'FILE'} (FID=".($ages{'FID'})."; IID=".($ages{'IID'})."; AGE=".($ages{'AGE'}).")\n"; print $LOG "Age file: $ages{'FILE'} (FID=".($ages{'FID'})."; IID=".($ages{'IID'})."; AGE=".($ages{'AGE'}).")\n";}

	if(!exists $sexes{'FILE'}){print "Sex file: none\n";}
	elsif(!$run_prePRIMUS && !-e $sexes{'FILE'}){print "Sex file: $sexes{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print "Sex file: $sexes{'FILE'} (FID=$sexes{'FID'}; IID=$sexes{'IID'}; SEX=$sexes{'SEX'}; MALE=$sexes{'MALE'}, FEMALE=$sexes{'FEMALE'})\n";}
	if(!exists $sexes{'FILE'}){print $LOG "Sex file: none\n";}
	elsif(!$run_prePRIMUS && !-e $sexes{'FILE'}){print $LOG "Sex file: $sexes{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print $LOG "Sex file: $sexes{'FILE'} (FID=$sexes{'FID'}; IID=$sexes{'IID'}; SEX=$sexes{'SEX'}; MALE=$sexes{'MALE'}, FEMALE=$sexes{'FEMALE'})\n";}

	print "Segment data file: $ersa_data\n" if $ersa_data ne "";
	print $LOG "Segment data file: $ersa_data\n" if $ersa_data ne "";

	print "Port number: $port_number\n" if $port_number ne "" && $port_number != 6000;
	print $LOG "Port number: $port_number\n" if $port_number ne "" && $port_number != 6000;

	# if ($ersa_data ne "") # checking if an argument was actually passed for ersa data
	# {
		#print "\nSegment file: $ersa_data\n";

	# get absolute path of compadre helper 
	my $libpath = $lib_dir;
	$libpath =~ s{/$}{};
	my $parent_dir = File::Spec->catdir(dirname($libpath));
	my $helper_path = File::Spec->catfile($parent_dir, 'compadre.py');

	my ($reader, $writer);

	# Check if ersa data is passed in at runtime. If it is, send that path, and if not, send 'NA'
	my $ersa_arg = ($ersa_data ne "") ? $ersa_data : "NA";
	if ($ersa_data ne "") {
		print "\nLaunching COMPADRE helper (with segment data) ...\n";
	}
	else {
		print "\nLaunching COMPADRE helper (no segment data) ...\n";
	}

	my $ersa_flags = "";
	foreach my $key (keys %ersa_settings) {
		my $value = $ersa_settings{$key};
		# Only include non-empty string values and non-zero numeric defaults where appropriate
		if ($value ne "") {
			$ersa_flags .= "$key|$value|";
		}
	}

	# Open socket
	#my $pid = open2($reader, $writer, "python3 $helper_path $ersa_arg $port_number \"$ersa_flags\" $output_dir");
	my $pid = open2($reader, $writer, 'python3', $helper_path, $ersa_arg, $port_number, $ersa_flags, $output_dir);

	if (!defined $pid) {
		die "Failed to launch COMPADRE helper: $!\n\n";
	}
	$compadre_pid = $pid;

	# Wait for the server to be ready and check for port changes
    my $actual_port = $port_number;  # Default to requested port

	while (my $line = <$reader>) {

        # Check for port change notification (comes via stderr which open2 merges)
        if ($line =~ /Port changed: (\d+)/) {
            $actual_port = $1;
			$port_number = $actual_port;

            print "\n[COMPADRE] Port $port_number was in use, using port $actual_port instead.\n";
            print $LOG "[COMPADRE] Port $port_number was in use, using port $actual_port instead.\n";
        }

        # Always parse the actual port from the "ready" line on stdout
		elsif ($line =~ /COMPADRE helper is ready on port\s+(\d+)/) {
			$actual_port = $1;  # <-- CRITICAL: take the authoritative port from stdout
			chomp $line;
			print "\n$line\n";

			# Overwrite the working port so all later steps use it
			$port_number = $actual_port;

			last;
		}

        elsif ($line =~ /ERROR|FAILED|Exception/) {
            die "COMPADRE error: $line\n";
        }
    }
	
	our $compadre_pid = $pid;

	print "\nReference file specification: $reference_pop\n\n" if $reference_pop ne "";
	print $LOG "Reference file specification: $reference_pop\n\n" if $reference_pop ne "";

	# Standard arguments 

	our $ersa_data_glob = $ersa_data;
	our $port_number_glob = $port_number;
	#print "\n[COMPADRE] Global port number: $port_number_glob\n";
	our $reference_pop_glob = $reference_pop;

	# ERSA additional arguments 

	our %ersa_settings_glob = %ersa_settings;

	#############################################################################

	if(!exists $affections{'FILE'}){print "Affection file: none\n";}
	elsif(!-e $affections{'FILE'}){print "Affection file: $affections{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print "Affection file: $affections{'FILE'} (FID=$affections{'FID'}; IID=$affections{'IID'}; AFFECTION=$affections{'AFFECTION'}; AFFECTION_VALUE=$affections{'AFFECTION_VALUE'})\n";}
	if(!exists $affections{'FILE'}){print $LOG "Affection file: none\n";}
	elsif(!-e $affections{'FILE'}){print $LOG "Affection file: $affections{'FILE'} does not exists\n"; $pod2usage->(2)}
	else{print $LOG "Affection file: $affections{'FILE'} (FID=$affections{'FID'}; IID=$affections{'IID'}; AFFECTION=$affections{'AFFECTION'}; AFFECTION_VALUE=$affections{'AFFECTION_VALUE'})\n";}

	print "Trait weighting:\n";
	print $LOG "Trait weighting:\n";
	foreach my $trait_file (@trait_order)
	{
		if($get_max_unrelated_set)
		{
			print "\t$trait_file ";
			print $LOG "\t$trait_file ";
			if(exists $traits{$trait_file}){print"($traits{$trait_file})"};
			if(exists $traits{$trait_file}){print $LOG "($traits{$trait_file})"};
			print"\n";
			print $LOG "\n";
		}
	}


	print "\nSETTINGS\n";
	print "Get PLINK IBD ESTIMATES with prePRIMUS: $run_prePRIMUS\n";
	print "Automatic reference population selection: " .!$no_automatic_IBD."\n";
	print "Verbose: $verbose\n";
	print "Relatedness threshold: $relatedness_threshold\n";
	print "Initial likelihood cutoff: $initial_likelihood_cutoff\n";
	print "Max generations: $max_generations\n";
	print "Max generational mating gap: $max_generation_gap\n";
	print "Get max unrelated set: $get_max_unrelated_set\n";
	print "Reconstruct pedigrees: $reconstruct_pedigrees\n";
	print $LOG "\nSETTINGS\n";
	print $LOG "Get PLINK IBD ESTIMATES with prePRIMUS: $run_prePRIMUS\n";
	print $LOG "Automatic reference population selection: " .!$no_automatic_IBD."\n";
	print $LOG "Verbose: $verbose\n";
	print $LOG "Relatedness threshold: $relatedness_threshold\n";
	print $LOG "Initial likelihood cutoff: $initial_likelihood_cutoff\n";
	print $LOG "Max generations: $max_generations\n";
	print $LOG "Max generational mating gap: $max_generation_gap\n";
	print $LOG "Get max unrelated set: $get_max_unrelated_set\n";
	print $LOG "Reconstruct pedigrees: $reconstruct_pedigrees\n";
}

sub get_sample_info
{
	my $file = shift;
	my $FID_col = shift;
	my $IID_col = shift;
	my $trait_col = shift;

	my %hash;
	
	open(IN,$file) or die "Can't open $file; $!";

	#print "file: $file\n";
	#print "$FID_col :: $IID_col\n";
	while(my $line = <IN>)
	{
		$line =~ s/^\s+//;
		my @temp = split(/\s+/,$line);
		my $FID = @temp[$FID_col-1];
		my $IID = @temp[$IID_col-1];

		#print "$FID :: $IID\n";

		if($FID =~ /\*\*/ || $FID =~ /,/)
		{
			print "File $file has invalid FID: $FID\n";
			print $LOG "File $file has invalid FID: $FID\n";
			die "FIDs and IIDs cannot contain '**' or ','\n";
		}
		if($IID =~ /\*\*/ || $IID =~ /,/)
		{
			print "File $file has invalid IID: $IID\n";
			print $LOG "File $file has invalid IID: $IID\n";
			die "FIDs and IIDs cannot contain '**' or ','\n";
		}

		$hash{"$IID"} = @temp[$trait_col-1];
	}
	
	return (\%hash);
}


sub apply_options {
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
		'public_html' => \$public_html_dir,

		# Settings
		"run_padre" => \$run_padre,
		"rel_threshold|threshold|t=f" => \$relatedness_threshold, 
		"degree_rel_cutoff|d=i" => \$degree_rel_cutoff,
		"max_memory|m=i" => \$max_memory,
		"max_gens=i" => \$max_generations,
		"max_gen_gap=i" => \$max_generation_gap,
		"int_likelihood_cutoff=f" => \$initial_likelihood_cutoff,
		"no_IMUS" => sub{$get_max_unrelated_set = 0},
		"no_PR" => sub{$reconstruct_pedigrees = 0},
    	"generate_likelihoods_only" => sub{$get_max_unrelated_set = 0;$reconstruct_pedigrees = 0;$generate_likelihood_vectors_only=1},
		"missing_val=f"=> \$missing_data_value,
		"genome" => sub{$run_prePRIMUS = 1},
		"no_automatic_IBD" => \$no_automatic_IBD,
		"remove_AIMs" => \$remove_AIMs,
		"keep_AIMs" => \$keep_AIMs,
		"MT_error_rate=f" => \$MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH,
		"Y_error_rate=f" => \$Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH,
		"no_mito" => \$no_mito,
		"no_y" => \$no_y,
		"use_mito_match" => \$use_mito_match,
		"use_y_match" => \$use_y_match,
		"internal_ref" => \$internal_ref,
		"alt_ref_stem|alt-ref|alt=s" => \$alt_ref_stem,
		"keep_inter_files" => \$keep_prePRIMUS_intermediate_files,
		"no_PCA_plot" => \$no_PCA_plot,
		"rerun" => \$rerun,
		"min_pihat_threshold=f" => \$min_pihat_threshold,
		"ref_pops=s" => sub 
		{
			@ref_pops = split(/,/,@_[1]);
			uc(@ref_pops);
			foreach my $pop (@ref_pops)
			{
				if($pop =~ /none/i){next}
				
				if(!grep ($_ eq $pop, @onekg_pops ))
				{
					print "\n\n[COMPADRE] Error: Invalid 1KG population: $pop\n";
					print "Must be a comma seperated list these: @onekg_pops\n";
					print "For example: --ref_pops CEU,TSI,YRI\n\n";
					print $LOG "\n\n[COMPADRE] Error: Invalid 1KG population: $pop\n";
					print $LOG "Must be a comma seperated list these: @onekg_pops\n";
					print $LOG "For example: --ref_pops CEU,TSI,YRI\n\n";
					$pod2usage->(2);
				}
			}
		},

		# Files and directories
		"bin=s" => \$bin_dir,
		"lib=s" => \$lib_dir,
		"plink_ex=s" => \$plink_path,
		"bfile=s" => sub  		{
			$data_stem=$_[1];
		},
		"file=s" => sub  
		{
			$data_stem=$_[1];
		},
		"ped=s" => sub  ## included for compatibility
		{
			$data_stem = $_[1] if $data_stem eq "";
			$data_stem =~ s/\.ped$//;
		}, 		
		"map=s" => sub  ## included for compatibility
		{
			$data_stem = $_[1] if $data_stem eq "";
			$data_stem =~ s/\.map$//;
		}, 		
		"plink_ibd|PLINK|p=s" => sub{$ibd_estimates{'FILE'} = $_[1]},
		"sex_file=s" => sub{$sexes{'FILE'} = $_[1]},
		"age_file=s" => sub{$ages{'FILE'} = $_[1]},
		"affection_file=s" => sub{$affections{'FILE'} = $_[1]},
		"output_dir|o=s" => \$output_dir,
		"segment_data|s:s" => \$ersa_data,
		"port_number:i" => \$port_number,
		"reference_pop:s" => \$reference_pop,
		"ages=s{1,4}" => sub
		{ 
			my @possible_keys = qw(FILE FID IID AGE);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$ages{$key} = $value;
			}
			else
			{
				print "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0] option: @_[1]\n";
				print $LOG "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0] option: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},
		"sexes=s{1,6}" => sub
		{ 
			my @possible_keys = qw(FILE FID IID SEX MALE FEMALE);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$sexes{$key} = $value;
			}
			else
			{
				print "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},
		"affections=s{1,5}" => sub
		{ 
			my @possible_keys = qw(FILE FID IID AFFECTION AFFECTION_VALUE);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$affections{$key} = $value;
			}
			else
			{
				print "\n\n[COMPADRE] Error: Invalid key before \"=\" for --affections option: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\n[COMPADRE] Error: Invalid key before \"=\" for --affections option: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},	
		"mito_data|mito_matches|mito=s{1,7}" => sub
		{ 
			my @possible_keys = qw(FILE FID1 IID1 FID2 IID2 MATCH MATCH_VAL);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$mito_matches{$key} = $value;
			}
			else
			{
				print "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},
		"y_data|y_matches|y=s{1,7}" => sub
		{ 
			my @possible_keys = qw(FILE FID1 IID1 FID2 IID2 MATCH MATCH_VAL);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$y_matches{$key} = $value;
			}
			else
			{
				print "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},
		"input|i=s{1,9}" => sub
		{ 
			my @possible_keys = qw(FILE FID1 IID1 FID2 IID2 IBD0 IBD1 IBD2 PI_HAT RELATEDNESS);
			my ($key,$value) = split(/=/,@_[1]);
			uc($key);
			if(grep ($_ eq $key, @possible_keys ))
			{
				$ibd_estimates{$key} = $value;
			}
			else
			{
				print "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print "Must be one of these: @possible_keys\n\n";
				print $LOG "\n\n[COMPADRE] Error: Invalid key before \"=\" for --@_[0]: @_[1]\n";
				print $LOG "Must be one of these: @possible_keys\n\n";
				$pod2usage->(2);
			}
		},

		## PRIMUS+ERSA options
		"ersa_model_output=s" => \$ersa_model_output,
		"ersa_results=s" => \$ersa_results,
		"project_summary_file=s" => \$project_summary_file,
		"PADRE_multiple_test_correct" => \$PADRE_multiple_test_correct,		

		## ERSA options (new)
		"min_cm=f" => \$ersa_settings{'min_cm'}, # 2.5
		"max_cm=f" => \$ersa_settings{"max_cm"}, # 10
		"max_meioses=i" => \$ersa_settings{"max_meioses"}, # 40
		"rec_per_meioses=f" => \$ersa_settings{"rec_per_meioses"}, # 35.2548101
		"ascertained_chromosome=s" => \$ersa_settings{"ascertained_chromosome"}, # 'no_ascertainment' (string)
		"ascertained_position=i" => \$ersa_settings{"ascertained_position"}, # -1
		"control_files=s" => \$ersa_settings{"control_files"}, # file specified or not 
		"control_sample_size=i" => \$ersa_settings{"control_sample_size"}, # none -- only used with ersa_control_files
		"exp_mean=f" => \$ersa_settings{"exp_mean"}, # 3.197036753
		"pois_mean=f" => \$ersa_settings{"pois_mean"}, # 13.73
		"number_of_ancestors=i" => \$ersa_settings{"number_of_ancestors"}, # default none
		"number_of_chromosomes=i" => \$ersa_settings{"number_of_chromosomes"}, # 22
		"parent_offspring_option=s" => \$ersa_settings{"parent_offspring_option"}, # true
		"parent_offspring_zscore=f" => \$ersa_settings{"parent_offspring_zscore"}, # 2.33
		"adjust_pop_dist=s" => \$ersa_settings{"adjust_pop_dist"}, # false
		"confidence_level=f" => \$ersa_settings{"confidence_level"}, # 0.95
		"mask_common_shared_regions=s" => \$ersa_settings{"mask_common_shared_regions"}, # false
		"mask_region_cross_length=i" => \$ersa_settings{"mask_region_cross_length"}, # 1000000
		"mask_region_file=s" => \$ersa_settings{"mask_region_file"}, # file specified or not
		"mask_region_threshold=f" => \$ersa_settings{"mask_region_threshold"}, # 4.0
		"mask_region_sim_count=i" => \$ersa_settings{"mask_region_sim_count"}, # 0
		"recombination_files=s" => \$ersa_settings{"recombination_files"}, # file specified or not
		"beagle_markers_files=s" => \$ersa_settings{"beagle_markers_files"}, # file specified or not

		## IMUS options
		"size:s" => sub{ push(@trait_order,"size");$traits{'size'} = "-size";},
		"high_btrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"low_btrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"high_qtrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"low_qtrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"mean_qtrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"tails_qtrait=s" => sub{ push(@trait_order,@_[1]);$traits{@_[1]} = @_[0];},
		"plink_ibd|PLINK|p=s" => sub{$ibd_estimates{'FILE'} = $_[1]},
    	"likelihood_vectors=s" => sub{$ibd_estimates{'FILE'} = $_[1]; $ibd_estimates{'likelihood_vectors'}=1},
		"sex_file=s" => sub{$sexes{'FILE'} = $_[1]},
		"age_file=s" => sub{$ages{'FILE'} = $_[1]},
		"affection_file=s" => sub{$affections{'FILE'} = $_[1]},
		"output_dir|o=s" => \$output_dir,
		#"ersa_data=s" => \$ersa_data,
		)
		or $pod2usage->(2);
    	}
    	if ( $man or $help ) {
		$pod2usage->(1) if $help;
		$pod2usage->(VERBOSE => 2) if $man;
    	}
	
	if($data_stem ne "" && $run_prePRIMUS != 1)
	{
		die "INVALID INPUT OPTIONS: --file input also requires the --genome flag to run.\n\n";
	}

	if($run_padre && $ersa_data eq "")
	{
		die "\n[COMPADRE] Error: --run_padre requires --segment_data to be specified.\n\n";
	}

	## PRIMUS + ERSA requires ERSA likelihood file and results file as well as a PRIMUS summary file.
	if($ersa_model_output ne "" || $ersa_results ne "")
	{
		if($ersa_model_output ne "" && $ersa_results ne "")
		{
			if($project_summary_file ne "")
			{
				$run_PRIMUS_plus_ERSA = 1;
			}
			elsif(!exists $ibd_estimates{'FILE'} && $data_stem eq "")
			{
				die "INVALID INPUTS: PRIMUS+ERSA requires a PRIMUS project summary file; either provide one on the commandline or generate one by running PRIMUS.\n\n";
			}
			else
			{
				$run_PRIMUS_plus_ERSA = 1;
			}
		}
		else
		{
			die "INVALID INPUTS: PRIMUS+ERSA requires but the model_output_fule and the results files output by ERSA.\n\n";
		}
	}

	if($run_prePRIMUS)
	{
		if(exists $ibd_estimates{'FILE'})
		{
			die "[COMPADRE] Error: If you already have IBD estimates, then you don't need to calculate new ones with --genome option. If you do need to calculate IBD estimates, don't input anything for the IBD estimates.\n\n";
		}
		$study_name = PRIMUS::prePRIMUS_pipeline_v7::get_file_name_from_stem($data_stem);
		$output_dir = "$data_stem\_PRIMUS" if $output_dir eq ""; 
		#print "out_dir: $output_dir\n" if $verbose > 0;
	}
	
	## Set output directory if not passed in
	if($output_dir eq "" && exists $ibd_estimates{'FILE'}){$output_dir = "$ibd_estimates{'FILE'}\_PRIMUS"}
	

	## Test that ibd_estimates were passed in
	if(!exists $ibd_estimates{'FILE'} && $data_stem eq "")
	{
		if(!$run_PRIMUS_plus_ERSA)
		{
			print "\n[COMPADRE] No arguments provided. Please run with the --help flag for details.\n\n";
			$pod2usage->(2);
		}
		else
		{
			$reconstruct_pedigrees = 0;
			$get_max_unrelated_set = 0;
		}
	}
	if(!-e $ibd_estimates{'FILE'} && $data_stem eq "")
	{
		if(!$run_prePRIMUS && !$run_PRIMUS_plus_ERSA)
		{
			print "[COMPADRE] Error: IBD INPUT FILE DOES NOT EXIST: $ibd_estimates{'FILE'}\n";
			$pod2usage->(2);
		}
	}


	## Set genome_file_name: if there is a path to the file, ignore all the path, and select the name of the file
	if($ibd_estimates{'FILE'} =~ /\/([^\/]+)$/){$dataset_name = $1}
	else{$dataset_name = $ibd_estimates{'FILE'} }

	## Add size as first trait if it is a qtrait
	if($traits{$trait_order[0]} =~ /qtrait/)
	{
		unshift(@trait_order,"size");
	}
	## Add size to traits if not already in it
	if(!grep ($_ eq "size", @trait_order ))
	{
		push (@trait_order,"size");
		$traits{'size'} = "size";
	}
	
	if(exists $ibd_estimates{'RELATEDNESS'})
	{
		$ibd_estimates{'PI_HAT'} = $ibd_estimates{'RELATEDNESS'};
	}
	
	## If the relatedness threshold or degree_rel_cutoff changed, then adjust them both here
	if($relatedness_threshold != .09375 || $degree_rel_cutoff != 3)
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
			die "INVALID value for -d|--degree_rel_cutoff; Valid options are are integers from 1 to 3, representing 1st through 3rd degree relatives\n\n";  
		}
	}

	$log_file = "$output_dir/COMPADRE_output.log";
}

sub do_arrays_match
{
	my $arr1_ref = shift;
	my $arr2_ref = shift;

	if(@$arr1_ref ne @$arr2_ref)
	{
		return 0;
	}
	foreach (@$arr1_ref)
	{
		my $val = $_;
		if(!grep($_ eq $val,@$arr2_ref))
		{
			return 0;
		}
	}
	foreach (@$arr2_ref)
	{
		my $val = $_;
		if(!grep($_ eq $val,@$arr1_ref))
		{
			return 0;
		}
	}
	return 1;
}

sub check_names
{
	## Load additional attributes files into simple hash references
	my $samples1_ref;
	my $samples2_ref;
	my $sexes_ref; 
	my $affections_ref; 
	my $ages_ref; 
	
	$samples1_ref = get_sample_info($ibd_estimates{'FILE'},$ibd_estimates{'FID1'},$ibd_estimates{'IID1'},$ibd_estimates{'PI_HAT'});
	$samples2_ref = get_sample_info($ibd_estimates{'FILE'},$ibd_estimates{'FID2'},$ibd_estimates{'IID2'},$ibd_estimates{'PI_HAT'});
	
	if(exists $sexes{'FILE'})
	{
		if(!-e $sexes{'FILE'})
		{
			warn "$sexes{'FILE'} does not exists; continuing without sex info.\n";
			delete $sexes{'FILE'};
		}
		else
		{
			$sexes_ref = get_sample_info($sexes{'FILE'},$sexes{'FID'},$sexes{'IID'},$sexes{'SEX'});
			foreach my $id (keys %$sexes_ref)
			{
				if($$sexes_ref{$id} eq $sexes{'MALE'}){$$sexes_ref{$id} = 1}
				elsif($$sexes_ref{$id} eq $sexes{'FEMALE'}){$$sexes_ref{$id} = 2}
				else{$$sexes_ref{$id} = 0}
				if($verbose > 1)
				{
					print "$id sex = $$sexes_ref{$id}\n";
					print $LOG "$id sex = $$sexes_ref{$id}\n";
				}
			}
			foreach my $sample (keys %$samples1_ref)
			{
				if(!exists $$sexes_ref{$sample})
				{
					warn "WARNING: sample $sample is not included in the sex file $sexes{'FILE'}\n";
				}
			}
			foreach my $sample (keys %$samples2_ref)
			{
				if(!exists $$sexes_ref{$sample})
				{
					warn "WARNING: sample $sample is not included in the sex file $sexes{'FILE'}\n";
				}
			}
		}
	}
	if(exists $affections{'FILE'})
	{
		$affections_ref = get_sample_info($affections{'FILE'},$affections{'FID'},$affections{'IID'},$affections{'AFFECTION'});
		foreach my $sample (keys %$samples1_ref)
		{
			if(!exists $$affections_ref{$sample})
			{
				warn "WARNING: sample $sample is not included in the affections file $affections{'FILE'}\n";
			}
		}
		foreach my $sample (keys %$samples2_ref)
		{
			if(!exists $$affections_ref{$sample})
			{
				warn "WARNING: sample $sample is not included in the affections file $affections{'FILE'}\n";
			}
		}

	}
	if(exists $ages{'FILE'})
	{
		$ages_ref = get_sample_info($ages{'FILE'},$ages{'FID'},$ages{'IID'},$ages{'AGE'});
		foreach my $sample (keys %$samples1_ref)
		{
			if(!exists $$ages_ref{$sample})
			{
				warn "WARNING: sample $sample is not included in the age file $ages{'FILE'}\n";
			}
		}
		foreach my $sample (keys %$samples2_ref)
		{
			if(!exists $$ages_ref{$sample})
			{
				warn "WARNING: sample $sample is not included in the age file $ages{'FILE'}\n";
			}
		}
	}
}

__END__

################ Documentation ################

=head1 NAME

run_COMPADRE.pl - Run an IBD file through COMPADRE

=head1 SYNOPSIS

B<run_COMPADRE.pl> [I<options>] -p F<file> | -i FILE=F<file> [I<options>] | -h 


=head1 DESCRIPTION

B<run_COMPADRE.pl> will read genome-wide IBD estimates and will identify a maximum unrelated set and/or all possible pedigrees for each family within the dataset. You may weight the selection of your unrelated set by loading in files containing values on which to weight. For pedigree reconstruction, you may also load information on sex, age, and affection status for each sample. PRIMUS is also able to generate IBD estimates for SNP data in PLINK input file formats (ped/map or bed/bim/fam) as long as the system also has the following software installed: R, PLINK.

=head1 OPTIONS

 For usage and documentation:
   -h, --help		Brief help message
 
 Required options (one of the following):
   -p, --plink_ibd	Specify path to a .genome IBD estimates file produced by PLINK
   -i, --input		Specify path to an IBD estimates file and additional column information
   
   (or both):
   --file 		Path to PLINK formatted data without the file extensions; behaves the same as in PLINK (requires --genome)
   --genome		Read in --file and calculate IBD estimates using PLINK

 COMPADRE options (new):
   --segment_data       Shared segments data in format readable by ERSA (see full documentation)
   --port_number        Port number for additional Python computation (default: 6000)
   --run_padre		Run PADRE after standard pedigree reconstruction is complete

 General options:
   -t, --rel_threshold	Set the minimum level of relatedness for two people to be considered related (default=0.1)
   --degree_rel_cutoff	Set the maximum degree of relatedness for two people to be considered related (default=3; i.e. 3rd degree relatives)
   -o, --output_dir	Specify path to the output directory for all results(default=[PATH_TO_IBD_FILE]_PRIMUS/)
   -v, --verbose	Set verbosity level (0=none; 1=default; 2=more; 3=max)

 prePRIMUS IBD estimation options:
   --file 		Path to PLINK formatted data without the file extensions; behaves the same as in PLINK (requires --genome)
   --genome		Read in the specified PLINK formatted data file and calculate IBD estimates using PLINK
   --plink_ex		Path to the plink executable file (searches environment variables by default)
   --ref_pops		Comma separated list of 1000 Genomes populations used for reference allele freqs (overrides default method) 
   --no_automatic	Turn off automatic selection of the HapMap3 populations for reference allele freqs (On by default)
   --remove_AIMs	Automatically remove ancestry informative markers (off by default)
   --keep_AIMs		Do not remove ancestry informative markers(off by default)
   --internal_ref	Use the dataset provided in --file to get reference allele frequencies
   --alt_ref_stem	Path to PLINK formatted data (no file extensions) used for allele frequencies
   --keep_inter_files	Keep intermediate files used to create the IBD estimates with prePRIMUS
   --min_pihat_threshold Set a minimum pi-hat threshold that will be used in the plink --genome calculation
   --max_memory		Specify amount of memory to be used in PLINK prePRIMUS commands (in MB)

 Identification of maximum unrelated set (IMUS) options:
   --no_IMUS		Don't identify a maximum unrelated set (runs IMUS by default)
   --missing_val	Set value that denotes missing data in IBD file
   -s, --size		Specify to weight on set size (Default unless a binary trait is specified first)
   --high_btrait	File with FID, IID, and binary trait to weight for the higher value
   --low_btrait		File with FID, IID, and binary trait to weight for the lower value
   --high_qtrait	File with FID, IID, and quantitative trait to weight for higher values
   --low_qtrait		File with FID, IID, and quantitative trait to weight for lower values
   --mean_qtrait	File with FID, IID, and quantitative trait to weight towards the mean value
   --tails_qtrait	File with FID, IID, and quantitative trait to weight against the middle values

 Pedigree reconstruction (PR) options:
   --no_PR		Don't reconstruct pedigrees (runs pedigree reconstruction by default)
   --max_gens		Max number of generations sampled in reconstructed pedigree (default = no limit)
   --max_gen_gap	Max number of generations between two people that have a child (default = 0)
   --age_file		Specify path to the file containing the age of each sample
   --ages		Like --age_file but requires FILE=[file], optional specification of file columns
   --sex_file		Specify path to the file containing the sex of each sample
   --sexes		Like --sex_file but requires FILE=[file], optional specification of file columns
   --mito_matches	Path to mito matching status for each pair of samples (requires FILE=[file])
   --y_matches		Path to y matching status for each pair of samples (requires FILE=[file])
   --MT_error_rate	Proportion of the MT sequence that must not match to be called a non-match
   --Y_error_rate	Proportion of the Y sequence that must not match to be called a non-match
   --affection_file	Specify path to the file containing the affection status of each sample
   --affections		Like --affection_file; need FILE=[file], optional specification of file columns
   --int_likelihood_cutoff	Initial minimum likelihood for a relationship to reconstruction (default = 0.1)

 ERSA options:
   --min_cm			minimum segment size to consider [default: 2.5]
   --max_cm			maximum segment size to consider for estimating the exponential distribution of 
   				segment sizes in the population [default: 10.0]
   --max_meioses		maximum number of meioses to consider [default: 40]
   --rec_per_meioses		expected number of recombination events per meioses [default: 35.2548101]	
   --ascertained_chromosome	chromosome number of ascertained disease locus (int)
   --ascertained_position	chromosomal position of ascertained disease locus (int)
   --control_files		GERMLINE or Beagle fibd output file(s) for population control
   --control_sample_size	Sample size of control population. Used with --control_files
   --exp_mean			Mean of exp distribution of shared segment size in population [default: 3.197036753]
   --pois_mean			Mean of the Poisson distribution of the number of segments shared between 
				a pair of individuals in the population [default: 13.73]
   --pair_file			Restrict pairwise comparisons to the ID pairs specified in this file
   --single_pair		Restrict pairwise comparisons to the pairs specified in this flag (id1:id2)
   --number_of_ancestors	Restrict relationships to [1] one parent (half-sibs/cousins), 
   				[2] two parents (full-sibs/cousins), or [0] (parent-offspring/grandparent-granchild).
   --number_of_chromosomes	Number of chromosomes [default: 22]
   --parent_offspring_option	Option to evaluate potential parent-offspring and sibling relationships 
   				based on total proportion of the genome that is shared ibd1 [default: true]
   --parent_offspring_zscore	Z-score for rejecting a sibling relationship in favor of a parent-offspring 
   				relationship [default: 2.33, alpha=0.01] 
   --adjust_pop_dist		Option to adjust the population distribution of shared segments downward 
   				for segments that could not be detected due to recent ancestry [default: false]
   --confidence_level		Confidence level for confidence interval around the estimated degree 
   				of relationship. [default: 0.95]
   --mask_common_shared_regions	excludes chromosomal regions that are commonly shared from evaluation.
   				Used only when the control_files or mask_region_file parameter is specified [default: false].
   --mask_region_cross_length	length in base pairs that a shared segment must extend past a masked segment 
   				in order to avoid truncation. Used only when mask_common_shared_regions parameter is specified 
				[default: 1000000].
   --mask_region_file		file containing chromosomal regions to exclude from evaluation.
   --mask_region_threshold	Threshold for the ratio of observed vs. expected segment sharing in controls 
   				before a region will be masked.
   --mask_region_sim_count	number of simulations performed of the null distribution of shared 
   				segment locations in controls; results written to output_file.sim
   --recombination_files	file containing genetic distances for all chromosomes. This parameter must be 
   				specified with Beagle fibd input files
   --beagle_markers_files	Beagle marker files (one file required for each chromosome, wildcards required, 
   				ex: chr*beagle.marker). Each filename must begin with the chromosome name followed by a period.
				This parameter must be specified with Beagle fibd input files


