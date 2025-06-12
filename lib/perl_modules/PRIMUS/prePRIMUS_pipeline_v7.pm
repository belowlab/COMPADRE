package PRIMUS::prePRIMUS_pipeline_v7;

use strict;
use Getopt::Long qw(GetOptionsFromArray);
use PRIMUS::IMUS qw(run_IMUS);
use File::Path qw(make_path);
use IO::Socket::INET; 
use IPC::Open2; 
use Cwd qw(abs_path);
use File::Find;
use File::Basename;

###################################################################################

# Eventual TODO: update the rest of the plink commands to plink2
# Currently, only the PCA command uses plink2

###################################################################################

# sub send_to_compadre_helper {

#     my ($data, $port) = @_;
#     $port //= 6000;  
# 	my $host = $ENV{COMPADRE_HOST} // 'localhost';

#     my $socket = IO::Socket::INET->new(
#         PeerAddr => $host,
#         PeerPort => $port,
#         Proto    => 'tcp',
#     ) or die "Cannot connect to Python server: $!\n";

#     $socket->print($data);
    
#     my $response = <$socket>;
#     if (defined $response) {
#         chomp $response;
#     } else {
#         $response = "No response";
#     }
    
#     close($socket);
#     return $response;
# }

sub send_to_compadre_helper {
    my ($data, $port) = @_;
    $port //= 6000;  
    my $host = $ENV{COMPADRE_HOST} // 'localhost';

    my $socket = IO::Socket::INET->new(
        PeerAddr => $host,
        PeerPort => $port,
        Proto    => 'tcp',
        Timeout  => 0,  # Set to 0 for no timeout (indefinite)
    ) or die "Cannot connect to Python server: $!\n";

    # Enable keep-alive on the socket
    $socket->sockopt(SO_KEEPALIVE, 1) if $socket->can('sockopt');
    $socket->blocking(1);

    my $bytes_sent = $socket->print($data);
    if (!defined $bytes_sent) {
        close($socket);
        die "Failed to send data to Python server: $!\n";
    }
    
    $socket->flush();
    
    # Read response with error checking
    my $response = <$socket>;
    if (defined $response) {
        chomp $response;
    } else {
        $response = "No response";
        warn "No response received from Python server\n";
    }
    
    close($socket);
    return $response;
}

# Get absolute location of this file and use that to define additional paths of interest
my $preprimus_script_path = abs_path($0);
my $pmloc = dirname(dirname($preprimus_script_path)); 

# These match what is currently installed via the Dockerfile
my $PLINK = "plink";
my $PLINK2 = "plink2";
my $R = "R";

#my $HM3_STEM = "$lib_dir/hapmap3/allhapmapUNREL_r2_b36_fwd.qc.poly";
my $lib_dir = "../lib";
my $onekg = "$pmloc/lib/1KG";
#print "\n\n1KG directory: $onekg\n\n";

#$HM3_STEM = "$lib_dir/hapmap3/allhapmapUNREL_r2_b36_fwd.qc.poly";
#my $onekg_STEM = "$onekg/all_unrelateds_NEW";
my $onekg_STEM = "$onekg/1KG_reference"; # new

my $public_html_dir;
open(my $LOG,"");

## Global paramenters
my $verbose = 3;
my $plink_silent = "";
my $test = 0;
my $rerun = 0;
my $MIND = .1;
my $GENO = .05;
my $MAF = .1;
my $MIN_POP_LIKELIHOOD_CUTOFF = 0.3;
my $ASSOC_P_VAL = 0.00001;
my $NO_MITO = 0;
my $NO_Y = 0;
my $MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH = 0.01; ## similar to genotyping/sequencing error rate for MT
my $Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH = 0.01; ## similar to genotyping/sequencing error rate for Y
my $MAX_PERCENT_UNKNOWN_FOR_MATCH = 0.05;
my $MIN_SAMPLES_WITHOUT_REF = 50;
my $THRESHOLD = .1;
my $EXCLUDE_VALUE = 0;

my $memory_flag = "";

## Maybe make this local and pass them around
my $INDIVIDUAL_ANCESTRY = 0;
my $ADMIXED = 0;
my $INTERNAL_REF = 0;
my $ALTERNATIVE_REF = 0;
my %intermediate_files;

#########################
my @onekg_pops = ("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI");

my %onekg_colors = (
	"ASW"=>"steelblue4", 		# african
	"LWK"=>"steelblue3",		# african
	"YRI"=>"steelblue2",		# african
	"ACB"=>"steelblue1",		# african
	"ESN"=>"slateblue1",		# african
	"GWD"=>"slateblue2",		# african
	"MSL"=>"slateblue3",		# african
	"CEU"=>"gray1", 			# european
	"TSI"=>"gray2",				# european
	"FIN"=>"gray3",				# european	
	"GBR"=>"gray4",				# european	
	"IBS"=>"gray5",				# european
	"CHB"=>"hotpink",			# asian
	"CHS"=>"hotpink1",			# asian
	"GIH"=>"hotpink2",			# asian
	"JPT"=>"hotpink3",			# asian	
	"STU"=>"hotpink4",			# asian
	"BEB"=>"indianred",			# asian
	"CDX"=>"indianred1",		# asian
	"ITU"=>"indianred2",		# asian
	"KHV"=>"indianred3",		# asian
	"PJL"=>"indianred4",		# asian
	"MXL"=>"orange",			# hispanic
	"CLM"=>"orange1",			# hispanic
	"PEL"=>"orange2",			# hispanic
	"PUR"=>"orange3",			# hispanic
);


## Need to change the test_paths() module to be correct, then uncomment next line
#test_paths();
#set_paths();

######
sub run_prePRIMUS_main {

	my $study_name;
	my $genome_file;
	my $sex_file = "";
	my $mt_file = "";
	my $y_file = "";
	my $output_dir = "./";
	my $data_stem;
	my $ped_file;
	my $map_file;
	my $rerun;
	my $alt_ref_stem = "";
	my $remove_AIMs = 0;
	my $keep_AIMs = 0;
	my $no_automatic_IBD = 0;
	my @ref_pops;
	my %ibd_estimates = ("FID1",1,"IID1",2,"FID2",3,"IID2",4,"IBD0",7,"IBD1",8,"IBD2",9,"PI_HAT",10);
	my %MT_estimates = ("FID1",1,"IID1",2,"FID2",3,"IID2",4,"MATCH",5,"MATCH_VAL",1);
	my %Y_estimates = ("FID1",1,"IID1",2,"FID2",3,"IID2",4,"MATCH",5,"MATCH_VAL",1);
	my $no_PCA_plot = 0;
	my $IBD0_vs_IBD1_plot;
	my $keep_intermediate_files = 1;
	my $max_memory = 0;

	GetOptionsFromArray(
		\@_,
		# Diagnostic options
		'verbose=i' => \$verbose,
		'test=i' => \$test,
		
		# General settings
		"study_name=s" => \$study_name, 
		"output_dir=s" => \$output_dir,
		"lib=s"=>\$lib_dir,
		"log_file_handle=s"=>\$LOG,

		# PLINK IBD PIPELINE settings
		"file=s" => \$data_stem,		## behaves the same as bfile
		"bfile=s" => \$data_stem,
		"rerun=i" => \$rerun,                 	## currently has no affect
		"internal_ref=i" => \$INTERNAL_REF,
		"alt_ref=s" => sub
		{
			$alt_ref_stem = $_[1];
			$ALTERNATIVE_REF = 1 if $alt_ref_stem ne "";
		},
		"ref_pops_ref=s" => sub
		{
			@ref_pops = @{$_[1]};
		},
		"remove_AIMs=i" => \$remove_AIMs,
		"keep_AIMs=i" => \$keep_AIMs,
		"no_automatic_IBD=i" => \$no_automatic_IBD,
		"no_PCA_plot=i" => \$no_PCA_plot,
		"keep_intermediate_files=i" => \$keep_intermediate_files,
		"MT_error_rate=f" => \$MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH,
		"Y_error_rate=f" => \$Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH,
		"no_mito=i" => \$NO_MITO,
		"no_y=i" => \$NO_Y,
		"max_memory=i" => \$max_memory,

		# IMUS Settings
		"rel_threshold=f" => \$THRESHOLD, 

	) or die "Failed to parse options for Pedigree Reconstruction\n";

	# Build PLINK --memory flag string if max_memory is greater than 0
	if ($max_memory > 0) {
		$memory_flag = "--memory $max_memory";
	} else {
		$memory_flag = "";
	}
	
	
	#### Check/set inputs
	test_paths();

	if($alt_ref_stem ne "" && @ref_pops > 0){die "can't use --alt_ref and --ref_pops_ref at the same time\n";}
	if($remove_AIMs == 1 && $keep_AIMs == 1){die "can't use --remove_AIMs and --keep_AIMs at the same time\n";}
	$study_name = get_file_name_from_stem($data_stem) if $study_name eq "";
	$plink_silent = "--silent" if $verbose < 2;
	make_path($output_dir) if !-d $output_dir;
	open($LOG,">$output_dir/$study_name.log") if($LOG eq "");

	print "\n\nUsing PLINK to calculate IBD estimates for $study_name\n" if($verbose > 0);
	print $LOG "\n\nUsing PLINK to calculate IBD estimates for $study_name\n" if($verbose > 0);
	print "Study_name: $study_name\n" if $test;
	print $LOG "Study_name: $study_name\n" if $test;

	$data_stem = make_binary_version($data_stem,"$output_dir/$study_name");
	
	## Set PID and MID to 0 so it doesn't mess up the founder allele frequency stuff 
        system("cp $data_stem.fam $data_stem.fam_temp"); 
        open(IN,"$data_stem.fam_temp"); 
        open(OUT,">$data_stem.fam"); 
        while(my $line = <IN>) 
        {       
                my ($FID, $IID, $PID,$MID,$SEX,$AFF) = split(/\s+/,$line); 
                print OUT "$FID\t$IID\t0\t0\t$SEX\t$AFF\n"; 
        } 
        close(OUT); 
        close(IN); 
	
	$mt_file = get_MT_estimates($data_stem) if !$NO_MITO;
	$y_file = get_Y_estimates($data_stem) if !$NO_Y;
	
	#exit;

	#### use internal_ref
	if($INTERNAL_REF)
	{
		print "\nRunning $data_stem using internal allele frequencies\n" if $verbose > 0;
		my ($no_dup_stem, @dup_SNPs) = remove_dups($data_stem);
		my $autosomal_no_dup_stem = remove_non_autosomal_SNPs($no_dup_stem);

		## Use unrelated samples
		my ($unrelated_stem, @unrelated_samples) = get_unrelateds($autosomal_no_dup_stem);
		my $PCA_plot = make_PCA_plot($autosomal_no_dup_stem,"",$study_name) if !$no_PCA_plot;
		my($AIMs_file,@AIMs) = get_AIMs($autosomal_no_dup_stem) if $remove_AIMs == 1;
		my $allele_freqs = get_allele_freqs($unrelated_stem);
		
		## Use all samples
		my $cleaned_all_samples_stem = remove_SNPs($no_dup_stem,\@AIMs,"$data_stem\_cleaned");
		$genome_file = calculate_IBD_estimates($cleaned_all_samples_stem,"",$allele_freqs);
		$IBD0_vs_IBD1_plot = make_IBD0_vs_IBD1_plot($genome_file,$study_name);
		call_callrate($cleaned_all_samples_stem);
		call_het($cleaned_all_samples_stem,$allele_freqs);
		$sex_file = call_sex($cleaned_all_samples_stem,$allele_freqs);
	}
	elsif($ALTERNATIVE_REF)
	{
		print "\nRunning $data_stem using $alt_ref_stem as a reference\n" if $verbose > 0;
		print $LOG "\nRunning $data_stem using $alt_ref_stem as a reference\n" if $verbose > 0;
		my ($no_dup_stem, @dup_SNPs) = remove_dups($data_stem);
		my $autosomal_no_dup_stem = remove_non_autosomal_SNPs($no_dup_stem);
		
		## Use unrelated samples
		my ($unrelated_stem, @unrelated_samples) = get_unrelateds($autosomal_no_dup_stem);
		my($merged_stem,$flipped_SNP_arr_ref,$remove_SNP_arr_ref) = merge($autosomal_no_dup_stem,$alt_ref_stem);
		my $PCA_plot = make_PCA_plot($merged_stem,$alt_ref_stem,$study_name) if !$no_PCA_plot;
		my($AIMs_file,@AIMs) = get_AIMs($merged_stem) if $remove_AIMs == 1;
		my $allele_freqs = get_allele_freqs($unrelated_stem);
		
		#### Use all non-reference samples
		## Remove the SNPs that were remove for the merge
		my $cleaned_all_samples_stem = remove_SNPs($no_dup_stem,$remove_SNP_arr_ref,"$data_stem\_cleaned");
		## Remove the AIMs
		$cleaned_all_samples_stem = remove_SNPs($cleaned_all_samples_stem,\@AIMs,$cleaned_all_samples_stem) if @AIMs > 0;
		## Flip the SNPs that were flipped for the merge
		#$cleaned_all_samples_stem = flip_SNPs($cleaned_all_samples_stem,$cleaned_all_samples_stem,@$flipped_SNP_arr_ref);
		## Remove homozyous SNPs so the .frq file from the merge data will work
		#$cleaned_all_samples_stem = remove_homozygous_SNPs($cleaned_all_samples_stem,$cleaned_all_samples_stem); ## This is necessary until --read-freq will allow homozygosity in the dataset (currently, 0/A in data fails if C/A in .frq);
		$genome_file = calculate_IBD_estimates($cleaned_all_samples_stem,"",$allele_freqs);
		$IBD0_vs_IBD1_plot = make_IBD0_vs_IBD1_plot($genome_file,$study_name);
		call_callrate($cleaned_all_samples_stem);
		call_het($cleaned_all_samples_stem,$allele_freqs);
		$sex_file = call_sex($cleaned_all_samples_stem,$allele_freqs);
	}
	else ## Use 1KG as a reference
	{
		print "\nRunning $data_stem using 1KG to find an appropriate reference\n" if $verbose > 0;
		print $LOG "\nRunning $data_stem using 1KG to find an appropriate reference\n" if $verbose > 0;
		my ($no_dup_stem, @dup_SNPs) = remove_dups($data_stem);
		my $autosomal_no_dup_stem = remove_non_autosomal_SNPs($no_dup_stem);
		
		my ($unrelated_stem, @unrelated_samples) = get_unrelateds($autosomal_no_dup_stem);
		my($merged_stem,$flipped_SNP_arr_ref,$remove_SNP_arr_ref);
		if(@ref_pops < 1) ## Select reference populations if not specified on commandline
		{
			($merged_stem,$flipped_SNP_arr_ref,$remove_SNP_arr_ref) = merge($autosomal_no_dup_stem,$onekg_STEM);
			#my $PCA_plot = make_PCA_plot($merged_stem,$onekg_STEM,$study_name,"onekg",1) if !$no_PCA_plot;
			#if($no_automatic_IBD){die "--no_automatic_IBD was selected. Please see $PCA_plot to select a reference population and rerun\n";}
			@ref_pops = pick_reference_populations($merged_stem,$study_name,1);
		}

		# Now we have ref_pops from pop_classifier code
		

		#### Proceeds down one of 4 paths depending on whether it is admixed and the number of unrelated samples
		my $allele_freqs;
		my $cleaned_all_samples_stem;

		if((is_admixed(@ref_pops) || $remove_AIMs ) && @unrelated_samples >= $MIN_SAMPLES_WITHOUT_REF && !$keep_AIMs)
		{
			print "\nREMOVING AIMS and >= $MIN_SAMPLES_WITHOUT_REF UNRELATED SAMPLES\n" if $verbose > 0;
			print $LOG "\nREMOVING AIMS and >= $MIN_SAMPLES_WITHOUT_REF UNRELATED SAMPLES\n" if $verbose > 0;
			my $PCA_plot = make_PCA_plot($autosomal_no_dup_stem,"",$study_name) if !$no_PCA_plot;
			my($AIMs_file,@AIMs) = get_AIMs($autosomal_no_dup_stem);
			$allele_freqs = get_allele_freqs($unrelated_stem);
			$cleaned_all_samples_stem = remove_SNPs($no_dup_stem,$remove_SNP_arr_ref,"$data_stem\_cleaned");
			$cleaned_all_samples_stem = remove_SNPs($cleaned_all_samples_stem,\@AIMs,$cleaned_all_samples_stem) if @AIMs > 0;
		}

		elsif((is_admixed(@ref_pops) || $remove_AIMs) && @unrelated_samples < $MIN_SAMPLES_WITHOUT_REF && !$keep_AIMs)
		{
			print "\nREMOVING AIMS and < $MIN_SAMPLES_WITHOUT_REF UNRELATED SAMPLES\n" if $verbose > 0;
			print $LOG "\nREMOVING AIMS and < $MIN_SAMPLES_WITHOUT_REF UNRELATED SAMPLES\n" if $verbose > 0;
			my ($merged_stem,$ref_stem,$flipped_SNP_arr_ref,$remove_SNP_arr_ref) = multiple_ref_pop_merge($unrelated_stem,"",@ref_pops);
			my $PCA_plot = make_PCA_plot($merged_stem,$ref_stem,$study_name) if !$no_PCA_plot;
			my($AIMs_file,@AIMs) = get_AIMs($merged_stem);
			my $merged_stem_no_AIMs = remove_SNPs($merged_stem,\@AIMs,"$merged_stem\_no_AIMs") if @AIMs > 0;
			$allele_freqs = get_allele_freqs($merged_stem);
			$cleaned_all_samples_stem = remove_SNPs($no_dup_stem,$remove_SNP_arr_ref,"$data_stem\_cleaned");
			$cleaned_all_samples_stem = remove_SNPs($cleaned_all_samples_stem,\@AIMs,$cleaned_all_samples_stem) if @AIMs > 0;
			## Flip cleaned_all_samples_stem so it matches the allele_freqs file
			$cleaned_all_samples_stem = flip_SNPs($cleaned_all_samples_stem,$cleaned_all_samples_stem,@$flipped_SNP_arr_ref);
			#$cleaned_all_samples_stem = remove_homozygous_SNPs($cleaned_all_samples_stem,$cleaned_all_samples_stem); ## This is necessary until --read-freq will allow homozygosity in the dataset (currently, 0/A in data fails if C/A in .frq);
		}

		elsif((!is_admixed(@ref_pops) || $keep_AIMs) && @unrelated_samples >= $MIN_SAMPLES_WITHOUT_REF)
		{
			print "\nNOT REMOVING AIMS and >= $MIN_SAMPLES_WITHOUT_REF UNRELATED SAMPLES\n" if $verbose > 0;
			print $LOG "\nNOT REMOVING AIMS and >= $MIN_SAMPLES_WITHOUT_REF UNRELATED SAMPLES\n" if $verbose > 0;
			my $PCA_plot = make_PCA_plot($unrelated_stem,"",$study_name) if !$no_PCA_plot;
			$allele_freqs = get_allele_freqs($unrelated_stem);
			$cleaned_all_samples_stem = remove_SNPs($no_dup_stem,$remove_SNP_arr_ref,"$data_stem\_cleaned");
		}

		elsif((!is_admixed(@ref_pops) || $keep_AIMs) && @unrelated_samples < $MIN_SAMPLES_WITHOUT_REF)
		{
			print "\nNOT REMOVING AIMS and < $MIN_SAMPLES_WITHOUT_REF UNRELATED SAMPLES\n" if $verbose > 0;
			print $LOG "\nNOT REMOVING AIMS and < $MIN_SAMPLES_WITHOUT_REF UNRELATED SAMPLES\n" if $verbose > 0;
			my ($merged_stem,$ref_stem,$flipped_SNP_arr_ref,$remove_SNP_arr_ref) = multiple_ref_pop_merge($unrelated_stem,"",@ref_pops);
			my $PCA_plot = make_PCA_plot($merged_stem,$ref_stem,$study_name) if !$no_PCA_plot;
			$allele_freqs = get_allele_freqs($merged_stem);
			
			#### Use all non-reference, clean data
			$cleaned_all_samples_stem = remove_SNPs($no_dup_stem,$remove_SNP_arr_ref,"$data_stem\_cleaned");
			$cleaned_all_samples_stem = flip_SNPs($cleaned_all_samples_stem,$cleaned_all_samples_stem,@$flipped_SNP_arr_ref);
			#$cleaned_all_samples_stem = remove_homozygous_SNPs($cleaned_all_samples_stem,$cleaned_all_samples_stem); ## This is necessary until --read-freq will allow homozygosity in the dataset (currently, 0/A in data fails if C/A in .frq);
		}

		else
		{
			die "INVALID INPUT OPTIONS\n";
		}

		# calculate_IBD_estimates subroutine must take the stem with the reference people merged in 

		$genome_file = calculate_IBD_estimates($cleaned_all_samples_stem,"",$allele_freqs);
		$IBD0_vs_IBD1_plot = make_IBD0_vs_IBD1_plot($genome_file,$study_name);


        #call_callrate($cleaned_all_samples_stem);
        #call_het($cleaned_all_samples_stem,$allele_freqs);
        #$sex_file = call_sex($cleaned_all_samples_stem,$allele_freqs);
	}

	remove_intermediate_files(11,$output_dir) if !$keep_intermediate_files;
	print "IBD estimates are in $genome_file\n" if $verbose > 0;
	print $LOG "IBD estimates are in $genome_file\n" if $verbose > 0;
	print "IBD0 vs IBD1 plot: $IBD0_vs_IBD1_plot\n" if $verbose > 0;
	print $LOG "IBD0 vs IBD1 plot: $IBD0_vs_IBD1_plot\n" if $verbose > 0;
	print "\n\nPREPRIMUS 1KG POPCLASSIFIER VERSION DONE!\n\n";
	return ($genome_file,$sex_file,$mt_file,$y_file);
}

#############################################################################
############## EXTERNALLY USEFUL METHODS ####################################
#############################################################################

sub get_MT_estimates
{
	my $data_stem = shift;
	my $new_stem = shift;
	$new_stem = "$data_stem\_MT_estimates" if $new_stem eq "";
	
	make_binary_version($data_stem);
	print "\nCalculating MT estimates from $data_stem => $new_stem.txt\n" if $verbose > 0;
	print $LOG "\nCalculating MT estimates from $data_stem => $new_stem.txt\n" if $verbose > 0;
	
	my $temp = system("$PLINK --noweb --bfile $data_stem --recode --mind 0.1 --geno 0.1 --chr 26 $plink_silent --out $data_stem\_chr26 $memory_flag");
	$intermediate_files{"$data_stem\_chr26"} = 5;

	return "" if !-e "$data_stem\_chr26.ped";

	my %data;
	open(IN,"$data_stem\_chr26.ped");
	while(my $line = <IN>)
	{
		my ($fid,$iid,$pid,$mid,$sex,$aff,@snps) = split(/\s+/,$line);
		$data{"$fid\__$iid"} = \@snps;
	}
	close(IN);

	## Calculate matches for each pair and write to file
	open(OUT,">$new_stem.txt");
	print OUT "FID1 IID1 FID2 IID2 MT_MATCH NUM_DIFF NUM_UNKNOWN LENGTH\n";
	foreach my $name1 (keys %data)
	{
		foreach my $name2 (keys %data)
		{
			last if $name1 =~ /^$name2$/;
			my ($match,$num_diff,$num_unknown,$length) = do_MT_sequences_match($data{$name1},$data{$name2});
			my ($fid1,$iid1) = split(/__/,$name1);
			my ($fid2,$iid2) = split(/__/,$name2);
            next if $match > 0;
			print OUT "$fid1 $iid1 $fid2 $iid2 $match $num_diff $num_unknown $length\n";
		}
	}
	close(OUT);
	return "$new_stem.txt";
}

sub get_Y_estimates
{
	my $data_stem = shift;
	my $new_stem = shift;
	$new_stem = "$data_stem\_Y_estimates" if $new_stem eq "";
	my %MT_estimates;
	
	make_binary_version($data_stem);
	print "\nCalculating Y estimates from $data_stem => $new_stem.txt\n" if $verbose > 0;
	print $LOG "\nCalculating Y estimates from $data_stem => $new_stem.txt\n" if $verbose > 0;
	
	## NEED TO REMOVE SNPs UP TO 2.65M BP, BECAUSE THAT IS THE PSEUDOAUTOSOMAL REGION
	#my $temp = system("$PLINK --noweb --bfile $data_stem --recode --mind 0.05 --geno 0.05 --chr 24 --from-bp 2650000 --to-bp 60000000 $plink_silent --out $data_stem\_chr24");
	my $temp = system("$PLINK --noweb --bfile $data_stem --recode --mind 0.05 --chr 24 --from-bp 2650000 --to-bp 60000000 $plink_silent --out $data_stem\_chr24 $memory_flag");
	$intermediate_files{"$data_stem\_chr24"} = 5;
	return "" if !-e "$data_stem\_chr24.ped";

	system("perl -pi -e 's/^24/22/g' $data_stem\_chr24.map");

	my $temp = system("$PLINK --noweb --file $data_stem\_chr24 --recode --geno 0.05 $plink_silent --out $data_stem\_chr2 $memory_flag");
	
	system("perl -pi -e 's/^22/24/g' $data_stem\_chr24.map");

	my %data;
	open(IN,"$data_stem\_chr24.ped");
	while(my $line = <IN>)
	{
		my ($fid,$iid,$pid,$mid,$sex,$aff,@snps) = split(/\s+/,$line);
		$data{"$fid\__$iid"} = \@snps;
	}
	close(IN);

	## Calculate matches for each pair and write to file
	open(OUT,">$new_stem.txt");
	print OUT "FID1 IID1 FID2 IID2 Y_MATCH NUM_DIFF NUM_UNKOWN LENGTH\n";
	foreach my $name1 (keys %data)
	{
		foreach my $name2 (keys %data)
		{
			last if $name1 =~ /^$name2$/;
			my ($match,$num_diff,$num_unknown,$length) = do_Y_sequences_match($data{$name1},$data{$name2});
			my ($fid1,$iid1) = split(/__/,$name1);
			my ($fid2,$iid2) = split(/__/,$name2);
            next if $match > 0;
			print OUT "$fid1 $iid1 $fid2 $iid2 $match $num_diff $num_unknown $length\n";
		}
	}
	close(OUT);
	return "$new_stem.txt";
}

sub do_MT_sequences_match
{
	my $seq1 = shift;
	my $seq2 = shift;
	
	my $num_differences = 0;
	my $num_matches = 0;
	my $num_unknown = 0;
	my $num_SNPs = @{$seq1}/2;
	for(my $pos = 0; $pos < $num_SNPs; $pos += 2) ## skip 2 each iteration because MT is haploid but is represented as diploid
	{
		my $snp1 = @{$seq1}[$pos];
		my $snp2 = @{$seq2}[$pos];
		if($snp1 =~ /[0Nn]/ || $snp2 =~ /[0Nn]/)
		{
			$num_unknown++;
		}
		else
		{
			$num_differences++ if $snp1 ne $snp2;
			$num_matches++ if $snp1 eq $snp2;
		}
	}
	#print "num_diff: $num_differences\n";

	my $percent_diff = $num_differences/($num_SNPs-$num_unknown); 

	if($num_unknown/$num_SNPs > $MAX_PERCENT_UNKNOWN_FOR_MATCH)
	{
		return (-1,$num_differences,$num_unknown,$num_SNPs);
	}
	elsif($percent_diff <= $MT_MAX_PERCENT_DIFFERENCE_FOR_MATCH)
	{
		return (1,$num_differences,$num_unknown,$num_SNPs);
	}
	else
	{
		return (0,$num_differences,$num_unknown,$num_SNPs);
	}

}

sub do_Y_sequences_match
{
	my $seq1 = shift;
	my $seq2 = shift;
	
	my $num_differences = 0;
	my $num_matches = 0;
	my $num_unknown = 0;
	my $num_SNPs = @{$seq1}/2;
	for(my $pos = 0; $pos < $num_SNPs; $pos += 2) ## skip 2 each iteration because Y is haploid but is represented as diploid
	{
		my $snp1 = @{$seq1}[$pos];
		my $snp2 = @{$seq2}[$pos];
		if($snp1 =~ /[0Nn]/ || $snp2 =~ /[0Nn]/)
		{
			$num_unknown++;
		}
		else
		{
			$num_differences++ if $snp1 ne $snp2;
			$num_matches++ if $snp1 eq $snp2;
		}
	}
	#print "num_diff: $num_differences\n";

	my $percent_diff = $num_differences/($num_SNPs-$num_unknown); 

	if($num_unknown/$num_SNPs > $MAX_PERCENT_UNKNOWN_FOR_MATCH)
	{
		return (-1,$num_differences,$num_unknown,$num_SNPs);
	}
	elsif($percent_diff <= $Y_MAX_PERCENT_DIFFERENCE_FOR_MATCH)
	{
		return (1,$num_differences,$num_unknown,$num_SNPs);
	}
	else
	{
		return (0,$num_differences,$num_unknown,$num_SNPs);
	}

}

sub remove_homozygous_SNPs
{
	my $data_stem = shift;
	my $new_stem = shift;
	$new_stem = "$data_stem\_non_homozygous" if $new_stem eq "";
	
	print "\nRemoving homozugous SNPs from $data_stem => $new_stem\n" if $verbose > 0;
	print $LOG "\nRemoving homozugous SNPs from $data_stem => $new_stem\n" if $verbose > 0;

	my $temp = system("$PLINK --noweb --bfile $new_stem --maf 0.001 --make-bed --indiv-sort 0 $plink_silent --out $new_stem $memory_flag");
	$intermediate_files{"$new_stem.bed"} = 5;
	$intermediate_files{"$new_stem.bim"} = 5;
	$intermediate_files{"$new_stem.fam"} = 5;
	$intermediate_files{"$new_stem.log"} = 5;

	return $new_stem;
}

sub make_IBD0_vs_IBD1_plot
{
	my $genome_file = shift;
	my $study_name = shift;

	print "\nMaking IBD0 vs IBD1 plot for $genome_file\n" if $verbose > 1;
	print $LOG "\nMaking IBD0 vs IBD1 plot for $genome_file\n" if $verbose > 1;

	## make R script
	open(R,">$genome_file\_IBD0_vs_IBD1.R");
	print R "data <-read.table(\"$genome_file\",header=T)\n";
	print R "jpeg(\"$genome_file\_IBD0_vs_IBD1.jpeg\", height=480, width=480)\n"; 
	print R "plot(data[,'Z0'],data[,'Z1'],xlab=\"IBD0\",xlim=c(0,1),ylim=c(0,1),ylab=\"IBD1\",main=\"IBD0 vs IBD1 for $study_name\")\n";
	print R "dev.off()\n";
	close(R);
	$intermediate_files{"$genome_file\_IBD0_vs_IBD1.R"} = 5;

	## Run R script
	my $temp = system("R --vanilla --slave < $genome_file\_IBD0_vs_IBD1.R > $genome_file\_IBD0_vs_IBD1.R_output");
	system("rm $genome_file\_IBD0_vs_IBD1.R_output");

	return "$genome_file\_IBD0_vs_IBD1.jpeg";
}

sub remove_non_autosomal_SNPs
{
	my $data_stem = shift; ## Your data that you want to flip
	my $new_stem = shift; ## The stem name of your filed data. If left blank, the name will be $data_stem\_flipped.
	$new_stem = "$data_stem\_autosomal" if $new_stem eq "";
	
	make_binary_version($data_stem);
	$intermediate_files{"$new_stem.bed"} = 5;
	$intermediate_files{"$new_stem.bim"} = 5;
	$intermediate_files{"$new_stem.fam"} = 5;

	## Get non-autosomal SNPs from bim file
	my @non_autosomal_SNPs;
	open(IN,"$data_stem.bim") or die "can't open $data_stem.bim; $!\n";
	while(my $line = <IN>)
	{
		my ($chr,$name) = split(/\s+/,$line);
		push(@non_autosomal_SNPs,$name) if($chr !~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)$/);
	}
	close(IN);
	close(OUT);
	
	remove_SNPs($data_stem,\@non_autosomal_SNPs,$new_stem);
	
	return $new_stem;
}

sub is_admixed
{
    my @ref_pops = @_;
    my $is_admixed = 0;

    print "Checking if @ref_pops are admixed\n" if $verbose > 1;
    print $LOG "Checking if @ref_pops are admixed\n" if $verbose > 1;
    
    # Known admixed populations
    # ASW - African Ancestry in Southwest US (African + European admixture)
    # ACB - African Caribbean (African + European admixture)
    # MXL - Mexican Ancestry in Los Angeles (Native American + European + African admixture)
    # CLM - Colombian in Medellin (Native American + European + African admixture)
    # PEL - Peruvian in Lima (Native American + European admixture)
    # PUR - Puerto Rican in Puerto Rico (European + African + Native American admixture)
    # GIH - Gujarati Indian in Houston (complex South Asian genetics)
    # BEB - Bengali in Bangladesh (complex South Asian genetics)
    # ITU - Indian Telugu in the UK (complex South Asian genetics)
    # PJL - Punjabi in Lahore, Pakistan (complex South Asian genetics)
    # STU - Sri Lankan Tamil in the UK (complex South Asian genetics)

	# Automatically set to admixed TRUE if it's one of the default admixed pops 

    if(grep {$_ eq "ASW" || $_ eq "ACB" || $_ eq "MXL" || $_ eq "CLM" || $_ eq "PEL" || $_ eq "PUR" || 
            $_ eq "GIH" || $_ eq "BEB" || $_ eq "ITU" || $_ eq "PJL" || $_ eq "STU"} @ref_pops) {
        $is_admixed = 1;
    }
    
    my $eur = 0;
    my $asn = 0;
    my $afr = 0;
    my $sas = 0;

    # European populations (EUR)
    if(grep {$_ eq "CEU" || $_ eq "TSI" || $_ eq "FIN" || $_ eq "GBR" || $_ eq "IBS"} @ref_pops) {
        $eur = 1;
    }
    
    # African populations (AFR)
    if(grep {$_ eq "YRI" || $_ eq "LWK" || $_ eq "GWD" || $_ eq "MSL" || $_ eq "ESN"} @ref_pops) {
        $afr = 1;
    }
    
    # East Asian populations (EAS)
    if(grep {$_ eq "CHB" || $_ eq "CHS" || $_ eq "JPT" || $_ eq "CDX" || $_ eq "KHV"} @ref_pops) {
        $asn = 1;
    }
    
	# All the SAS subpops are de facto admixed so they're handled above for now

    # If more than one continental ancestry is involved, consider it admixed
    if(($eur + $asn + $afr + $sas) > 1) {
        $is_admixed = 1;
    }

    $ADMIXED = $is_admixed;
    return $is_admixed;
}

sub pick_reference_populations {
	
	my $data_stem = shift;
	my $study_name = shift;
	my $rerun_pca= shift;

	print "\nSelecting reference population(s) for $data_stem\n" if $verbose > 0;
	print $LOG "\nSelecting reference population(s) for $data_stem\n" if $verbose > 0;
	$study_name = get_file_name_from_stem($data_stem) if $study_name eq "";

	#if (!-e "$data_stem.eigenvec" || $rerun_pca == 1) { print "\nNeed to run PCA\n"; }
	run_pca($data_stem, 1) if (!-e "$data_stem.eigenvec" || $rerun_pca == 1);
	
	my %samples; 
	my %KDEs;
	my @PC1;
	my @PC2;
	my @ref_pops;

	my $port_number = $main::port_number_glob;
	my $onekg_idfile = "$onekg/1KG_pop_classifier_ids.txt";
	my $socket_data = "$data_stem.eigenvec|$onekg_idfile|pop_classifier";

	# Run population classifier and return the top populations as @ref_pops
	my $ref_pops_str = send_to_compadre_helper($socket_data, $port_number);

	if($ref_pops_str eq 'No response')
	{
		die "Error. Run population classifier script in isolation with PLINK's PCA *.eigenvec file for more information.  You can also manually select reference sub-populations (ex., CEU) and rerun COMPADRE with the --ref_pops [POP] option\n";
	}

	@ref_pops = split(/\|/, $ref_pops_str);

	print "Ref pops: @ref_pops\n" if $verbose > 0;
	print $LOG "Ref pops: @ref_pops\n" if $verbose > 0;


	return @ref_pops;
}

sub multiple_ref_pop_merge
{
	my $data_stem = shift;
	my $new_stem = shift;
	my @ref_pops = @_;
	
	$new_stem = "$data_stem\_multi_merged" if $new_stem eq "";
	
	print "\nMerging $data_stem with @ref_pops => $new_stem\n" if $verbose > 0;
	print $LOG "\nMerging $data_stem with @ref_pops => $new_stem\n" if $verbose > 0;
	
	my $ref_samples_to_keep_stem = "$new_stem\_1KG_reference_samples";
    #my $files = "";
    #foreach my $pop (@ref_pops)
    #{
    #	$files .= "$lib_dir/hapmap3/hapmap3_r2_b36_fwd.$pop.qc.poly.genome_maximum_independent_set ";
    #}
    my $grepping_term = join('|',@ref_pops); 
    system("grep -E \"$grepping_term\" $onekg_STEM.fam > $ref_samples_to_keep_stem.txt");
    #system("cat $files > $ref_samples_to_keep_stem.txt");
	$intermediate_files{"$ref_samples_to_keep_stem.txt"} = 1;

	keep_samples($onekg_STEM,"$new_stem\_1KG_reference_samples","$ref_samples_to_keep_stem.txt");
	$intermediate_files{"$new_stem\_1KG_reference_samples.bed"} = 10;
	$intermediate_files{"$new_stem\_1KG_reference_samples.bim"} = 10;
	$intermediate_files{"$new_stem\_1KG_reference_samples.fam"} = 10;
	
	my ($new_stem,$flipped_SNP_arr_ref,$remove_SNP_arr_ref) = merge($data_stem,$ref_samples_to_keep_stem,$new_stem);
	return ($new_stem,$ref_samples_to_keep_stem,$flipped_SNP_arr_ref,$remove_SNP_arr_ref);
}

sub merge {

	my $data_stem = shift;
	my $ref_stem = shift;
	my $new_stem = shift;
	$new_stem = "$data_stem\_merged" if $new_stem eq "";
	
	print "\nMerging $data_stem and $ref_stem => $new_stem\n" if $verbose > 0;
	print $LOG "\nMerging $data_stem and $ref_stem => $new_stem\n" if $verbose > 0;
	
	make_binary_version($data_stem);
	make_binary_version($ref_stem);

	## Remove Dups
	remove_dups($data_stem,$data_stem); 
	
	## Flip $data_stem to match $ref_stem
	my($flipped_stem,$flipped_SNP_arr_ref,$remove_SNP_arr_ref) = flip($data_stem,$ref_stem);
	
	## Merge the $data_stem and $ref_stem
	my $temp = system("$PLINK --allow-no-sex --noweb --bfile $flipped_stem --bmerge $ref_stem --extract $flipped_stem.bim --geno $GENO --maf $MAF --make-bed --indiv-sort 0 $plink_silent --out $new_stem $memory_flag");
	if($temp > 0){die "ERROR!!! Plink failed to merge ped files; Check $new_stem.log\n";}
	
	## Finish selecting the intersection of SNPs between the two
	my $temp = system("$PLINK --allow-no-sex --noweb --bfile $new_stem --extract $ref_stem.bim --make-bed --indiv-sort 0 $plink_silent --out $new_stem $memory_flag");
	$intermediate_files{"$new_stem.bed"} = 10;
	$intermediate_files{"$new_stem.bim"} = 10;
	$intermediate_files{"$new_stem.fam"} = 10;

	return ($new_stem,$flipped_SNP_arr_ref,$remove_SNP_arr_ref);
}


## Flips SNPs to match a reference file and removes bad and ambiguous SNPs
## Reads in your data stem, the reference data stem, and an optional output stem. 
## Will also generate several other files that include ambigious SNPs, bad SNPs (triallelics), and a list of those that need to be flipped.
## It will then flip the SNPs that need to be flipped and remove the ambigious and bad SNPs write out the results to the $new_stem_name.
sub flip
{
	my $data_stem = shift; ## Your data that you want to flip
	my $ref_stem = shift; ## The reference data you want your data to be flipped to match 
	my $new_stem = shift; ## The stem name of your filed data. If left blank, the name will be $data_stem\_flipped.
	$new_stem = "$data_stem\_flipped" if $new_stem eq "";
	$ref_stem = "$onekg_STEM" if $ref_stem eq "";
	
	print "\nFlipping $data_stem to match $ref_stem => $new_stem\n" if $verbose > 0;
	print $LOG "\nFlipping $data_stem to match $ref_stem => $new_stem\n" if $verbose > 0;
	
	make_binary_version($data_stem);
	make_binary_version($ref_stem);
	
	## Load data
	my %data;
	my %ref_data;
	my $data_has_rsIDs = 0;
	my $ref_data_has_rsIDs = 0;
	my $data_genotypes_numeric = 0;
	my $ref_data_genotypes_numeric = 0;
	
	open(IN,"$data_stem.bim");
	while(my $line = <IN>)
	{
		my ($chr,$snpID,$genetic_pos,$pos,$a1,$a2) = split(/\s+/,$line);
		#$data{"$chr\__$pos"}{$a1}=1;
		#$data{"$chr\__$pos"}{$a2}=1;
		$data{$snpID}{$a1}=1 if $a1 ne 0;
		$data{$snpID}{$a2}=1 if $a2 ne 0;
		$data_has_rsIDs = 1 if $snpID =~ /^rs/;
		$data_genotypes_numeric = 1 if $a1 =~ /[1|2]/;
	}
	close(IN);
	open(IN,"$ref_stem.bim");
	while(my $line = <IN>)
	{
		my ($chr,$snpID,$genetic_pos,$pos,$a1,$a2) = split(/\s+/,$line);
		#$ref_data{"$chr\__$pos"}{$a1}=1;
		#$ref_data{"$chr\__$pos"}{$a2}=1;
		$ref_data{$snpID}{$a1}=1 if $a1 ne 0;
		$ref_data{$snpID}{$a2}=1 if $a2 ne 0;
		$ref_data_has_rsIDs = 1 if $snpID =~ /^rs/;
		$ref_data_genotypes_numeric = 1 if $a1 =~ /[1|2]/;
	}
	close(IN);

	## Check compatability between two dataset
	if($data_has_rsIDs != $ref_data_has_rsIDs)
	{
		warn "$ref_stem and $data_stem are not compatible: one has rsIDs and the other does not\n";
	}
	if($data_genotypes_numeric != $ref_data_genotypes_numeric)
	{
		warn "$ref_stem and $data_stem are not compatible: one has numeric genotypes and the other does not\n";
	}
	elsif($data_genotypes_numeric && $ref_data_genotypes_numeric)
	{
		warn "Can't flip numeric genotypes; attempting to merge without flipping\n";
		return $new_stem; ## I am not returning reference hashes like below, so that might cause a problems downstream
	}

	my @overlap_SNPs;
	my @non_overlap_SNPs;
	my @triallelic;
	my @ambigious;
	my @SNPs_to_be_flipped;
	my @weird;
	foreach my $snp (keys %data)
	{
		if(!exists $ref_data{$snp})
		{
			push(@non_overlap_SNPs,$snp);
			next;	
		}
		push(@overlap_SNPs,$snp);
		
		my %alleles_hash;
		my @a1 = keys %{$data{$snp} };
		my @a2 = keys %{$ref_data{$snp} };
		foreach(@a1){$alleles_hash{$_}=1;}
		foreach(@a2){$alleles_hash{$_}=1;}
		my @alleles = keys %alleles_hash;
		
		foreach my $allele(@alleles)
		{
			if($allele !~ /^[a|c|g|t|A|C|G|T]$/)
			{
				push(@weird,$snp);
			}
		}


		#print "$snp: a1(@a1) a2(@a2) alleles(@alleles)\n";
		if(@a1 > 2 || @a2 > 2)
		{
			#print "tri1\n";
			push(@triallelic,$snp);
			next;
		}
		if(is_ambiguous(@a1) || is_ambiguous(@a2))
		{
			#print "amb1\n";
			push(@ambigious,$snp);
			next;
		}
		
		if(@alleles > 4) ## should never happen since .bim files only have two alleles each
		{
			push(@weird,$snp);
		}
		elsif(@alleles == 4)
		{
			#print "flip1\n";
			push(@SNPs_to_be_flipped,$snp);
		}
		elsif(@alleles == 3)
		{
			if(@a1 == 2 && @a2 == 2)
			{
				#print "tri2\n";
				push(@triallelic,$snp);
			}
			else
			{
				#print "flip2\n";
				#if (@a1 == 1 && @a2 == 2) || (@a1 == 2 && @a2 == 1) then we need to flip to match
				push(@SNPs_to_be_flipped,$snp);
			}
	
		}
		elsif(@alleles == 2)
		{
			if(is_ambiguous(@alleles))
			{
				#print "amb3\n";
				push(@ambigious,$snp);
			}
			else
			{
				## if(@a1 == 2 && @a2 == 2) then no need to flip
				## if(@a1 == 1 && @a2 == 2) then @a1 is in @a2; no need to flip, but is this weird?
				## if(@a1 == 2 && @a2 == 1) then @a2 is in @a1; no need to flip, but is this weird?
				## if(@a1 == 1 && @a2 == 1) then @a1 and @a2 have no overlap and are not ambiguous; no need to flip, but is this weird?
			}
		}
		else
		{
			## no need to do anything because this "SNP" appears to be homozygous in the population.
		}
	}

	## Combine SNPs and remove
	my @SNPs_to_remove = (@triallelic,@ambigious,@weird,@non_overlap_SNPs);
	print "# of SNPs overlap between datasets: " . @overlap_SNPs . "\n" if $verbose > 1;
	print $LOG "# of SNPs overlap between datasets: " . @overlap_SNPs . "\n" if $verbose > 1;
	print "# of SNPs non-overlap between datasets: " . @non_overlap_SNPs . "\n" if $verbose > 1;
	print $LOG "# of SNPs non-overlap between datasets: " . @non_overlap_SNPs . "\n" if $verbose > 1;
	print "# of SNPs to remove before merging: " . @SNPs_to_remove . "\n" if $verbose > 1;
	print $LOG "# of SNPs to remove before merging: " . @SNPs_to_remove . "\n" if $verbose > 1;
	print "# of SNPs to flip before merging: " . @SNPs_to_be_flipped . "\n" if $verbose > 1;
	print $LOG "# of SNPs to flip before merging: " . @SNPs_to_be_flipped . "\n" if $verbose > 1;
	
	remove_SNPs($data_stem,\@SNPs_to_remove,$new_stem);
	
	$new_stem = flip_SNPs($new_stem,$new_stem,@SNPs_to_be_flipped);

	return ($new_stem,\@SNPs_to_be_flipped,\@SNPs_to_remove);
}

sub flip_SNPs
{
	my $data_stem = shift;
	my $new_stem = shift;
	my @SNPs_to_flip = @_;
	$new_stem = "$data_stem\_flipped" if $new_stem eq "";
	
	my $flip_file = "$data_stem\_flip_SNPs.txt";
	## Write file of SNPs that need to be flipped
	open(OUT, ">$flip_file");
	foreach(@SNPs_to_flip){print OUT "$_\n";}
	close(OUT);
	$intermediate_files{$flip_file} = 10;

	## Run plink to flip the strands of SNPs in $file_ped
	my $temp = system("$PLINK --allow-no-sex --noweb --bfile $data_stem --flip $flip_file --make-bed $plink_silent --indiv-sort 0 --out $new_stem $memory_flag");
	$intermediate_files{"$new_stem.bed"} = 5;
	$intermediate_files{"$new_stem.bim"} = 5;
	$intermediate_files{"$new_stem.fam"} = 5;
	
	return $new_stem;
}

sub call_callrate
{
	my $stem_name = shift;
	print "\nCalling callrate for $stem_name\n" if $verbose > 0;
	print $LOG "\nCalling callrate for $stem_name\n" if $verbose > 0;
	make_binary_version($stem_name);
	system("$PLINK --allow-no-sex --noweb --bfile $stem_name --maf $MAF --geno $GENO --missing $plink_silent --out $stem_name $memory_flag");
	$intermediate_files{"$stem_name.log"} = 1;
	return "$stem_name.imiss";
}

sub call_het
{
	my $data_stem = shift;
	my $ref_freq_file = shift;
	print "\nCalling het rate for $data_stem\n" if $verbose > 0;
	print $LOG "\nCalling het rate for $data_stem\n" if $verbose > 0;
	make_binary_version($data_stem);
	system("$PLINK --noweb --bfile $data_stem --maf $MAF --geno $GENO --read-freq $ref_freq_file --het $plink_silent --out $data_stem $memory_flag");
	$intermediate_files{"$data_stem.log"} = 1;
	return "$data_stem.het";
}

sub call_sex
{
	my $data_stem = shift;
	my $ref_freq_file = shift;
	print "\nCalling sex for $data_stem\n" if $verbose > 0;
	print $LOG "\nCalling sex for $data_stem\n" if $verbose > 0;
	make_binary_version($data_stem);
	system("$PLINK --noweb --bfile $data_stem --mind $MIND --maf $MAF --geno $GENO --read-freq $ref_freq_file --check-sex $plink_silent --out $data_stem $memory_flag");
	$intermediate_files{"$data_stem.log"} = 1;
	return "$data_stem.sexcheck";
}

sub get_allele_freqs
{
	my $stem_name = shift;
	print "\nGetting allele frequencies for $stem_name\n" if $verbose > 0;
	print $LOG "\nGetting allele frequencies for $stem_name\n" if $verbose > 0;
	make_binary_version($stem_name);
	system("$PLINK --noweb --bfile $stem_name --nonfounders --freq $plink_silent --out $stem_name $memory_flag");
	$intermediate_files{"$stem_name.frq"} = 10;
	$intermediate_files{"$stem_name.log"} = 1;
	return "$stem_name.frq";
}

## Identifies Ancestry Informative Markers (AIMs) = SNPs that genome wide significantly associate with either PCV1 or PCV2
## INPUT1: stem_name to the data
## INPUT2: the .eigenvec file generated from the stem_name data (optional)
## INPUT3: You can force a run of pca even if the .eigenvec file already exists (optional)
## OUtPUT1: Path to the file that has the list of AIMs
## OUTPUT2: An array of the AIMs

sub get_AIMs
{
	my $stem_name = shift;
	my $eigenvec_file = shift;
	my $rerun_pca= shift;
	
	make_binary_version($stem_name);

	$eigenvec_file = "$stem_name.eigenvec" if $eigenvec_file eq "";
	
	print "\nGetting AIMs for $stem_name\n" if $verbose > 0;
	print $LOG "\nGetting AIMs for $stem_name\n" if $verbose > 0;

	## run_pca (e.g. EIGENSTRAT)
	if(!-e $eigenvec_file || $rerun_pca == 1)
	{
		run_pca($stem_name);
	}

	## Convert .eigenvec file into format for PLINK for association tests
	open(IN,"$eigenvec_file") or die "Can't open $eigenvec_file; $!\n";
	open(OUT,">$eigenvec_file\_temp") or die "Can't write to $eigenvec_file\_temp; $!\n";
	my $header = <IN>;
	my %pops;
	while(my $line = <IN>)
	{
		$line =~ s/^\s+//;
		my ($fid,$iid,@pcv) = split(/\s+/,$line);
		print OUT "$fid $iid @pcv\n"; 
	}
	close(IN);
	close(OUT);
	$intermediate_files{"$eigenvec_file\_temp"} = 1;

	## Run Association test for eavh PCV
	system("$PLINK --noweb --bfile $stem_name --pheno $eigenvec_file\_temp --allow-no-sex --mpheno 1 --assoc $plink_silent --out $stem_name\_PCV1 $memory_flag");
	$intermediate_files{"$stem_name\_PCV1.log"} = 1;
	$intermediate_files{"$stem_name\_PCV1.qassoc"} = 1;
	system("$PLINK --noweb --bfile $stem_name --pheno $eigenvec_file\_temp --allow-no-sex --mpheno 2 --assoc $plink_silent --out $stem_name\_PCV2 $memory_flag");
	$intermediate_files{"$stem_name\_PCV2.log"} = 1;
	$intermediate_files{"$stem_name\_PCV2.qassoc"} = 1;
	
	## Pull AIMs from .
	my %AIMs;	
	for(1..2)
	{
		print "\nPulling AIMs from $stem_name\_PCV$_.qassoc\n" if $verbose > 1;
		print $LOG "\nPulling AIMs from $stem_name\_PCV$_.qassoc\n" if $verbose > 1;
		open(IN,"$stem_name\_PCV$_.qassoc");
		my $header = <IN>;
		while(my $line = <IN>)
		{
			$line =~ s/^\s+//;
			my ($chr,$snpID,$pos,$nmiss,$beta,$se,$r2,$T,$p_val) = split(/\s+/,$line);
			if($p_val =~ 'NA'){next;}
			if($p_val <= $ASSOC_P_VAL)
			{
				$AIMs{$snpID} = $p_val;
			}
		}
		close(IN);
	}
	
	## Write out AIMs to file
	open(OUT,">$stem_name\_AIMs.txt");
	foreach my $snp (keys %AIMs)
	{
		print OUT "$snp\n";
	}
	close(OUT);
	$intermediate_files{"$stem_name\_AIMs.txt"} = 10;
	
	print "# of AIMs found: ". (keys %AIMs). "\n" if $verbose > 0;
	print $LOG "# of AIMs found: ". (keys %AIMs). "\n" if $verbose > 0;

	return ("$stem_name\_AIMs.txt",(keys %AIMs));
}

## Makes a PCV1 vs PCV2 plot; will run smart pca to get .eigenvec file if the .eigenvec file is not provided.
## Input1: read in the stem name of the data or the stem name of the .eigenvec file.
## Input2: If you want the points ot be color coded, provide a hash reference with the name=>color of the samples in the data/eigenvec file
## Input3: For the subroutine to rerun pca, if the stem_name.evev file already exists, then specify 1;
## Output: Path to the PFD for the plot.

sub make_PCA_plot {

	## Set file paths
	my $stem_name = shift;
	my $ref_stem = shift;
	my $study_name = shift;
	my $ref_name = shift;
    my $project_onto_1KG = shift;
	my $rerun_pca= shift;

	print "\nMaking PCA plot for $stem_name\n" if $verbose > 0;
	print $LOG "\nMaking PCA plot for $stem_name\n" if $verbose > 0;
	$study_name = "Sample data" if $study_name eq "";
	$ref_name = get_file_name_from_stem($ref_stem) if $ref_name eq "" ;

  #my %population_color_codes_and_sample_population_hash_ref = get_population_assignments_and_colors($stem_name,$ref_stem,$study_name,$ref_name);
	
	if($rerun_pca == 1 || !-e "$stem_name.eigenvec")
    {
        run_pca($stem_name,$project_onto_1KG);
    }
	
	## Convert .eigenvec file and color codes into new input file for R
	open(IN,"$stem_name.eigenvec") or die "Can't open $stem_name.eigenvec; $!\n";
	my $header = <IN>;

	my $max_pop_name_length = 1;
	my $name_length_limit = 30; # use?
	my @lines;
	my %pops;
	my @pch;
	while(my $line = <IN>)
	{
		$line =~ s/^\s+//;
		my ($fid,$iid,@pcv) = split(/\s+/,$line);
        my $color = "blue";
        my $pop = $fid;
        if(exists $onekg_colors{$pop})
        {
            $color = $onekg_colors{$pop};
        }
        $pops{$pop} = $color;
        $max_pop_name_length = length $pop  if length $pop > $max_pop_name_length; ## need this to correctly set the legend
        
        if(!exists $onekg_colors{$pop})
        {
            #print "non ref $pop\n";
            push(@lines,"$fid $iid \"$pop\" \"$color\" @pcv\n"); ## Add the non-reference to back of array to be plotted last 
            push(@pch,13);
        }
        else
        {
            unshift(@lines,"$fid $iid \"$pop\" \"$color\" @pcv\n"); ## Add the reference samples to the front of the array to be plotted first
            unshift(@pch,1);
        }
    }
    close(IN);


	## Write out R input
	open(OUT,">$stem_name\_PCV1vPCV2.txt");
	foreach (@lines){print OUT $_};
	close(OUT);
	$intermediate_files{"$stem_name\_PCV1vPCV2.txt"} = 1;
	
	## Make legend input
	my @legend_pops;
	foreach my $pop (sort keys %pops)
	{
		if($pop =~ /^$study_name$/) 
		{
			unshift(@legend_pops,$pop); ## Add the non-reference to front of array to be at top of legend
		}
		else
		{
			push(@legend_pops,$pop); ## Add the reference samples to the front of the array to be plotted first
		}
	}
	my @legend_colors;
	foreach(@legend_pops)
	{
		push(@legend_colors,$pops{$_});
	}

	## Make R script to plot the PCVs and colors
	my $cex_main = 2 - ((length $study_name) - 10)*.028;
	$cex_main = .4 if $cex_main < .4;
	$cex_main = 2 if $cex_main > 2;
	my $legend_cex = 1 - ($max_pop_name_length-30)*.015;
	$legend_cex = 1 if $legend_cex > 1;
	$legend_cex = .4 if $legend_cex < .4;
	
	open(R,">$stem_name\_PCV1vPCV2.R");
	print R "data <-read.table(\"$stem_name\_PCV1vPCV2.txt\",stringsAsFactors=FALSE,header=F)\n";

	print R "pdf(\"$stem_name\_PCV1vPCV2.pdf\", height=8, width=10.5,family=\"mono\")\n";
	print R "par(mai = c(1,1,1,4), oma = c(1,1,1,1),xpd=T)\n";
	print R "plot(data[,5],data[,6],col=data[,4],xlab=\"PCV1\",ylab=\"PCV2\",pch=c(".join(',',@pch)."),main=\"PCV1 vs PCV2 for $study_name\",cex.main=$cex_main)\n";
	print R "legend(x=\"left\",inset=1.01,legend=c(\"".join('","',@legend_pops)."\"),text.col=c(\"".join('","',@legend_colors)."\"),fill=c(\"".join('","',@legend_colors)."\"),cex=$legend_cex)\n";
	print R "dev.off()\n";
	close(R);
	$intermediate_files{"$stem_name\_PCV1vPCV2.R"} = 1;

	## Run R script
	my $temp = system("R --vanilla --slave < $stem_name\_PCV1vPCV2.R > $stem_name\_R_output");
	system("rm $stem_name\_R_output");
	if($temp > 0)
	{
		die "ERROR!!! Failed to draw eigenstrat results.\n";
	}
	
	print "PCA plot: $stem_name\_PCV1vPCV2.pdf\n" if $verbose > 0;
	print $LOG "PCA plot: $stem_name\_PCV1vPCV2.pdf\n" if $verbose > 0;
	return ("$stem_name\_PCV1vPCV2.pdf");
}

sub run_pca
{	

	my $stem_name = shift;
	my $project_onto_1KG = shift;

	print "\nRunning plink's PCA (new) on $stem_name => $stem_name.eigenvec\n" if $verbose > 0;
	print $LOG "\nRunning plink's PCA (new) on $stem_name => $stem_name.eigenvec\n" if $verbose > 0;

	make_binary_version($stem_name);

	
	#LD PRUNING STEP 
	#-----------------------------
	#inputs: window size (10kb), step size (10 or 100 variants), and r2 threshold (0.2)
	#outputs: <input>.prune.in and <input>.prune.out files 
	#add --exclude/--extract to skip them in the subsequent PCA-approx step
	

	print "\nRunning LD pruning\n" if $verbose > 0;
	print $LOG "\nRunning LD pruning\n" if $verbose > 0;

	my $temp = system("$PLINK --bfile $stem_name --indep-pairwise 10 10 0.2 --out $stem_name\_pruned $memory_flag");
	#my $temp = system("plink --bfile $stem_name --indep-pairwise 50 5 0.1 --out $stem_name\_pruned");

	################################################################################################################

    ## RUN PLINK PCA
	my $cluster_names;
    if($project_onto_1KG eq 1)
    {
        $cluster_names = "--pca-cluster-names ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI";
    }


	# added --extract flag from LD pruning
	#my $temp = system("$PLINK --allow-no-sex --bfile $stem_name --family --pca $cluster_names --maf $MAF --geno $GENO --out $stem_name");
    my $temp = system("$PLINK2 --allow-no-sex --bfile $stem_name --pca approx --extract $stem_name\_pruned.prune.in --maf $MAF --geno $GENO --out $stem_name $memory_flag");


	if($temp > 0)
	{
		die "ERROR!!! PLINK's PCA failed (plink2, pca approx step).\n";
	}
	$intermediate_files{"$stem_name.log"} = 1;
	return "$stem_name.eigenvec";
}	

## Run get rough IBDs with PLINK, then run IMUS to get unrelated set, then --keep the unrelated set from the input data and return
## Input: 1. the stem to the plink data files in either ped/map or bed/bim/fam (required)
## Input: 2. the stem to the new plink data files in either ped/map or bed/bim/fam (optional)
sub get_unrelateds
{
	## Set file paths
	my $stem_name = shift;
	my $new_stem_name = shift;
	$new_stem_name = "$stem_name\_unrelateds" if $new_stem_name eq "";
	make_binary_version($stem_name);

	print "\nGet Unrelated set for $stem_name => $new_stem_name\n" if $verbose > 0;
	print $LOG "\nGet Unrelated set for $stem_name => $new_stem_name\n" if $verbose > 0;

    my $temp = system("$PLINK --allow-no-sex --bfile $stem_name --maf $MAF --geno $GENO --thin-count 10000 --rel-cutoff 0.09375 --mind --out $stem_name $memory_flag");
	if($temp > 0)
	{
		die "ERROR!!! PLINK's rel-cutoff failed.\n";
	}
	$intermediate_files{"$stem_name.log"} = 1;
    

	## Set rough IBD estimates
    #my $rough_IBDs_genome = calculate_IBD_estimates($stem_name,"$stem_name\_naive");
    #$intermediate_files{"$rough_IBDs_genome"} = 3;

    #my %ibd_estimates;
    #$ibd_estimates{'FILE'}= $rough_IBDs_genome;
    #$ibd_estimates{'FID1'}= 1;
    #$ibd_estimates{'IID1'}= 2;
    #$ibd_estimates{'FID2'}= 3;
    #$ibd_estimates{'IID2'}= 4;
    #$ibd_estimates{'PI_HAT'}= 10;
	
	## RUN IMUS to get number of unrelated samples
    #print "\nRun IMUS on IBDs for all original samples\n" if $verbose > 1;
    #print $LOG "\nRun IMUS on IBDs for all original samples\n" if $verbose > 1;
    #my @IMUS_commands = ("--do_IMUS",1,"--do_PR",0,"--ibd_estimates",\%ibd_estimates,"--verbose",$verbose,"--output_dir","$stem_name\_IMUS","--lib",$lib_dir, "--rel_threshold",$THRESHOLD,"--log_file_handle",$LOG);
    #my ($unrelated_file, @unrelated_samples) = PRIMUS::IMUS::run_IMUS(@IMUS_commands);

	my $unrelated_file = "$stem_name.rel.id";
    my @unrelated_samples = `cat $unrelated_file`;
    chomp(@unrelated_samples);

    if(!-e $unrelated_file || @unrelated_samples < 1){die "IMUS failed to return an unrelated set\n"};
	$intermediate_files{$unrelated_file} = 3;

	my $num_unrelated_samples = @unrelated_samples; ## Decrement one for the header
	print "\nNum unrelated samples: $num_unrelated_samples\n" if $verbose > 1;
	print $LOG "\nNum unrelated samples: $num_unrelated_samples\n" if $verbose > 1;

	## Make an unrelated version of the data
	keep_samples($stem_name,$new_stem_name,$unrelated_file);
	$intermediate_files{"$new_stem_name.bed"} = 10;
	$intermediate_files{"$new_stem_name.bim"} = 10;
	$intermediate_files{"$new_stem_name.fam"} = 10;
	
	return ($new_stem_name,@unrelated_samples);
}

sub calculate_IBD_estimates {

	# Generates a .genome file with the cleaned_all_samples_stem files 

	my $stem_name = shift;
	my $new_stem_name = shift;
	my $freq_file = shift;

	$new_stem_name = "$stem_name" if $new_stem_name eq "";
	print "\nCalculating IBDs for $stem_name (.frq = $freq_file) => $new_stem_name.genome\n" if $verbose > 0;
	print $LOG "\nCalculating IBDs for $stem_name (.frq = $freq_file) => $new_stem_name.genome\n" if $verbose > 0;
	
	make_binary_version($stem_name);
	
	my $read_freq = "--read-freq $freq_file";
	$read_freq = "" if $freq_file eq "";

	# Maybe roll back this change? 
	
	system("$PLINK --noweb --bfile $stem_name $read_freq --maf $MAF --geno $GENO $plink_silent --make-bed --out $stem_name\_temp $memory_flag");
	my $temp = system("$PLINK --noweb --bfile $stem_name\_temp $read_freq --genome --maf $MAF --geno $GENO --mind $MIND $plink_silent --out $new_stem_name --min 0 $memory_flag");
	if($temp > 0)
	{
		die "ERROR!!! PLINK failed to calculate IBD estimates; check log file: $new_stem_name.log\n";
	}

	return "$new_stem_name.genome";
}

## Removes duplicate SNP names or duplicate SNP positions from PLINK data files
## Input: 1. the stem to the plink data files in either ped/map or bed/bim/fam (required)
## Input: 2. the stem to the new plink data files in either ped/map or bed/bim/fam (optional)
## Output: 1. A duplicate free version of the data provided saved at either the output stem provided or at with the
## Output: 2. An array containing the names of all the SNPs that were removed.
sub remove_dups
{
	## Set file paths
	my $stem_name = shift;
	my $no_dups_name = shift;
	$no_dups_name = "$stem_name\_noDups" if $no_dups_name eq "";
	make_binary_version($stem_name);

	print "\nRemoving Dups from $stem_name => $no_dups_name\n" if ($verbose > 0);
	print $LOG "\nRemoving Dups from $stem_name => $no_dups_name\n" if ($verbose > 0);
	
	## Load SNP data
	my %snp_names;
	my %snp_poss;
	my %pos_to_name;
	open(IN,"$stem_name.bim") or die "Cannot open $stem_name.bim: $!\n";
	while(my $line = <IN>)
	{
		my ($chr,$rsID,$genetic_pos,$pos,$a1,$a2) = split(/\s+/,$line);
		$snp_names{$rsID}++;
		$snp_poss{"$chr\__$pos"}++;
		push(@{$pos_to_name{"$chr\__$pos"} },$rsID);
	}
	close IN;
	
	## Find the duplicated rsIDs and duplicate positions
	my @bad_SNPs;
	foreach my $snp (keys %snp_names)
	{
		push(@bad_SNPs,$snp)  if($snp_names{$snp} > 1);
	}
	foreach my $pos (keys %snp_poss)
	{
		if($snp_poss{$pos} > 1)
		{
			push(@bad_SNPs,@{$pos_to_name{$pos} }) if($snp_poss{$pos} > 1);
		}
	}
	
	remove_SNPs($stem_name,\@bad_SNPs,$no_dups_name);
	$intermediate_files{"$no_dups_name.bed"} = 10;
	$intermediate_files{"$no_dups_name.bim"} = 10;
	$intermediate_files{"$no_dups_name.fam"} = 10;

	print "# of Dup_SNPs: " . @bad_SNPs. "\n" if $verbose > 1;
	print $LOG "# of Dup_SNPs: " . @bad_SNPs. "\n" if $verbose > 1;
	return ($no_dups_name,@bad_SNPs);
}

#############################################################################
############## INTERNALLY USEFUL METHODS ####################################
#############################################################################

sub print_intermediate_files
{
	print "INTERMEDIATE FILES:\n";
	foreach my $file (keys %intermediate_files)
	{
		print "$file $intermediate_files{$file}\n";
	}
}

sub remove_intermediate_files
{
	my $max_priority_to_remove = shift;
	my $output_dir = shift;
	my @dirs_to_remove;
	foreach my $file (keys %intermediate_files)
	{
		if(-d $file)
		{
			system("rm -r $file/*") if $intermediate_files{$file} <= $max_priority_to_remove;
			push(@dirs_to_remove,$file) if $intermediate_files{$file} <= $max_priority_to_remove;;
		}
		else
		{
			system("rm -r $file") if $intermediate_files{$file} <= $max_priority_to_remove && -e $file;
		}
	}
	system("rm $output_dir/*.log") if(-e "$output_dir/*.log");
	system("rm $output_dir/*~") if(-e "$output_dir/*~");
	system("rm $output_dir/*.hh") if(-e "$output_dir/*.hh");
	foreach(@dirs_to_remove)
	{
		#system("rmdir $_");
	}

}

sub get_population_assignments_and_colors
{
	my $data_stem = shift;
	my $ref_stem = shift;
	my $study_name = shift;
	my $ref_name = shift;
	make_binary_version($data_stem);
	make_binary_version($ref_stem) if $ref_stem ne "";
	
	$study_name = get_file_name_from_stem($data_stem) if $study_name eq "";
	$ref_name = get_file_name_from_stem($ref_stem) if $ref_name eq "";

	my %IDs_to_pops;
	open(IN,"$data_stem.fam");
	while(my $line = <IN>)
	{
		my ($FID,$IID,@rest) = split(/\s+/,$line);
    ##$IDs_to_pops{"$FID\__$IID"} = $study_name;
		$IDs_to_pops{"$IID"} = $study_name;
	}
	close(IN);

	# my @hm3_pops = ("ASW","CEU","CHB","CHD","GIH","JPT","LWK","MEX","MKK","TSI","YRI");
	# my %hm3_colors = ("ASW","steelblue4","CEU","red","CHB","darkred","CHD","green","GIH","darkgreen","JPT","slateblue1","LWK","slateblue4","MEX","sienna1","MKK","sienna4","TSI","black","YRI","darkgoldenrod3");
	if($ref_name =~ /onekg/i)
	{
		open(IN,"$onekg/1KG_sample_populations.txt");
		while(my $line = <IN>)
		{
			chomp($line);
			my ($FID,$IID,$pop) = split(/\s+/,$line);
      #$IDs_to_pops{"$FID\__$IID"} = $pop;
			$IDs_to_pops{"$IID"} = $pop;
		}
		close(IN);
		foreach my $pop (keys %onekg_colors)
		{
			$IDs_to_pops{$pop} = $onekg_colors{$pop};
		}
		$IDs_to_pops{$study_name} = "blue";
	}
	else
	{
		open(IN,"$ref_stem.fam");
		while(my $line = <IN>)
		{
			my ($FID,$IID,@rest) = split(/\s+/,$line);
      #$IDs_to_pops{"$FID\__$IID"} = $ref_name;
			$IDs_to_pops{"$IID"} = $ref_name;
		}
		close(IN);
		$IDs_to_pops{$study_name} = "blue"; ## R color used to plot non-ref data
		$IDs_to_pops{$ref_name} = "red"; ## R color used to plot the alt-ref data;
	}
	return %IDs_to_pops;
}

sub is_ambiguous
{
	my @alleles = @_;
	my @amb1 = qw(g c);
	my @amb2 = qw(a t);
	my @amb3 = qw(G C);
	my @amb4 = qw(A T);
	if(do_arrays_match(\@alleles,\@amb1) || do_arrays_match(\@alleles,\@amb2) || do_arrays_match(\@alleles,\@amb3) || do_arrays_match(\@alleles,\@amb4))
	{
		#print "yes @alleles\n";
		return 1;
	}
	else
	{
		#print "no @alleles\n";
		return 0;
	}
}

sub keep_samples
{
	my $stem_name = shift;
	my $new_stem_name = shift;
	my $samples_to_keep = shift;
	make_binary_version($stem_name);

	if($stem_name eq "")
	{
		$new_stem_name = "$stem_name\_unrelated";
	}

	print "\nKeep $samples_to_keep samples from $stem_name => $new_stem_name\n" if $verbose > 1;
	print $LOG "\nKeep $samples_to_keep samples from $stem_name => $new_stem_name\n" if $verbose > 1;

	## Make an unrelated version of the data
	my $temp = system("$PLINK --allow-no-sex --noweb --geno $GENO --keep $samples_to_keep --bfile $stem_name --indiv-sort 0 --make-bed $plink_silent --out $new_stem_name $memory_flag");
	if($temp > 0){die "ERROR!!! PLINK failed to do --keep; check log file: $new_stem_name.log\n";}
	my $temp = system("$PLINK --allow-no-sex --noweb --geno $GENO --mind $MIND --bfile $new_stem_name --make-bed $plink_silent --out $new_stem_name $memory_flag");
	if($temp > 0){die "ERROR!!! PLINK failed to do --keep; check log file: $new_stem_name.log\n";}
	return $new_stem_name;
}

sub remove_SNPs
{
	my $stem_name = shift;
	my $SNPs_to_remove_array_ref = shift;
	my $new_stem_name = shift;
	$new_stem_name = "$stem_name" if $new_stem_name eq "";

	if($SNPs_to_remove_array_ref eq "")
	{
		my @temp = ();
		$SNPs_to_remove_array_ref = \@temp;
	}

	if(@$SNPs_to_remove_array_ref < 1 && $stem_name =~ /^$new_stem_name$/)
	{
		return $new_stem_name;
	}
	print "\nRemoving ". @$SNPs_to_remove_array_ref ." SNPs from $stem_name => $new_stem_name\n" if $verbose > 1;
	print $LOG "\nRemoving ". @$SNPs_to_remove_array_ref ." SNPs from $stem_name => $new_stem_name\n" if $verbose > 1;

	$intermediate_files{"$new_stem_name.bed"} = 10;
	$intermediate_files{"$new_stem_name.bim"} = 10;
	$intermediate_files{"$new_stem_name.fam"} = 10;

	my %bad_SNPs = map { $_ => 1 } @{$SNPs_to_remove_array_ref};
	make_binary_version($stem_name);

	my %snp_names;
	open(IN,"$stem_name.bim") or die "Cannot open $stem_name.bim: $!\n";
	while(my $line = <IN>)
	{
		my ($chr,$rsID,$genetic_pos,$pos,$a1,$a2) = split(/\s+/,$line);
		$snp_names{$rsID}++;
	}
	close IN;

	## Write file with SNPs to keep
	open(OUT,">$stem_name.SNPs_to_keep.txt");
	foreach my $snp (keys %snp_names)
	{
		print OUT "$snp\n" if !exists $bad_SNPs{$snp};
	}
	close OUT;

	## Run plink to generate the no_dups version of stem_file
	system("$PLINK --allow-no-sex --noweb --bfile $stem_name --extract $stem_name.SNPs_to_keep.txt --indiv-sort 0 --make-bed $plink_silent --out $new_stem_name $memory_flag");
	system("rm $stem_name.SNPs_to_keep.txt") if $test == 0;
	return ($new_stem_name);
}

sub make_non_binary_version
{
	my $stem_name = shift;
	my $new_stem_name = shift;
	my $overwrite = shift;
	$new_stem_name = $stem_name if $new_stem_name eq "";
	
	## Don't run if new_stem_name already exists and if overwrite is not 1
	if(!-e "$new_stem_name.ped" || !-e "$new_stem_name.map" || $overwrite == 1){
		print "\nMake_non_binary_version $stem_name => $new_stem_name\n" if $verbose > 1;
		print $LOG "\nMake_non_binary_version $stem_name => $new_stem_name\n" if $verbose > 1;
		if(-e "$stem_name.map" && -e "$stem_name.ped"){
			system("$PLINK --file $stem_name --recode $plink_silent --out $new_stem_name $memory_flag");
		}
		elsif(-e "$stem_name.bed" && -e "$stem_name.bim" && -e "$stem_name.fam"){
			print "here: $PLINK --bfile $stem_name --recode $plink_silent --out $new_stem_name $memory_flag\n";
			system("$PLINK --bfile $stem_name --recode $plink_silent --out $new_stem_name $memory_flag");
		}
		else
		{
			die "$stem_name is not a valid PLINK --file or --bfile input name";
		}
	}
	else
	{
		return $stem_name;
	}
	return $new_stem_name;
}

sub make_binary_version
{
	my $stem_name = shift;
	my $new_stem_name = shift;
	my $overwrite = shift;
	$new_stem_name = $stem_name if $new_stem_name eq "";
	
	## Don't run if new_stem_name already exists and if overwrite is not 1
	if(!-e "$new_stem_name.bed" || !-e "$new_stem_name.bim" || !-e "$new_stem_name.fam" || $overwrite == 1){
		print "\nMake binary version of $stem_name => $new_stem_name\n" if $test == 1;
		print $LOG "\nMake binary version of $stem_name => $new_stem_name\n" if $test == 1;
		if(-e "$stem_name.bed" && -e "$stem_name.bim" && -e "$stem_name.fam"){
			system("$PLINK --bfile $stem_name --indiv-sort 0 --make-bed $plink_silent --out $new_stem_name $memory_flag");
		}
		elsif(-e "$stem_name.map" && -e "$stem_name.ped"){
			system("$PLINK --file $stem_name --indiv-sort 0 --make-bed $plink_silent --out $new_stem_name $memory_flag");
		}
		else
		{
			die "$stem_name is not a valid PLINK --file or --bfile input name";
		}
		$intermediate_files{"$new_stem_name.bed"} = 10;
		$intermediate_files{"$new_stem_name.bim"} = 10;
		$intermediate_files{"$new_stem_name.fam"} = 10;
	}
	return $new_stem_name;
}

sub get_file_name_from_stem
{
	my $file_stem = shift;
	$file_stem =~ /([^\/]+)$/;
	my $file_name = $1;
	return $file_name;
}

sub do_arrays_match
{
        my $arr1_ref = shift;
        my $arr2_ref = shift;

	#print "@$arr1_ref @$arr2_ref\n";
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

sub normalize_vector
{
	my $ref = shift;
	my $sum;
	foreach my $val (@$ref)
	{
		$sum += $val;
	}
	for(my $i = 0; $i < @$ref; $i++)
	{
		@$ref[$i] = @$ref[$i]/$sum;
	}
	return $sum;
}

sub set_paths
{
	my $hash_ref = shift;
	foreach my $name (keys %$hash_ref)
	{
		if($name eq "PLINK"){$PLINK=$$hash_ref{$name} }
		elsif($name eq "LIB"){$lib_dir = $$hash_ref{$name} }
		else{die "INVALID SETTING FOR prePRIMUS_pipeline_v2.pm\n";}
	}
	test_paths();
}

sub test_paths
{
	return; ## This is broken right now and needs to be fixes, but for now I am not running it because everyone thinks that the code is broken
}


