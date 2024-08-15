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


my $useage = "\n\nUSAGE: $0\t  -input [IBD_file]  -output_dir [output_dir]  -threshold [num; default = 0.1]  -[high|low|mean|tails]/[b|q]trait [trait_file]

\n\nType '$0 -help' for explanation of input\n\n";

my $do_IMUS = 1;
my $do_PR = 1;
my $verbose;
my $lib_dir;
my $LOG;

my %arg;
my $relatedness_file;
my $relatedness_file_name;
my $ped_file = "none";
my $missingness_file = "none";
my $output_dir;
my $ersa_data;
my $network_ctr; # the next network # to add to %networks hash
my $LOWEST_MAX_NETWORK_SIZE = 60;
my $THRESHOLD = .1; # default is .1
my $MIN_LIKELIHOOD = .1;
my $EXCLUDE_VALUE = 0;
my $OUTFILE_HEADER;
my $IBD_file_ref;
my $RELATEDNESS_COLUMN = -1;
my $ID1_COLUMN = -1;
my $ID2_COLUMN = -1;
my $FID1_COLUMN = -1;
my $FID2_COLUMN = -1;
my %TRAIT_FID_COLUMNS = -1;
my %TRAIT_ID_COLUMNS = -1;
my %TRAIT_DATA_COLUMNS = -1;
my $PRINT_ALTERNATE_RESULTS = 1;


### Here 'networks' could also be considered related groups or families; 
my %networks; ## Hash of arrays: Key = network; Value = array (list) of IIDs in that network
my %id_network; ## Regular hash; Key = IID1; Value = network that the individual belongs to
my %id_id_scores; ## Regular hash; Key = IID1_IID2; Value = PI_HAT value
my %id_id_all_info;
my @trait_refs;
my @trait_order;
my %trait_files; ## Key = file name; Value = trait type (quantitative or binary)
my %child_parents;

sub run_IMUS
{

	print "\nRunning IMUS ...\n\n";
	#parseCommandLine(@_);
	reset_values();
	set_values2(@_);

	open($LOG,">$output_dir/$relatedness_file_name.log") if($LOG eq "");

	print "Relatedness_file: $relatedness_file\n" if $verbose > 0;
	print "Threshold: $THRESHOLD\n" if $verbose > 0;
	print "Selection criteria are based on the following:\n" if $verbose > 0;
	print $LOG "Relatedness_file: $relatedness_file\n" if $verbose > 0;
	print $LOG "Threshold: $THRESHOLD\n" if $verbose > 0;
	print $LOG "Selection criteria are based on the following:\n" if $verbose > 0;
	foreach(@trait_order)
	{
		print "\t$_ ($trait_files{$_})\n" if $verbose > 0;
		print $LOG "\t$_ ($trait_files{$_})\n" if $verbose > 0;
	}
	
	if($verbose >= 1){print "\nIDENTIFYING FAMILY NETWORKS IN DATA\n";}
	if($verbose >= 1){print "Writing network files to $output_dir/\n";}
	if($verbose >= 1){print $LOG "\nIDENTIFYING FAMILY NETWORKS IN DATA\n";}
	if($verbose >= 1){print $LOG "Writing network files to $output_dir/\n";}

	if($verbose >= 1){print "Loading data...\n";}
	if($verbose >= 1){print $LOG "Loading data...\n";}
	load_data($relatedness_file);
	load_trait_data();
	if($verbose >= 1){print "done.\n";}
	if($verbose >= 1){print $LOG "done.\n";}

	if($verbose >= 1){print $LOG "colapsing networks...\n";}
	colapse_networks();
	if($verbose >= 1){print "done.\n";}
	if($verbose >= 1){print $LOG "done.\n";}

	foreach my $key (keys %networks)
	{
		#print "key: $key; @{ $networks{$key} }\n";
	}
	
	write_out_networks($relatedness_file_name);
	
	if(!$do_IMUS){return 1;}

	if($verbose > 0){print "\nIDENTIFYING A MAXIMUM UNRELATED SET\n";}
	if($verbose > 0){print $LOG "\nIDENTIFYING A MAXIMUM UNRELATED SET\n";}
	if($verbose eq 1){print "Checking for large networks...\n";}
	if($verbose eq 1){print $LOG "Checking for large networks...\n";}
	breakup_large_networks();
	if($verbose > 0){print "done.\n";}
	if($verbose > 0){print $LOG "done.\n";}

	my $num_networks = keys %networks;

	if($verbose > 0){print "# of family networks: $num_networks\n";}
	if($verbose > 0){print $LOG "# of family networks: $num_networks\n";}
	if($verbose >= 1){print "Writing out unrelated set\n";}
	if($verbose >= 1){print $LOG "Writing out unrelated set\n";}
	my %PRIMUS_unrelated_set = write_out_independent_set($relatedness_file_name);
	if($verbose >= 1){print "done.\n";}
	if($verbose >= 1){print $LOG "done.\n";}

	if($verbose >= 1){print "Testing alternative methods...\n";}
	if($verbose >= 1){print $LOG "Testing alternative methods...\n";}
	compare_alternative_methods($relatedness_file_name,\%PRIMUS_unrelated_set,%networks);
	if($verbose >= 1){print "done.\n";}
	if($verbose >= 1){print $LOG "done.\n";}

	print "unrelated_file: $relatedness_file_name\_maximum_independent_set\n" if $verbose > 0;
	print "unrelated_set size: " . (keys %PRIMUS_unrelated_set) . "\n" if $verbose > 0;
	print $LOG "unrelated_file: $relatedness_file_name\_maximum_independent_set\n" if $verbose > 0;
	print $LOG "unrelated_set size: " . (keys %PRIMUS_unrelated_set) . "\n" if $verbose > 0;

	return ("$output_dir/$relatedness_file_name\_maximum_independent_set",(keys %PRIMUS_unrelated_set));
}

#############################################################################################
#### SUBROUTINES ############################################################################
#############################################################################################
sub reset_values
{
	$do_IMUS = 1;
	$do_PR = 1;
	$verbose = 1;
	$lib_dir = "";

	%arg = ();
	$relatedness_file = "";
	$relatedness_file_name = "";
	$ped_file = "none";
	$missingness_file = "none";
	$output_dir = "";
	$network_ctr = 0; # the next network # to add to %networks hash
	$LOWEST_MAX_NETWORK_SIZE = 60;
	$THRESHOLD = .1; # default is .1
	$MIN_LIKELIHOOD = .1;
	$EXCLUDE_VALUE = 0;
	$OUTFILE_HEADER = "";
	$IBD_file_ref = "";
	$RELATEDNESS_COLUMN = -1;
	$ID1_COLUMN = -1;
	$ID2_COLUMN = -1;
	$FID1_COLUMN = -1;
	$FID2_COLUMN = -1;
	%TRAIT_FID_COLUMNS = -1;
	%TRAIT_ID_COLUMNS = -1;
	%TRAIT_DATA_COLUMNS = -1;
	$PRINT_ALTERNATE_RESULTS = 1;


	### Here 'networks' could also be considered related groups or families; 
	%networks = (); ## Hash of arrays: Key = network; Value = array (list) of IIDs in that network
	%id_network = (); ## Regular hash; Key = IID1; Value = network that the individual belongs to
	%id_id_scores = (); ## Regular hash; Key = IID1_IID2; Value = PI_HAT value
	%id_id_all_info = ();
	@trait_refs = ();
	@trait_order = ();
	%trait_files = (); ## Key = file name; Value = trait type (quantitative or binary)
	%child_parents = ();
	
}

sub help
{
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

sub parseCommandLine 
{
	my $trait_ctr = 1;
        for (my $i = 0; $i <= $#ARGV; $i++)
        {
                if ($ARGV[$i] =~ /^-/)
                {
                        if($ARGV[$i] =~ /trait/)
			{
				if($ARGV[$i+1] ne "NA" && $ARGV[$i+1] !~ /^-/)
				{
					## Check if user is trying to weight on a qtrait first; if so, default to weighting on size, then qtrait
					if($ARGV[$i] =~ /qtrait/ && @trait_order eq 0)
					{
						push(@trait_order,"size");
					}
					push(@trait_order,$ARGV[$i+1]);
					$trait_files{$ARGV[$i+1]} = $ARGV[$i]; ## Key is file and value is trait type
				}
			}
			$arg{$ARGV[$i]} = $ARGV[$i+1];
                }
	}
	
	## Append the size trait to the end of traits if it is not already in it.
	my @arr = %trait_files;
	if(!grep(/size/i,@trait_order))
	{
		push (@trait_order,"size");
	}
	$trait_files{"size"} = "-size"; ## Key is supposed to be file and value is trait type

	if(exists $arg{-help})
	{
		help();
		exit;
	}
        die("\n\nInput file required $useage") if (!($arg{-input}));

}

sub set_values2
{
	
	GetOptionsFromArray(
		\@_,
		# Diagnostic options
		'verbose=i' => \$verbose,
		'help|?'  => sub{help()},
		
		# Settings
		"rel_threshold=f" => \$THRESHOLD, 
		"int_likelihood_cutoff=f" => \$MIN_LIKELIHOOD, 
		"do_IMUS=i" => \$do_IMUS,
		"do_PR=i" => \$do_PR,
		"missing_val=f"=> \$EXCLUDE_VALUE,
		"output_dir=s" => \$output_dir,
		#"ersa_data=s" => \$ersa_data,
		"trait_order=s" => sub{ @trait_order = @{$_[1]} },
		"lib=s"=>\$lib_dir,
		"traits=s" => sub{ %trait_files = %{$_[1]} },
		"log_file_handle=s"=>\$LOG,
		"ibd_estimates=s" => sub
		{
			my %ibds = %{@_[1]};
			$IBD_file_ref = \%ibds;
			$relatedness_file = $ibds{'FILE'};
			if($relatedness_file =~ /\/([^\/]+)$/)  ## if there is a path to the file, ignore all the path, and select the name of the file
			{	
				$relatedness_file_name = $1;
			}
			else
			{
				$relatedness_file_name = $relatedness_file;
			}
			$FID1_COLUMN = $ibds{'FID1'} -1;
			$ID1_COLUMN = $ibds{'IID1'} -1;
			$FID2_COLUMN = $ibds{'FID2'} -1;
			$ID2_COLUMN = $ibds{'IID2'} -1;
			$RELATEDNESS_COLUMN = $ibds{'PI_HAT'} -1;


		},


	) or die "Failed to parse options for IMUS\n";

	if(!-d $output_dir)
	{
		mkdir("$output_dir") or die "Can't make $output_dir; $!\n";
	}

	## If traits were empty, then add size
	if (!exists $trait_files{'size'})
	{
		$trait_files{'size'} = 'size';
		push(@trait_order,'size');
	}
	foreach my $file (keys %trait_files)
	{
		#print "Trait file: $file\n";
		## Check that trait files exist
		if($file eq "NA" || $file =~ /^-/)
		{
			delete $trait_files{$file};
			next;
		}
		elsif($file =~ /size/i)
		{
			next;
		}
		elsif(!-e $file)
		{
			die "Trait file $file does not exist\n";
		}

		my $trait_type = $trait_files{$file};
		
	
		open(IN,$file) or die "Trait file $file failed to open; $!\n";
		my $header = <IN>; chomp($header);
		close(IN);
		$header =~ s/^\s+//;
		
		## If there are only two columns, assume columns "ID  Trait_value"
		my @temp = split(/\s+/,$header);
		if(@temp == 3)
		{
			$TRAIT_FID_COLUMNS{$file} = 0;
			$TRAIT_ID_COLUMNS{$file} = 1;
			$TRAIT_DATA_COLUMNS{$file} = 2;
		}
		elsif($file =~ /.ped$/i)
		{
			$TRAIT_FID_COLUMNS{$file} = 0;
			$TRAIT_ID_COLUMNS{$file} = 1;
			$TRAIT_DATA_COLUMNS{$file} = 5;
		}
		else ## get columns
		{
			print "For $file:\n";
			$TRAIT_FID_COLUMNS{$file} = get_correct_column("Enter the column # containing the FIDs:  ",$header);
			$TRAIT_ID_COLUMNS{$file} = get_correct_column("Enter the column # containing the IDs:  ",$header);
			$TRAIT_DATA_COLUMNS{$file} = get_correct_column("Enter the column # containing the $trait_type data:  ",$header);
		}
	}
}

sub set_values
{
	$relatedness_file = $arg{-input};
	if($relatedness_file =~ /\/([^\/]+)$/)  ## if there is a path to the file, ignore all the path, and select the name of the file
	{	
		$relatedness_file_name = $1;
	}
	else
	{
		$relatedness_file_name = $relatedness_file;
	}
	if(!-e $relatedness_file)
	{
		die "input_file $relatedness_file does not exist\n";
	}

	if(exists $arg{-BKcutoff})
	{
		print "ERROR! -BKcutoff invalid with versions 2.7 and above\n";
		exit;
		#$LOWEST_MAX_NETWORK_SIZE = $arg{-BKcutoff};
	}

	if(exists $arg{-threshold})
	{
		$THRESHOLD = $arg{-threshold};
	}

	if(exists $arg{-missing_value})
	{
		$EXCLUDE_VALUE = $arg{-missing_value};
	}

	if(exists $arg{-output_dir}){
		$output_dir = $arg{-output_dir};
	}
	else{
		$output_dir = "$relatedness_file\_results";
	}
	if(!-d $output_dir)
	{
		mkdir("$output_dir") or die "Can't make $output_dir; $!\n";
	}

	## Determine the type of intput file or column containing relatedness
	open(IN, $relatedness_file);
	my $header = <IN>; chomp($header);
	$header =~ s/^\s+//;
	close IN;
	my @h_elements = split(/\s+/,$header);

	for(my $col = 0; $col < @h_elements; $col++)
	{
		if(@h_elements[$col] =~ /PI_HAT/i || @h_elements[$col] =~ /Kinship/i)
		{
			$RELATEDNESS_COLUMN = $col;
		}
		if(@h_elements[$col] =~ /IID1/i || (@h_elements[$col] =~ /ID1/i && @h_elements[$col] =~ /FID1/i))
		{
			$ID1_COLUMN = $col;
		}
		if(@h_elements[$col] =~ /IID2/i || (@h_elements[$col] =~ /ID2/i && @h_elements[$col] =~ /FID2/i))
		{
			$ID2_COLUMN = $col;
		}
		if(@h_elements[$col] =~ /FID1/i || @h_elements[$col] =~ /FID$/i)
		{
			$FID1_COLUMN = $col;
		}
		if(@h_elements[$col] =~ /FID2/i || @h_elements[$col] =~ /FID$/i)
		{
			$FID2_COLUMN = $col;
		}
	}
	if($RELATEDNESS_COLUMN eq -1)
	{
		print "Unable to identify the relatedness column in input file\n";
		$RELATEDNESS_COLUMN = get_correct_column("Enter the column # containing the relatedness data:  ",$header);
	}
	if($ID1_COLUMN eq -1)
	{
		print "Unable to identify the ID1 column in input file\n";
		$ID1_COLUMN = get_correct_column("Enter the column # containing ID1:  ",$header);
	}
	if($ID2_COLUMN eq -1)
	{
		print "Unable to identify the ID2 column in input file\n";
		$ID2_COLUMN = get_correct_column("Enter the column # containing ID2:  ",$header);
	}
	if($FID1_COLUMN eq -1)
	{
		print "Unable to identify the FID1 column in input file\n";
		$FID1_COLUMN = get_correct_column("Enter the column # containing FID1:  ",$header);
	}
	if($FID2_COLUMN eq -1)
	{
		print "Unable to identify the FID2 column in input file\n";
		$FID2_COLUMN = get_correct_column("Enter the column # containing FID2:  ",$header);
	}

	foreach my $file (keys %trait_files)
	{
		print "Trait file: $file\n";
		## Check that trait files exist
		if($file eq "NA" || $file =~ /^-/)
		{
			delete $trait_files{$file};
			next;
		}
		elsif($file =~ /size/i)
		{
			next;
		}
		elsif(!-e $file)
		{
			die "Trait file $file does not exist\n";
		}

		my $trait_type = $trait_files{$file};
		
	
		open(IN,$file) or die "Trait file $file failed to open; $!\n";
		my $header = <IN>; chomp($header);
		close(IN);
		$header =~ s/^\s+//;
		
		## If there are only two columns, assume columns "ID  Trait_value"
		my @temp = split(/\s+/,$header);
		if(@temp == 3)
		{
			$TRAIT_FID_COLUMNS{$file} = 0;
			$TRAIT_ID_COLUMNS{$file} = 1;
			$TRAIT_DATA_COLUMNS{$file} = 2;
		}
		elsif($file =~ /.ped$/i)
		{
			$TRAIT_FID_COLUMNS{$file} = 0;
			$TRAIT_ID_COLUMNS{$file} = 1;
			$TRAIT_DATA_COLUMNS{$file} = 5;
		}
		else ## get columns
		{
			print "For $file:\n";
			$TRAIT_FID_COLUMNS{$file} = get_correct_column("Enter the column # containing the FIDs:  ",$header);
			$TRAIT_ID_COLUMNS{$file} = get_correct_column("Enter the column # containing the IDs:  ",$header);
			$TRAIT_DATA_COLUMNS{$file} = get_correct_column("Enter the column # containing the $trait_type data:  ",$header);
		}
	}
}

## Method used for setting the data column of interest
sub get_correct_column
{
	my $question = shift;
	my $header = shift; chomp($header);
	my $column = -1;

	my @h_elements = split(/\s+/,$header);

	while($column eq -1)
	{
		for(my $col = 1; $col <= @h_elements; $col++)
		{
			print "$col: @h_elements[$col-1]\n";
		}
		print "$question\n";
		## READ IN RESPONSE, and check if valid. set it as the column	
		my $answer = <STDIN>; chomp $answer;
		
		if($answer >= 1 && $answer <= @h_elements)
		{
			$column = $answer;
		}
		else
		{
			print "Invalid response\n\n";
		}
	}

	return $column - 1; # needs to be zero based, not one based like printed out to the user.
}

sub load_data
{
	my $file = shift;
	open(IN,$file) or die "ERROR!!! Relatedness input file $file cannot be read in; $!\n";
	$OUTFILE_HEADER = <IN>; ## skip header
	$network_ctr = 0;
	while(my $line = <IN>)
	{
		$line =~ s/^\s+//;
		chomp($line);
		if($line =~ /^FID/){$OUTFILE_HEADER = $line; next;} ## skip header
		my @temp = split(/\s+/,$line);
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

		my $PI_HAT = @temp[$RELATEDNESS_COLUMN];
		if($PI_HAT > $THRESHOLD)
		{
			$id_id_scores{"$name1\;$name2"} = $PI_HAT;
			$id_id_scores{"$name2\;$name1"} = $PI_HAT;
		}
		$id_id_all_info{"$name1\;$name2"} = $line;
		$id_id_all_info{"$name2\;$name1"} = $line;

		#$iid_to_fid{$IID1} = $FID1;
		#$iid_to_fid{$IID2} = $FID2;
		
		if(!exists $id_network{$name1})
		{
			$id_network{$name1} = $network_ctr;
			push @{ $networks{$network_ctr} }, "$name1";
			$network_ctr++
		}
		if(!exists $id_network{$name2})
		{
			$id_network{$name2} = $network_ctr;
			push @{ $networks{$network_ctr} }, "$name2";
			$network_ctr++
		}
	}
	close(IN);
}

sub load_trait_data
{
	foreach my $file (@trait_order)
	{
		#print "Loading trait data: $file\n";
		my $trait_type = $trait_files{$file};
		my %trait_hash;
		
		## input the trait size values of 1 for each ID from the input file
		if($file =~ /size/i)
		{
			foreach my $ID (keys %id_network)
			{
				$trait_hash{$ID} = 1;
			}
			push(@trait_refs, \%trait_hash); ## The hash references are loaded onto this array in the order of selection
			next;
		}
		open(IN,$file) or die "Trait file $file failed to open; $!\n";
		while(my $line = <IN>)
		{
			$line =~ s/^\s+//; chomp($line);
			my @temp = split(/\s+/,$line);

			my $FID_COLUMN = $TRAIT_FID_COLUMNS{$file};
			my $IID_COLUMN = $TRAIT_ID_COLUMNS{$file};
			my $FID = @temp[$FID_COLUMN];
			my $IID = @temp[$IID_COLUMN];
      #my $ID = "$FID**$IID";
			my $ID = "$IID";
			my $DATA_COLUMN = $TRAIT_DATA_COLUMNS{$file};
			my $trait_val = @temp[$DATA_COLUMN];

			if($trait_val == $EXCLUDE_VALUE)
			{
				$trait_val = "NA";
			}
			elsif($trait_type =~ /btrait/i)
			{
				if($trait_val == 2)
				{
					$trait_val = 1;
				}
				elsif($trait_val == 1)
				{
					$trait_val = 0;
				}
				else
				{
					die "Binary trait file $file contains a binary value other than 1 or 2 in the data column $DATA_COLUMN:\n\t$line\n";
				}
			}
			$trait_hash{$ID} = $trait_val;
		}
		close(IN);

		if($trait_type =~ /(mean|tail|\d+)/)
		{
			if($trait_type =~ /btrait/i)
			{
				die "INVALID TRAIT OPTION: $trait_type for file $file.\nOnly options 'high' and 'low' are valid with binary traits\n\n";
			}

			my $selection_val = fold_trait_data($trait_type,\%trait_hash);
			print "Selection value for $file: ". $selection_val ."\n";
		}
		push(@trait_refs, \%trait_hash); ## The hash references are loaded onto this array in the order of selection
	}
}

## Method used for doing weighting of tails, mean, or user specified value. It will adjust all the trait values to be the distance from the mean or the user specified value. This is NOT the most effective way to do the tail weighting of skewed traits where the mean is not actually halfwaye between the highest and lowest value. A better way would be to find the highest and lowest value, and for each value change it to the min(max-val,val-min). This will set it to the distance from the max or min whichever it is closer to.
sub fold_trait_data
{
	my $trait_type = shift;
	my $hash_ref = shift;

	my $fold_value = 0;
	if($trait_type =~ /(mean|tail)/)
	{
		## set fold value to the mean of all trait values
		my $ctr = 0;
		my $sum = 0;
		foreach(keys %$hash_ref)
		{
			$ctr++;
			$sum = $sum + $$hash_ref{$_};
		}
		$fold_value = $sum/$ctr;
	}
	if($trait_type =~ /\d+/)
	{
		## up out the fold_value from the trait type
		$trait_type =~ /-?(-?\.?\d+)/;
		$fold_value = $1;
	}


	## subtract $fold_value from each trait value and then take the absolute value, resulting in possitive distance from the fold value
	foreach(keys %$hash_ref)
	{
		my $old_val = $$hash_ref{$_};
		my $new_val = abs($old_val - $fold_value);
		$$hash_ref{$_} = $new_val;
	}
	return $fold_value;
}

## Network processing subroutines
sub colapse_networks
{
	# print all unique entries
	#foreach my $network (sort {$a <=> $b} keys %networks)
	#{
	#	my @temp = @{ $networks{$network} };
	#}
	my $relationships_ref;
	if($do_PR)
	{
		eval
		{
			($relationships_ref) = PRIMUS::predict_relationships_2D::get_relationship_likelihood_vectors($IBD_file_ref,$MIN_LIKELIHOOD,0,$lib_dir,$output_dir);
		};
		if($@)
		{
			print "ERROR: IBD ESTIMATES are not good enough. $@\n";
			die;
		};
	}
	
	foreach my $id_pair (keys %id_id_scores)
	{
		my ($id1,$id2) = split (/\;/,$id_pair);
		my $score = $id_id_scores{$id_pair};
		my $id1_network = $id_network{$id1};
		my $id2_network = $id_network{$id2};
	
		#if unrelated, do nothing
		if(!$do_PR)
		{
			if($score <= $THRESHOLD)
			{
				next;
			}
		}
		else
		{
      #$id1 =~ s/\*\*/__/;
      #$id2 =~ s/\*\*/__/;
			foreach(keys %$relationships_ref)
			{
				#print "key: $_\n";
			}
			if(!exists $$relationships_ref{$id1}{$id2})
			{
				my $temp = $id1;
				$id1 = $id2;
				$id2 = $temp;
			}
			#print "id1: $id1; id2: $id2\n";
			my @vector = @{$$relationships_ref{$id1}{$id2} };
			my @possibilities = PRIMUS::predict_relationships_2D::predict_relationship(@vector);
			## Check that there is at least one.
			if(@possibilities < 1)
			{
				
			}
			## If that one is UN then treat as unrelated
			if(@possibilities == 1)
			{
				if($possibilities[0] eq "UN")
				{
					next;
				}
			}
		}
		
		#if already in the same network, do nothing
		if($id1_network eq $id2_network)
		{
			next;
		}

		# combine $id1_network and $id2_network
		my @network1 = @{ $networks{$id1_network} };
		my @network2 = @{ $networks{$id2_network} };
		my @new_array = (@network1,@network2);
		$networks{$id1_network} = [@new_array];
		delete $networks{$id2_network};

		foreach my $id (@network2)
		{
			my $old_network = $id_network{$id};
			if($old_network ne $id2_network)
			{
				print "ERROR!!! $id claims to be in $old_network; actually in $id2_network\n";
				exit;
			}
			$id_network{$id} = $id1_network;
		}
	}
}

sub write_out_networks
{
	my $num_networks = keys %networks;
	my $file = shift;
	
	open(NETWORKS_OUT, ">$output_dir/$file\_networks") or die "Can't write to $output_dir/$file\_networks; $!\n";
	open(UNREL_OUT, ">$output_dir/$file\_unrelated_samples.txt") or die "Can't write to $output_dir/$file\_unrelateds; $!\n";
	print NETWORKS_OUT "Network\t$OUTFILE_HEADER";
	print UNREL_OUT "FID\tIID\n";
	my $network_ctr = 1;
	foreach my $network (sort {$a <=> $b} keys %networks)
	{	
		
		my @temp = @{ $networks{$network} };
		if(@temp < 2)
		{
			my $name = @temp[0];
			$name =~ s/\*\*/\t/;
			print UNREL_OUT "$name\n";
			next;
		}
		if(@temp > 4)
		{
			write_out_dot_file($file,$network,$network_ctr);
		}
		## This will write out the .genome file for each network (if a .genome file was read in)
		open(GENOMES_OUT, ">$output_dir/$file\_network$network_ctr.genome") or die "Can't write to $output_dir/$file\_networks; $!\n";
		print GENOMES_OUT "$OUTFILE_HEADER";

		for my $i (0 .. @temp-1)
		{
			for(my $j = 0; $j < $i; $j++)
			{
				my $id1 = $temp[$i];
				my $id2 = $temp[$j];
				my $info = $id_id_all_info{"$id1\;$id2"};
				if($info ne "")
				{
					print GENOMES_OUT "$info\n";
				}
				my $score = $id_id_scores{"$id1\;$id2"};
				if($score > $THRESHOLD)
				{
					print NETWORKS_OUT "$network_ctr\t$info\n";
				}
			}
		}
		$network_ctr++;
	}
	close(NETWORKS_OUT);
	close(UNREL_OUT);

}

sub get_connectedness
{
	my $network_ref = shift;
	my $num_connections;
	my %visited;
	my @ids = keys %$network_ref;
	my $num_ids = keys %$network_ref;
	my $max_connections = 0;
	for(my $i = 0; $i < $num_ids; $i++)
	{
		for(my $j = $i + 1; $j < $num_ids; $j++)
		{
			my $key = @ids[$i];
			my $key2 = @ids[$j];
			$max_connections++;
			my $PI_HAT = $id_id_scores{"$key\;$key2"};
			if($PI_HAT > $THRESHOLD)
			{
				$num_connections++;
			}
		}
	}
	if($max_connections eq 0)
	{
		return 1;
	}
	my $connectedness = $num_connections/$max_connections;
	return $connectedness;
}

sub breakup_large_networks
{
	foreach my $network (sort {$a <=> $b} keys %networks)
	{
		my @temp = @{ $networks{$network} };
		my %P = map { $_ => 1} @temp;	

		if(keys %P > $LOWEST_MAX_NETWORK_SIZE)
		{
			my $connectedness = get_connectedness(\%P);
			my $size = keys %P;
			my $MAX_NETWORK_SIZE = get_max_network_size($connectedness,$size);
			if($size > $MAX_NETWORK_SIZE)
			{
				print "WARNING!!! ". @temp . " NODES WITH CONNECTIVITY OF $connectedness WILL TAKE TOO LONG TO RUN; USING NEXT BEST SOLUTION.\n";
				breakup_large_network(\%P,$MAX_NETWORK_SIZE);
				@{ $networks{$network} } = keys %P;
			}
		}
	}
	

}

sub breakup_large_network
{
	my $P_ref = shift;
	my $MAX_NETWORK_SIZE = shift;
	my %degrees;
	my %neighbors;

	load_degrees_and_neighbors($P_ref,\%degrees,\%neighbors);
	
	while (keys %$P_ref > $MAX_NETWORK_SIZE)
	{
		my ($node_to_remove,$degree) = get_highest_degree_node(\%degrees);
		if($degree > 0)
		{
			delete $$P_ref{$node_to_remove};
			delete $degrees{$node_to_remove};
			reduce_neighbors($node_to_remove,\%neighbors,\%degrees);
			delete $neighbors{$node_to_remove};
			# Get connected components of new graph
			my @component_refs = get_connected_components($P_ref,\%neighbors);
			%{ $P_ref} = %{$component_refs[0]};
			shift(@component_refs);
			# Add smaller connected components to end of $networks{$network_ctr}
			#foreach my $hash_ref (sort {keys %{ $component_refs[$a] } <=> keys %{ $component_refs[$b]} } @component_refs)
			foreach my $hash_ref (@component_refs)
			{
				my @temp = keys %$hash_ref;
				foreach(@temp)
				{
					delete $degrees{$_};
					delete $neighbors{$_};
				}
				
				## if still too big, call this routine recursively
				if(keys %$hash_ref > $LOWEST_MAX_NETWORK_SIZE)
				{
					my $connectedness = get_connectedness($hash_ref);
					my $LOCAL_MAX_NETWORK_SIZE = get_max_network_size($connectedness);
					print "Connectedness: $connectedness\n";
					print "LOCAL_MAX_NETWORK_SIZE: $LOCAL_MAX_NETWORK_SIZE\n";
					if(keys %$hash_ref > $LOCAL_MAX_NETWORK_SIZE)
					{
						## DECEND RECURSIVELY
						breakup_large_network($hash_ref,$LOCAL_MAX_NETWORK_SIZE);
					}
				}
				$network_ctr++;
				@{ $networks{$network_ctr} } = keys %$hash_ref;
			}
			# proceed with larger network %P until it is smaller than its $MAX_NETWORK_SIZE for its connectivity
		}
		else
		{
			# It should never get here
			die "ERROR!!! PRUNING NODE WITHOUT RELATIVES!!!\n";
		}
		my $connectedness = get_connectedness($P_ref);
		$MAX_NETWORK_SIZE = get_max_network_size($connectedness);
		#print "Connectedness: $connectedness\n";
		#print "MAX_NETWORK_SIZE: $MAX_NETWORK_SIZE\n";
	}
}

sub get_connected_components
{
	my $network_ref = shift;
	my $neighbors_ref = shift;
	my @component_refs;
	my %remaining_ids = %$network_ref;

	## This loop will cycle through all connected components;
	while (keys %remaining_ids > 0)
	{
		## For this node, identify the connected component, and remove all nodes in the connected component from %remaining_ids
		my @temp = keys %remaining_ids;	
		my $node = @temp[0];
		my %connected_component;
		my %waiting_to_visit;
		$waiting_to_visit{$node} = 1;

		while(keys %waiting_to_visit > 0)
		{
			foreach my $node (keys %waiting_to_visit)
			{
				$connected_component{$node} = 1;
				delete $waiting_to_visit{$node};
				delete $remaining_ids{$node};
				my @neighbors = keys %{ $neighbors_ref->{$node}};
				foreach(@neighbors)
				{
					if(!exists $connected_component{$_})
					{
						$waiting_to_visit{$_}=1;
					}
				}
			}
		}
		push(@component_refs, \%connected_component);
	}
	return @component_refs;
}

sub write_out_independent_set
{
	my $file = shift;
	my %unrelated_set;
	open (UNIQUE_OUT, ">$output_dir/$file\_maximum_independent_set_PRIMUS");
	print UNIQUE_OUT "FID\tIID\n";
	## get most unrelateds in each network with Bron-Kerbosch Algorithm
	foreach my $network (sort {$a <=> $b} keys %networks)
	{
		my @temp = @{ $networks{$network} };
		
		my %P = map { $_ => 1} @temp;	
		my %R;
		my %X;
		my @maximal_cliques;
		my %degrees;
		my %neighbors;
	
		my $nodes_visited = 0;
		if(@temp > 45)
		{
			print "Running BronKerbosh for network $network (size = ".@temp.")\n";
		}
		BronKerbosh(\@maximal_cliques,\%R,\%P,\%X,\$nodes_visited);
		
		my $maximum_clique = get_maximum_clique(@maximal_cliques);
		my @maximum_ids = keys %{ $maximal_cliques[$maximum_clique] };
		write_out_maximum_clique_ids(@maximum_ids);
		foreach(@maximum_ids){$unrelated_set{$_}=1;}
	}
	close(UNIQUE_OUT);
	return %unrelated_set;
}

sub load_degrees_and_neighbors
{
	my $network_ref =shift;
	my $degree_ref = shift;
	my $neighbors_ref = shift;
	
	foreach my $node (keys %$network_ref)
	{
		my %neighbors = get_actual_neighbors($node,$network_ref);
		
		$$neighbors_ref{$node} = \%neighbors; #get_actual_neighbors($node,$network_ref);
		my @temp = keys %{ $neighbors_ref->{$node}}; # %neighbors;
		my $degree = @temp;
		$$degree_ref{$node} = $degree;
	}
}

sub reduce_neighbors
{
	my $node = shift;
	my $neighbors_ref = shift;
	my $degree_ref = shift;
	foreach my $neighbor (keys %{ $$neighbors_ref{$node}})
	{
		$$degree_ref{$neighbor}--;
		delete $$neighbors_ref{$neighbor}{$node};
	}
}

sub get_highest_degree_node
{
	my $degrees_ref = shift;
	my $max_degree = -1;
	my $max_node;
	my @max_trait_values = qw(NA);

	foreach my $node (keys %$degrees_ref)
	{
		my $degree = $$degrees_ref{$node};
		
		if($degree > $max_degree)
		{
			$max_degree = $degree;
			$max_node = $node;
		}
		elsif($degree == $max_degree)
		{
			my @temp_values;
			for(my $i = 0; $i < @trait_order; $i++)
			{
				@temp_values[$i] = $trait_refs[$i]{$node};
			}
			
			my $max = weighted_comparison(\@max_trait_values, \@temp_values); ## return 1 for max, and 2 for temp
			## We want the lower of the two to remove; therefore if max has higher trait values, switch max to temp
			if($max == 1 || @max_trait_values[0] eq "NA")
			{
				$max_degree = $degree;
				$max_node = $node;
				@max_trait_values = @temp_values;
			}
		}
	}
	return ($max_node,$max_degree);
}

sub get_maximum_clique
{
	my @maximal_cliques = @_;
	my @max_values = qw(NA);
	my $maximum_clique;
	
	for(my $i = 0; $i < @maximal_cliques; $i++)
	{
		my @temp_values;
		my $ctr = 0;
		foreach my $id (keys %{ $maximal_cliques[$i] })
		{
			$ctr++;
			for(my $i = 0; $i < @trait_order; $i++)
			{
				my $trait_type = $trait_files{$trait_order[$i]};
				my $trait_val = $trait_refs[$i]{$id};
				if(exists $trait_refs[$i]{$id})
				{
					@temp_values[$i] += $trait_val;
				}
			}
		}
		
		## Average  non-binary trait values
		for(my $i = 0; $i < @trait_order; $i++)
		{
			my $trait_type = $trait_files{$trait_order[$i]};
			if($trait_type =~ /qtrait/i && $ctr ne 0)
			{
				@temp_values[$i] = @temp_values[$i]/$ctr;
			}
		}

		## Compare
		my $max = weighted_comparison(\@max_values, \@temp_values); ## returns 1 for max  and 2 for temp
		
		## Update max if necessary
		if($max eq 2)
		{
			$maximum_clique = $i;
			@max_values = @temp_values;
		}
	}
	return $maximum_clique;	
}

sub weighted_comparison
{
	my $a_ref = shift;
	my $b_ref = shift;

	if ($$a_ref[0] eq "NA")
	{
		return 2;
	}

	## This will loop through the traits in order and returns the proper number of array once a difference is found. 
	##If no differences, the first array # is returned.
	for(my $i = 0; $i < @$a_ref; $i++)
	{
		my $trait_type = $trait_files{$trait_order[$i]};
		if($trait_type =~ /high/i || $trait_type =~ /ped/i || $trait_type =~ /size/i || $trait_type =~ /tail/i)
		{
			if($$a_ref[$i] > $$b_ref[$i])
			{
				return 1;
			}
			if($$a_ref[$i] < $$b_ref[$i])
			{
				return 2;
			}
		}
		elsif($trait_type =~ /low/i || $trait_type =~ /mean/i || $trait_type =~ /\d+/i)
		{
			if($$a_ref[$i] < $$b_ref[$i])
			{
				return 1;	
			}
			if($$a_ref[$i] > $$b_ref[$i])
			{
				return 2;
			}
		}
	}
	return 1;
}

sub write_out_maximum_clique_ids
{
	my @clique = @_;

	foreach my $ID (@clique)
	{
		my ($FID,$IID) = split(/\*\*/,$ID);
		print UNIQUE_OUT "$FID\t$IID\n";
	}

}

### BronKerbosh Algorithm subroutines
sub BronKerbosh
{
	my $maximal_cliques_ref = shift;
	my $R_ref = shift;
	my $P_ref = shift;
	my $X_ref = shift;
	my $num_visited_ref = shift;
	$$num_visited_ref++;
	my $nodes_visited = 1;
	
	if(keys %$P_ref == 0 and %$X_ref == 0)
	{
		push @{ $maximal_cliques_ref}, { %$R_ref };
		return;
	}

	## Select pivot node from P union X that maximizes the cardinality of P intersection N(u);
	my ($u, %u_neighbors) = select_pivot($P_ref,$X_ref);
	
	foreach my $v(keys %$P_ref)
	{
		## Use the pivot node: continue only if $v is not a neighbor of the pivot $u
		if(exists $u_neighbors{$v}){next;}

		my %temp_R = %$R_ref;
		$temp_R{$v} = $$P_ref{$v}; # Load temp_R = R union v

		my %temp_P;
		get_inverse_neighbors($v,$P_ref,\%temp_P); # Load temp_P = P intersection N(v)
		my @temp_arr = keys %temp_P;

		my %temp_X;
		get_inverse_neighbors($v,$X_ref,\%temp_X); # Load temp_X = X intersect N(v)

		BronKerbosh($maximal_cliques_ref,\%temp_R,\%temp_P,\%temp_X,$num_visited_ref);
		$$X_ref{$v} = $$P_ref{$v};
		delete($$P_ref{$v});
	}
	return $nodes_visited;
}

## Select pivot node from P union X that maximizes the cardinality of P intersection N(u);
sub select_pivot
{
	my $P_ref = shift;
	my $X_ref = shift;
	my $u;
	my $max_size = -1;
	my %u_neighbors;
	
	foreach my $key(keys %$P_ref)
	{
		my %neighbors; 
		get_inverse_neighbors($key,$P_ref,\%neighbors);
		if ($max_size < keys %neighbors)
		{
			$max_size = keys %neighbors;
			$u = $key;
			%u_neighbors = %neighbors;
		}
	}
	return ($u,%u_neighbors);
}

sub get_actual_neighbors
{
	my $v = shift;
	my $hash_ref = shift;
	my %neighbors;

	foreach my $n_v (keys %$hash_ref)
	{
		if($n_v eq $v)
		{
			next;
		}
		my $score = $id_id_scores{"$v\;$n_v"};
		if($score > $THRESHOLD)
		{
			$neighbors{$n_v} = $$hash_ref{$n_v};
		}
	}
	
	return %neighbors;
}

sub get_inverse_neighbors
{
	my $v = shift;
	my $hash_ref = shift;
	my $neighbors = shift;

	foreach my $n_v(keys %$hash_ref)
	{
		
		if($n_v eq $v)
		{
			next;
		}
		my $score = $id_id_scores{"$v\;$n_v"};
		if($score <= $THRESHOLD)
		{
			$$neighbors{$n_v} = $$hash_ref{$n_v};
		}
	}
	return $neighbors;
}


### Write the .dot file that can be read into a graph visualization program like Graphviz
sub write_out_dot_file
{
	my $file = shift;
	my $network = shift;
	my $network_ctr = shift;

	open(GRAPH_OUT,">$output_dir/$file\_network$network_ctr.dot");
	print GRAPH_OUT "graph network$network_ctr {\n";
	print GRAPH_OUT "\tnode [shape=circle];\n\n";
	
	my @temp = @{ $networks{$network} };
	
	## Write out all connections
	for my $i (0 .. @temp-1)
	{
		for(my $j = 0; $j < $i; $j++)
		{
			my $id1 = $temp[$i];
			my $id2 = $temp[$j];
			my $info = $id_id_all_info{"$id1\;$id2"};
			my $score = $id_id_scores{"$id1\;$id2"};
			
			## Change the ** delimiter between FID and IID to and _
      #$id1 =~ s/\*\*/_/;
      #$id2 =~ s/\*\*/_/;
			
			if($score > $THRESHOLD)
			{
				my $color = "black";
				print GRAPH_OUT "\t\"$id1\" -- \"$id2\" [color=$color];\n";
			}
			elsif($network eq "all")
			{
				print GRAPH_OUT "\tnode \"$id1\";\n";
				print GRAPH_OUT "\tnode \"$id2\";\n";
			}
		}
	}
	print GRAPH_OUT "}";
	close(GRAPH_OUT);
}

### Subroutine is necessary if processing more than one unrelated file in a single run. Call this inbetween runs.
sub get_max_network_size
{
	my $connectedness = shift;
	if($connectedness <= .15){return 60;}
	if($connectedness <= .2){return 65;}
	if($connectedness <= .25){return 70;}
	if($connectedness <= .35){return 100;}
	if($connectedness <= .45){return 130;}
	if($connectedness <= .55){return 170;}
	if($connectedness <= .65){return 230;}
	if($connectedness <= .75){return 330;}
	if($connectedness <= .85){return 500;}
	if($connectedness <= .95){return 500;}
	return 500;
}



############ ALTERNATIVE METHOD SUBROUTINES #################################
sub compare_alternative_methods
{
	### Run alternative methods
	my $file = shift;
	my $PRIMUS_unrelated_set = shift;
	my %networks_alt = @_;
	my %KING_network = write_out_independent_set_KING($file,%networks_alt);

	breakup_large_networks_PLINK(\%networks_alt);
	
	my %PLINK_network = write_out_independent_set_PLINK($file,%networks_alt);
	

	### Compare KING, PLINK, and PRIMUS
	my $val = get_maximum_clique(\%PLINK_network,\%KING_network,$PRIMUS_unrelated_set);
	
	if($val eq 0)
	{
		## PLINK performed the best
		system("cp $output_dir/$file\_maximum_independent_set_PLINK $output_dir/$file\_maximum_independent_set");
	}
	elsif($val eq 1)
	{
		## KING performed the best
		system("cp $output_dir/$file\_maximum_independent_set_KING $output_dir/$file\_maximum_independent_set");
	}
	elsif($val eq 2)
	{
		## PRIMUS performed the best
		system("cp $output_dir/$file\_maximum_independent_set_PRIMUS $output_dir/$file\_maximum_independent_set");
	}
	else
	{
		die "ERROR!!!!!\n";
	}

	if($PRINT_ALTERNATE_RESULTS eq 0)
	{
		system("rm $output_dir/$file\_maximum_independent_set_PLINK");
		system("rm $output_dir/$file\_maximum_independent_set_KING");
		system("rm $output_dir/$file\_maximum_independent_set_PRIMUS");
	}
}



sub write_out_independent_set_KING
{
	my $file = shift;
	my %networks_king = @_;
	my %unrelated_set;
	open (UNIQUE_OUT, ">$output_dir/$file\_maximum_independent_set_KING");
	print UNIQUE_OUT "FID\tIID\ta_status\n";
	
	## get most unrelateds in each network
	foreach my $network (sort {$a <=> $b} keys %networks_king)
	{
		my @temp = @{ $networks_king{$network} };
		my %P = map { $_ => 1} @temp;
		my @maximal_cliques;
	
		King_method(\@maximal_cliques,\%P);

		my $maximum_clique = get_maximum_clique(@maximal_cliques);
		my @maximum_ids = keys %{ $maximal_cliques[$maximum_clique] };
		write_out_maximum_clique_ids(@maximum_ids);
		foreach(@maximum_ids){$unrelated_set{$_}=1;}
	}
	close(UNIQUE_OUT);
	return %unrelated_set;
}

sub King_method
{
	my $maximal_cliques_ref = shift;
	my $network_ref = shift;
	my %actual_degrees;
	my %actual_neighbors;
	my %inverse_degrees;
	my %inverse_neighbors;
	my %i_set = ();
	load_inverse_degrees_and_neighbors($network_ref,\%inverse_degrees,\%inverse_neighbors);
	load_degrees_and_neighbors($network_ref,\%actual_degrees,\%actual_neighbors);
	
	my @sorted_ids = (sort {$inverse_degrees{$b} cmp $inverse_degrees{$a} } keys %inverse_degrees);
	
	while(keys %inverse_degrees)
	{
		my ($id,$degree) = get_highest_degree_node(\%inverse_degrees);
		delete $inverse_degrees{$id};
		my $add = 1;
		foreach my $test_id (keys %i_set)
		{
			if(exists $actual_neighbors{$id}{$test_id})
			{
				$add = 0;
				last;
			}
		}
		if($add eq 1)
		{
			$i_set{$id} = 1;
		}
	}
	push(@$maximal_cliques_ref,\%i_set);
}

sub load_inverse_degrees_and_neighbors
{
	my $network_ref =shift;
	my $degree_ref = shift;
	my $neighbors_ref = shift;
	
	foreach my $node (keys %$network_ref)
	{
		my $neighbors = get_inverse_neighbors($node,$network_ref);
		
		$$neighbors_ref{$node} = $neighbors; 
		my @temp = keys %{ $neighbors_ref->{$node}}; # %neighbors;
		my $degree = @temp;
		$$degree_ref{$node} = $degree;
	}
}


sub write_out_independent_set_PLINK
{
	my $file = shift;
	my %networks = @_;
	my %unrelated_set;
	open (UNIQUE_OUT, ">$output_dir/$file\_maximum_independent_set_PLINK");
	print UNIQUE_OUT "FID\tIID\n";
	
	## get most unrelateds in each network
	foreach my $network (sort {$a <=> $b} keys %networks)
	{
		my @node = @{ $networks{$network} };
		
		if(@node > 1)
		{
			die "ERROR!!! ". @node . " NODES. Should only be 1\n";
		}
		write_out_maximum_clique_ids(@node);
		foreach(@node){$unrelated_set{$_} = 1;}
	}
	close(UNIQUE_OUT);
	return %unrelated_set;
}

sub breakup_large_networks_PLINK
{
	my $networks_ref = shift;
	foreach my $network (sort {$a <=> $b} keys %$networks_ref)
	{
		
		my @temp = @{ $$networks_ref{$network} };
		#print "\n\nNEXT NETWORK $network: @temp\n";
		
		my %P = map { $_ => 1} @temp;	
	
		if(keys %P > 0)
		{
			breakup_large_network_PLINK($networks_ref,\%P);
			@{ $$networks_ref{$network} } = keys %P;
		}
	}
}

sub breakup_large_network_PLINK
{
	my $networks_ref = shift;
	my $P_ref = shift;
	my %degrees;
	my %neighbors;

	load_degrees_and_neighbors($P_ref,\%degrees,\%neighbors);
	
	while (keys %$P_ref > 1)
	{
		my @P_nodes = keys %$P_ref;
		my @degrees = keys %degrees;
		my @neighbors = keys %neighbors;
		
		my ($node_to_remove,$degree) = get_node_to_remove(\%degrees,\%neighbors);
		if($degree > 0)
		{
			delete $$P_ref{$node_to_remove};
			reduce_neighbors($node_to_remove,\%neighbors,\%degrees);
			delete $degrees{$node_to_remove};
			delete $neighbors{$node_to_remove};
			# Get connected components of new graph
			my @component_refs = get_connected_components($P_ref,\%neighbors);
			%{ $P_ref} = %{$component_refs[0]};
			shift(@component_refs); #remove the first element from the $hash_ref
			
			foreach my $hash_ref (@component_refs)
			{
				my @temp = keys %$hash_ref;
				foreach(@temp)
				{
					delete $degrees{$_};
					delete $neighbors{$_};
				}
				
				## if still too big, call this routine recursively
				if(keys %$hash_ref > 1)
				{
					breakup_large_network_PLINK($networks_ref,$hash_ref);
				}
				$network_ctr++;
				@{ $$networks_ref{$network_ctr} } = keys %$hash_ref;
			}
			# proceed with network %P until it is smaller than $MAX_NETWORK_SIZE
		}
		else
		{
			# I should never get here
			die "ERROR!!! PRUNING NODE WITHOUT RELATIVES!!!\n";
		}
		#exit;
	}
	# Replace $networks{$network} with keys %P
}

sub get_node_to_remove
{
	my $degrees_ref = shift;
	my $neighbors_ref = shift;
	my @nodes = keys %$degrees_ref;
	my $node = @nodes[0];
	my $degree = $$degrees_ref{$node};
	if($degree == 0)
	{
		die "ERROR REMOVING A NODE WITHOUT NEIGHBORS $node\n";
	}
	my @neighbors = keys %{ $neighbors_ref->{$node}}; # %neighbors;
	my $neighbor = @neighbors[0];
	my $neighbor_degree = $$degrees_ref{$neighbor};
	my @node_values;
	my @neighbor_values;
	
	for(my $i = 0; $i < @trait_order; $i++)
	{
		@node_values[$i] = $trait_refs[$i]{$node};
	}
	for(my $i = 0; $i < @trait_order; $i++)
	{
		@neighbor_values[$i] = $trait_refs[$i]{$neighbor};
	}
	
	my $max = weighted_comparison(\@node_values, \@neighbor_values); ## return 1 for node, and 2 for neighbor
	if($max == 2)
	{
		return ($node,$degree);
	}
	return ($neighbor,$neighbor_degree);
}

1;
