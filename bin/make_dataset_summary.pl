#!/usr/bin/env perl

use strict;

my $results_dir = shift; ## ../data/JHS/justJHSARICnoAIMS
my $dataset_name = shift; ## justJHSARICnoAIMS
my @lines;

if($results_dir eq "")
{
	print "\n\n $0 [dataset, e.g. JHS, mexam, hapmap3]\n\n";
	exit;
}

#my $results_dir = "JHS";

if($results_dir eq "reichenberger")
{
	#print "doing JHS\n";
	$results_dir = "../data/Reichenberger4/";
	$dataset_name = "justReichenberger.genome";
}
if($results_dir eq "yale")
{
	my $t = $dataset_name;
	#print "doing JHS\n";
	$results_dir = "../data/yale/justYale011412_v2_$t";
	$dataset_name = "justYale011412.genome";
}
if($results_dir eq "jarvick")
{
	my $t = $dataset_name;
	#print "doing JHS\n";
	$results_dir = "../data/jarvik/justJarvick_v2_$t";
	$dataset_name = "justJarvick.genome";
}
if($results_dir eq "JHS")
{
	my $t = $dataset_name;
	#print "doing JHS\n";
	$results_dir = "../data/JHS/justJHSARICnoAIMS_v7_$t";
	$dataset_name = "justJHSARICnoAIMS.genome";
}
if($results_dir eq "mexam")
{
	#print "doing mexam\n";
	$results_dir = "../data/mexam/mexamCC_noaims_.3";
	$dataset_name = "mexamCC_noaims.genome";
}
if($results_dir eq "hapmap3")
{
	my $t = "$dataset_name";
	
	my @hapmap_pops = ('ASW','CEU','CHB','CHD','GIH','JPT','LWK','MEX','MKK','TSI','YRI');
	#my @hapmap_pops = ('ASW');
	foreach(@hapmap_pops)
	{
		$results_dir = "../data/hapmap3_v8_$t/$_";
		$dataset_name = "1kLDpruned_hm3pop$_.genome";
		print "res: $results_dir; name: $dataset_name\n";
		system("./make_dataset_summary.pl $results_dir $dataset_name");
	}

	#system("cp ../data/hapmap3/@hapmap_pops[0]/Summary_@hapmap_pops[0].txt ../data/hapmap3/Summary_Hapmap3.txt");
	my $files = "";
	
	foreach(@hapmap_pops)
	{
		$files .= "../data/hapmap3_v8_$t/$_/Summary_hapmap3_r2_b36_fwd.$_.qc.poly.genome.txt ";
	}
	print "$files\n";
	system("cat $files > ../data/hapmap3_v8_$t/Summary_Hapmap3.txt");
}

my $continue = 1;
my $curr_network = 0;
my $last_network = 0;
my $num_networks = 0;

while($continue)
{
	if($curr_network > $last_network + 1000)
	{
		$continue = 0;
		next;
	}
	my $dir = "$results_dir/$dataset_name\_network$curr_network";
	my $file = "$dir/Summary_$dataset_name\_network$curr_network.txt";
	if(!-d $dir)
	{
    #print "$file\n";
		$curr_network++;
		next;
	}
	if(!-e $file)
	{
		my $new_line = "$curr_network\t$dataset_name\_network$curr_network\t0\t?\tNA\tNA\t?";
		push(@lines,$new_line);
		$curr_network++;
		next;
	}
	$num_networks++;
	$last_network = $curr_network;
	my $new_line = "$curr_network";
	open(IN,$file) or die "can't open $file; $!\n";
	while(my $line = <IN>)
	{
		chomp($line);
		if($line !~ /^##/)
		{
			last;
		}
		else
		{
			my @temp = split(/\s+/,$line);
			shift(@temp);
			shift(@temp);
			$new_line .= "\t@temp";
		}
	}
	push(@lines,$new_line);
	close(IN);
	$curr_network++;
}
open(OUT,">$results_dir/Summary_$dataset_name.txt");
print OUT "family_name: $dataset_name\n";
print OUT "num_networks: $num_networks\n";
print OUT "\n";
print OUT "Network_num\tnetwork_directory\tnum_possible_pedigrees\tnum_unflagged_pedigrees\tnum_samples\tnum_at_max_score\tSamples_IIDs\tError_messages\n";
foreach my $line (@lines)
{
	print OUT "$line\n";
}
print OUT "";
print OUT "";
print OUT "";
print OUT "\n\n\n\n\n\n";
close(OUT);
