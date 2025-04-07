#! /usr/bin/perl

use strict;

use Cwd 'abs_path';

#print "run_PADRE.pl\n";

my $lib_dir;
my $bin_dir;

## Get full paths to PADRE directories so PADRE can be run from anywhere
my $run_PADRE_remote_path = abs_path($0);

if($run_PADRE_remote_path =~ m/(.*)\/bin\/run_PADRE.pl/)
{
	$lib_dir = "$1/lib";
	$bin_dir = "$1/bin";
}
else
{
	die "Change this script name back to \"run_PADRE.pl\" or change this script so it will work with the new name\n";
}

## Set environmental variables
#print "env: $ENV{'PERL5LIB'}\n";
$ENV{'PERL5LIB'} = "$ENV{'PERL5LIB'}:$lib_dir/perl_modules:$lib_dir/perl_modules/PADRE";
#print "env: $ENV{'PERL5LIB'}\n";

$ENV{'PERL5LIB'} =~ s/5\.14\.2/5\.10\.1/g;
#print "env: $ENV{'PERL5LIB'}\n";

## RUN PADRE SCRIPT
#print "$bin_dir/padre_kickoff.pl @ARGV --bin $bin_dir --lib $lib_dir\n";
system("$bin_dir/padre_kickoff.pl @ARGV --bin $bin_dir --lib $lib_dir");


