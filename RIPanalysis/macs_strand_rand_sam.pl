#!opt/bin/perl -w
use strict;

# v2 contains a counter to ouptput # of reads assigned to each strand

#2/18/22 to randomize strand of samfile see if that improves MACS
# performance on RIP peak calls

my ($in, $out)=@ARGV;

open (IN, "$in") or die "cant open $in in\n";
open (OUT, ">$out") or die "cant open $out in\n";

my ($pos, $neg);

while (my $line=<IN>) {
    chomp $line; #  remove the last trailing newline from the input string

    if ($line=~/^@/) {
	print OUT "$line\n"; # print header
    } else {
	my ($id, $strand, @array)=split(/\t/, $line);
	
	my $rand=rand(2);
	
	my $nstrand;
	if ($rand>1) {
	    $nstrand="16";
	    $neg++;
	} else {
	    $nstrand="0";
	    $pos++;
	}
	
	my $nl="$id\t$nstrand\t";
	foreach my $val (@array) {
	    $nl.="$val\t";
	}
	chop $nl;
	print OUT "$nl\n";
    }
}

print "$pos reads assigned to positive strand, $neg to negative\n";
