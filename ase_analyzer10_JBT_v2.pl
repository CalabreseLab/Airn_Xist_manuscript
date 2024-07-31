#!opt/bin/perl -w
use strict;

#2/7/23 JBT ("v2"):
#fixed a typo causing the column names of the the final output table to be out of the intended order
#used to be: print OUT "id\tchr\tstrand\tstart\tend\tB6_exon\tB6_ejc\tB6_intron\tB6_total\tcast_exon\tcast_ejc\tcast_intron\tcast_total\n"
#changed to: print OUT "id\tchr\tstart\tend\tstrand\tB6_exon\tB6_ejc\tB6_intron\tB6_total\tcast_exon\tcast_ejc\tcast_intron\tcast_total\n" 

#1/20/22
# v9 was 129/cast specific, now B6/cast specific
# also accepts single and pe
# also fixed error in pe intersection from v9
# use with snp_reads_final output file from intersect_reads_snps16.pl

#6/18/21
# substantial reworking -- designed for paired end reads, to count allelic
# splicejunction and non-junction overlapping reads
# takes output of intersect_reads_snps15.pl and gencode basic gtf file


#1/15/15
# modified from v6
# strand info still counts
# v8 tracks the total number of unique snpids counted per gene

#7/27/14
# no longer output rpm (from ase4)
# UDP method returns reverse strand of reality, take this into account when parsing.

#1/15 modified to output rpm in addition to reads
#10/29 modified to input my collapsed list of gene symbol coords
#10/26 to take strand into account
#10/24/09
#take output of intersect_reads_snps2.pl (list of snp-overlapping reads) 
# and intersect it with a list of genomic coordinates (like ref-seq genes)

my($in_gtf, $in_snp_reads, $peflag, $outpref)=@ARGV;

#pe flag must be y or n
if (($peflag ne "y") && ($peflag ne "n")) {
    die "pe flag must be y or n\n";
}

my ($genehoa, $exonhoa, $ejchoa, $genedata)=readgtf($in_gtf);
my ($genecounts)=intersect($genehoa, $exonhoa, $ejchoa, $in_snp_reads, $peflag);
printer($genecounts, $genedata, $outpref);

### SUBS ###

#take snp overlapping readand parse into 
# broken up into 10kb windows
sub readgtf {
    my ($in)=@_;
    open (IN, "$in") or die "Cant open infile $in\n";
    
    my %coords;
    my %exons;
    my %ejcs;
    my %genedata;
    
    <IN>;
    <IN>;
    <IN>;
    <IN>;
    <IN>;

    #string to store upstream exon end.
    my $eend;

    #string to store gene name for exons
    my $gene_name;
    
    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	my $chr=$array[0];
	my $type=$array[2];
	my $start=$array[3];
	my $end=$array[4];
	my $strand=$array[6];
	my $ginfo=$array[8];

	if ($type=~/gene/) {
	    my @garray=split(/;/,$ginfo);
	    my $lgn=$garray[2];
	    $lgn=~/"(\S+)"/;
	    my $gn="$1"."_"."$chr"."_"."$start"."_"."$end"."_"."($strand)";
	    $gene_name=$gn;
	    
	    #bins to add gene to
	    my $fbin=int ($start/5000);
	    my $ebin=int ($end/5000);
	    for (my $i=$fbin; $i<=$ebin; $i++) {
		@{$coords{$chr}}[$i].="$start\t$end\t$strand\t$gene_name\n";
	    }

	    #reset stored data
	    $eend=undef();

	    #print "$type:$start\t$end\t$strand\t$gene_name\n";

	    #store gene info for final file printing
	    $genedata{$gene_name}="$chr\t$start\t$end\t$strand";
	    
	} elsif ($type=~/transcript/) {
	    
	    #reset stored data
	    $eend=undef();

	} elsif ($type=~/exon/) {

	    #store pure exon coords
	    my $fbin=int ($start/5000);
	    my $ebin=int ($end/5000);
	    for (my $i=$fbin; $i<=$ebin; $i++) {
		@{$exons{$chr}}[$i].="$start\t$end\t$strand\t$gene_name\n";
	    }

	    #print "$type:$start\t$end\t$strand\t$gene_name\n";

	    #order of storing things differs depending on strand
	    if ($strand eq '+') {
		
		if (!defined $eend) {
		    #then this is a first exon, store the end
		    $eend=$end;
		} else {
		    #then this is at least a second exon and there exists a junction
		    my $estart=$start;

		    #store the junction
		    my $fbin=int ($eend/5000);
		    my $ebin=int ($estart/5000);
		    for (my $i=$fbin; $i<=$ebin; $i++) {
			@{$ejcs{$chr}}[$i].="$eend\t$estart\t$strand\t$gene_name\n";
		    }

		    #print "ejc pos:$eend\t$estart\t$strand\t$gene_name\n";
		    
		    #overwrite the eend for the next junction
		    $eend=$end;
		    
		}
		
	    } else {
		
		#this is for negative stranded transcrips, whose exons are listed
		# in the opposite direction -- last exon first in gtf
		
		if (!defined $eend) {
		    #then this is a first exon, store the end
		    $eend=$start;
		} else {
		    #then this is at least a second exon and there exists a junction
		    my $estart=$end;

		    #store the junction
		    my $fbin=int ($estart/5000);
		    my $ebin=int ($eend/5000);
		    for (my $i=$fbin; $i<=$ebin; $i++) {
			@{$ejcs{$chr}}[$i].="$estart\t$eend\t$strand\t$gene_name\n";
		    }

		    #print "ejc neg:$estart\t$eend\t$strand\t$gene_name\n";
		    
		    #overwrite the eend for the next junction
		    $eend=$start;
		    
		}
	    }
	}	
    }
    return \%coords, \%exons, \%ejcs, \%genedata;
    close IN;
}






# assign snp reads to genes
sub intersect {
    my ($genes, $exons, $ejcs, $in, $pe)=@_;
    
    open (IN, "$in") or die "Cant open infile $in\n";  

    my %counted; #hash to store reads that have been counted
    my %genecounts; #geneid->$gen->(intron|exon|ejc)=count

    my $tr; #all reads
    my $gr; #all reads that overlap genes
    
    <IN>;
    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	my $rid=$array[0];
	my $chr=$array[2];
	my $rstrand=$array[3];
	my ($rstart, $rend);

	$tr++;

	if ($pe=~/y/) {
	    
	    if ($rstrand eq '+') {
		$rstart=$array[4];
		$rend=$array[11];
	    } else {
		$rstart=$array[10];
		$rend=$array[5];
	    }
	} else {
	    $rstart=$array[4];
	    $rend=$array[5];
	}
	
	my $gen=$array[15];
	
	#define bin
	my $rsbin=int ($rstart/5000);
	my $rebin=int ($rend/5000);
	#define intervals to search
	my $binstart=$rsbin-1;
	my $binend=$rebin+1;
	
	#define array of potential interacting features -- genes, exons, ejcs
	my @igenes;
	my @iexons;
	my @iejcs;

	#if there is data in any of the bins surrounding the gene,
	# genes
	for (my $i=$binstart; $i<=$binend; $i++) {
	    if (defined $genes->{$chr}[$i]) {

		my @bininfo;

		#split the info in the bins by\n
		@bininfo=split(/\n/, $genes->{$chr}[$i]);
		push(@igenes, @bininfo);

	    }
	}
	# exons
	for (my $i=$binstart; $i<=$binend; $i++) {
	    if (defined $exons->{$chr}[$i]) {

		my @bininfo;

		#split the info in the bins by\n
		@bininfo=split(/\n/, $exons->{$chr}[$i]);
		push(@iexons, @bininfo);

	    }
	}
	# ejcs
	for (my $i=$binstart; $i<=$binend; $i++) {
	    if (defined $ejcs->{$chr}[$i]) {

		my @bininfo;

		#split the info in the bins by\n
		@bininfo=split(/\n/, $ejcs->{$chr}[$i]);
		push(@iejcs, @bininfo);

	    }
	}

	
	#now search if genes overlap read -- if yes, then need to decide on
	# exon vs ejc vs intron
	foreach my $geneline (@igenes) {

	    my @geneinfo=split(/\t/, $geneline);
	    my $gstart=$geneinfo[0];
	    my $gend=$geneinfo[1];
	    my $gstrand=$geneinfo[2];
	    my $gname=$geneinfo[3];

	    if (!defined $counted{$rid}) {
		#$counted{$rid}=1; #define the read as counted
			
		if ("$rstrand" ne "$gstrand") { #data are reverse stranded
		    
		    #determine whether it overlaps gene -- if so, then examine exons and ejcs
		    if (($rend>=$gstart && $rend<=$gend) ||
			($rstart>=$gstart && $rstart<=$gend)) {
			
			#then the read overlaps the gene in some capacity
			# but I only want to count each read once, so it will hit the first
			# overlapping gene/exon/ejc/intron. %counted		    
			

			#print "$geneline\n\t$rstart, $rend, $rstrand, $chr\n";
			
			#now search exons
			foreach my $exonline (@iexons) {
			    
			    my @exoninfo=split(/\t/, $exonline);
			    my $exonstart=$exoninfo[0];
			    my $exonend=$exoninfo[1];
			    my $exonstrand=$exoninfo[2];
			    my $gname=$exoninfo[3];
			    
			    if (!defined $counted{$rid}) {
				if ("$rstrand" ne "$exonstrand") { #data are reverse stranded
				    #determine whether it falls completely within the exon
				    if ($rstart>=$exonstart && $rend<=$exonend) {
					#print "\texon overlap: $chr:$rstart-$rend, $rstrand, $exonstart, $exonend, $exonstrand\n";
					$genecounts{$gname}{$gen}{'exon'}+=1;
					$counted{$rid}=1; #define the read as counted
					$gr++;
				    }
				}
			    }
			}
			
			#now if it didnt get defined, search ejcs
			if (!defined $counted{$rid}) {
			    foreach my $ejcline (@iejcs) {
				
				my @ejcinfo=split(/\t/, $ejcline);
				my $ejcstart=$ejcinfo[0];
				my $ejcend=$ejcinfo[1];
				my $ejcstrand=$ejcinfo[2];
				my $gname=$ejcinfo[3];
				
				if (!defined $counted{$rid}) {
				    if ("$rstrand" ne "$ejcstrand") { #data are reverse stranded
					#determine whether it completely encompasses the junction
					if ($rstart<=$ejcstart && $rend>=$ejcend) {
					    #print "\tejc overlap! $chr:$rstart-$rend, $rstrand, $ejcstart, $ejcend, $ejcstrand\n";
					    $genecounts{$gname}{$gen}{'ejc'}+=1;
					    $counted{$rid}=1; #define the read as counted
					    $gr++;
					}
				    }
				}
			    }
			}

			#now if it still wasnt counted, it falls in an intron!
			if (!defined $counted{$rid}) {
			    #print "\tits in an intron! $chr:$rstart-$rend, $rstrand\n";
			    $genecounts{$gname}{$gen}{'intron'}+=1;			    
			    $counted{$rid}=1; #define the read as counted
			    $gr++;
			}
		    }   					    
		}
	    }
	}
    }
    print "$gr reads overlap genes from $tr reads total\n";
    return \%genecounts;
}

sub printer{
    my ($gcounts, $geneinfo, $outname)=@_;
    open (OUT, ">$outname") or die "cant open outfile name $outname\n";
    
    print OUT "id\tchr\tstart\tend\tstrand\tB6_exon\tB6_ejc\tB6_intron\tB6_total\tcast_exon\tcast_ejc\tcast_intron\tcast_total\n";

    foreach my $gene (keys %{$geneinfo}) {

	my $b6exon="0";
	my $b6ejc="0";
	my $b6intron="0";
	my $b6tot="0";
	my $cexon="0";
	my $cejc="0";
	my $cintron="0";
	my $ctot="0";
	
	if (defined $gcounts->{$gene}) {
	    if (defined $gcounts->{$gene}{'B6'}{'exon'}) {
		$b6exon=$gcounts->{$gene}{'B6'}{'exon'};
		$b6tot+=$b6exon;
		
	    }
	    if (defined $gcounts->{$gene}{'B6'}{'ejc'}) {
		$b6ejc=$gcounts->{$gene}{'B6'}{'ejc'};
		$b6tot+=$b6ejc;

	    }
	    if (defined $gcounts->{$gene}{'B6'}{'intron'}) {
		$b6intron=$gcounts->{$gene}{'B6'}{'intron'};
		$b6tot+=$b6intron;

	    }
	    if (defined $gcounts->{$gene}{'CAST'}{'exon'}) {
		$cexon=$gcounts->{$gene}{'CAST'}{'exon'};
		$ctot+=$cexon;
		
	    }
	    if (defined $gcounts->{$gene}{'CAST'}{'ejc'}) {
		$cejc=$gcounts->{$gene}{'CAST'}{'ejc'};
		$ctot+=$cejc;
		
	    }
	    if (defined $gcounts->{$gene}{'CAST'}{'intron'}) {
		$cintron=$gcounts->{$gene}{'CAST'}{'intron'};
		$ctot+=$cintron;
		
	    }
	}

	print OUT "$gene\t$geneinfo->{$gene}\t$b6exon\t$b6ejc\t$b6intron\t$b6tot\t$cexon\t$cejc\t$cintron\t$ctot\n";	
	    
    }
}
