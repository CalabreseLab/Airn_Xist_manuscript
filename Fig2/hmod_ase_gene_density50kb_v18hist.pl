#!opt/bin/perl -w
use strict;
use POSIX;

#9/18/24
# v18 requires total read num of dataset just for good measure

#9/11/24
# v17 prints out number of genes and does
# simple normalization by snp count in print out as well

#4/14/24
# v16 reads in input from extract_alignable_tss.pl (multiple tss per gene)
# and prints row by row data for each tss along with tss # per gene,
# so that counts per tss can be divided by the number of tss per gene

#4/8/24
# update of v9, to print out snp and read data per gene per bin rather than sum
# metagenes

# sam file cannot have a header
# snp file is output of parse_alignable_snps.pl
# tsslist is output of find_bestguess_tss4.pl
# or extract_alignable_tss.pl

# count of snps is based on alignablity score from parse_alignable_snps.pl

my ($in_snps, $in_sam, $in_tsslist, $total_readnum, $outname)=@ARGV;

my ($snplochoa, $snpidhash) =readsnps($in_snps);
my $readhoa=readsam($in_sam, $snpidhash);
my ($tsshoa, $tsscount)=readtss($in_tsslist);

my $tssisnps= snpintersect($snplochoa, $tsshoa);
my $tssireads= readintersect($readhoa, $tsshoa);

printhist($tsshoa, $tsscount, $tssisnps, $tssireads, $total_readnum, $outname);



### SUBS ###

sub readsnps {
    print "reading snps...\n";
    my ($in)=@_;
    open (IN, "$in") or die "Cant open infile $in\n";
    
    my %snplocs;
    my %snpids;
    my $i;

    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	
	my $snpid=$array[0];
	my $chr=$array[1];
	my $start=$array[2];
	my $al=$array[10];

	#print "$snpid, $start, $al\n";

	#store id
	$snpids{$snpid}=1;
	
	#store position and alignablity score  by bin
	#store 3rd value to be compatible with sam storage format
	my $pos=int ($start/50000);
	@{$snplocs{$chr}{'snps'}}[$pos].="$start,$al,NA\t";
	$i++;
	
    }
    print "$i snps counted\n";
    return \%snplocs, \%snpids;

}

    
sub readsam {
    print "reading reads...\n";
    my ($in, $snps)=@_;
    open (IN, "$in") or die "Cant open infile $in\n";
    
    my %reads;
    my $i;
    
    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);
	
	my $chr=$array[2];
	my $strand=$array[3];
	my $start=$array[4];
	my $end=$array[5];
	my $lsnpid=$array[7];
	$lsnpid=~/(snpid\S+)\|/;
	my $snpid=$1;
	my $genmatch=$array[15];

	#print "$chr, $strand, $start, $end, $snpid, $genmatch\n";

	#store a value of 0.5 for start end end of  each read to make the file compatible
	# with snp storage format
	
	my $pos=int ($start/50000);
	@{$reads{$chr}{$genmatch}}[$pos].="$start,0.5,$strand\t";
	@{$reads{$chr}{$genmatch}}[$pos].="$end,0.5,$strand\t";
	$i++;
	 
    }
    print "$i reads counted\n";
    return \%reads;
}


sub readtss {
    print "reading tss list...\n";
    my $in = shift;
    open (IN, "$in") or die "no infile $in given\n";

    my %tss; #store tss per gene
    my %tcount; #store # of tss per gene
    my $i;
    
    while (my $line=<IN>) {
	chomp $line;
	my @array=split(/\t/, $line);

	my $lid=$array[0];
	my $tss=$array[1];
	$lid=~/_(chr\S+)_\d+_\d+_\((\S)\)/;
	my $chr=$1;
	my $strand=$2;

	#print "$lid, $tss, $chr, $strand\n";

	#store unique tss 
	$tss{$chr}{$lid}{$tss}="$strand";
	$tcount{$lid}+=1;
	$i++;
    }
    print "$i tss counted\n";
    return \%tss, \%tcount;

}


sub snpintersect {
    my ($snps, $tsss)=@_;

    #genes that intersect reads
    my ($i, $j, $k, $l);

    #keep running tab of intersecting snps per gene
    my $snphisto;

    #mm10 chr lengths     
    my %lengths = (					
	'chr1' => '195471971',
	'chr2' => '182113224',
	'chr3' => '160039680',
	'chr4' => '156508116',
	'chr5' => '151834684',
	'chr6' => '149736546',
	'chr7' => '145441459',
	'chr8' => '129401213',
	'chr9' => '124595110',
	'chr10' => '130694993',
	'chr11' => '122082543',
	'chr12' => '120129022',
	'chr13' => '120421639',
	'chr14' => '124902244',
	'chr15' => '104043685',
	'chr16' => '98207768',
	'chr17' => '94987271',
	'chr18' => '90702639',
	'chr19' => '61431566',
	'chrX' => '171031299',
	#'chrY' => '91744698',
    );

    foreach my $chr (sort keys %{$tsss}) {

	foreach my $gene (sort keys %{$tsss->{$chr}}) {

	    foreach my $tss (sort keys %{$tsss->{$chr}{$gene}}) {
		
		my $tstrand=$tsss->{$chr}{$gene}{$tss};
	    
		#print "$gene, $chr, $tss, $tstrand\n";

		#retrieve the read list for array positions surrounding gene
		my $ref=floor $tss/50000;
		
		#populate arrays of nearby snps
		#array of potential interacting seqs
		my @nearbysnps;
		
		#define intervals to search
		my $a1pos=$ref-1;
		my $a2pos=$ref;
		my $a3pos=$ref+1;
		
		#if there are snps in any of the 3 bins surrounding the tss,
		# store the snp positions in a temporary list @nearbysnps
		if (defined $snps->{$chr}{'snps'}[$a1pos]  || 
		    defined $snps->{$chr}{'snps'}[$a2pos]  || 
		    defined $snps->{$chr}{'snps'}[$a3pos]) {
		    
		    my (@a1, @a2, @a3);
		    
		    if (defined $snps->{$chr}{'snps'}[$a1pos]) {
			@a1=split(/\t/, $snps->{$chr}{'snps'}[$a1pos]);
			push(@nearbysnps, @a1);
		    }
		    
		    if (defined $snps->{$chr}{'snps'}[$a2pos]) {
			@a2=split(/\t/, $snps->{$chr}{'snps'}[$a2pos]);
			push(@nearbysnps, @a2);
		    }
		    
		    if (defined $snps->{$chr}{'snps'}[$a3pos]) {
			@a3=split(/\t/, $snps->{$chr}{'snps'}[$a3pos]);
			push(@nearbysnps, @a3);
		    }
		    
		}
		
		#if $len >0,
		#calculate the distance of each snp from each tss and populate 
		# an array.
		my $len=@nearbysnps;
		
		#for snps, there is no b6/cast so define "genomematch" as snp
		my $genmatch="snps";
		
		if ($len>0) {
		    $i++;
		    ($snphisto)=histogram (\@nearbysnps, 
					   $tss, 
					   $tstrand,
					   $snphisto, 
					   $gene,
					   $genmatch);
		}
	    }
	}
    }
    
    return ($snphisto);
}

sub readintersect {
    my ($reads, $tsss)=@_;

    #genes that intersect reads
    my ($i, $j, $k, $l);

    #keep running tab of intersecting reads per gene
    my $readhisto;

    #mm10 chr lengths     
    my %lengths = (					
	'chr1' => '195471971',
	'chr2' => '182113224',
	'chr3' => '160039680',
	'chr4' => '156508116',
	'chr5' => '151834684',
	'chr6' => '149736546',
	'chr7' => '145441459',
	'chr8' => '129401213',
	'chr9' => '124595110',
	'chr10' => '130694993',
	'chr11' => '122082543',
	'chr12' => '120129022',
	'chr13' => '120421639',
	'chr14' => '124902244',
	'chr15' => '104043685',
	'chr16' => '98207768',
	'chr17' => '94987271',
	'chr18' => '90702639',
	'chr19' => '61431566',
	'chrX' => '171031299',
	#'chrY' => '91744698',
    );

    foreach my $chr (sort keys %{$tsss}) {

	foreach my $gene (sort keys %{$tsss->{$chr}}) {

	    foreach my $tss (sort keys %{$tsss->{$chr}{$gene}}) {
		
		my $tstrand=$tsss->{$chr}{$gene}{$tss};
	    	    
		#print "$gene, $chr, $tss, $tstrand\n";
		
		#retrieve the read list for array positions surrounding gene
		my $ref=floor $tss/50000;
		
		#define intervals to search
		my $a1pos=$ref-1;
		my $a2pos=$ref;
		my $a3pos=$ref+1;
		
		#store data one genome at a time
		my @genomes=('B6', 'CAST');
		foreach my $genmatch (@genomes) {
		    
		    #populate arrays of nearby snps
		    #array of potential interacting seqs for genmatch of interest
		    my @nearbyreads;
		    
		    #if there are snps in any of the 3 bins surrounding the tss,
		    # store the snp positions in a temporary list @nearbysnps
		    if (defined $reads->{$chr}{$genmatch}[$a1pos]  || 
			defined $reads->{$chr}{$genmatch}[$a2pos]  || 
			defined $reads->{$chr}{$genmatch}[$a3pos]) {
			
			my (@a1, @a2, @a3);
			
			if (defined $reads->{$chr}{$genmatch}[$a1pos]) {
			    @a1=split(/\t/, $reads->{$chr}{$genmatch}[$a1pos]);
			    push(@nearbyreads, @a1);
			}
			
			if (defined $reads->{$chr}{$genmatch}[$a2pos]) {
			    @a2=split(/\t/, $reads->{$chr}{$genmatch}[$a2pos]);
			    push(@nearbyreads, @a2);
			}
			
			if (defined $reads->{$chr}{$genmatch}[$a3pos]) {
			    @a3=split(/\t/, $reads->{$chr}{$genmatch}[$a3pos]);
			    push(@nearbyreads, @a3);
			}
			
			#if $len >0 (ie matching sequences in genmatch of interest),
			#calculate the distance of each snp from each tss and populate 
			# an array.
			my $len=@nearbyreads;
			
			if ($len>0) {
			    $i++;
			    ($readhisto)=histogram (\@nearbyreads, 
						    $tss, 
						    $tstrand,
						    $readhisto, 
						    $gene,
						    $genmatch);
			}
		    }
		}
	    }
	}
    }
        
    return ($readhisto);
}

#sub to populate histogram of snps and reads
sub histogram {
    my ($reads, $tss, $tstrand, $readhist, $gene, $genome)=@_;

    #go through each potential seq for an interection
    foreach my $hitinfo (@{$reads}) {

	my ($hit,$score,$rstrand)=split(/,/, $hitinfo);

	#print "hit $hit score $score strand $rstrand\n";
	#print "$hitinfo\n";
	
	#count read/snp if it falls within 50kb of tss
	my $flag=abs($tss-$hit);
	if ($flag<=50000) {
	    
	    #now define position in readhist array based on strand of gene/tss
	    if ($tstrand=~/\+/) {
		#if the read is upstream of start, position will be smaller than
		# start
		
		#find distance relative to start of $tss (50kb before gene)
		my $arraystart=$tss-50000;
		my $arraydist=$hit-$arraystart;
		my $tssbin=floor($arraydist/500);

		${$readhist}{$gene}{$tss}{$genome}[$tssbin]+=$score;

	    } else {
		#now upstream reads will have data that is greater than start
		my $arraystart=$tss+50000;
		my $arraydist=$arraystart-$hit;
		my $tssbin=floor($arraydist/500);

		${$readhist}{$gene}{$tss}{$genome}[$tssbin]+=$score;

	    }
	}
    }

    return ($readhist);
}


#sub to print out histogram info
sub printhist {
    my ($tsss, $tcounts, $snpdata, $readdata, $trn, $file)=@_;
    open (OUT, ">$file") or die "no outprefix given\n";

    #hist format row #, pos rel to tss, snpsum, b6-f/r, cast-f/r

    #declare and populate arrays for printing
    my (@snps, @b6, @cast, $genesum);

    foreach my $chr (sort keys %{$tsss}) {
	
	foreach my $gene (sort keys %{$tsss->{$chr}}) {

	    #keep tabs of total number of genes
	    $genesum++;
	    
	    #divide snps b6 cast counts by # of tss in each gene
	    # so that genes with multiple tss's do not contribute
	    # disproportionately to histogram
	    my $tcount=$tcounts->{$gene};

	    foreach my $tss (sort keys %{$tsss->{$chr}{$gene}}) {
		
		my $tstrand=$tsss->{$chr}{$gene}{$tss};

		#count snps, divide by tss num
		for (my $i='0'; $i<200; $i++) {
		    
		    my $snpcount=$snpdata->{$gene}{$tss}{'snps'}[$i];	
		    if (!defined $snpcount) {
			$snpcount="0";
		    }
		    
		    $snps[$i]+=$snpcount/$tcount;
		}
		
		#count b6, divide by tss num
		for (my $i='0'; $i<200; $i++) {
		    
		    my $readcount=$readdata->{$gene}{$tss}{'B6'}[$i];	
		    if (!defined $readcount) {
			$readcount="0";
		    }
		    $b6[$i]+=$readcount/$tcount;
		}
		
		#count cast, divide by tss num
		for (my $i='0'; $i<200; $i++) {
		    
		    my $readcount=$readdata->{$gene}{$tss}{'CAST'}[$i];	
		    if (!defined $readcount) {
			$readcount="0";
		    }
		    $cast[$i]+=$readcount/$tcount;
		}
	    }
	}	
    }

    #now go through arrays and print out histogram
    print OUT "row_id\trelpos\tgene_count\tsnp_count\tb6\tcast\tb6_by_snp_rpm\tcast_by_snp_rpm\n";
    
    for (my $i='0'; $i<200; $i++) {
	
	my $relpos=($i-100)*500;

	my $b6snl=$b6[$i]/$snps[$i]/$trn*1000000;
	my $b6sn=sprintf("%.4f", $b6snl);
	my $castsnl=$cast[$i]/$snps[$i]/$trn*1000000;
	my $castsn=sprintf("%.4f", $castsnl);
	
	print OUT "$i\t$relpos\t$genesum\t$snps[$i]\t$b6[$i]\t$cast[$i]\t$b6sn\t$castsn\n";
    }

    close OUT;
}
