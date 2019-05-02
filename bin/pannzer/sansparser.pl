# perl SANSparser.pl qpid spid qcov scov bitscore < sans.raw > sans.tab
# read raw SANS output
# write requested columns
#
# SANS tags: 
#	QUERY nid= vote_cutoff= LSEQ= 
#	SBJCT VOTE= TUPS= PIDE= LALI= BITS= EVALUE= DIAG= LSEQ= 
# computed:
# 	qpid, spid
# 	qcov,scov = LALI/LSEQ
# 	desc, species, gene
# 	rank
# 	qseq, seq

use strict;
use Switch;

if($#ARGV<0) { die "USAGE: $0 COLNAMES < raw-SANS-output\n\nValid COLNAMES: qpid spid qcov scov desc species gene rank qseq seq vote tups pide lali qfrom qto sfrom sto evalue diag lseq\n"; }

# desired colnames
my(@col)=@ARGV;

# raw SANS data read from STDIN
# query blocks separated by < /QUERY >
my @lines;
my @sbjct_list; # list of hashes of sbjct keys/values
my %query; # hash of query keys/values
while(<STDIN>) {
	chomp;
	push(@lines,$_);
	if(/<\/QUERY/) { 
		&table(@lines); 
		undef(@lines);
	}
}

###############################################################################

sub table {
	my(@lines)=@_;
	my($ierr)=&parse_result_keywords(@lines);
	last if($ierr != 0);
	next if($#sbjct_list < 0); # no hits
	# select key to sort sbjctlist 
	my $sortkey='BITS';
	if(!defined($sbjct_list[0]->{$sortkey})) { $sortkey='VOTE'; }
	# parse query protein identifier from header
	my $qpid=get_pid($query{'hdr'});
	my $qseqlen=$query{'LSEQ'};
	# print one line per sbjct
	my $rank=0;
	foreach my $rec (sort { $b->{$sortkey} <=> $a->{$sortkey} } @sbjct_list) {
		$rank++;
		my @row;
		my $lali=$rec->{'LALI'};
		foreach my $col (@col) {
			switch(lc($col)) {
				case 'qpid' { push(@row,$qpid); }
				case 'spid' { push(@row,get_pid($rec->{'hdr'})); }
				case 'qcov' { push(@row,$lali/$qseqlen); }
				case 'scov' { push(@row,$lali/$rec->{'LSEQ'});}
				case 'desc' { push(@row,get_desc($rec->{'hdr'})); }
				case 'species' { push(@row,get_species($rec->{'hdr'})); }
				case 'gene' { push(@row,get_gene($rec->{'hdr'})); }
				case 'rank' { push(@row,$rank); }
				case 'vote' { push(@row,$rec->{'VOTE'}); }
				case 'tups' { push(@row,$rec->{'TUPS'}); }
				case 'pide' { push(@row,$rec->{'PIDE'}); }
				case 'lali' { push(@row,$lali); }
				case 'bits' { push(@row,$rec->{'BITS'}); }
				case 'evalue' { push(@row,$rec->{'EVALUE'}); }
				case 'diag' { push(@row,$rec->{'DIAG'}); }
				case 'lseq' { push(@row,$rec->{'LSEQ'}); }
				case 'seq' { push(@row,$rec->{'seq'}); }
				case 'qseq' { push(@row,$query{'seq'}); }
				case 'qfrom' { push(@row,$rec->{'QFROM'}); }
				case 'qto' { push(@row,$rec->{'QTO'}); }
				case 'sfrom' { push(@row,$rec->{'SFROM'}); }
				case 'sto' { push(@row,$rec->{'STO'}); }
			}
		}
		# print tab-separated fields
		print join("\t",@row),"\n";
	}
}	

sub get_pid {
	$_=shift(@_);
	my($pid)=/^\s*\>(\S+)/;
	return $pid;
}

sub get_species {
	$_=shift(@_);
	my($species)='';
        if(/OS=(\S.+)\s\w{2}\=/) { $species=$1; $species =~ s/\w{2}\=.*$//; }
	return $species;
}

sub get_gene {
	$_=shift(@_);
	my $gene='';
        if(/GN=(\S+)/) { $gene=$1; }
	return $gene;
}

sub get_desc {
	$_=shift(@_);
        s/^\S+\s+//;
        my($desc)=split(/\s\w{2}\=\S/,$_);
	return($desc);
}

###############################################################################

sub parse_result_keywords {
        # returns %query, @sbjct_list == array of hashes with keywords
        my $inquery=0;
        my $insbjct=0;
        my $hdr='';
        @sbjct_list=();
        my $nsbjct=0;
        foreach(@_) {
                next if(/^\s*#/);
                if(/^\s*ERROR/) { print $_; return(1); }  # echo error messages
                if(/<QUERY/) {
                        $inquery=1;
			s/^\s*<QUERY\s+//;
                        my(@x)=split(/[= ]+/);
			while($#x>0) {
				my($key)=shift(@x);
				my($value)=shift(@x);
				$query{$key}=$value;
			}
                } elsif(/<\/QUERY>/) {
                        $inquery=0;
                } elsif(/<SBJCT/) {
                        $insbjct=1;
                        s/^.SBJCT\s*//;
                        my(@x)=split(/[= ]+/);
                        my $rec={};
                        while($#x>0) {
                                my($key)=shift(@x);
                                my($value)=shift(@x);
                                $rec->{$key}=$value; # keywords: VOTE, PIDE, LALI, LSEQ, TUPS, BITS, EVALUE
                        }
                        push(@sbjct_list, $rec);
                        $nsbjct=$#sbjct_list;
                } elsif(/<\/SBJCT/) {
                        $insbjct=0;
                } elsif(/^>/) {
                        if($insbjct) { $sbjct_list[$nsbjct]{'hdr'}=$_; } elsif($inquery) { $query{'hdr'}=$_; }
                } else {
                        if($insbjct) { $sbjct_list[$nsbjct]{'seq'}=$_; } elsif($inquery) { $query{'seq'}=$_; }
                }
        }
        return(0);
}

###############################################################################
