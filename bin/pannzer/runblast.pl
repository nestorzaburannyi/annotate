# run blast and parse format to emulate SANS.tab

use strict;

$|=1;

if($#ARGV<2) { die "USAGE: $0 fastafile dbname num_threads\n"; }
my($fastafile,$db,$num_threads)=@ARGV;
my $printquery=0;

my $BLAST_EXE="/data/szchandr/program/ncbi-blast-2.2.31+/bin/blastp";
my $cmd="$BLAST_EXE -num_threads $num_threads -query $fastafile -db $db -outfmt \"7 qseqid sseqid qlen slen length pident evalue bitscore stitle\"";
my(@lines)=`$cmd`;

# memorize query sequences
my %seqs;
my %desc;
my $pid;
my $desc;
my $seq;
open(IN,"<$fastafile");
while(<IN>) {
	chomp;
	if(/^>(\S+)\s*(.*)$/) {
		if($pid) { $seq=~s/\W//g; $seqs{$pid}=$seq; $desc{$pid}=$desc; }
		$pid=$1;
		$desc=$2;
		$seq='';
	} else {
		$seq.=$_;
	}
}
close(IN);

# convert to SANS-tabular-format: pident->pide stitle->desc,species,genename, isquery, qseq
print join("\t","nid","isquery","qpid","spid","qcov","scov","bits","pide","lali","desc","species","qseq","genename"),"\n";
open(IN,"$cmd |");
my $oldqpid='';
my $nid=0;
while(<IN>) {
	next if(/^#/);
	chomp;
	my($qpid,$spid,$qlen,$slen,$lali,$pident,$evalue,$bits,$stitle)=split(/\t/);
	if($qpid ne $oldqpid) { # output query row
		$nid++;
		$oldqpid=$qpid;
		if ($printquery>0) {
			my($desc,$species,$genename)=&parse_header($desc{$qpid});
			print join("\t",$nid,1,$qpid,$spid,1.0,1.0,0.0,1.0,$qlen,$desc,$species,$seqs{$qpid},$genename),"\n";
		}
	} else { # output sbjct row
		my($desc,$species,$genename)=&parse_header($stitle);
		my $qcov=$lali/$qlen;
		my $scov=$lali/$slen;
		my $pide=$pident/100;
		print join("\t",$nid,0,$qpid,$spid,$qcov,$scov,$bits,$pide,$lali,$desc,$species,"n.d.",$genename),"\n";
	}
}
close(IN);

sub parse_header {
	chomp;
	my($stitle)=@_;
	$stitle=~s/^\S+\s+//; # remove pid
	my $desc=$stitle;
	$desc=~s/ \w{2}\=.*$//;
	$_=$stitle;
	my ($species)=/OS=(.*)$/;
	$species=~s/ \w{2}\=.*$//;
	$_=$stitle;
	my ($genename)=/GN=(.*)$/;
	$genename=~s/ \w{2}\=.*$//;
	return($desc,$species,$genename);
}
