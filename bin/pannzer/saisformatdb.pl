# index database for SANS
# input: uniprot uniprot_sprot.fasta.gz uniprot_trembl.fasta.gz
# outputs: name.[i].psq name.[i].pin name.[i].phr name.[i].SAP name.[o].sres
# algo:
#       1. split headers, sequence, pointers (phr, psq, pin)
#       2. split database into 2 Gaa chunks
#       3. foreach chunk
#       	3.1 suftest(8 bytes)
#       	3.2 create ISA (8 bytes)
#       	3.3 . create SAP, sres (4 + 2 bytes)

use strict;
use lib "/annotate/lib";

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check the SANS_HOME env setting or you modify it here directly
#___________________________________________________________________________
my $SANS_HOME="/annotate/bin/pannzer/";
if(exists $ENV{'SANS_HOME'}){
    $SANS_HOME=$ENV{'SANS_HOME'};
    if($SANS_HOME=~/\/ *$/){ chop($SANS_HOME) }
}else{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
    # The SANS home dir if unset in environment: set the variable
    # $HOMEDIR in the Makefile, which will replace the below default location
    # during `make install`.
    # Comment the second line and uncomment the first line to have
    # SANS_HOME defaulting to the current directory.
    #________________________________________________________________________
    #$SANS_HOME=$ENV{'PWD'};
    $SANS_HOME="/annotate/bin/pannzer/";
}

if($#ARGV<1) { die "USAGE: $0  db-name fastafiles\n"; }
my($project,@fastafiles)=@ARGV;
my $SUFTEST_EXE="$SANS_HOME\/sais8"; # 8 bytes SA
my $ISA_EXE="$SANS_HOME\/isa8"; # 8 bytes SA, 8 bytes ISA
my $SAP_EXE="$SANS_HOME\/sap8"; # 8 bytes ISA, 4 bytes SAP, 2 bytes sres
my $removegaps=0;
my $CHUNKSIZE=2000000000; # 2 Gb (not effective in splitter(1)
my $MAXRES=64000; # truncate longer

# create database.{phr|psq|pin}
my($n)=&splitter($project,1,$removegaps,@fastafiles);
warn "# splitter(1) returned $n\n";
# run sais,generate isa, sap,sres

my $cmd="$SUFTEST_EXE $project";
warn "$cmd\n";
&system_cmd($cmd);
my $cmd="$ISA_EXE $project";
warn "$cmd\n";
&system_cmd($cmd);
my $cmd="$SAP_EXE $project";
warn "$cmd\n";
&system_cmd($cmd);

exit();


sub splitter {
        my($project,$nosplit,$removegaps,@fastafiles)=@_;
        my $section=1;
        #&openfiles($project,$section,$nosplit);
        my $x="$project";
        if($nosplit<1) { $x.="\.$section"; }
        open(PHR,"\>$x\.phr");
        open(PSQ,"\>$x\.psq");
        open(PIN,"\>$x\.pin");

        my $ptr=0;
        my $hdrptr=0;
        my $seq='';
        my $header='';

        my $compressed=0;
        foreach(@fastafiles) {
                if(/\.gz$/) {
                        $compressed=1;
                } else {
                        if($compressed) { die "Can't mix compressed and uncompressed inputs!\n"; }
                }
        }
        my $x="cat ".join(' ',@fastafiles)." |";
        if($compressed) { $x="gzip -dc ".join(' ',@fastafiles)." |"; }
        warn "#to open: $x\n";
        open(IN,"$x");
        while(<IN>) {
                chomp;
                if(/^>(.*)$/) {
                        my $lseq=length($seq);
			if($lseq>$MAXRES) {
				warn "# Truncating sequence of $lseq aa (longer than $MAXRES)\n";
				$seq=substr($seq,$[,$MAXRES);
				$lseq=$MAXRES;
			}
                        my $lhdr=length($header);
                        if(($nosplit<1) && ($ptr+$lseq>$CHUNKSIZE)) { # split kb.psq
                                warn "# split $section\: text at $ptr characters\n";
                                close(PIN);
                                close(PSQ);
                                close(PHR);
                                $section++;
                                #&openfiles($project,$section,$nosplit);
                                my $x="$project";
                                if($nosplit<1) { $x.="\.$section"; }
                                open(PHR,"\>$x\.phr");
                                open(PSQ,"\>$x\.psq");
                                open(PIN,"\>$x\.pin");

                                $ptr=0;
                                $hdrptr=0;
                        }
                        if($lseq > 0) {
                                $seq=~tr/[a-z]/[A-Z]/; # uppercase
                                print PHR $header,"\n";
                                print PSQ $seq,"\n";
                                print PIN join("\t",$ptr,$lseq,$hdrptr,$lhdr),"\n";
                                $hdrptr+=$lhdr+1; # newline
                                $ptr+=$lseq+1; # newline
                        }
                        $header=$1;
                        $seq='';
                } else {
                        if($removegaps>0) { s/\s//g; }
                        $seq.=$_;
                }
        }
        if($seq ne '') {
                $seq=~tr/[a-z]/[A-Z]/; # uppercase
                print PHR $header,"\n";
                print PSQ $seq,"\n";
                my $lseq=length($seq);
                print PIN join("\t",$ptr,$lseq),"\n";
        }

        close(PIN);
        close(PSQ);
        close(PHR);
        close(IN);
        return($section);
}

sub openfiles {
        my($project,$section,$nosplit)=@_;
        my $x="$project";
        if($nosplit<1) { $x.="\.$section"; }
        open(PHR,"\>$x\.phr");
        open(PSQ,"\>$x\.psq");
        open(PIN,"\>$x\.pin");
}

#######################################################################
# my batch queue
# bsub(npara,pollinterval,@cmdlist) waits until all processes finished
#######################################################################

sub bsub {
        my($NPARA,$POLLINTERVAL,@cmdlist)=@_;
        # fork subprocesses, wait for finish, so that max NPARA process run in parallel
        my @submitted;
        my $ncmd=$#cmdlist;
        my $icmd=-1;
        my $nrunning=0;
        while($icmd<$ncmd) {
                # how many processes are currently running?
                $nrunning=&howmanyrunning(@submitted);
                # submit more jobs
                while($nrunning<$NPARA && $icmd<$ncmd) {
                        $icmd++;
                        my $cmd=$cmdlist[$icmd];
                        my $pid=fork();
                        if($pid==0) {
                                &system_cmd($cmd);
                                exit 0;
                        }
                        push(@submitted,$pid);
                        $nrunning++;
                }
                # sleepy
                sleep($POLLINTERVAL);
        }
        # wait until all processes have finished
        while($nrunning>0) {
                sleep($POLLINTERVAL);
                $nrunning=&howmanyrunning(@submitted);
        }
}

sub howmanyrunning {
        my(@pidlist)=@_;
        if($#pidlist<0) { return(0); }
        my $cmd="ps -p @pidlist | grep -v defunct";
        my(@x)=`$cmd`;
        my $nrunning=$#x;
        return($nrunning);
}

#######################################################################
#
#######################################################################


sub system_cmd {
        my($cmd)=@_;
        my $status=system($cmd);
        warn "# system returned status: $status cmd: $cmd\n";
        if($status != 0) { die "FATAL ERROR\n"; }
}

#######################################################################
#
#######################################################################

