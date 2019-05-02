# input: scaffold sequeces (fasta)
# output: three-frame translations

use strict;

my $hdr="";
my $seq="";
my $MINORFLEN=80;
my($x)=@ARGV;
if($x>0) { $MINORFLEN=$x; }

while(<STDIN>) {
	chomp;
	if(/^>(\S+)/) {
		if($hdr ne "") { &output($hdr,$seq); }
		$hdr=$1;
		$seq="";
	} else {
		$seq.=$_;
	} 	
}
# last entry
if($hdr ne "") { &output($hdr,$seq); }

sub output { # all orfs, filter by length
        my($hdr,$seq)=@_;
	$hdr=~s/^>//;
        $seq=~s/\W//g;
       my $rev=reverse $seq;
       $rev =~ tr/[ACGT]/[TGCA]/;
	&output_orfs($hdr,&translate($seq));
	&output_orfs($hdr,&translate($rev));
}

sub output_orfs {
	my($hdr,$seq)=@_;
	$seq=~s/X+/X/g; # remove masked segments so SANS doesn't match runs of Xs!
	my(@x)=split(/[\*X]+/,$seq);
	my $i=0;
	foreach(@x) {
		my $len=length($_);
		next if($len<$MINORFLEN);
		$i++;
		print ">orf$i\_$hdr\n$_\n";
	}
}

sub output_translation {
	my($hdr,$seq)=@_;
	$seq=~s/\W//g; 
	my $rev=reverse $seq;
	foreach my $frame (1..3) {
                print "$hdr frame $frame\n";
                print &translate($seq),"\n";
                substr($seq,0,1)=""; # remove first character to shift frame
        }
	# reverse direction
	$rev =~ tr/[ACGT]/[TGCA]/;
        foreach my $frame (-3..-1) {
                print "$hdr frame $frame\n";
                print &translate($rev),"\n";
                substr($rev,0,1)=""; # remove first character to shift frame
        }
}

sub translate {
        my($seq)=@_;
        # translate ORFs
        my $protein="";
        my(@x)=split(//,$seq);
        while($#x > 0) {
                my($first)=shift(@x);
                my($second)=shift(@x);
                my($third)=shift(@x);
                if($first eq "T") {
                        if($second eq "T") {
                                if($third eq "T" || $third eq "C") {
                                        $protein.="F";
                                } else {
                                        $protein.="L";
                                }
                        } elsif($second eq "C") {
                                $protein.="S";
                        } elsif($second eq "A") {
                                if($third eq "T" || $third eq "C") {
                                        $protein.="Y";
                                } else {
                                        $protein.="*"; # stop
                                }
                        } elsif($second eq "G") {
                                if($third eq "A") {
                                        $protein.="*"; # stop
                                } elsif($third eq "G") {
                                        $protein.="W";
                                } else {
                                        $protein.="C";
                                }
                        }
                } elsif ($first eq "C") {
                        if($second eq "T") {
                                $protein.="L";
                        } elsif($second eq "C") {
                                $protein.="P";
                        } elsif($second eq "A") {
                                if($third eq "T" or $third eq "C") {
                                        $protein.="H";
                                } else {
                                        $protein.="Q";
                                }
                        } else {
                                $protein.="R";
                        }
                } elsif ($first eq "A") {
                        if($second eq "T") {
                                if($third eq "G") {
                                        $protein.="M";
                                } else {
                                        $protein.="I";
                                }
                        } elsif($second eq "C") {
                                $protein.="T";
                        } elsif($second eq "A") {
                                if($third eq "T" || $third eq "C") {
                                        $protein.="N";
                                } else {
                                        $protein.="K";
                                }
                        } else {
                                if($third eq "T" || $third eq "C") {
                                        $protein.="S";
                                } else {
                                        $protein.="R";
                                }
                        }
                } elsif($first eq "G") {
                        if($second eq "T") {
                                $protein.="V";
                        } elsif($second eq "C") {
                                $protein.="A";
                        } elsif($second eq "A") {
                                if($third eq "T" || $third eq "C") {
                                        $protein.="D";
                                } else {
                                        $protein.="E";
                                }
                        } else {
                                $protein.="G";
                        }
                } else {
                        $protein.="X";
                }
        }
        #warn "# protein: $protein\n";
        return($protein);
}

