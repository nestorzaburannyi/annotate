#!/usr/bin/perl
# convert annotation table to HTML
# USAGE: perl anno2html.pl species.anno
# write to ARG1.html

use strict;
#use CGI;
use Switch;

my($title,$pr)=@ARGV;
my $predictor="";
if($pr eq "argot" || $pr eq "jac" || $pr eq "hyge") { $predictor=$pr; }

my $MAXPRED=10;
my $SANSLINK="http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/sans.cgi?mode=table&seq=";

#my(@annotab)=<STDIN>;
print "<TITLE>$title</TITLE>\n";
print &anno2html($predictor); ##,@annotab);
#print "<PRE><H1>INPUTS</H1>\n"; print join("\n",@annotab),"\n</PRE>\n";
exit;

###############################################################################

#######################################################################
sub anno2html {
my($predictor)=@_;
# print html header
my $tableheaderrow="<TR><TH>".join("</TH><TH>","Query header","GN","Description<BR>Estimated PPV, description","Biological process<BR>Estimated PPV, GO-id, description","Molecular function<BR>Estimated PPV, GO-id, description","Cellular component<BR>Estimated PPV, GO-id, description","Inverse ec2go, kegg2go")."</TH></TR>";
my @lines;
push(@lines, &rgbstyles);
# read sansparser output from STDIN
my @result;
my $oldqpid='';
my $nid=1;
while(<STDIN>) {
#foreach (@_) {
        # skip comments
        next if (/^\#/);
	# collect block
        chomp;
        my($row)=$_;
        my $qpid;
        if (/^(\S+)/) { $qpid=$1; } else { next; }
        if($qpid ne $oldqpid) {
                # print annotation summary row of previous query
                if($oldqpid ne '') { push(@lines,&nested_tables($predictor,@result)); }
                $oldqpid=$qpid;
                @result=();
                # table break
                if($nid % 10 == 0 || $nid==1) { push(@lines, "</TR>\n</TABLE>\n<TABLE BORDER>\n$tableheaderrow\n"); }
                $nid++;
        }
        push(@result,$row);
}
# last block
if($oldqpid ne '') { push(@lines, &nested_tables($predictor,@result)); }

# print html footer
push(@lines, "</TR>\n</TABLE>\n");
return(@lines);
}

###############################################################################

sub nested_tables {
	my($predictor)=shift;
#        shift(@_); # skip header row
        my $orig="";
        my @bp; #("<TR><TH>RM3</TH><TH>GOid</TH><TH>desc</TH></TR>\n");
        my @mf; #("<TR><TH>RM3</TH><TH>GOid</TH><TH>desc</TH></TR>\n");
        my @cc; #("<TR><TH>RM3</TH><TH>GOid</TH><TH>desc</TH></TR>\n");
        my @ec;
	my @kegg;
	my $qpid='';
        my @de; #("<TR><TH>RM2</TH><TH>DE</TH></TR>\n");
        # qpid    type    score   id      description
        # limit number of predicted classes to MAXPRED
        my $nde=0;
        my $nbp=0;
        my $nmf=0;
        my $ncc=0;
	my $nec=0;
	my $nkegg=0;
	my $gn="";

        foreach (@_) {
                chomp;
                my($qpid,$type,$s,$ppv,$id,$desc)=split(/\t/);
                next if($qpid eq '');
                $qpid=~s/^>//;
                my $score=$ppv;
                if($ppv !~ /^[\d\.]+$/) { $score=0.0; } else { $score=sprintf("%4.2f",$ppv); }
                my $color=int($score*100); if($color<0) { $color=0; } if ($color>100) { $color=99; }
                # hyperlink EC and KEGG to dbget.jp
                if($id=~/KEGG:/) { $id="<A HREF=http://www.genome.jp/dbget-bin/www_bget?$id>$id</A>"; } 
		if($id=~/EC:(\S+)/) {
			my $s=$1;
			$_=$s;
			my $url="<A HREF=http://enzyme.expasy.org/EC";
			if(/\d+\.\d+\.\d+\.\d+/) { $id="$url\/$s\>$id</A>"; }
			elsif(/\d+\.\d+\.\d+/) { $id="$url\/$s\.\-\>$id</A>"; }
                        elsif(/\d+\.\d+/) { $id="$url\/$s\.\-\.\-\>$id</A>"; }
                        else {$id="$url\/$s\.\-\.\-\.\-\>$id</A>"; }
		}
		switch($type) {
                        case 'original_DE' { $orig="$qpid<BR>$desc"; }
                        case 'qseq' { $orig.="<BR><A HREF=$SANSLINK$desc&hdr=$qpid target=_blank>Search</A>"; }
                        case 'DE' { $nde++; if($nde<=$MAXPRED) { push(@de,"<TR><TD><span class=color$color>$score</span></TD><TD>$desc</TD></TR>\n"); } }
			case 'GN' { $gn=$desc; }
		}
		if($predictor eq '') {
		  switch($type) {
                        case 'BP_RM3' { $nbp++; if($nbp<=$MAXPRED) { push(@bp,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'MF_RM3' { $nmf++; if($nmf<=$MAXPRED) { push(@mf,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'CC_RM3' { $ncc++; if($ncc<=$MAXPRED) { push(@cc,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
			case 'EC_RM3' { $nec++; if($nec==1) { push(@ec,"<TR><TD><span class=color$color>$score</span></TD><TD>$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'KEGG_RM3' { $nkegg++; if($nkegg==1) { push(@kegg,"<TR><TD><span class=color$color>$score</span></TD><TD>$id</TD><TD>$desc</TD></TR>\n"); } }
		  } 
                }elsif($predictor eq 'argot') {
                  switch($type) {
                        case 'BP_ARGOT' { $nbp++; if($nbp<=$MAXPRED) { push(@bp,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'MF_ARGOT' { $nmf++; if($nmf<=$MAXPRED) { push(@mf,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'CC_ARGOT' { $ncc++; if($ncc<=$MAXPRED) { push(@cc,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'EC_ARGOT' { $nec++; if($nec==1) { push(@ec,"<TR><TD><span class=color$color>$score</span></TD><TD>$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'KEGG_ARGOT' { $nkegg++; if($nkegg==1) { push(@kegg,"<TR><TD><span class=color$color>$score</span></TD><TD>$id</TD><TD>$desc</TD></TR>\n"); } }
                  }
                }elsif($predictor eq 'jac') {
                  switch($type) {
                        case 'BP_JAC' { $nbp++; if($nbp<=$MAXPRED) { push(@bp,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'MF_JAC' { $nmf++; if($nmf<=$MAXPRED) { push(@mf,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'CC_JAC' { $ncc++; if($ncc<=$MAXPRED) { push(@cc,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'EC_JAC' { $nec++; if($nec==1) { push(@ec,"<TR><TD><span class=color$color>$score</span></TD><TD>$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'KEGG_JAC' { $nkegg++; if($nkegg==1) { push(@kegg,"<TR><TD><span class=color$color>$score</span></TD><TD>$id</TD><TD>$desc</TD></TR>\n"); } }
                  }
                }elsif($predictor eq 'hyge') {
                  switch($type) {
                        case 'BP_HYGE' { $nbp++; if($nbp<=$MAXPRED) { push(@bp,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'MF_HYGE' { $nmf++; if($nmf<=$MAXPRED) { push(@mf,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'CC_HYGE' { $ncc++; if($ncc<=$MAXPRED) { push(@cc,"<TR><TD><span class=color$color>$score</span></TD><TD>GO:$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'EC_HYGE' { $nec++; if($nec==1) { push(@ec,"<TR><TD><span class=color$color>$score</span></TD><TD>$id</TD><TD>$desc</TD></TR>\n"); } }
                        case 'KEGG_HYGE' { $nkegg++; if($nkegg==1) { push(@kegg,"<TR><TD><span class=color$color>$score</span></TD><TD>$id</TD><TD>$desc</TD></TR>\n"); } }
                  }
		}
        }
        # one row = one query's annotations
        my $html="<TR><TD>$orig</TD><TD>$gn</TD>\n<TD><TABLE style=\"width:300 px\">\n@de</TABLE></TD>\n<TD><TABLE style=\"width:300 px\">\n@bp</TABLE></TD>\n<TD><TABLE style=\"width:300 px\">\n@mf</TABLE></TD>\n<TD><TABLE style=\"width:300 px\">\n@cc</TABLE></TD><TD><TABLE style=\"width:300 px\">\n@ec\n@kegg</TABLE></TD>\n</TR>";
        return($html);
}

sub rgbstyles {
        my $style="<style type='text/css'>\ntable,td,th{border-collapse:collapse;padding:2px 10px;}\nth{cursor:n-resize;background-color:lightgrey;}\n</style>\n<style>\n";
# define colors in steps of 0.1
        my $i=0;
        my $x=0.0;
        while($x<=1) {
                $x+=0.01;
                my($r,$g,$b)=&rgb_colormap($x);
                $style.="span.color$i { background-color: rgb($r,$g,$b); }\n";
                $i++;
        }
        $style.="</style>\n";
        return($style);
}

sub rgb_colormap {
        my($x)=@_;
        my $blue=0;
        my $green=1;
        my $red=0;
        if($x<0.5) { $red=1; $green=2*$x; } else { $green=1; $red=2*(1-$x); }
# return RGB
        return(int(255*$red),int(255*$green),int(255*$blue));
}

sub rgb_colormap_rainbow {
        my ($x)=@_;
# scale values 0..1 to blue-green-red gradient
        my $blue=0;
        my $green=1;
        my $red=0;
# blue
        if($x<0.25) { $blue=1; } elsif($x<0.5) { $blue=1-4*($x-0.25); }
# green
        if($x<0.25) { $green=4*$x; } elsif($x>0.75) { $green=1-4*($x-0.75); }
# red
        if($x>0.75) { $red=1; } elsif($x>0.5) { $red=4*($x-0.5); }
# return RGB
        return(int(255*$red),int(255*$green),int(255*$blue));
}

###############################################################################
