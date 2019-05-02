# convert raw SANS output to XML

use strict;
use Switch;

my $programtag='';
my $databasetag='';
my $paramtag='';
my %query;
my @errors;
my @timeinfo;
my @sbjct_list;
my @lines=();

print "<?xml version=\"1.0\"?>\n";
my $inquery=0;
while(<STDIN>) {
	if(/ERROR/) {
		print "<EBIApplicationError>$_</EBIApplicationError>\n"
	}
	if(/<QUERY/) { 
		$inquery=1; 
	}
	if($inquery) { 
		push(@lines,$_); 
	}	
	if(/<\/QUERY/) { 
		$inquery=0; 
		# returns %query, @errors, @params, @timeinfo, $programtag, @sbjct_list == array of hashes with keywords
		&parse_result_keywords(@lines);
		&xml();
		@lines=();
	}
}

exit();

###############################################################################

sub xml {
	# collect header items
	my $nhits=$#sbjct_list+1-$[;
	# parse database tag
	my $dbname='';
	my $dbdate='';
	my $H='';
	my $naa='';
	my $nprot='';
	my $evalue_cutoff='';
	my $searchtime='';
	my $timestamp_start='';
	my $timestamp_end='';
	#<DATABASE= swiss.Jun2015 letters=              49059460 sequences=       137147 >
	$_=$databasetag;
	if(/DATABASE=\s*(\S+)\.(\S+)/) {
		$dbname=$1;
		$dbdate=$2;
	}
	if(/letters=\s*(\d+)/) { $naa=$1; }
	if(/sequences=\s*(\d+)/) { $nprot=$1; }
	# parse SANS parameters
	# <PARAM H=          10  EVALUE_CUTOFF=   1.0000000     >
	$_=$paramtag;
	if(/H=\s*(\d+)/) { $H=$1; }
	if(/EVALUE_CUTOFF=\s*(\S+)/) { $evalue_cutoff=$1; }
	# parse time stamps
	foreach(@timeinfo) {
		if(/search=\s*(\S+)/) { $searchtime=$1; }
		if(/start=\"(.*)\"/) { $timestamp_start=$1; }
		if(/end=\"(.*)\"/) { $timestamp_end=$1; }
	}
	# XML header
	print<<EOB;
<EBIApplicationResult xsi:noNamespaceSchemaLocation="http://www.ebi.ac.uk/schema/AppplicationResult.xsd" xmlns:xsi="http://www.w3/org/2001/XMLSchema-instance" xmlns="http://www.ebi.ac.uk/schema">
 - <Header>
	$programtag
     -  <parameters>
         - <sequences total="1">
	 	<sequence name="$query{'hdr'}" length="$query{'LSEQ'}" type="p" number="1"/>
	   </sequences>
	 - <databases total="1" letter="$naa" sequences="$nprot">
	 	<database name="$dbname" type="p" number="1" created="$dbdate"/>
	   </databases>
	   <scores>$H</scores>
	   <alignments>0</alignments>
	   <matrix>BLOSUM62</matrix>
	   <expectationUpper>$evalue_cutoff</expectationUpper>
	</parameters>
	<timeInfo search="$searchtime" end="$timestamp_end" start="$timestamp_start"/>
    </Header>
 - <SequenceSimilaritySearchResult>
    - <hits total="$nhits">
EOB
        # print one Hits element per sbjct
        my $rank=0;
	my $sortkey='BITS';
        if(!defined($sbjct_list[0]->{$sortkey})) { $sortkey='VOTE'; }
        foreach my $rec (sort { $b->{$sortkey} <=> $a->{$sortkey} } @sbjct_list) {
                $rank++;
                my $lali=$rec->{'LALI'};
		my($desc)=get_desc($rec->{'hdr'});
		$_=get_pid($rec->{'hdr'});
		my($db,$acc,$id)=split(/\|/); # check exceptions: Uniref, pdb
		if(/Uniref/) { ($db,$acc)=split(/\_/); $_=$desc; ($id)=/RepId=(\S+)/; }
		elsif(!/\|/) { $db='pdb'; $id=$_; }
		my $qcov=int(100*$lali/(0.1+$query{'LSEQ'}));
		my $scov=int(100*$lali/(0.1+$rec->{'LSEQ'}));
		my $ide=100*$rec->{'PIDE'};
		#>> check exception: VOTE (verifast)
		print<<EOB;
        - <hit length=\"$rec->{'LSEQ'}\" number=\"$rank\" description=\"$desc\" ac=\"$acc\" id=\"$id\" database=\"$db\">
	   - <alignments total="1">
	      - <alignment number="1">
	          <initn>$rec->{'VOTE'}</initn>
		  <bits>$rec->{'BITS'}</bits>
		  <expectation>$rec->{'EVALUE'}</expectation>
		  <identity>$ide</identity>
		  <querySeq start=\"$rec->{'QFROM'}\" end=\"$rec->{'QTO'}\">
		  <matchSeq start=\"$rec->{'SFROM'}\" end=\"$rec->{'STO'}\">
		  <overlap>$lali</overlap>
		  <queryMatch>$qcov</queryMatch>
		  <match>$scov</match>
	        </alignment>
	     </alignments>
	  </hit>
EOB
	}
	# close XML tags
	print "       </hits>\n    </SequenceSimilaritySearchResult>\n</EBIApplicationResult>\n";
}

###############################################################################

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
        # returns %query, $paramtag, @timeinfo, $programtag,$databasetag, @sbjct_list == array of hashes with keywords
	undef %query;
	@timeinfo=();
	$programtag='';
	$databasetag='';
	$paramtag='';
        @sbjct_list=();
        my $inquery=0;
        my $insbjct=0;
        my $hdr='';
        my $nsbjct=0;
        foreach(@_) {
		chomp;
                next if(/^\s*#/);
                if(/^\s*ERROR/) { 
			push(@errors,$_);  # echo error messages
		} elsif(/<program/) {
			$programtag=$_;
		} elsif(/timeInfo/) {
			push(@timeinfo,$_);
		} elsif(/<PARAM/) {
			$paramtag=$_;
		} elsif(/<DATABASE/) {
			$databasetag=$_;
		} elsif(/<QUERY/) {
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
__END__

<QUERY nid=     1 vote_cutoff=          55        LSEQ=         154     >
<DATABASE= swiss.Jun2015 letters=              49059460 sequences=       137147 >
<program name=SANSparallel version=2.1 citation=PMID:25855811>
 <PARAM H=          10  EVALUE_CUTOFF=   1.0000000     >
 <timeInfo start="Date 06/20/2015; time 16:34:47"/>
>101mA   MYOGLOBIN
MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRVKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGNFGADAQGAMNKALELFRKDIAAKYKELGYQG
<SBJCT VOTE=    3086 TUPS=     781 PIDE=  0.99 LALI=   153 BITS=   331.6 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   153 SFROM=     1 STO=   153 >
>sp|P02185|MYG_PHYCD Myoglobin OS=Physeter catodon GN=MB PE=1 SV=2
MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG
</SBJCT>
<SBJCT VOTE=    2988 TUPS=     760 PIDE=  0.95 LALI=   153 BITS=   322.8 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   153 SFROM=     1 STO=   153 >
>sp|P02184|MYG_KOGSI Myoglobin OS=Kogia sima GN=MB PE=1 SV=2
MVLSEGEWQLVLHVWAKVEADIAGHGQDILIRLFKHHPETLEKFDRFKHLKSEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPADFGADAQGAMSKALELFRKDIAAKYKELGYQG
</SBJCT>
<SBJCT VOTE=    2797 TUPS=     720 PIDE=  0.90 LALI=   152 BITS=   306.0 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   152 SFROM=     1 STO=   152 >
>sp|P68279|MYG_TURTR Myoglobin OS=Tursiops truncatus GN=MB PE=1 SV=2
MGLSDGEWQLVLNVWGKVEADLAGHGQDVLIRLFKGHPETLEKFDKFKHLKTEADMKASEDLKKHGNTVLTALGAILKKKGHHDAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPAEFGADAQGAMNKALELFRKDIAAKYKELGFHG
</SBJCT>
<SBJCT VOTE=    2769 TUPS=     718 PIDE=  0.90 LALI=   152 BITS=   305.2 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   152 SFROM=     1 STO=   152 >
>sp|P68277|MYG_PHODA Myoglobin OS=Phocoenoides dalli dalli GN=MB PE=1 SV=2
MGLSEGEWQLVLNVWGKVEADLAGHGQDVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHGNTVLTALGGILKKKGHHDAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPAEFGADAQGAMNKALELFRKDIATKYKELGFHG
</SBJCT>
<SBJCT VOTE=    2835 TUPS=     727 PIDE=  0.91 LALI=   152 BITS=   308.9 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   152 SFROM=     1 STO=   152 >
>sp|Q0KIY3|MYG_PENEL Myoglobin OS=Peponocephala electra GN=MB PE=2 SV=3
MGLSDGEWQLVLNVWGKVEADLAGHGQDILIRLFKGHPETLEKFDKFKHLKTEADMKASEDLKKHGITVLTALGAILKKKGHHDAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPAEFGADAQGAMNKALELFRKDIAAKYKELGFHG
</SBJCT>
<SBJCT VOTE=    2797 TUPS=     720 PIDE=  0.90 LALI=   152 BITS=   306.0 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   152 SFROM=     1 STO=   152 >
>sp|P68276|MYG_DELDE Myoglobin OS=Delphinus delphis GN=MB PE=1 SV=2
MGLSDGEWQLVLNVWGKVEADLAGHGQDVLIRLFKGHPETLEKFDKFKHLKTEADMKASEDLKKHGNTVLTALGAILKKKGHHDAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPAEFGADAQGAMNKALELFRKDIAAKYKELGFHG
</SBJCT>
<SBJCT VOTE=    2769 TUPS=     718 PIDE=  0.90 LALI=   152 BITS=   305.2 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   152 SFROM=     1 STO=   152 >
>sp|P68278|MYG_PHOPH Myoglobin OS=Phocoenoides phocoena GN=MB PE=1 SV=2
MGLSEGEWQLVLNVWGKVEADLAGHGQDVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHGNTVLTALGGILKKKGHHDAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPAEFGADAQGAMNKALELFRKDIATKYKELGFHG
</SBJCT>
<SBJCT VOTE=    2608 TUPS=     699 PIDE=  0.86 LALI=   153 BITS=   297.2 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   153 SFROM=     1 STO=   153 >
>sp|P14396|MYG_CASFI Myoglobin OS=Castor fiber GN=MB PE=1 SV=2
MGLSDGEWQLVLHVWGKVEADLAGHGQEVLIRLFKGHPETLEKFNKFKHIKSEDEMKASEDLKKHGVTVLTALGGVLKKKGHHEAEIKPLAQSHATKHKIPIKYLEFISEAIIHVLQSKHPGBFGADABGAMNKALELFRKDIAAKYKELGFQG
</SBJCT>
<SBJCT VOTE=    2988 TUPS=     760 PIDE=  0.95 LALI=   153 BITS=   322.8 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   153 SFROM=     1 STO=   153 >
>sp|Q0KIY5|MYG_KOGBR Myoglobin OS=Kogia breviceps GN=MB PE=2 SV=3
MVLSEGEWQLVLHVWAKVEADIAGHGQDILIRLFKHHPETLEKFDRFKHLKSEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPADFGADAQGAMSKALELFRKDIAAKYKELGYQG
</SBJCT>
<SBJCT VOTE=    2813 TUPS=     724 PIDE=  0.91 LALI=   152 BITS=   307.7 EVALUE=  0.00000000E+00 DIAG=      -1 LSEQ=     154 QFROM=     1 QTO=   152 SFROM=     1 STO=   152 >
>sp|P02174|MYG_GLOME Myoglobin OS=Globicephala melas GN=MB PE=1 SV=2
MGLSDGEWQLVLNVWGKVEADLAGHGQDILIRLFKGHPETLEKFDKFKHLKTEADMKASEDLKKHGNTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPAEFGADAQGAMNKALELFRKDIAAKYKELGFHG
</SBJCT>
 <timeInfo search=   3.96220684051513672E-002  >
 <timeInfo end="Date 06/20/2015; time 16:34:47"/>
</QUERY>

