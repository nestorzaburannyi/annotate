use strict;
use List::Util qw(first);

# scatterplot
my %marker; $marker{'Eukaryota'}='diamond'; $marker{'Bacteria'}='circle'; $marker{'Archaea'}='cross'; $marker{'other'}='square';
my %kingdom;
my %species;
my %median;
my %average;
my %matched;
my $MAXGENUS=8;
my $ngenus=0;
my @genuslist;
my %other;
#my $MINIDENTITY=0.5;
# barchart
my @headers; 
my @traces; 
#my $H=20;
my $n=0;
# wsum    count   kingdom genus   matched_fraction        average_identity        median  pide_bins       lineage species
my $species_col=9;
my $pide_col=5;
my $median_col=6;
my $core_col=4;
my $genus_col=3;
my $king_col=2;
my $lineage_col=8;
my $vector_col=7;
my $multiplicity_col=-1;

if($#ARGV < 8) { die "USAGE: $0 TITLE onesided|bidirectional HIDFOLDER MINIDENTITY H MIN_X MAX_X MIN_Y MAX_Y \n"; }
my($TITLE,$METRIC,$HIDFOLDER,$MINIDENTITY,$H,$MIN_X,$MAX_X,$MIN_Y,$MAX_Y)=@ARGV; # scatterplot, barchart parameters
warn "arguments: ($TITLE,$METRIC,$HIDFOLDER,$MINIDENTITY,$H,$MIN_X,$MAX_X,$MIN_Y,$MAX_Y)\n";
my $check2='';
my $check1='checked';
my $taxtabfile="$HIDFOLDER\/tax.tab";
if($METRIC eq 'bidirectional') { $check2='checked'; $check1=''; $taxtabfile="$HIDFOLDER\/tax2.tab"; }
print<<EOB;
<html>
<head>
<title>AAI-profiler</title>
<head>
<body>
<A HREF=\"/AAI/\">AAI-profiler home</A>
<h1>$TITLE</h1>
<head>
<script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js"></script>
<script>
\$(document).ready(function () {

    \$(\"show\").hide();

    \$(\"#btn1\").click(function(){
         \$(\"#p1\").toggle();
    });
});
</script>
</head>
<body>
<div id=\"btn1\" class=\"button\"><a href=\"#\">Toggle parameters</a></div>
<show id=\"p1\">
<form action="/cgi-bin/sans/tax2scatter.cgi" method="POST">
<input type="hidden" name="TITLE" value="$TITLE">
<input type="hidden" name="HIDFOLDER" value="$HIDFOLDER">
 <fieldset>
    <legend>AAI method</legend>
    <input type="radio" name="metric" value="onesided" $check1> Onesided  
    <input type="radio" name="metric" value="bidirectional" $check2> Bidirectional
 </fieldset>
 <fieldset>
    <legend>Scatterplot</legend>
    Exclude data points with AAI less than <input type="text" name="MINIDENTITY" value=$MINIDENTITY size=5><br>
  </fieldset>
  <fieldset>
    <legend>AAI histograms</legend>
    Show top <input type="text" name="H" value="$H" size=5> species with AAI values between 
    <input type="text" name="MIN_X" value="$MIN_X" size=5> and <input type="text" name="MAX_X" value="$MAX_X" size=5>
    and matched fraction between <input type="text" name="MIN_Y" value="$MIN_Y" size=5> and <input type="text" name="MAX_Y" value="$MAX_Y" size=5> 
  </fieldset>
  <input type="submit" value="Redraw plots">
</form>

</show>
<h2>AAI of Uniprot species</h2>
Species are grouped and coloured by genus. Literature suggests &gt 95 \% identity at the species boundary. Matched fraction is maximal in completely sequenced genomes.
Use Plotly option (at top right of plot area) to show closest data on hover and to autoscale or zoom.
<div id='myDiv'></div>
<script>
EOB
my $firstline=1;
my $size=18;
open(IN,"<$taxtabfile") || die "Can't open $taxtabfile\n"; 
#warn "#reading $taxtabfile\n";
while(<IN>) {
        chomp;
        my(@data)=split(/\t/);
	# get column indeces from header line
        if($firstline) {
		$species_col = first { $data[$_] eq 'species' } 0..$#data;
                $pide_col = first { $data[$_] eq 'average_identity' } 0..$#data;
                $median_col = first { $data[$_] eq 'median' } 0..$#data;
                $core_col = first { $data[$_] eq 'matched_fraction' } 0..$#data;
                $genus_col = first { $data[$_] eq 'genus' } 0..$#data;
                $king_col = first { $data[$_] eq 'kingdom' } 0..$#data;
                $lineage_col = first { $data[$_] eq 'lineage' } 0..$#data;
                $vector_col = first { $data[$_] eq 'pide_bins' } 0..$#data;
		$multiplicity_col = first { $data[$_] eq 'multiplicity' } 0..$#data;
		$firstline=0;
		#warn "# @data gave species=$species_col average_idenity=$pide_col median=$median_col core=$core_col genus=$genus_col kingdom=$king_col lineage=$lineage_col vector=$vector_col\n";
		next;
	}
	# data
        my $species=$data[$species_col];
        my $average_identity=sprintf("%5.3f",$data[$pide_col]);
        my $median=sprintf("%4.2f",$data[$median_col]);
        my $matched_fraction=sprintf("%5.3f",$data[$core_col]);
	my $multiplicity=sprintf("%6.1f",$data[$multiplicity_col]);
	# scatterplot
        if($data[$pide_col]>=$MINIDENTITY) { # filter too many data points
	        my $genus=$data[$genus_col];
	        if(!defined($species{$genus})) { # top ten genus colored differently
	                $ngenus++;
			my $kingdom=$data[$king_col]; 
			if( ($kingdom ne 'Eukaryota') && ($kingdom ne 'Bacteria') && ($kingdom ne 'Archaea') ) { $kingdom='other'; }
	                if($ngenus>=$MAXGENUS) {
	                        $genus=$kingdom;
				$other{$kingdom}=1;
	                } else {
	                        push(@genuslist,$genus);
	                }
                        $kingdom{$genus}=$kingdom;
		}
	        push(@{$species{$genus}},'"'.$species.'"');
	        push(@{$median{$genus}},$median);
	        push(@{$average{$genus}},$average_identity);
	        push(@{$matched{$genus}},$matched_fraction);
	}
	# barcharts
	next if($average_identity < $MIN_X || $average_identity > $MAX_X);
	next if($matched_fraction < $MIN_Y || $matched_fraction > $MAX_Y);
	if($n<$H) {
		$n++;
	        push(@headers,"<h3>$n\. $species</h3>\nAverage amino acid identity = $average_identity, Median amino acid identity = $median, Matched fraction = $matched_fraction, Multiplicity = $multiplicity\n<br>Lineage: $data[$lineage_col]\n<div id=\'myDiv$n\'></div>\n");
	        push(@traces,"var trace$n = [\{ y: \[$data[$vector_col]\], type: \'bar\'\}];\nvar layout$n = \{ legend: { font: { size: $size } }, xaxis: \{ title: 'Amino-acid identity', titlefont: { size: $size }, tickfont: { size: $size}, showticklabels: true\}, yaxis: \{ title: 'Number of proteins', titlefont: { size: $size }, tickfont: { size: $size }, showticklabels: true\}, height: 400\};\n");
	}
}
close(IN);

# scatterplot: color genus differently in data series
my @traces1;
foreach my $genus (@genuslist,(sort keys %other)) {
	next if(!defined($species{$genus}));
	my $x=$genus.'_1'; $x=~s/\W/_/g; # $x=~s/[\/\(\)\[\]\-\s]/\_/g;
	warn "genus $genus transformed to trace $x\n";
	my $marker=$marker{$kingdom{$genus}}; 
        print "var $x = { name: '$genus', \nmarker: { symbol: \"$marker\", size: $size },\n y: [",join(",\n",@{$matched{$genus}}),"], x: [",join(",\n",@{$average{$genus}}),"], mode: 'markers', type: 'scatter', text: [",join(",\n",@{$species{$genus}}),"] };\n";
        push(@traces1,"$x");
}
print "var layout = { legend: { font: { size: $size } }, height: 800, yaxis: {  title: 'Matched fraction', titlefont: { size: $size }, tickfont: { size: $size } }, xaxis: { title: 'Average Amino-acid Identity', titlefont: { size: $size }, tickfont: { size: $size } } };\n";
print "var data = [",join(",",@traces1),"];\n";
print "Plotly.newPlot('myDiv',data,layout);\n";
print "</script>\n";
# barcharts
print "<h2>AAI histograms of Uniprot species</h2>\n";
print @headers;
print "<script>\n";
print @traces;
foreach my $i (1..$n) { print "Plotly.newPlot(\'myDiv$i\',trace$i, layout$i );\n"; }
print "</script>\n</body>\n";
print "</html>\n";


