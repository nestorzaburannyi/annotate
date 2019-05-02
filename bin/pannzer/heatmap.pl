use strict;

# input: kt_self.txt, matrix.txt
# output: plotly HTML

my $H=50;
my $n=0;
my @x; # labels
my @y; # labels
my @z; # csv row-vectors
open(IN,"<kt_self.txt");
while(<IN>) {
	next if(/^#/);
	$n++;
	last if($n>$H);
	chomp;
	my(@tmp)=split(/\t/);
        my($species)=pop(@tmp);
	push(@x,$species);
}
close(IN);
open(IN,"<matrix.txt");
while(<IN>) {
	next if(/qpid/);
	chomp;
	my($qpid,$csv)=split(/\t/);
	push(@y,$qpid);
	push(@z,$csv);
}
close(IN);
warn "# data sizes: x $#x y $#y z $#z\n";
# output HTML 
my $zstring=join('],[',@z);
my $ystring=join('","',@y);
my $xstring=join('","',@x);
my $height=5*($#y+1);
print<<EOB;
<head>
  <!-- Plotly.js -->
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>

<body>
Browser compatibility: Chrome recommended. Windows Edge does not display Plotly heatmap correctly.  
  <div id="myDiv"><!-- Plotly chart will be drawn inside this DIV --></div>
  <script>
    <!-- JAVASCRIPT CODE GOES HERE -->
var layout = {
  autosize: false,
  width: 1500,
  height: $height,
  margin: {
    l: 160,
    r: 100,
    b: 400,
    t: 100,
    pad: 4
  },
}
var data = [
  {
z: [ [$zstring] ],
y: [ "$ystring" ],
x: [ "$xstring" ],
colorscale: 'Picnic',
type: 'heatmap'
  }
];

Plotly.plot('myDiv', data, layout);

  </script>
</body>
EOB

