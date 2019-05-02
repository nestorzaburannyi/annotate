# create inverse dictionary of external2go tables
use strict;

my %inverse;
while(<STDIN>) {
	next if(/^!/);
	chomp;
	my(@x)=split(/ [>;] /);
	my($goid)=$x[2];
	($goid)=~ s/GO://;
	$inverse{$goid}.="$x[0] ";	
}
# output tabular: if multiple links then map to common ancestor
foreach my $goid (keys %inverse) {
	print join("\t",$goid,$inverse{$goid}),"\n";
}
