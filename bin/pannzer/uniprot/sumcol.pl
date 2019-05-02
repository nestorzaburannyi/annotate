use strict;

my($col)=@ARGV;

my $summa=0;
while(<STDIN>) {
	my(@x)=split(/\s+/);
	$summa+=$x[$col];
}
print "$summa\n";
