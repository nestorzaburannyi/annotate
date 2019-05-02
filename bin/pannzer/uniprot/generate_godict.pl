#!/usr/bin/perl

# inputs: go_data = goid, parentlist; gene_goa from STDIN: accession_number, GO:goid
# outputs: godict.txt (STDOUT), GO.IEA.counts (STDERR)

use strict;

if($#ARGV<2) { die "USAGE: perl $0 go_data obo_with_ic.tab eg2go.tab
kegg2go.tab > mergeGO.out\n"; }
my($GO_DATA,$IC_TAB,$EC_TAB,$KEGG_TAB)=@ARGV;

my %gocount;
my %ontology;
my %name;
my %plist;

# hash external2go resoureces
my %ec;
my %kegg;
my %ic;
open(IN,"<$EC_TAB") || warn "No ec2go mapping provided\n";
while(<IN>) {
         chomp;
         my($goid,$list)=split(/\t/);
         $ec{$goid}=$list;
}
close(IN);
open(IN,"<$KEGG_TAB") || warn "No kegg2go mapping provided\n";
while(<IN>) {
         chomp;
         my($goid,$list)=split(/\t/);
         $kegg{$goid}=$list;
}
close(IN);
open(IN,"<$IC_TAB") || warn "No IC data provided\n";
while(<IN>) {
         chomp;
         next if(/^goid/);
my($goid,$ic,$namespace,$children,$parentlist,$gocount,$parentcount,$name)=split(/\t/);
         next if($namespace eq '');
         $ontology{$goid}=$namespace;
         $name{$goid}=$name;
         $gocount{$goid}=$gocount;
         $ic{$goid}=$ic;
}
close(IN);
# propagated parents from go_data
open(IN,"<$GO_DATA") || warn "No parent lists provided\n";
while(<IN>) {
         chomp;
         my($goid,$namespace,$name,$parentlist)=split(/\t/);
         $plist{$goid}=join(",",split(/,/,$parentlist)); # sort parent list
}
close(IN);

# output  gocounts to STDERR
foreach my $goid (keys %gocount) {
         my $ec=$ec{$goid};
         my $kegg=$kegg{$goid};
         my $p=$plist{$goid};
         my $ic=$ic{$goid};
         print
join("\t",$gocount{$goid},$ontology{$goid},$goid,$name{$goid},$p,$ec,$kegg,$ic),"\n";
}
