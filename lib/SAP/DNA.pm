package SAP::DNA;

use strict;
use warnings;
use SAP::Common;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(dna_prediction);

sub dna_prediction {
    my ( $o ) = @_;
    print_log( $o, "Starting DNA feature predictions..." );
    if ( $o->{"dna-tandem"} and $o->{"dna-tandem-program"} eq "trf" ) {
        run_trf ( $o );
    }

    if ( $o->{"dna-tandem"} ) {
        parse_tandem_prediction ( $o );
    }
}

sub run_trf {
    my ( $o ) = @_;
    print_log( $o, "Running Tandem Repeat Finder ".$o->{"trf-version"}."..." );
    my ( $c ) = get_command ( $o, "trf" );
    run_program ( $o, $c );
}


sub parse_tandem_prediction {
  my ( $o, $s ) = @_;
  print_log( $o, "Parsing tandem repeat predictions..." );
  foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"dna-tandem-program"}, "\\s+", $o->{"dna-tandem-program"} ) ) {
    # program-independent block
    # skip low scores if set by the user
    next if ( $l->{"score"} < $o->{"dna-tandem-score"} );
    # set the primary tag
    $l->{"primary_tag"} = "repeat_region";

    # create repeat_region sequence feature candidate
    my $feature = create_feature ( $o, $l );

    # program-specific block
    if ( $o->{"dna-tandem-program"} eq "trf" ) {
      # rpt_type tag
      $feature->add_tag_value ("rpt_type", "tandem");
      # inference tag
      $feature->add_tag_value ("inference", "COORDINATES:ab initio prediction:".$o->{"dna-tandem-program"}.":".$o->{$o->{"dna-tandem-program"}."-version"} );
    }

    # check and store repeat_region sequence feature
    check_and_store_feature ( $o, $feature );
  }
}

1;