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
    if ( $o->{"dna-c"} and $o->{"dna-c-program"} eq "pilercr" ) {
        run_pilercr ( $o );
    }
    if ( $o->{"dna-c"} and $o->{"dna-c-program"} eq "minced" ) {
        run_minced ( $o );
    }
    if ( $o->{"dna-t"} and $o->{"dna-t-program"} eq "trf" ) {
        run_trf ( $o );
    }

    if ( $o->{"dna-c"} ) {
        parse_crispr_prediction ( $o );
    }
    if ( $o->{"dna-t"} ) {
        parse_tandem_prediction ( $o );
    }
}

sub run_pilercr {
    my ( $o ) = @_;
    print_log( $o, "Running PILER-CR ".$o->{"pilercr-version"}."..." );
    my ( $c ) = get_command ( $o, "pilercr" );
    #run_program ( $o, $c );
}

sub run_minced {
    my ( $o ) = @_;
    print_log( $o, "Running MinCED ".$o->{"minced-version"}."..." );
    my ( $c ) = get_command ( $o, "minced" );
    run_program ( $o, $c ) if ! $o->{"test"};
    system("cp /annotate/databases/temp/minced ".$o->{"job"}."/minced") if $o->{"test"};
}

sub run_trf {
    my ( $o ) = @_;
    print_log( $o, "Running Tandem Repeat Finder ".$o->{"trf-version"}."..." );
    my ( $c ) = get_command ( $o, "trf" );
    #run_program ( $o, $c );
}

sub parse_crispr_prediction {
  my ( $o ) = @_;
  print_log( $o, "Parsing CRISPR predictions..." );
  foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"dna-c-program"}, "\\t", $o->{"dna-c-program"} ) ) {
    # program-independent block
    # skip low scores if set by the user
    next if ( $l->{"score"} < $o->{"dna-c-score"} );
    # create repeat_region sequence feature candidate
    my $feature = create_feature ( $o, $l );

    # program-specific block
    if ( $o->{"dna-c-program"} eq "pilercr" ) {
      # inference tag
      $feature->add_tag_value ("inference", "COORDINATES:ab initio prediction:".$o->{"dna-c-program"}.":".$o->{$o->{"dna-c-program"}."-version"} );
    }
    elsif ( $o->{"dna-c-program"} eq "minced" ) {
      # inference tag
      $feature->add_tag_value ("inference", "COORDINATES:ab initio prediction:".$o->{"dna-c-program"}.":".$o->{$o->{"dna-c-program"}."-version"} );
    }

    # check and store repeat_region sequence feature
    check_and_store_feature ( $o, $feature );
  }
}

sub parse_tandem_prediction {
  my ( $o, $s ) = @_;
  print_log( $o, "Parsing tandem repeat predictions..." );
  foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"dna-t-program"}, "\\s+", $o->{"dna-t-program"} ) ) {
    # program-independent block
    # skip low scores if set by the user
    next if ( $l->{"score"} < $o->{"dna-t-score"} );
    # set the primary tag
    $l->{"primary_tag"} = "repeat_region";

    # create repeat_region sequence feature candidate
    my $feature = create_feature ( $o, $l );

    # program-specific block
    if ( $o->{"dna-t-program"} eq "trf" ) {
      # rpt_type tag
      $feature->add_tag_value ("rpt_type", "tandem");
      # inference tag
      $feature->add_tag_value ("inference", "COORDINATES:ab initio prediction:".$o->{"dna-t-program"}.":".$o->{$o->{"dna-t-program"}."-version"} );
    }

    # check and store repeat_region sequence feature
    check_and_store_feature ( $o, $feature );
  }
}

1;
