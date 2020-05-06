package SAP::CDS;

use strict;
use warnings;
use SAP::Common;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(cds_prediction);

sub cds_prediction {
  my ( $o ) = @_;
  # cds-h
  print_log( $o, "Starting homology CDS predictions..." ) if ( $o->{"cds-h"} );
  prepare_input ( $o, "sequences" ) if ( $o->{"cds-h"} ) ;
  run_blast ( $o ) if ( $o->{"cds-h"} and $o->{"cds-h-program"} eq "blast" );
  run_diamond ( $o ) if ( $o->{"cds-h"} and $o->{"cds-h-program"} eq "diamond" );
  parse_homology ( $o ) if ( $o->{"cds-h"} ) ;
  # cds-i
  print_log( $o, "Starting ab initio CDS predictions..." ) if ( $o->{"cds-i"} );
  prepare_input ( $o, "sequences" ) if ( $o->{"cds-i"} ) ;
  run_glimmer ( $o ) if ( $o->{"cds-i"} and $o->{"cds-i-program"} eq "glimmer" );
  run_prodigal ( $o ) if ( $o->{"cds-i"} and $o->{"cds-i-program"} eq "prodigal" );
  run_genemarks ( $o ) if ( $o->{"cds-i"} and $o->{"cds-i-program"} eq "genemarks" );
  parse_ab_initio ( $o ) if ( $o->{"cds-i"} );
  # cds-a
  prepare_input ( $o, "annotation" ) if ( $o->{"cds-a"} );
  print_log( $o, "Starting protein annotation..." ) if ( $o->{"cds-a"} );
  run_pannzer ( $o ) if ( $o->{"cds-a"} and $o->{"cds-a-program"} eq "pannzer" );
  run_emapper ( $o ) if ( $o->{"cds-a"} and $o->{"cds-a-program"} eq "emapper" );
  parse_annotation ( $o ) if ( $o->{"cds-a"} );
}

sub run_glimmer {
  my ( $o ) = @_;
    print_log( $o, "Running Glimmer ".$o->{"glimmer-version"}."..." );
  my ( $c ) = get_command ( $o, "glimmer1" );
  run_program ( $o, $c ) if ! $o->{"test"};
  ( $c ) = get_command ( $o, "glimmer2" );
  run_program ( $o, $c ) if ! $o->{"test"};
  ( $c ) = get_command ( $o, "glimmer3" );
  run_program ( $o, $c ) if ! $o->{"test"};
}

sub run_prodigal {
  my ( $o ) = @_;
  print_log( $o, "Running Prodigal ".$o->{"prodigal-version"}."..." );
  my ( $c ) = get_command ( $o, "prodigal" );
    run_program ( $o, $c );
}

sub run_genemarks {
  my ( $o ) = @_;
    print_log( $o, "Running GeneMarkS ".$o->{"genemarks-version"}."..." );
    my ( $c ) = get_command ( $o, "genemarks" );
    run_program ( $o, $c );
}

sub run_blast {
  my ( $o ) = @_;
    print_log( $o, "Running BLASTX ".$o->{"blast-version"}."..." );
    my ( $c ) = get_command ( $o, "blast" );
    run_program ( $o, $c );
}

sub run_diamond {
  my ( $o ) = @_;
  print_log( $o, "Running DIAMOND ".$o->{"diamond-version"}."..." );
  my ( $c ) = get_command ( $o, "diamond" );
    run_program ( $o, $c );
}

sub run_pannzer {
  my ( $o ) = @_;
  print_log( $o, "Running PANNZER ".$o->{"pannzer-version"}."..." );
  my ( $c ) = get_command ( $o, "pannzer" );
    run_program ( $o, $c );
}

sub run_emapper {
  my ( $o ) = @_;
  print_log( $o, "Running emapper ".$o->{"emapper-version"}."..." );
  my ( $c ) = get_command ( $o, "emapper" );
    run_program ( $o, $c );
}

sub parse_ab_initio {
  my ( $o ) = @_;
  print_log( $o, "Annotating ab initio CDS predictions..." );
  foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"cds-i-program"}, "\\s+", $o->{"cds-i-program"} ) ) {
    # program-independent block
    # skip low scores if set by the user
    next if ( $l->{"score"} < $o->{"cds-i-score"} );
    # set the primary tag
    $l->{"primary_tag"} = "CDS";
    # create CDS sequence feature candidate
    my $feature = create_feature ( $o, $l );

    # program-specific block
    if ( $o->{"cds-i-program"} eq "prodigal" ) {
      # inference tag
      $feature->add_tag_value ("inference", "COORDINATES:ab initio prediction:".$o->{"cds-i-program"}.":".$o->{$o->{"cds-i-program"}."-version"} );
    }
    elsif ( $o->{"cds-i-program"} eq "glimmer" ) {
      # inference tag
      $feature->add_tag_value ("inference", "COORDINATES:ab initio prediction:".$o->{"cds-i-program"}.":".$o->{$o->{"cds-i-program"}."-version"} );
    }
    elsif ( $o->{"cds-i-program"} eq "genemarks" ) {
      # inference tag
      $feature->add_tag_value ("inference", "COORDINATES:ab initio prediction:".$o->{"cds-i-program"}.":".$o->{$o->{"cds-i-program"}."-version"} );
    }

    # check and store CDS sequence feature
    check_and_store_feature ( $o, $feature );
  }
}

sub parse_homology {
  my ( $o ) = @_;
  print_log( $o, "Annotating homology CDS predictions..." );
  foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"cds-h-program"}, "\\t", $o->{"cds-h-program"} ) ) {
    # program-independent block
    # skip low scores if set by the user
    next if ( $l->{"score"} < $o->{"cds-h-score"} );
    # set the primary tag
    $l->{"primary_tag"} = "CDS";
    # create CDS sequence feature candidate
    my $feature = create_feature ( $o, $l );

    # program-specific block
    if ( $o->{"cds-h-program"} eq "blast" or $o->{"cds-h-program"} eq "diamond" ) {
      # product tag
      $feature->add_tag_value ( "product", $l->{"product"} );
      # gene tag
      $feature->add_tag_value ( "gene", $l->{"gene_name"} );
      # inference tag
      $feature->add_tag_value ( "inference", "COORDINATES:alignment:swissprot:".( $l->{"hit_id"} =~ m/^\S+\|(\S+_\S+)$/ )[0] );
      # check and store CDS sequence feature
      check_and_store_feature ( $o, $feature );
    }
  }
}

sub parse_annotation {
    my ( $o ) = @_;
    print_log( $o, "Annotating annotation results..." );
    foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"cds-a-program"}, "\\t", $o->{"cds-a-program"} ) ) {
      # program-independent block
      # skip low scores if set by the user
      next if ( $l->{"score"} < $o->{"cds-a-score"} );
      # get the CDS feature
      my $feature = $o->{"feature_by_id"}->{ $l->{"seq_id"} };
      # remove all previous product tags
      $feature->remove_tag( "product" ) if $feature->has_tag("product");
      # product tag
      $feature->add_tag_value ( "product", clean_up_description( $l->{"product"} ) );
      # gene tag
      $feature->add_tag_value ( "gene", $l->{"gene_name"} );
      print_verbose( $o, "Found product ".$l->{"product"}." for query".$l->{"seq_id"} );
    }
}



1;
