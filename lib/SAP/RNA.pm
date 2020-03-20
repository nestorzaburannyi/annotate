package SAP::RNA;

use strict;
use warnings;
use SAP::Common;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(rna_prediction);

sub rna_prediction {
    my ( $o ) = @_;
    print_log( $o, "Starting RNA gene predictions..." );
    if ( $o->{"rna-r"} and $o->{"rna-r-program"} eq "rnammer" ) {
        run_rnammer ( $o );
    }
    if ( $o->{"rna-t"} and $o->{"rna-t-program"} eq "trnascanse" ) {
        run_trnascanse ( $o ) ;
    }
    if ( $o->{"rna-t"} and $o->{"rna-t-program"} eq "aragorn"
    or $o->{"rna-tm"} and $o->{"rna-tm-program"} eq "aragorn" ) {
        run_aragorn ( $o );
    }
    if ( $o->{"rna-r"} and $o->{"rna-r-program"} eq "infernal"
    or $o->{"rna-t"} and $o->{"rna-t-program"} eq "infernal"
    or $o->{"rna-tm"} and $o->{"rna-tm-program"} eq "infernal"
    or $o->{"rna-nc"} and $o->{"rna-nc-program"} eq "infernal" ) {
        run_infernal ( $o );
    }
    if ( $o->{"rna-r"} ) {
        parse_rrna_prediction ( $o );
    }
    if ( $o->{"rna-tm"} ) {
        parse_tmrna_prediction ( $o )
    }
    if ( $o->{"rna-t"} ) {
        parse_trna_prediction ( $o )
    }
    if ( $o->{"rna-nc"} ) {
        parse_ncrna_prediction ( $o )
    }
}

sub run_rnammer {
    my ( $o ) = @_;
    print_log( $o, "Running RNAmmer ".$o->{"rnammer-version"}."..." );
    my ( $c ) = get_command ( $o, "rnammer" );
    run_program ( $o, $c );
}

sub run_trnascanse {
    my ( $o, $s ) = @_;
    print_log( $o, "Running tRNAscan-SE ".$o->{"trnascanse-version"}."..." );
    my ( $c ) = get_command ( $o, "trnascanse" );
    run_program ( $o, $c );
}

sub run_aragorn {
    my ( $o, $s ) = @_;
    print_log( $o, "Running ARAGORN ".$o->{"aragorn-version"}."..." );
    my ( $c ) = get_command ( $o, "aragorn" );
    run_program ( $o, $c );
}

sub run_infernal {
    my ( $o, $s ) = @_;
    print_log( $o, "Running Infernal ".$o->{"infernal-version"}."..." );
    my ( $c ) = get_command ( $o, "infernal" );
    run_program ( $o, $c );
}

sub parse_rrna_prediction {
  my ( $o ) = @_;
  print_log( $o, "Annotating rRNA predictions..." );
  foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"rna-r-program"}, "\\s+", $o->{"rna-r-program"} ) ) {
    # program-independent block
    # skip low scores if set by the user
    next if ( $l->{"score"} < $o->{"rna-r-score"} );
    # set the primary tag
    $l->{"primary_tag"} = "rRNA";
    # create rRNA sequence feature candidate
    my $feature = create_feature ( $o, $l );

    # program-specific block
    if ( $o->{"rna-r-program"} eq "rnammer" ) {
      # skip non-rRNA rnammer predictons
      next if not ( $l->{"product"} =~ m/S subunit ribosomal rRNA$/ );
      # product tag
      $feature->add_tag_value ( "product", $l->{"product"} );
      # inference tag
      $feature->add_tag_value ( "inference", "profile:".$o->{"rna-r-program"}.":".$o->{$o->{"rna-r-program"}."-version"} );
    }
    elsif ( $o->{"rna-r-program"} eq "infernal" ) {
      # skip non-rRNA infernal predictons
      next if not ( $l->{"type"} =~ m/^Gene; rRNA;$/ );
      # product tag
      $feature->add_tag_value ( "product", $l->{"product"} );
      # inference tag
      $feature->add_tag_value ( "inference", "profile:".$o->{"rna-r-program"}.":".$o->{$o->{"rna-r-program"}."-version"}.":rfam:".$l->{"accession"} );
    }

    # check and store rRNA sequence feature
    check_and_store_feature ( $o, $feature );
  }
}

sub parse_trna_prediction {
  my ( $o ) = @_;
  print_log( $o, "Annotating tRNA predictions..." );
  foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"rna-t-program"}, "\\s+", $o->{"rna-t-program"} ) ) {
    ## program-independent block
    # skip low scores if set by the user
    next if ( $l->{"score"} < $o->{"rna-t-score"} );
    # set the primary tag
    $l->{"primary_tag"} = "tRNA";
    # create tRNA sequence feature candidate
    my $feature = create_feature ( $o, $l );

    # program-specific block
    if ( $o->{"rna-t-program"} eq "trnascanse" ) {
      # skip non-tRNA trnascanse predictons
      next if not ( $l->{"product"} =~ m/^tRNA/ );
      # product tag
      $feature->add_tag_value ( "product", $l->{"product"} );
      # inference tag
      $feature->add_tag_value ( "inference", "profile:".$o->{"rna-t-program"}.":".$o->{$o->{"rna-t-program"}."-version"} );
    }
    elsif ( $o->{"rna-t-program"} eq "aragorn" ) {
      # skip non-tRNA aragorn predictons
      next if not ( $l->{"product"} =~ m/^tRNA/ );
      # product tag
      $feature->add_tag_value ( "product", $l->{"product"} );
      # inference tag
      $feature->add_tag_value ( "inference", "profile:".$o->{"rna-t-program"}.":".$o->{$o->{"rna-t-program"}."-version"} );
    }
    elsif ( $o->{"rna-t-program"} eq "infernal" ) {
      # skip non-tRNA infernal predictons
      next if not ( $l->{"type"} =~ m/^Gene; tRNA;$/ );
      # product tag
      $feature->add_tag_value ( "product", $l->{"product"} );
      # inference tag
      $feature->add_tag_value ( "inference", "profile:".$o->{"rna-t-program"}.":".$o->{$o->{"rna-t-program"}."-version"}.":rfam:".$l->{"accession"} );
    }

    # check and store tRNA sequence feature
    check_and_store_feature ( $o, $feature );
  }
}

sub parse_tmrna_prediction {
  my ( $o ) = @_;
  print_log( $o, "Annotating tmRNA predictions..." );
  foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"rna-tm-program"}, "\\s+", $o->{"rna-tm-program"} ) ) {
    # program-independent block
    # skip low scores if set by the user
    next if ( $l->{"score"} < $o->{"rna-tm-score"} );
    # set the primary tag
    $l->{"primary_tag"} = "tmRNA";
    # create tmRNA sequence feature candidate
    my $feature = create_feature ( $o, $l );

    # program-specific block
            # skip low scores if set by the user
            next if ( $score < $o->{"rna-tm-score"} );
            # skip non-tmRNA predictons
            next if not ( $product =~ m/^tmRNA/ );
        }
        # create tmRNA sequence feature
        my $feature = create_feature ( "tmRNA", $seq_id, $start, "EXACT", $end, "EXACT", $strand, $score );
        $feature->add_tag_value ( "product", $product );
        # generate the inference tag
        my $inference = "profile:".$o->{"rna-tm-program"}.":".$o->{$o->{"rna-tm-program"}."-version"};
        $feature->add_tag_value ("inference", $inference );
        #store tmRNA sequence feature
        store_feature ( $o, $feature );
    }
}

sub parse_ncrna_prediction {
  my ( $o ) = @_;
  print_log( $o, "Annotating ncRNA predictions..." );
  foreach my $l ( parse_file( $o, $o->{"job"}."/".$o->{"rna-nc-program"}, "\\s+", $o->{"rna-nc-program"} ) ) {
    # program-independent block
    # skip low scores if set by the user
    next if ( $l->{"score"} < $o->{"rna-nc-score"} );
    # set the primary tag
    $l->{"primary_tag"} = "ncRNA";
    # create ncRNA sequence feature candidate
    my $feature = create_feature ( $o, $l );

            # create ncRNA sequence feature
            my $feature = create_feature ( "ncRNA", $seq_id, $start, "EXACT", $end, "EXACT", $strand, $score );
            # these are "true" ncRNA
            if ( $l->[17] eq "Gene;" ) {
                if ( $l->[18] eq "antisense;" ) {
                    $feature->add_tag_value ( "ncRNA_class", "antisense" );
                    $feature->add_tag_value ( "product", $product );
                }
                elsif ( $l->[18] eq "antitoxin;" ) {
                    $feature->add_tag_value ( "ncRNA_class", "antitoxin" );
                    $feature->add_tag_value ( "product", $product );
                }
                elsif ( $l->[18] eq "CRISPR;" ) {
                    $feature->add_tag_value ( "ncRNA_class", "CRISPR" );
                    $feature->add_tag_value ( "product", $product );
                }
                elsif ( $l->[18] eq "lncRNA;" ) {
                    $feature->add_tag_value ( "ncRNA_class", "lncRNA" );
                    $feature->add_tag_value ( "product", $product );
                }
                elsif ( $l->[18] eq "miRNA;" ) {
                    $feature->add_tag_value ( "ncRNA_class", "miRNA" );
                    $feature->add_tag_value ( "product", $product );
                }
                elsif ( $l->[18] eq "ribozyme;" ) {
                    # special cases
                    if ( $accession eq "RF00010" ) {
                        $feature->add_tag_value ( "ncRNA_class", "RNase_P_RNA" );
                        $feature->add_tag_value ( "product", $product );
                    }
                    else {
                        $feature->add_tag_value ( "ncRNA_class", "ribozyme" );
                        $feature->add_tag_value ( "product", $product );
                    }
                }
                elsif ( $l->[18] eq "snRNA;" ) {
                    $feature->add_tag_value ( "ncRNA_class", "snRNA" );
                    $feature->add_tag_value ( "product", $product );
                }
                elsif ( $l->[18] eq "sRNA;" ) {
                    $feature->add_tag_value ( "ncRNA_class", "sRNA" );
                    $feature->add_tag_value ( "product", $product );
                }
                else {
                    # special cases
                    if ( $accession eq "RF00169" ) {
                        $feature->add_tag_value ( "ncRNA_class", "Bacteria_small_SRP" );
                        $feature->add_tag_value ( "product", $product );
                    }
                    else {
                        $feature->add_tag_value ( "ncRNA_class", "other" );
                        $feature->add_tag_value ( "product", $product );
                    }
                }
            }
            elsif ( $l->[17] eq "intron;" ) {
                $feature->add_tag_value ( "ncRNA_class", "autocatalytically_spliced_intron" );
                $feature->add_tag_value ( "product", $product );
            }
            # these are "false" ncRNA features, so we change the primary_tag, but with some exceptions
            elsif ( $l->[17] eq "Cis-reg;" ) {
                if ( $l->[18] eq "frameshift_element" ) {
                    $feature->primary_tag("misc_feature" );
                    $feature->add_tag_value ( "note", $product );
                }
                elsif ( $l->[18] eq "IRES;" ) {
                    $feature->primary_tag("misc_feature" );
                    $feature->add_tag_value ( "note", $product );
                }
                elsif ( $l->[18] eq "leader;" ) {
                    $feature->primary_tag("misc_feature" );
                    $feature->add_tag_value ( "note", $product );
                }
                elsif ( $l->[18] eq "riboswitch;" ) {
                    # special cases
                    if ( $accession eq "RF00234" ) {
                        $feature->add_tag_value ( "ncRNA_class", "ribozyme" );
                        $feature->add_tag_value ( "product", $product );
                        $feature->add_tag_value ( "gene", "glmS" );
                    }
                    # Riboswitches should be annotated as misc_binding feature when there is a known bound moiety and they aren't defined as ribozymes.
                    elsif ( $accession eq "RF00050" ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "flavin mononucleotide" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( $accession eq "RF00059" ) {

                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "thiamine pyrophosphate" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( ( $accession eq "RF01482") or ( $product eq "RF01689") or ( $product eq "RF00174") ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "adenosylcobalamin" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( $product eq "RF00504" ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "glycine" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( ( $product eq "RF00634") or ( $product eq "RF01826") or ( $product eq "RF01725") or ( $product eq "RF00162") or ( $product eq "RF01767") or ( $product eq "RF00521") ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "S-adenosylmethionine" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( $product eq "RF01727" ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "S-adenosylmethionine and/or S-adenosylhomocysteine" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( $product eq "RF01057" ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "S-adenosylhomocysteine" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( $product eq "RF00167" ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "guanine and/or adenine" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( $product eq "RF00168" ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "lysine" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( ( $product eq "RF00522") or ( $product eq "RF01054") ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "pre-queuosine1" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( $product eq "RF01831" ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "tetrahydrofolate" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( $product eq "RF01055" ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "molybdenum cofactor" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    elsif ( $product eq "RF01786" ) {
                        $feature->primary_tag("misc_binding" );
                        $feature->add_tag_value ( "bound_moiety", "cyclic di-GMP" );
                        $feature->add_tag_value ( "note", $product );
                    }
                    # Riboswitches with unknown functions or bound moieties are annotated as misc_features with the name in a note
                    else {
                        $feature->primary_tag("misc_feature" );
                        $feature->add_tag_value ( "note", $product );
                    }
                }
                elsif ( $l->[18] eq "thermoregulator;" ) {
                    $feature->primary_tag("misc_feature" );
                    $feature->add_tag_value ( "note", $product );
                }
                else {
                    $feature->primary_tag("misc_feature" );
                    $feature->add_tag_value ( "note", $product );
                }
            }
            else {
                next;
            }
        # generate the inference tag
        my $inference = "profile:".$o->{"rna-nc-program"}.":".$o->{$o->{"rna-nc-program"}."-version"}.":rfam:$l->[1]";
        $feature->add_tag_value ( "inference", $inference );
        # store ncRNA sequence feature
        store_feature ( $o, $feature );
        }
    }
}

1;
