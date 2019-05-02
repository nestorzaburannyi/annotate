package SAP::RNA;

use strict;
use warnings;
use SAP::Common;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(rna_prediction);

sub rna_prediction {
    my ( $o, $s ) = @_;
    print_log( $o, "Starting RNA gene predictions..." );
    if ( $o->{"rna-r"} and $o->{"rna-r-program"} eq "rnammer" ) {
        run_rnammer ( $o, $s );
    }
    if ( $o->{"rna-t"} and $o->{"rna-t-program"} eq "trnascanse" ) {
        run_trnascanse ( $o, $s ) ;
    }
    if ( $o->{"rna-t"} and $o->{"rna-t-program"} eq "aragorn"
    or $o->{"rna-tm"} and $o->{"rna-tm-program"} eq "aragorn" ) {
        run_aragorn ( $o, $s );
    }
    if ( $o->{"rna-r"} and $o->{"rna-r-program"} eq "infernal"
    or $o->{"rna-t"} and $o->{"rna-t-program"} eq "infernal"
    or $o->{"rna-tm"} and $o->{"rna-tm-program"} eq "infernal"
    or $o->{"rna-nc"} and $o->{"rna-nc-program"} eq "infernal" ) {
        run_infernal ( $o, $s );
    }
    if ( $o->{"rna-r"} ) {
        parse_rrna_prediction ( $o, $s );
    }
    if ( $o->{"rna-t"} ) {
        parse_trna_prediction ( $o, $s )
    }
    if ( $o->{"rna-tm"} ) {
        parse_tmrna_prediction ( $o, $s )
    }
    if ( $o->{"rna-nc"} ) {
        parse_ncrna_prediction ( $o, $s )
    }
}

sub run_rnammer {
    my ( $o, $s ) = @_;
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
    my ( $o, $s ) = @_;
    print_log( $o, "Parsing rRNA predictions..." );
    while ( my $l = parse_file( $o, $o->{"o"}."/".$o->{"rna-r-program"}, "line", "\\s+", $o->{"rna-r-program"} ) ) {
        my ( $seq_id, $start, $end, $strand, $product, $score );
        #program-specific part
        if ( $o->{"rna-r-program"} eq "rnammer" ) {
            ( $seq_id, $start, $end, $strand, $product, $score ) = ( $l->[0], $l->[3], $l->[4], $l->[6] eq "+" ? 1 : -1, $l->[8], $l->[5] );
            # skip low scores if set by the user
            next if ( $score < $o->{"rna-r-score"} );
            # skip non-rRNA predictons
            next if not ( $product =~ m/s_rRNA$/ );
            $product =~ s/s_r/S ribosomal /i;
        }
        elsif ( $o->{"rna-r-program"} eq "infernal" ) {
            ( $seq_id, $start, $end, $strand, $product, $score ) = ( $l->[2], $l->[9] eq "+" ? $l->[7] : $l->[8], $l->[9] eq "+" ? $l->[8] : $l->[7], $l->[9] eq "+" ? 1 : -1, $l->[0], $l->[14] );
            # skip low scores if set by the user
            next if ( $score < $o->{"rna-r-score"} );
            # skip non-rRNA predictons
            next if not ( $product =~ m/SU_rRNA_bacteria|SU_rRNA_archaea|5S_rRNA/ );
            $product =~ s/LSU_rRNA_bacteria/23S ribosomal RNA/;
            $product =~ s/LSU_rRNA_archaea/23S ribosomal RNA/;
            $product =~ s/SSU_rRNA_bacteria/16S ribosomal RNA/;
            $product =~ s/SSU_rRNA_archaea/16S ribosomal RNA/;
            $product =~ s/5S_rRNA/5S ribosomal RNA/;
        }
        #create rRNA sequence feature
        my $feature = create_feature ( "rRNA", $seq_id, $start, "EXACT", $end, "EXACT", $strand, $score );
        $feature->add_tag_value ( "product", $product );
        # generate the inference tag
        my $inference = "profile:".$o->{"rna-r-program"}.":".$o->{$o->{"rna-r-program"}."-version"};
        $feature->add_tag_value ( "inference", $inference );
        #store rRNA sequence feature
        store_and_check_feature ( $o, $feature );
    }
}

sub parse_trna_prediction {
    my ( $o, $s ) = @_;
    print_log( $o, "Parsing tRNA predictions..." );
    while ( my $l = parse_file( $o, $o->{"o"}."/".$o->{"rna-t-program"}, "line", "\\s+", $o->{"rna-t-program"} ) ) {
        my ( $seq_id, $start, $end, $strand, $product, $score );
        # program-specific part
        if ( $o->{"rna-t-program"} eq "trnascanse" ) {
            ( $seq_id, $start, $end, $strand, $product, $score ) = ( $l->[0], $l->[2] < $l->[3] ? $l->[2] : $l->[3], $l->[2] < $l->[3] ? $l->[3] : $l->[2], $l->[2] < $l->[3] ? 1 : -1, $l->[4] eq "Pseudo" ? "tRNA-Xxx" : "tRNA-$l->[4]", $l->[8] );
            # skip low scores if set by the user
            next if ( $score < $o->{"rna-t-score"} );
            # skip non-tRNA predictons
            next if not ( $product =~ m/^tRNA/ );
        }
        elsif ( $o->{"rna-t-program"} eq "aragorn" ) {
            ( $seq_id, $start, $end, $strand, $product, $score ) = ( $l->[-1], ( $l->[2] =~ m/^(c?)\[(\d+),(\d+)\]$/ )[1], ( $l->[2] =~ m/^(c?)(\[\d+),(\d+)\]$/ )[2], ( $l->[2] =~ m/^(c?)(\[\d+),(\d+)\]$/ )[0] ne "c" ? 1 : -1, $l->[1], 'inf' );
            # skip low scores if set by the user
            next if ( $score < $o->{"rna-t-score"} );
            # skip non-tRNA predictons
            next if not ( $product =~ m/^tRNA/ );
        }
        elsif ( $o->{"rna-t-program"} eq "infernal" ) {
            ( $seq_id, $start, $end, $strand, $product, $score ) = ( $l->[2], $l->[9] eq "+" ? $l->[7] : $l->[8], $l->[9] eq "+" ? $l->[8] : $l->[7], $l->[9] eq "+" ? 1 : -1, $l->[0], $l->[14] );
            # skip low scores if set by the user
            next if ( $score < $o->{"rna-t-score"} );
            # skip non-tRNA predictons
            next if not ( $product =~ m/^tRNA/ );
        }
        # create tRNA sequence feature
        my $feature = create_feature ( "tRNA", $seq_id, $start, "EXACT", $end, "EXACT", $strand, $score );
        $feature->add_tag_value( "product", $product );
        $feature->add_tag_value( "pseudo", "_no_value" ) if ( $product =~ m/^tRNA-Xxx$/ );
        # generate the inference tag
        my $inference = "profile:".$o->{"rna-t-program"}.":".$o->{$o->{"rna-t-program"}."-version"};
        $feature->add_tag_value ( "inference", $inference );
        # store tRNA sequence feature
        store_and_check_feature ( $o, $feature );
    }
}

sub parse_tmrna_prediction {
    my ( $o, $s ) = @_;
    print_log( $o, "Parsing tmRNA predictions..." );
    while ( my $l = parse_file( $o, $o->{"o"}."/".$o->{"rna-tm-program"}, "line", "\\s+", $o->{"rna-tm-program"} ) ) {
        my ( $seq_id, $start, $end, $strand, $product, $score );
        #program-specific part
        if ( $o->{"rna-tm-program"} eq "aragorn" ) {
            ( $seq_id, $start, $end, $strand, $product, $score ) = ( $l->[-1], ( $l->[2] =~ m/^(c?)\[(\d+),(\d+)\]$/ )[1], ( $l->[2] =~ m/^(c?)(\[\d+),(\d+)\]$/ )[2], ( $l->[2] =~ m/^(c?)(\[\d+),(\d+)\]$/ )[0] ne "c" ? 1 : -1, $l->[1] =~ m/\*/ ? "$l->[1] (permuted)" : $l->[1], 'inf' );
            # skip low scores if set by the user
            next if ( $score < $o->{"rna-tm-score"} );
            # skip non-tmRNA predictons
            next if not ( $product =~ m/^tmRNA/ );
        }
        elsif ( $o->{"rna-tm-program"} eq "infernal" ) {
            ( $seq_id, $start, $end, $strand, $product, $score ) = ( $l->[2], $l->[9] eq "+" ? $l->[7] : $l->[8], $l->[9] eq "+" ? $l->[8] : $l->[7], $l->[9] eq "+" ? 1 : -1, $l->[0], $l->[14] );
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
        store_and_check_feature ( $o, $feature );
    }
}

sub parse_ncrna_prediction {
    my ( $o, $s ) = @_;
    print_log( $o, "Parsing ncRNA predictions..." );
    while ( my $l = parse_file( $o, $o->{"o"}."/".$o->{"rna-nc-program"}, "line", "\\s+", $o->{"rna-nc-program"} ) ) {
        my ( $seq_id, $start, $end, $strand, $product, $score );
        # program-specific part
        if ( $o->{"rna-nc-program"} eq "infernal" ) {                                                                                                                                                # taking trailing fields (17+) and removing subfields ending with ";", then joining the rest with spaces -> this goes to $product
            ( $seq_id, $start, $end, $strand, $product, $score ) = ( $l->[2], $l->[9] eq "+" ? $l->[7] : $l->[8], $l->[9] eq "+" ? $l->[8] : $l->[7], $l->[9] eq "+" ? 1 : -1, ( join " ", grep { not $_ =~ m/;$/ } @$l[17..$#$l] ), $l->[14] );
            # skip low scores if set by the user
            next if ( $score < $o->{"rna-nc-score"} );
            my $accession = $l->[1];
            # skip rRNA/tRNA/tmRNA genes - these are annotated from the same file elsewhere
            next if ( ( $l->[17] eq "Gene;" ) and ( $l->[18] eq "rRNA;" ) );
            next if ( ( $l->[17] eq "Gene;" ) and ( $l->[18] eq "tRNA;" ) );
            next if ( ( $l->[17] eq "Gene;" ) and ( $product =~ m/transfer-messenger RNA$/ ) );
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
        store_and_check_feature ( $o, $feature );
        }
    }
}

1;
