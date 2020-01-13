package SAP::CDS;

use strict;
use warnings;
use SAP::Common;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(cds_prediction);

sub cds_prediction {
    my ( $o, $s ) = @_;
    print_log( $o, "Starting ab initio CDS predictions..." ) if ( $o->{"cds-i"} );
    run_glimmer ( $o, $s ) if ( $o->{"cds-i"} and $o->{"cds-i-program"} eq "glimmer" );
    run_prodigal ( $o, $s ) if ( $o->{"cds-i"} and $o->{"cds-i-program"} eq "prodigal" );
    run_genemarks ( $o, $s ) if ( $o->{"cds-i"} and $o->{"cds-i-program"} eq "genemarks" );
    parse_ab_initio ( $o, $s ) if ( $o->{"cds-i"} );
    while ( $o->{"cds-h"} ) {
        print_log( $o, "Starting homology CDS predictions, iteration ".$o->{"cds-h"}."..." ) if ( $o->{"cds-h"} );
        prepare_input ( $o, $s, "homology" ) if ( $o->{"cds-h"} ) ;
        run_blast ( $o, $s ) if ( $o->{"cds-h"} and $o->{"cds-h-program"} eq "blast" );
        run_diamond ( $o, $s ) if ( $o->{"cds-h"} and $o->{"cds-h-program"} eq "diamond" );
        parse_homology ( $o, $s ) if ( $o->{"cds-h"} ) ;
    }
    prepare_input ( $o, $s, "annotation" ) if ( $o->{"cds-a"} );
    print_log( $o, "Starting protein annotation..." ) if ( $o->{"cds-a"} );
    run_pannzer ( $o, $s ) if ( $o->{"cds-a"} and $o->{"cds-a-program"} eq "pannzer" );
    run_emapper ( $o, $s ) if ( $o->{"cds-a"} and $o->{"cds-a-program"} eq "emapper" );
    parse_annotation ( $o, $s ) if ( $o->{"cds-a"} );
}

sub run_glimmer {
    my ( $o, $s ) = @_;
    print_log( $o, "Running Glimmer ".$o->{"glimmer-version"}."..." );
    my ( $c ) = get_command ( $o, "glimmer" );
    run_program ( $o, $c );
}

sub run_prodigal {
    my ( $o, $s ) = @_;
    print_log( $o, "Running Prodigal ".$o->{"prodigal-version"}."..." );
    my ( $c ) = get_command ( $o, "prodigal" );
    run_program ( $o, $c );
}

sub run_genemarks {
    my ( $o, $s ) = @_;
    print_log( $o, "Running GeneMarkS ".$o->{"genemarks-version"}."..." );
    my ( $c ) = get_command ( $o, "genemarks" );
    run_program ( $o, $c );
}

sub run_blast {
    my ( $o, $s ) = @_;
    print_log( $o, "Running BLASTX ".$o->{"blast-version"}."..." );
    my ( $c ) = get_command ( $o, "blast" );
    run_program ( $o, $c );
}

sub run_diamond {
    my ( $o, $s ) = @_;
    print_log( $o, "Running DIAMOND ".$o->{"diamond-version"}."..." );
    my ( $c ) = get_command ( $o, "diamond" );
    run_program ( $o, $c );
}

sub run_pannzer {
    my ( $o, $s ) = @_;
    print_log( $o, "Running PANNZER ".$o->{"pannzer-version"}."..." );
    my ( $c ) = get_command ( $o, "pannzer" );
    run_program ( $o, $c );
}

sub run_emapper {
    my ( $o, $s ) = @_;
    print_log( $o, "Running emapper ".$o->{"emapper-version"}."..." );
    my ( $c ) = get_command ( $o, "emapper" );
    run_program ( $o, $c );
}

sub parse_ab_initio {
    my ( $o, $s ) = @_;
    print_log( $o, "Parsing ab initio CDS predictions..." );
    my ( $c ) = get_command ( $o, $o->{"cds-i-program"} );
    while ( my $l = parse_file( $o, $o->{"o"}."/".$o->{"cds-i-program"}, "line", "\\s+", $o->{"cds-i-program"} ) ) {
        my ( $seq_id, $start, $end, $strand, $score, $start_before, $end_after );
        #program-specific part
        if ( $o->{"cds-i-program"} eq "prodigal" ) {
            ( $seq_id, $start, $end, $strand, $score ) = ( $l->[0], $l->[3], $l->[4], $l->[6] eq "+" ? 1 : -1, $l->[5] );
            ( $start_before, $end_after, my $start_type, my $stop_type ) = $l->[8] =~ m/;partial=([0|1])([0|1]);start_type=(.*);stop_type=(.*);rbs_motif=/;
            $start_before = "1" if ( ( ! $stop_type ) and ( $strand eq "-1" ) );
            $end_after = "1" if ( ( ! $start_type ) and ( $strand eq "1" ) );
            # skip low scores if set by the user
            next if ( $score < $o->{"cds-i-score"} );
        }
        elsif ( $o->{"cds-i-program"} eq "glimmer" ) {
            ( $seq_id, $start, $end, $strand, $score ) = ( $l->[-1], $l->[3] =~ m/\+/ ? $l->[1] : $l->[2], $l->[3] =~ m/\+/ ? $l->[2] : $l->[1], $l->[3] =~ m/\+/ ? 1 : -1, $l->[4] );
            # glimmer can extend beyong sequence, both sides
            $start = $start+3 if ( $start < 1 );
            $end = $end-3 if ( $end > $s->{$seq_id}->length );
            # do not skip
            # next if ( $score < $o->{"cds-i-score"} );
        }
        elsif ( $o->{"cds-i-program"} eq "genemarks" ) {
            ( $seq_id, $start, $end, $strand, $score ) = ( $l->[0], $l->[3], $l->[4], $l->[6] eq "+" ? 1 : -1, $l->[5] );
            # skip low scores if set by the user
            # next if ( $score < $o->{"cds-i-score"} );
        }
        # create CDS sequence feature
        my $feature = create_feature ( "CDS", $seq_id, $start, $start_before ? "BEFORE" : "EXACT", $end, $end_after ? "AFTER" : "EXACT", $strand, $score );
        # WАRNING: valid [SEQ_FEAT.PartialProblem] PartialLocation: 3' partial is not at stop AND is not at consensus splice site FEATURE: CDS: Alginate lyase [lcl|sequence_1:c985-<2] [lcl|sequence_1: raw, dna len= 76377] -> [lcl|sequence_1_1]
        # NОTE: valid [SEQ_FEAT.PartialProblem] PartialLocation: Stop does not include first/last residue of sequence (but is at consensus splice site) FEATURE: CDS: Type IV secretion protein Rhs [lcl|sequence_1:71735->76375] [lcl|sequence_1: raw, dna len= 76377] -> [
        # if $start_before and $strand = 1, or $end_after and $strand = -1, check the difference between sequence start / end and apply codon_start, if necessary
        # if ( $start_before and $strand eq 1 and $start > 1 ) {
        #   $feature->add_tag_value ("codon_start", $start );
        #   $start = 1;
        # }
        # if ( $end_after and $strand eq -1 and $end < $s->{$seq_id}->length ) {
        #   $feature->add_tag_value ("codon_start", $s->{$seq_id}->length );
        #   $end = $s->{$seq_id}->length;
        # }
        # generate the inference tag
        my $inference = "COORDINATES:ab initio prediction:".$o->{"cds-i-program"}.":".$o->{$o->{"cds-i-program"}."-version"};
        $feature->add_tag_value ("inference", $inference );
        #store CDS sequence feature
        check_and_store_feature ( $o, $feature );
    }
}

sub parse_homology {
    my ( $o, $s ) = @_;
    print_log( $o, "Parsing homology CDS predictions..." );
    my $success;
    while ( my $l = parse_file( $o, $o->{"o"}."/output_homology_".$o->{"cds-h"}, "line", "\\s+", $o->{"cds-h-program"} ) ) {
        # program-specific part
        if ( ( $o->{"cds-h-program"} eq "blast" ) or ( $o->{"cds-h-program"} eq "diamond" ) ) {
            #qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score slen cached
            my ( $qseqid, $sseqid, $qstart, $qend, $qstrand, $sstart, $score, $identity, $cached ) = ( $l->[0], $l->[1], $l->[6] < $l->[7] ? $l->[6] : $l->[7], $l->[6] < $l->[7] ? $l->[7] : $l->[6], $l->[6] < $l->[7] ? 1 : -1, $l->[8], $l->[11], $l->[2], $l->[13] );
            my ( $seq_id, $qoffset, $source_strand ) = ( $l->[0] =~ m/^(\d+)_offset(\S+)_strand(.?\d)_.*$/ );
            # set an "exists" flag so we will not repeat this search
            $o->{"homology"}->{$qseqid} = 1;
            # skip low scores if set by the user
            next if ( $score < $o->{"cds-h-score"} );
            # skip intraCDS hits unless we are quite sure we deal with the very same protein (>= 90% identity by default)
            next if ( ( abs($source_strand) eq 1 ) and ( $identity < $o->{"cds-h-identity"} ) );
            ################### DETECT ###################
            my $feature = find_appropriate_cds_feature ( $o, $s, $seq_id, $qstart + $qoffset - ( $qstrand eq 1 ? ($sstart-1)*3 : 0 ), $qend + $qoffset + ( $qstrand eq 1 ? 0 : ($sstart-1)*3 ), $qstrand, $score, $sseqid ) or next;

            # convert ORF to CDS only if they have been confirmed in the second round
            if ( ( $source_strand eq 2 ) and ( get_overlapped_features ( $o, $feature, "equals", ["ORF"] ) ) ) {
                my $new_qseqid = $qseqid;
                my $new_feature_strand = $feature->strand;
                $new_qseqid =~ s/_strand2/_strand$new_feature_strand/;
                $o->{"homology"}->{$new_qseqid} = 1;
            }
            elsif ( $source_strand eq 2 ) {
                next;
            }
            # generate the inference tag
            my $inference = "COORDINATES:alignment:swissprot:$sseqid";
            $feature->add_tag_value ("inference", $inference);
            #store CDS sequence feature
            $success++ if check_and_store_feature ( $o, $feature );
        }
    }
    # last iteration when no more success
    $o->{"cds-h"} = 0 if not $success;
    # else next iteration
    $o->{"cds-h"}++ if $success;
}

sub parse_annotation {
    my ( $o, $s ) = @_;
    print_log( $o, "Parsing annotation results..." );
    my $feature;
    # program-specific part
    if ( $o->{"cds-a-program"} eq "pannzer" ) {
        while ( my $l = parse_file( $o, $o->{"o"}."/output.desc", "line", "\t", $o->{"cds-a-program"} ) ) {
            # qpid    cluster_GSZ     cluster_RM1sum  cluster_size    cluster_desccount       RM2     val_avg jac_avg desc    genename
            my ( $query_id, $description, $symbol ) = ( $l->[0], clean_up_description($l->[8]), $l->[9] );
            if ( ! $feature or $feature->primary_id ne $query_id ) {
              store_feature ( $o, $feature ) if $feature;
              $feature = $o->{"dbh_uuid"}->fetch( $query_id );
              # remove all previous product tags, add the best hit from PANNZER
              $feature->remove_tag( "product" ) if $feature->has_tag("product");
              $feature->add_tag_value ( "product", $description );
            }
            print_verbose( $o, "Found 'product: $description' for query '$query_id'." );
        }

        while ( my $l = parse_file( $o, $o->{"o"}."/output.go", "line", "\t", $o->{"cds-a-program"} ) ) {
            # program-specific part
            if ( $o->{"cds-a-program"} eq "pannzer" ) {
                # qpid    ontology        goid    desc    ARGOT_score     ARGOT_PPV       ARGOT_rank
                my ( $query_id, $type, $id, $description ) = ( $l->[0], $l->[1], $l->[2], clean_up_description($l->[3]) );
                if ( ! $feature or $feature->primary_id ne $query_id ) {
                  store_feature ( $o, $feature ) if $feature;
                  $feature = $o->{"dbh_uuid"}->fetch( $query_id );
                }
                if ( $type eq "MF" ) {
                  $feature->add_tag_value ( "db_xref", "GO:$id" );
                  print_verbose( $o, "Found 'GO:$id, $description' for query '$query_id'." );
                }
                elsif ( $type eq "CC" ) {
                  $feature->add_tag_value ( "db_xref", "GO:$id" );
                  print_verbose( $o, "Found 'GO:$id, $description' for query '$query_id'." );

                }
                elsif ( $type eq "BP" ) {
                  $feature->add_tag_value ( "db_xref", "GO:$id" );
                  print_verbose( $o, "Found 'GO:$id, $description' for query '$query_id'." );
                }
            }
        }
    }

    elsif ( $o->{"cds-a-program"} eq "emapper" ) {
      while ( my $l = parse_file( $o, $o->{"o"}."/output.emapper.annotations", "line", "\t", $o->{"cds-a-program"} ) ) {
          my ( $query_id, $description ) = ( $l->[0], clean_up_description($l->[12]) );
          if ( ! $feature or $feature->primary_id ne $query_id ) {
            store_feature ( $o, $feature ) if $feature;
            $feature = $o->{"dbh_uuid"}->fetch( $query_id );
            # remove all previous product tags, add the best hit from emapper
            $feature->remove_tag( "product" ) if $feature->has_tag("product");
            $feature->add_tag_value ( "product", $description );
          }
          print_verbose( $o, "Found 'product: $description' for query '$query_id'." );
      }
    }
}



1;
