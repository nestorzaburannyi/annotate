package SAP::Output;

use strict;
use warnings;
use SAP::Common;
use Bio::SeqIO;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(output_and_validate);

sub output_and_validate {
    my ( $o, $s ) = @_;
    add_general_information ( $o, $s );
    add_source_features ( $o, $s );
    add_locus_numbering ( $o, $s ) if $o->{"prefix"};
    add_gene_features ( $o, $s );
    move_features_to_sequence_hash ( $o, $s );
    count_features ( $o, $s );
    # no more modification of features after this point, only writing sequences
    write_cmt_output ( $o, $s );
    write_fasta_output ( $o, $s );
    write_table_output ( $o, $s );
    write_gbf_and_sqn_output ( $o, $s );
    run_antismash ( $o, $s ) if $o->{"antismash"};
}

sub run_antismash {
    my ( $o, $s ) = @_;
    print_log ( $o, "Running optional antismash tool..." );
    run_program ( $o, "cp ".$o->{"o"}."/output.gbf ".$o->{"o"}."/output.gbk" );
    run_program ( $o, "PATH=\$PATH:".$o->{"cwd"}."/bin/antismash_deps ".$o->{"cwd"}."/bin/antismash/run_antismash.py --genefinding-tool prodigal --cb-general --cb-subclusters --cb-knownclusters --minlength 1 -c 1 --output-dir ".$o->{"o"}."/antismash ".$o->{"o"}."/output.gbk" );
}

sub add_general_information {
    my ( $o, $s ) = @_;
    foreach my $seq_id (sort keys %$s) {
        # set id to sequence + number
        $s->{$seq_id}->id( "sequence_$seq_id" );
    }
    # set the ORGANISM/SOURCE lines if it was specified
    if ( $o->{"taxid"} ) {
        my $node = $o->{"dbh_taxonomy"}->get_taxon( -taxonid => $o->{"taxid"} );
        $o->{"organism"} = $node->node_name;
        $o->{$node->rank} = $node->node_name;
        while ( $node = $node->ancestor ) {
            # skip unknown ranks
            next if ( $node->rank eq "no rank" );
            $o->{"classification"} = $node->node_name . ( $o->{"classification"} ? "; ". $o->{"classification"} : "" );
            $o->{$node->rank} = $node->node_name;
        }
    }
}

sub add_source_features {
    my ( $o, $s ) = @_;
    foreach my $seq_id (sort keys %$s) {
        my $feature = create_feature ( "source", $seq_id, 1, "EXACT", $s->{$seq_id}->length, "EXACT", 0, 0 );
        $feature->add_tag_value ( "organism", $o->{"organism"} // "unidentified" );
        $feature->add_tag_value ( "strain", $o->{"strain"} ) if ( $o->{"strain"} );
        # $feature->add_tag_value ( "db_xref", "taxon:".$o->{"taxid"} ) if ( $o->{"taxid"} );
        #store source feature
        check_and_store_feature ( $o, $feature );
    }
}

sub add_locus_numbering {
    my ( $o, $s ) = @_;
    print_log ( $o, "Assigning NCBI-compatible locus numbering..." );
    # counter
    my $last_locus_tag = 0;
    foreach my $seq_id (sort keys %$s) {
        # add locus numbers only for for CDS, ncRNA, rRNA, tmRNA and tRNA features
        my @features = $o->{"dbh_uuid"}->features(-seq_id => $seq_id, -type => ["CDS", "ncRNA", "rRNA", "tmRNA", "tRNA"]);
        my $num_features = ( $#features + 1 ) * $o->{"step"};
                            #we have to sort them to keep proper track of locus_tags
                                     # we sort first by start (smaller first), if starts are equal we sort by end (smaller first)
        foreach my $feature ( sort { $a->start<=>$b->start || $a->end<=>$b->end } @features ) {
            # do not add locus_tag if there is one already (e.g. derived from features in transparent mode)
            if ( $feature->has_tag("locus_tag") ) {
                # but make note of the number if possible
                my $digit_part = (( $feature->get_tag_values("locus_tag"))[0] =~ m/(\d+)(\D*)$/)[0];
                if ( $digit_part =~ m/^\d+$/) {
                    $last_locus_tag = $digit_part;
                }
                next;
            }
            # if there is no locus_tag
            else {
                $last_locus_tag = $last_locus_tag + $o->{"step"};
                my $new_locus_tag = $o->{"prefix"}."_".add_leading_zeros( $last_locus_tag, $o->{"locus_digit"} );
                # seach for new allowed locus tag until successful
                while ( exists $o->{"forbidden_locus_tags"}->{$new_locus_tag} ) {
                    $num_features = $num_features + $o->{"step"};
                    $new_locus_tag = $o->{"prefix"}."_".add_leading_zeros( $num_features, $o->{"locus_digit"} );
                }
                # forbid the use of the same locus_tag again
                $o->{"forbidden_locus_tags"}->{$new_locus_tag} = 1;
                $feature->add_tag_value ( "locus_tag", $new_locus_tag);
                $feature->add_tag_value ( "protein_id", "gnl|".$o->{"bioproject"}."|$new_locus_tag" );
                #update feature
                store_feature ( $o, $feature );
            }
        }
    }
}

sub add_gene_features {
    my ( $o, $s ) = @_;
    print_log ( $o, "Adding gene features..." );
    foreach my $feature_to_copy ( $o->{"dbh_uuid"}->features(-type => ["CDS", "ncRNA", "rRNA", "tmRNA", "tRNA"]) ) {
        # if there is a "pseudo" tag, we do not copy. Instead, we change the primary tag to "gene"
        if ( $feature_to_copy->has_tag("pseudo") ) {
            # change the primary_tag to "gene"
            $feature_to_copy->primary_tag("gene" );
            # change the "product" tag(s?) to "note"
            if ( $feature_to_copy->has_tag("product") ) {
                foreach my $product_tag ( $feature_to_copy->get_tagset_values("product") ) {
                    $feature_to_copy->add_tag_value( "note", $product_tag." pseudogene" );
                }
                $feature_to_copy->remove_tag( "product" );
            }
            # store gene feature
            check_and_store_feature ( $o, $feature_to_copy );
        }
        else {
            # we do not want to modify already present features, so we copy it
            my $feature = clone_feature( $o, $feature_to_copy );
            foreach my $tag ( $feature->get_all_tags ) {
                # keep gene tags
                next if $tag eq "gene";
                # keep locus tags
                next if $tag eq "locus_tag";
                # remove the rest, unwanted tags
                $feature->remove_tag( $tag );
            }
            # change the primary_tag to "gene"
            $feature->primary_tag( "gene" );
            # store gene feature
            check_and_store_feature ( $o, $feature );
        }
    }
}

sub move_features_to_sequence_hash {
    my ( $o, $s ) = @_;
    print_log ( $o, "Populating sequences with features..." );
    foreach my $seq_id ( sort keys %$s ) {
        # we have to sort the feature, otherwise they won't be in ascending numerical order in the output files
        foreach my $feature ( sort { $a->start<=>$b->start || $b->end<=>$a->end } $o->{"dbh_uuid"}->features(-seq_id => $seq_id ) ) {
            # removing the internal database ID tags
            $feature->remove_tag( "ID" );
            # removing the score tags, not needed any more
            $feature->remove_tag( "score" ) if $feature->has_tag( "score" );
            # make sure all the CDS features have at least some product tag
            $feature->add_tag_value ( "product", "Hypothetical protein" ) if ( ( $feature->primary_tag eq "CDS" ) and not ( $feature->has_tag( "product" ) ) );
            print_verbose ( $o, "Adding feature on sequence ".$feature->seq_id.", start ".$feature->start.", end ".$feature->end.", strand ".$feature->strand );
            $s->{$seq_id}->add_SeqFeature( $feature );
        }
    }
}

sub count_features {
    my ( $o, $s ) = @_;
    print_log ( $o, "Counting features..." );
    foreach my $feature_to_count ( $o->{"dbh_uuid"}->features ) {
        ++$o->{"feature_count"}->{$feature_to_count->primary_tag};
        # if there is a "pseudo" tag in "gene feature", increment pseudo as well
        if ( ( $feature_to_count->primary_tag eq "gene" ) and ( $feature_to_count->has_tag('pseudo') ) ) {
            ++$o->{"pseudo_count"}->{"pseudo"};
        }
    }
}

sub write_cmt_output {
    my ( $o, $s ) = @_;
    print_log ( $o, "Writing statistics .cmt output file..." );
    my $all_feature_types = join "; ", keys %{ $o->{"feature_count"} };
    my $current_date_and_time = current_date_and_time;
    open my $output_filehandle, ">", $o->{"o"}."/output.cmt" or die "Could not open ".$o->{"o"}."/output.cmt for writing - $!";
    print {$output_filehandle} "##Genome-Annotation-Data-START##
Annotation Provider          :: HZI/HIPS
Annotation Date              :: ".$current_date_and_time."
Annotation Pipeline          :: ".$o->{"program-name"}."
Annotation Software revision :: ".$o->{"program-version"}."
Feature Types Annotated      :: ".$all_feature_types."
rRNA                         :: ".($o->{"feature_count"}->{"rRNA"} // 0)."
tRNA                         :: ".($o->{"feature_count"}->{"tRNA"} // 0)."
tmRNA                        :: ".($o->{"feature_count"}->{"tmRNA"} // 0)."
ncRNA                        :: ".($o->{"feature_count"}->{"ncRNA"} // 0)."
CDS                          :: ".($o->{"feature_count"}->{"CDS"} // 0)."
Total gene                   :: ".($o->{"feature_count"}->{"gene"} // 0)."
##Genome-Annotation-Data-END##";
    close( $output_filehandle);
}

sub write_fasta_output {
    my ( $o, $s ) = @_;
    print_log ( $o, "Writing .fsa output file..." );
    my $output_filehandle = Bio::SeqIO->new(-file => ">".$o->{"o"}."/output.fsa", -format => "fasta" );
    my $counter = 1;
    foreach my $seq_id (sort keys %$s) {
        #add some data that is needed later by tbl2asn and also if *.fsa is required for the submission (BankIt)
        $s->{$seq_id}->description(
                                            ("[organism=".($o->{"organism"} ? $o->{"organism"} : "unidentified")."] ").
                                            "[gcode=11] ".
                                            "[division=BCT] ".
                                            ($o->{"classification"} ? "[lineage=".$o->{"classification"}."] " : "").
                                            (( $o->{"circular"}->{$counter} // 0 ) ? "[topology=circular] " : "[topology=linear] ").
                                            ($o->{"complete"} ? "[completeness=complete] " : "").
                                            # reuse the description without changes in transparent mode, if it exists
                                            # ( $o->{"transparent"} and $s->{$seq_id}->description ? $s->{$seq_id}->description :
                                              (
                                                  ($o->{"organism"} ? $o->{"organism"} : "Unidentified organism").
                                                  ($o->{"strain"} ? " strain ".$o->{"strain"} : "").
                                                  ", ".
                                                  ($o->{"complete"} ? "complete genome" : "draft genome")
                                              )
                                            # )

                                            );
        $output_filehandle->write_seq( $s->{$seq_id} );
    }
}

sub write_table_output {
    my ( $o, $s ) = @_;
    print_log ( $o, "Writing Sequin/BankIt .tbl output file..." );
    open my $output_filehandle, ">", $o->{"o"}."/output.tbl" or die "Could not open ".$o->{"o"}."/output.tbl for writing - $!";
    foreach my $seq_id (sort keys %$s) {
        print {$output_filehandle} ">Feature ".$s->{$seq_id}->id."\n";
        foreach my $feature ( $s->{$seq_id}->get_SeqFeatures ) {
            # $end is first in case of -1 strand, in all other cases "$start" should be first
            print {$output_filehandle} ( $feature->strand eq -1 ? ( $feature->location->end_pos_type eq "AFTER" ? "<".$feature->end : $feature->end ) : ( $feature->location->start_pos_type eq "BEFORE" ? "<".$feature->start : $feature->start ) )."\t".( $feature->strand eq -1 ? ( $feature->location->start_pos_type eq "BEFORE" ? ">".$feature->start : $feature->start ) : ( $feature->location->end_pos_type eq "AFTER" ? ">".$feature->end : $feature->end ) )."\t".$feature->primary_tag."\n";
            foreach my $tag_type ( $feature->get_all_tags ) {
                foreach my $tag_instance ( $feature->get_tag_values( $tag_type )) {
                    next if ( ! $tag_instance );
                    print {$output_filehandle} "\t\t\t$tag_type".( $tag_instance eq "_no_value" ? "\n" : "\t$tag_instance\n" );
                }
            }
        }
    }
    close( $output_filehandle );
}

sub write_gbf_and_sqn_output {
    my ( $o, $s ) = @_;
    print_log ( $o, "Writing Sequin .sqn output file..." );
    run_program ( $o, $o->{"cwd"}."/bin/tbl2asn/tbl2asn -T T -i ".$o->{"o"}."/output.fsa -w ".$o->{"o"}."/output.cmt -o ".$o->{"o"}."/output.sqn -V vb -s T -Z ".$o->{"o"}."/output.dis" );
}

1;
