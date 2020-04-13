package SAP::Output;

use strict;
use warnings;
no warnings "experimental";
use SAP::Common;
use Bio::SeqIO;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(output_and_validate);

sub output_and_validate {
    my ( $o ) = @_;
    add_general_information ( $o );
    add_source_features ( $o );
    add_locus_numbering ( $o ) if $o->{"prefix"};
    add_gene_features ( $o );
    count_features ( $o );
    mark_ambiguities ( $o );
    # no more modification of features after this point, only writing sequences
    comment_output ( $o );
    write_fasta_output ( $o );
    write_table_output ( $o );
    write_gbf_and_sqn_output ( $o );
    merge_transparent_features ( $o ) if $o->{"transparent"};
    copy_gbf_to_gbk ( $o ) if ! $o->{"transparent"};
    run_antismash ( $o ) if $o->{"antismash"};
}

sub merge_transparent_features {
  my ( $o ) = @_;
  print_log ( $o, "Merging transparent features..." );
  my $input_filehandle = Bio::SeqIO->new(-file => $o->{"job"}."/output.gbf", -format => "genbank" );
  my $output_filehandle = Bio::SeqIO->new(-file => ">".$o->{"job"}."/output.gbk", -format => "genbank" );
  foreach my $record ( $input_filehandle->next_seq() ) {
    # combine "new" and transparent features
    my @features = ( $record->get_SeqFeatures(), @{ $o->{"r"}->{$record->id}->{"transparent"} } );
    # remove features once again
    $record->remove_SeqFeatures();
    # loop all features, has to be sorted because "new" and transparent are probably out of their natural order...
    foreach my $feature ( sort { $a->start<=>$b->start || $b->end<=>$a->end } @features ) {
      $record->add_SeqFeature( $feature );
    }
    # write the record
    $output_filehandle->write_seq( $record );
  }
}

sub copy_gbf_to_gbk {
  my ( $o ) = @_;
  # antismash does not understand gbf file extensions
  run_program ( $o, "cp ".$o->{"job"}."/output.gbf ".$o->{"job"}."/output.gbk" );
}

sub run_antismash {
    my ( $o ) = @_;
    print_log( $o, "Running antiSMASH ".$o->{"antismash-version"}."..." );
    run_program ( $o, "PATH=\$PATH:".$o->{"cwd"}."/bin/antismash_deps ".$o->{"cwd"}."/bin/antismash/run_antismash.py --genefinding-tool prodigal --cb-general --cb-subclusters --cb-knownclusters --minlength 1 --output-dir ".$o->{"job"}."/antismash ".$o->{"job"}."/output.gbk" );
}

sub add_general_information {
    my ( $o ) = @_;
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
  my ( $o ) = @_;
  foreach my $seq_id ( sort keys %{$o->{"r"}} ) {
    my $feature = create_feature ( $o, { primary_tag => "source",
                                         seq_id => $seq_id,
                                         start => 1,
                                         start_type => "EXACT",
                                         end => $o->{"r"}->{$seq_id}->length,
                                         end_type => "EXACT",
                                         strand => 0,
                                       } );
    $feature->add_tag_value ( "organism", $o->{"organism"} // "unidentified" );
    $feature->add_tag_value ( "strain", $o->{"strain"} ) if ( $o->{"strain"} );
    # check and store source feature candidate
    check_and_store_feature ( $o, $feature );
  }
}

sub add_locus_numbering {
    my ( $o ) = @_;
    print_log ( $o, "Assigning NCBI-compatible locus numbering..." );
    # counter
    my $last_locus_tag = 0;
    foreach my $seq_id ( sort keys %{$o->{"r"}} ) {
        # add locus numbers only for for CDS, ncRNA, rRNA, tmRNA and tRNA features
        my @features = grep { $_->primary_tag ~~ ["CDS", "ncRNA", "rRNA", "tmRNA", "tRNA"] } sort { $a->start<=>$b->start || $b->end<=>$a->end } $o->{"r"}->{$seq_id}->get_SeqFeatures();
        my $num_features = ( $#features + 1 ) * $o->{"step"};
                            #we have to sort them to keep proper track of locus_tags
                                     # we sort first by start (smaller first), if starts are equal we sort by end (smaller first)
        foreach my $feature ( sort { $a->start<=>$b->start || $a->end<=>$b->end } @features ) {

            # make sure all the CDS features have at least some product tag
            $feature->add_tag_value ( "product", "Hypothetical protein" ) if ( ( $feature->primary_tag eq "CDS" ) and not ( $feature->has_tag( "product" ) ) );

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
                # protein_id tag is only valid in conjunction with CDS features
                # see http://www.insdc.org/documents/feature-table#7.3.1
                if ( $feature->primary_tag eq "CDS" and $o->{"bioproject"} ) {
                  $feature->add_tag_value ( "protein_id", "gnl|".$o->{"bioproject"}."|$new_locus_tag" );
                }
            }
        }
    }
}

sub add_gene_features {
  my ( $o ) = @_;
  print_log ( $o, "Adding gene features..." );
  foreach my $seq_id ( sort keys %{$o->{"r"}} ) {
    foreach my $feature ( grep { $_->primary_tag ~~ ["CDS", "ncRNA", "rRNA", "tmRNA", "tRNA"] } sort { $a->start<=>$b->start || $b->end<=>$a->end } $o->{"r"}->{$seq_id}->get_SeqFeatures() ) {
      my $gene_feature = create_feature ( $o, { primary_tag => "gene",
                                                seq_id => $seq_id,
                                                location => $feature->location(),
                                               } );
      foreach my $tag ( grep { $_ ~~ ["pseudo", "locus_tag", "gene"] } $feature->get_all_tags() ) {
        foreach my $value ( $feature->get_tagset_values($tag) ) {
          $gene_feature->add_tag_value( $tag, $value );
        }
      }

      # check and store gene feature candidate
      check_and_store_feature ( $o, $gene_feature );
    }
  }
}

sub mark_ambiguities {
  my ( $o ) = @_;
  print_log ( $o, "Marking ambiguities..." );
  foreach my $seq_id ( sort keys %{$o->{"r"}} ) {
    my $seq = $o->{"r"}->{$seq_id}->seq();
    while ( $seq =~ /[^ATGCN]{1,}/g ) {
      # create a feature first
      my $feature = create_feature ( $o, {
                                           primary_tag => "unsure",
                                           seq_id => $seq_id,
                                           start => $-[0] + 1,
                                           start_type => "EXACT",
                                           end => $+[0],
                                           end_type => "EXACT",
                                           strand => 0,
                                         } );
      # check and store gap feature candidate
      check_and_store_feature ( $o, $feature );
    }

  }
}


sub count_features {
  my ( $o ) = @_;
  print_log ( $o, "Counting features..." );
  foreach my $seq_id ( sort keys %{$o->{"r"}} ) {
    foreach my $feature ( sort { $a->start<=>$b->start || $b->end<=>$a->end } $o->{"r"}->{$seq_id}->get_SeqFeatures() ) {
        ++$o->{"feature_count"}->{$feature->primary_tag};
        # if there is a "pseudo" tag in "gene feature", increment pseudo as well
        if ( ( $feature->primary_tag eq "gene" ) and ( $feature->has_tag('pseudo') ) ) {
            ++$o->{"pseudo_count"}->{"pseudo"};
        }
    }
  }
}

sub comment_output {
  my ( $o ) = @_;
  print_log ( $o, "Writing comments..." );
  # structured comment file
  open my $output_filehandle, ">", $o->{"job"}."/output.cmt" or die "Could not open ".$o->{"job"}."/output.cmt for writing - $!";

  my $comment = "StructuredCommentPrefix    \t         Genome-Annotation-Data".
              "\nAnnotation Provider        \t         HZI/HIPS".
              "\nWebsite                    \t         https://prosnap.helmholtz-hzi.de".
              "\nAnnotation Date            \t       ".current_date_and_time().
              "\nAnnotation Pipeline        \t       ".$o->{"program-name"}.
              "\nSoftware Version           \t       ".$o->{"program-version"}.
              "\nCDS                        \t       ".( $o->{"cds"} ? "yes" : "no" ).
              "\nCDS ab initio              \t       ".( $o->{"cds-i"} ? $o->{"cds-i-program"} : "no" ).
              "\nCDS by reference           \t       ".( $o->{"cds-h"} ? $o->{"cds-h-program"} : "no" ).
              "\nCDS annotation             \t       ".( $o->{"cds-a"} ? $o->{"cds-a-program"} : "no" ).
              "\nRNA                        \t       ".( $o->{"rna"} ? "yes" : "no" ).
              "\nAnnotation UUID            \t       ".$o->{"job_uuid"}.
              "\nAnnotated Features         \t       ".( join "; ", keys %{ $o->{"feature_count"} } ).
              "\nrRNA                       \t       ".($o->{"feature_count"}->{"rRNA"} // 0).
              "\ntRNA                       \t       ".($o->{"feature_count"}->{"tRNA"} // 0).
              "\ntmRNA                      \t       ".($o->{"feature_count"}->{"tmRNA"} // 0).
              "\nncRNA                      \t       ".($o->{"feature_count"}->{"ncRNA"} // 0).
              "\nCDS                        \t       ".($o->{"feature_count"}->{"CDS"} // 0).
              "\nTotal genes                \t       ".($o->{"feature_count"}->{"gene"} // 0).
              "\nregulatory                 \t       ".($o->{"feature_count"}->{"regulatory"} // 0).
              "\nmisc_feature               \t       ".($o->{"feature_count"}->{"misc_feature"} // 0).
              "\nStructuredCommentSuffix    \t         Genome-Annotation-Data";
  print {$output_filehandle} $comment;
}

sub write_fasta_output {
    my ( $o ) = @_;
    print_log ( $o, "Writing .fsa output file..." );
    my $output_filehandle = Bio::SeqIO->new(-file => ">".$o->{"job"}."/output.fsa", -format => "fasta" );
    foreach my $seq_id ( sort keys %{$o->{"r"}} ) {
        #add some data that is needed later by tbl2asn and also if *.fsa is required for the submission (BankIt)
        $o->{"r"}->{$seq_id}->description(
                                            ("[organism=".($o->{"organism"} ? $o->{"organism"} : "unidentified")."] ").
                                          "\n[gcode=11] ".
                                          "\n[division=BCT] ".
                                            ($o->{"classification"} ? "[lineage=".$o->{"classification"}."] " : "").
                                            ($o->{"r"}->{$seq_id}->is_circular ? "[topology=circular] " : "[topology=linear] ").
                                            ($o->{"complete"} ? "[completeness=complete] " : "").
                                            # reuse the description without changes in transparent mode, if it exists
                                            # ( $o->{"transparent"} and $o->{"r"}->{$seq_id}->description ? $o->{"r"}->{$seq_id}->description :
                                              (
                                                  ($o->{"organism"} ? $o->{"organism"} : "Unidentified organism").
                                                  ($o->{"strain"} ? " strain ".$o->{"strain"} : "").
                                                "\n, ".
                                                  ($o->{"complete"} ? "complete genome" : "draft genome")
                                              )
                                            # )

                                            );
        $output_filehandle->write_seq( $o->{"r"}->{$seq_id} );
    }
}

sub write_table_output {
    my ( $o ) = @_;
    print_log ( $o, "Writing Sequin/BankIt .tbl output file..." );
    open my $output_filehandle, ">", $o->{"job"}."/output.tbl" or die "Could not open ".$o->{"job"}."/output.tbl for writing - $!";
    foreach my $seq_id ( sort keys %{$o->{"r"}} ) {
        print {$output_filehandle} ">Feature ".$o->{"r"}->{$seq_id}->id."\n";
        foreach my $feature ( $o->{"r"}->{$seq_id}->get_SeqFeatures ) {
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
    my ( $o ) = @_;
    print_log ( $o, "Writing Sequin .sqn output file..." );
    run_program ( $o, $o->{"cwd"}."/bin/tbl2asn/tbl2asn -i ".$o->{"job"}."/output.fsa -w ".$o->{"job"}."/output.cmt -o ".$o->{"job"}."/output.sqn -V vb -s T -a r10k" );
}

1;
