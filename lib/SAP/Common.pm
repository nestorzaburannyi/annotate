package SAP::Common;
use strict;
use warnings;
use Bio::DB::SeqFeature::Store;
use Bio::DB::Taxonomy;
use Bio::Location::Fuzzy;
use Bio::Tools::CodonTable;
use Bio::Tools::GuessSeqFormat;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Data::UUID;
use File::Path qw(make_path);
use DBI;
use DBD::SQLite::Constants qw/:file_open/;
use Getopt::Long;
use Digest::MD5  qw(md5_hex);

# suppress warnings from bioperl
$SIG{"__WARN__"} = sub { warn $_[0] unless ( caller eq "Bio::Root::RootI" ) };

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
                get_command
                run_program
                get_nucleotide_sequence_of_feature
                get_protein_sequence_of_feature
                current_date_and_time
                initialize_options
                initialize_pipeline
                print_log
                print_verbose
                exit_program
                create_seq_hash
                prepare_input
                parse_file
                create_feature
                clone_feature
                store_feature
                check_and_store_feature
                add_leading_zeros
                find_appropriate_cds_feature
                get_overlapped_features
                clean_up_description
            );

######################################################## PARAMETERS PART ########################################################
sub initialize_options {
    my ( $cwd ) = @_;
    #setting default values
    my $o = {
        #modules
            #submodules
                #parameters
        #RNA gene prediction
        "rna" => 1,
            #rRNA
            "rna-r" => 1,
                "rna-r-program" => "rnammer",
                "rna-r-score" => -inf,
            #tRNA
            "rna-t" => 1,
                "rna-t-program" => "trnascanse",
                "rna-t-score" => -inf,
            #tmRNA
            "rna-tm" => 1,
                "rna-tm-program" => "aragorn",
                "rna-tm-score" => -inf,
            #ncRNA
            "rna-nc" => 1,
                "rna-nc-program" => "infernal",
                "rna-nc-score" => -inf,
        #CDS gene prediction
        "cds" => 1,
            #de novo
            "cds-i" => 1,
                "cds-i-program" => "prodigal",
                "cds-i-score" => 20,
            #homology
            "cds-h" => 1,
                "cds-h-program" => "diamond",
                "cds-h-score" => 40,
                "cds-h-identity" => 90,
            #annotation
            "cds-a" => 1,
                "cds-a-program" => "pannzer",
        #additional parameters
        "program-name" => "Sequence Annotation Pipeline",
        "program-version" => "20200106",
        "rnammer-version" => "1.2",
        "trnascanse-version" => "2.0.5",
        "aragorn-version" => "1.2.38",
        "infernal-version" => "1.1.3",
        "glimmer-version" => "3.02",
        "prodigal-version" => "2.6.3",
        "genemarks-version" => "4.30",
        "blast-version" => "2.6.0+",
        "diamond-version" => "0.8.36",
        "pannzer-version" => "2.0.0",
        "emapper-version" => "1.0.3",
        "pilercr-version" => "1.06",
        "minced-version" => "0.4.2",
        "antismash-version" => "5.1.0",
        "trf-version" => "4.09",
        "locus_count" => "0",
        "locus_digit" => "5",
        "step" => "10",
        "complete" => "0",
	  };
        #getting options
    GetOptions (
        #modules
            #submodules
                #parameters
        #mandatory
        "input|i=s"                                         => \$o->{"i"},
        #update databases
        "update"                                            => \$o->{"update"},
        #RNA gene prediction
        "rna!"                                              => \$o->{"rna"},
            #rRNA
            "rna-r!"                                        => \$o->{"rna-r"},
                "rna-r-program=s"                           => \$o->{"rna-r-program"},
                "rna-r-score=s"                             => \$o->{"rna-r-score"},
            #tRNA
            "rna-t!"                                        => \$o->{"rna-t"},
                "rna-t-program=s"                           => \$o->{"rna-t-program"},
                "rna-t-score=s"                             => \$o->{"rna-t-score"},
            #tmRNA
            "rna-tm!"                                       => \$o->{"rna-tm"},
                "rna-tm-program=s"                          => \$o->{"rna-tm-program"},
                "rna-tm-score=s"                            => \$o->{"rna-tm-score"},
            #ncRNA
            "rna-nc!"                                       => \$o->{"rna-nc"},
                "rna-nc-program=s"                          => \$o->{"rna-nc-program"},
                "rna-nc-score=s"                            => \$o->{"rna-nc-score"},
        #CDS gene prediction
        "cds!"                                              => \$o->{"cds"},
            #ab initio
            "cds-i!"                                        => \$o->{"cds-i"},
                "cds-i-program=s"                           => \$o->{"cds-i-program"},
                "cds-i-score=s"                             => \$o->{"cds-i-score"},
            #homology
            "cds-h!"                                        => \$o->{"cds-h"},
                "cds-h-program=s"                           => \$o->{"cds-h-program"},
                "cds-h-score=s"                             => \$o->{"cds-h-score"},
                "cds-h-identity=s"                          => \$o->{"cds-h-identity"},
            #annotation
            "cds-a!"                                        => \$o->{"cds-a"},
                "cds-a-program=s"                           => \$o->{"cds-a-program"},
                "cds-a-score=s"                             => \$o->{"cds-a-score"},
        #meta data
        "taxid=s"                                           => \$o->{"taxid"},
        "strain=s"                                          => \$o->{"strain"},
        "prefix=s"                                          => \$o->{"prefix"},
        "bioproject=s"                                      => \$o->{"bioproject"},
        "uuid=s"                                            => \$o->{"uuid"},
        "complete"                                          => \$o->{"complete"},
        "circular=s"                                        => \$o->{"circular"},
        "step=i"                                            => \$o->{"step"},
        #troubleshooting and general options
        "help|h"                                            => \$o->{"help"},
        "version"                                           => \$o->{"version"},
        "quiet"                                             => \$o->{"quiet"},
        "verbose"                                           => \$o->{"verbose"},
        "transparent"                                       => \$o->{"transparent"},
        "antismash"                                         => \$o->{"antismash"},
        "output|o=s"                                        => \$o->{"output"},
    ) or exit;

sub initialize_pipeline {
    my ( $o, $cwd ) = @_;
    # set the option to cwd
    $o->{"cwd"} = $cwd;
    # print help and exit
    print_help($o) if $o->{"help"};
    # print version and exit
    print_version($o) if $o->{"version"};

    # generate unique uuid for temp job folder, or reuse the specified one
    $o->{"job_uuid"} = $o->{"uuid"} || Data::UUID->new->create_str;

    # the the base for output folder for files
    $o->{"job"} = $o->{"cwd"}."/public/jobs/".$o->{"job_uuid"};
    # create folder for the job if does not exits
    make_path($o->{"job"});
    # change the current working folder to the created one, as some tools (e.g. GeneMarks, rnammer) will spam and create temporary files in the cwd
    chdir $o->{"job"};
    # no input file specifed
    exit_program($o, "Input file specified via -input parameter is required.") if not $o->{"i"};
    # input file does not exist
    exit_program($o, "Required file ".$o->{"i"}." does not exist.") if not -e $o->{"i"};
    # guess format or exit
    $o->{"input_format"} = Bio::Tools::GuessSeqFormat->new(-file => $o->{"i"})->guess;
    exit_program($o, "Submitted file format is not supported. Please submit nucleotide FASTA or nucleotide-containing GenBank") if not ( $o->{"input_format"} =~ m/^genbank|fasta$/ );
    #print some starting info about the run
    print_log($o, "Starting ".$o->{"program-name"}.", build ".$o->{"program-version"} );
    print_log($o, "Job ID is ".$o->{"job_uuid"} );

    # auto-generate prefix if not specified
    if ( not $o->{"prefix"} ) {
      $o->{"prefix"} = "BOGUS";
      print_log($o, "No locus tag prefix specified, using ".$o->{"prefix"} );
    }
    else {
      print_log($o, "Locus tag prefix is ".$o->{"prefix"} );
    }

    print_log($o, "Strain name is ".$o->{"strain"} ) if $o->{"strain"};
    print_log($o, "BioProject is ".$o->{"bioproject"} ) if $o->{"bioproject"};

    # set the circular molecules
    $o->{"circular"} = { map { $_ => 1 } split ",", $o->{"circular"} } if $o->{"circular"};
    # generating codon_table object for reuse later
    $o->{"codon_table"} = Bio::Tools::CodonTable->new( -id => $o->{"gcode"} );
    $o->{"default_codon_table_id"} = $o->{"codon_table"}->add_table("prokaryotic",
                                                                                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                                                                                "---M-------------------------------M---------------M------------");
    print_log($o, "Using the ".$o->{"codon_table"}->name." genetic code (transl_table=".$o->{"codon_table"}->id."):" );
    $o->{"dbh_taxonomy"} = Bio::DB::Taxonomy->new( -source => "flatfile", -directory => $o->{"cwd"}."/databases/taxonomy", -nodesfile => $o->{"cwd"}."/databases/taxonomy/nodes.dmp", -namesfile => $o->{"cwd"}."/databases/taxonomy/names.dmp" );
    $o->{"forbidden_locus_tags"} = ();
    return $o;

}

sub print_help {
    my ( $o ) = @_;
    print '
Usage: annotate -input <input_file> [OPTIONS]
Example: annotate -input my_genome.fasta -taxid 13

        Parameters              Type            Default                                 Description

Mandatory parameters:
        -input                  string          -                                       input file name (FASTA or GenBank format)

Optional parameters:
        -rna                    boolean         true                                    prediction of RNA genes
        -rna-r                  boolean         true                                    prediction of ribosomal RNA genes
        -rna-t                  boolean         true                                    prediction of transport RNA genes
        -rna-tm                 boolean         true                                    prediction of transport-messenger RNA genes
        -rna-nc                 boolean         true                                    prediction of noncoding RNA genes
        -cds                    boolean         true                                    prediction of CDS genes
        -cds-i                  boolean         true                                    prediction of CDS genes (ab initio)
        -cds-h                  boolean         true                                    prediction of CDS genes (homology-based)
        -cds-a                  boolean         true                                    annotation of CDS genes
        -taxid                  integer         -                                       NCBI taxonomy id of the organism (e.g. 52)
        -strain                 string          -                                       strain name
        -prefix                 string          -                                       prefix for locus numbering (NCBI-compatibility)
        -bioproject             string          -                                       a BioProject identifier (NCBI-compatibility)
        -step                   integer         10                                      increment of locus numbering (NCBI-compatibility)
        -sbt                    string          -                                       file with sbt formular (NCBI-compatibility)
        -transparent            boolean         false                                   transparent mode (keep already annotated features)
        -help                   boolean         false                                   print this help message
        -quiet                  boolean         false                                   suppress messages to the console
        -verbose                boolean         false                                   print more messages to the console
        -output                 string          current working directory               directory name for output files
        -uuid                   string          output                                  filename prefix for output files

Advanced parameters:
        -rna-r-program          string          rnammer                                 rnammer, infernal
        -rna-t-program          string          trnascanse                              trnascanse, aragorn, infernal
        -rna-tm-program         string          aragorn                                 aragorn, infernal
        -rna-nc-program         string          infernal                                infernal
        -cds-i-program          string          prodigal                                prodigal, glimmer, genemarks
        -cds-i-score            string          5                                       score threshold for ab initio predictions
        -cds-h-program          string          blast                                   blast, diamond
        -cds-h-evalue           string          1e-5                                    e-value threshold for homology searches
        -cds-a-program          string          pannzer                                 pannzer, emapper

Note: you can disable modules or submodules of the program by prepending "-no" to the respective option: "-no-cds" will cause the pipeline not to predict any CDS
';
    exit (1);
}

########################################################## IO TERM ###########################################################

sub print_message {
  my ( $o, $message ) = @_;
  return "[".current_date_and_time()."] $message\n";
}

sub print_log {
    # turn off buffering
    local $| = 1;
    my ( $o, $message ) = @_;
    # print to console only if not quiet
    if ( ! $o->{"quiet"} ) {
      print print_message( $o, $message );
    }
    # print to log
    open my $output_filehandle, ">>", $o->{"job"}."/log" or die "Could not open ".$o->{"job"}."/log for writing - $!";
    print {$output_filehandle} print_message( $o, $message );
      close $output_filehandle;
      # print to verbose log
      open $output_filehandle, ">>", $o->{"job"}."/verbose" or die "Could not open ".$o->{"job"}."/verbose for writing - $!";
      print {$output_filehandle} print_message( $o, $message );
      close $output_filehandle;
    }

sub print_verbose {
    my ( $o, $message ) = @_;
    # print to console only in verbose mode
    if ( $o->{"verbose"} ) {
      print print_message( $o, $message );
    }
      open my $output_filehandle, ">>", $o->{"job"}."/verbose" or die "Could not open ".$o->{"job"}."/verbose for writing - $!";
      print {$output_filehandle} print_message( $o, $message );
      close $output_filehandle;
    }

sub run_program {
    my ( $o, $command ) = @_;
    # show full command verbose
    print_verbose( $o, "Command used: $command" );
    # this is the only system call allowed, always wrap calls in run_program()
    # tee everything to verbose log
    system ( $command . "| tee -a ".$o->{"job"}."/verbose 1>/dev/null 2>/dev/null" );
}

sub exit_program {
    my ( $o, $message) = @_;
    print_log( $o, $message );
    exit ( 1 );
}

sub current_date_and_time {
    my @t = map { sprintf("%02d", $_) } localtime(time);
    return ( $t[5] + 1900)."-". sprintf( "%02d", ( $t[4] + 1 ) ) ."-$t[3] $t[2]:$t[1]:$t[0]";
}

sub add_leading_zeros {
    my ( $number, $digits ) = @_;
    return sprintf("%0${digits}d", $number );
}

################################################### SEQUENCE AND FEATURE ##################################################
sub create_feature {
    my ( $o, $l ) = @_;

    my $feature = Bio::SeqFeature::Generic->new(-location => $l->{"location"} // Bio::Location::Fuzzy->new( -start => $l->{"start"},
                                                                                                            -end => $l->{"end"},
                                                                                                            -strand => $l->{"strand"},
                                                                                                            -start_fuz => $l->{"start_type"},
                                                                                                            -end_fuz => $l->{"end_type"}
                                                                                                            ),
                                                -seq_id => $l->{"seq_id"},
                                                -primary => $l->{"primary_tag"},
                                                );

    # attach sequence to a feature
    $feature->attach_seq( $o->{"r"}->{$l->{"seq_id"}} );

    # add custom annotations
    $feature->annotation( {
                          "gene_id" => $l->{"gene_id"}, # to track duplicate hits to the same gene
                          "hit_id" => $l->{"hit_id"}, # to track special feature exceptions
                          "exception" => $l->{"exception"} // undef,
                          "score" => $l->{"score"},
                          "type" => $l->{"type"} // undef,
                          "method" => $l->{"method"} // undef,
                          "unique_id" => ++$o->{"unique_id"}, # just an unique counter
    } );

    # set initial tags provided by the generator
    $feature->set_attributes( -tag => $l->{"tags"} );

    return $feature;

}

sub check_and_store_feature {
  my ( $o, $feature ) = @_;
  # perform feature checks before any feature storage
  return if feature_rules ( $o, $feature );
  # perform overlap checks before any feature storage
  return if overlap_rules ( $o, $feature );

sub delete_feature {
  # removes the feature from the feature pool
  # accepts: a Bio::SeqFeature::Generic object to be deleted
  # returns: undef
  my ( $o, $outgoing_feature ) = @_;
  if ( ! $outgoing_feature->isa("Bio::SeqFeature::Generic") ) {
    exit_program( $o, "outgoing feature is not a Bio::SeqFeature::Generic object" );
  }
  if ( $outgoing_feature->has_tag("locus_tag") ) {
    foreach my $locus_tag_value ( $outgoing_feature->get_tag_values("locus_tag") ) {
      if ( exists $o->{"forbidden_locus_tags"}->{$locus_tag_value} ) {
        delete $o->{"forbidden_locus_tags"}->{$locus_tag_value};
      }
    }
  }
  $o->{"dbh_uuid"}->delete( $outgoing_feature ) or die "Could not delete feature";
}

sub clone_feature {
  # clones a feature, that is creates an independently mutable copy
  # accepts: a Bio::SeqFeature::Generic object to be cloned
  # returns: a cloned Bio::SeqFeature::Generic object
  my ( $o, $source_feature ) = @_;
  $source_feature->isa("Bio::SeqFeature::Generic") or die "source feature is not a Bio::SeqFeature::Generic object";
  # create a feature wuth the same location parameters as the source feature
  my $new_feature = create_feature ( $source_feature->primary_tag,
                                     $source_feature->seq_id,
                                     $source_feature->start,
                                     $source_feature->location->start_pos_type,
                                     $source_feature->end,
                                     $source_feature->location->end_pos_type,
                                     $source_feature->strand,
                                     $source_feature->score );
  # transfer all the tags from the source feature to its clone
  foreach my $source_tag_type ( $source_feature->get_all_tags ) {
    foreach my $source_tag_value ( $source_feature->get_tag_values( $source_tag_type ) ) {
      $new_feature->add_tag_value( $source_tag_type, $source_tag_value );
    }
  }
  $new_feature->isa("Bio::SeqFeature::Generic") or die "new feature is not a Bio::SeqFeature::Generic object";
  return $new_feature;
}

sub get_overlapped_features {
  # checks if a feature overlaps / contains certain feature types
  # accepts: $new_feature - a feature to check for overlaps
  #          $rule - rule, according to which the check is being performed
  #          $source_features - a list of one or more feature types to include in checking
  # returns: a list of features that overlap the feature in question
  my ( $o, $new_feature, $rule, $source_features ) = @_;
  my $source_feature_iterator = $o->{"dbh_uuid"}->get_seq_stream( -seq_id => $new_feature->seq_id,
                                                                  -start => $new_feature->start,
                                                                  -end => $new_feature->end,
                                                                  -type => $source_features );
  my @return_features;
  while (my $source_feature = $source_feature_iterator->next_seq) {
      if ( $new_feature->$rule( $source_feature, "ignore") ) {
          push @return_features, $source_feature;
      }
  }
  return @return_features;
}



sub create_seq_hash {
    my ( $o ) = @_;
    # to decide on sequence numbering below
    my $tmp_filehandle = Bio::SeqIO->new(-file => $o->{"i"}, -format => $o->{"input_format"} );
    while ( my $input_record = $tmp_filehandle->next_seq ) {
        $o->{"input_count"}++;
        $o->{"input_size"} = ( $o->{"input_size"} ? $o->{"input_size"} : 0 ) + $input_record->length;
    }
    $o->{"sequence_digit"} = length( $o->{"input_count"} );
    print_log( $o, "Input has ".$o->{"input_count"}." sequences, setting sequence numbering to ".$o->{"sequence_digit"}."-digit" );
    # the real part
    my $input_filehandle = Bio::SeqIO->new(-file => $o->{"i"}, -format => $o->{"input_format"} );
    my %s;
    while ( my $input_record = $input_filehandle->next_seq ) {
      my $counter = ( keys %s ) + 1;
        # we have to change the name in the very beginning, otherwise transparency won"t work
        $input_record->id(add_leading_zeros($counter, $o->{"sequence_digit"}));

        # set the sequence to circular if specified
        $input_record->is_circular( $o->{"circular"}->{$counter} );

        # transparent
        if ( $o->{"transparent"} ) {
            foreach my $feature ( $input_record->get_SeqFeatures ) {
                # set organism taxonomy based on db_xref tag of source feature
                if ( $feature->primary_tag eq "source" ) {
                    if ( not $o->{"taxid"} ) {
                        foreach my $db_xref ( $feature->get_tagset_values( "db_xref" ) ) {
                            next if not ( $db_xref =~ m/taxon:(\d+)/ );
                            $o->{"taxid"} = $1;
                        }
                    }
                    if ( not $o->{"strain"} ) {
                        $o->{"strain"} = ( $feature->get_tagset_values( "strain" ) )[0];
                    }
                }
                print_verbose ( $o, "Transparent mode, keeping ".$feature->primary_tag." (".$feature->start."..".$feature->end.")" );
                my $transparent_feature = clone_feature( $o, $feature );
                $transparent_feature->seq_id( $input_record->id );
                store_feature ( $o, $transparent_feature );
            }
        }
        # remove features anyway, transparent or not -> they will be stored in the sqlite DB only
        $input_record->remove_SeqFeatures;
        # now store naked sequence in the database, to be able to query for id:start..end features
        $o->{"dbh_uuid"}->insert_sequence( $input_record->id, $input_record->seq );
        $s{$input_record->id} = $input_record;
    }

    $o->{"dbh_uuid"}->commit;
    $o->{"locus_digit"} = int (log( 100000 )/log(10) - 2);
    print_log( $o, "Input is ".$o->{"input_size"}." bp long, setting locus numbering to ".$o->{"locus_digit"}."-digit" );
    return wantarray ? %s : \%s;
}

sub overlap_rules {
    my ( $o, $check_feature) = @_;

    # catch-all rule for duplicates, no need to even mention it in the log
    foreach my $overlapped_feature ( get_overlapped_features ( $o, $check_feature, "equals", [$check_feature->primary_tag] ) ) {
        return 1;
    }

    if ( $check_feature->primary_tag eq "source" ) {
        # source features are not to be skipped
    }

    elsif ( $check_feature->primary_tag eq "gene" ) {
        # genes features are copied exclusively from other feature types (rRNA, tRNA, tmRNA, CDS),
        # so if those features pass the overlap detection, we should not check it again during storage
    }

    elsif ( $check_feature->primary_tag eq "rRNA" ) {
        # rRNA features are not to be skipped
    }

    elsif ( $check_feature->primary_tag eq "tRNA" ) {
        # tRNA features are not to be skipped
    }

    elsif ( $check_feature->primary_tag eq "tmRNA" ) {
        # tmRNA features are not to be skipped
    }

    elsif ( $check_feature->primary_tag eq "ncRNA" ) {
        # ncRNA features are not to be skipped
    }

    elsif ( $check_feature->primary_tag eq "misc_feature" ) {
        # misc_feature features are not to be skipped
    }

    elsif ( $check_feature->primary_tag eq "misc_binding" ) {
        # misc_binding features are not to be skipped
    }

    elsif ( $check_feature->primary_tag eq "ORF" ) {
        # ORF features are not to be skipped
    }

    # CDS features are most difficult
    elsif ( $check_feature->primary_tag eq "CDS" ) {

        #########################CONTAINS#########################
        # CDS features should not contain any important features, we remove the longer
        #--------------------------------------------------------------------------------------> CDS
        #                               <--------------- tRNA
        foreach my $overlapped_feature ( get_overlapped_features ( $o, $check_feature, "contains", ["tRNA", "tmRNA"] ) ) {
            print_verbose ( $o, "Skipping ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.") due to overlap with ".$overlapped_feature->primary_tag." (".$overlapped_feature->start."..".$overlapped_feature->end.")" );
            return 1;
        }
        # repeat_region inside a new CDS -> we remove the repeat_region, likely mispredicted
        #--------------------------------------------------------------------------------------> CDS
        #        ---------------> repeat_region
        foreach my $overlapped_feature ( get_overlapped_features ( $o, $check_feature, "contains", ["repeat_region"] ) ) {
            print_verbose ( $o, "Removing previously annotated ".$overlapped_feature->primary_tag." (".$overlapped_feature->start."..".$overlapped_feature->end.") fully contained within ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.")" );
            delete_feature ( $o, $overlapped_feature);
        }

        # old CDS is inside a larger new CDS, different starts/ends # we remove different starts/ends first if the score is greater
        #--------------------------------------------------------------------------------------> CDS1
        #        ---------------> CDS2
        # or
        #--------------------------------------------------------------------------------------> CDS1
        #                                                                       ---------------> CDS2
        foreach my $overlapped_feature ( get_overlapped_features ( $o, $check_feature, "contains", ["CDS"] ) ) {
            print_verbose ( $o, "Removing previously annotated ".$overlapped_feature->primary_tag." (".$overlapped_feature->start."..".$overlapped_feature->end.") fully contained within ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.")" );
            delete_feature ( $o, $overlapped_feature);
        }

        #########################OVERLAPS#########################
        # a CDS overlaping a rRNA is likely a mistake, we skip it
        #--------------------------------------------------------------------------------------> CDS
        #                                                                                   <------------------------------- rRNA
        #                                                                            ---------------> repeat_reagion
        foreach my $overlapped_feature ( get_overlapped_features ( $o, $check_feature, "overlaps", ["rRNA"] ) ) {
            print_verbose ( $o, "Skipping ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.") due to overlap with ".$overlapped_feature->primary_tag." (".$overlapped_feature->start."..".$overlapped_feature->end.")" );
            return 1;
        }
        # among ncRNA, for now only features with ncRNA_class=ribozyme and =RNase_P_RNA are not allowed to overlap with CDS
        # a CDS overlaping these types of ncRNA is likely a mistake, we skip it
        #--------------------------------------------------------------------------------------> CDS
        #                                                                                   <------------------------------- ncRNA, nc_RNA_class=ribozyme
        foreach my $overlapped_feature ( get_overlapped_features ( $o, $check_feature, "overlaps", ["ncRNA"] ) ) {
            if ( grep { $_ =~ m/^ribozyme$|^RNase_P_RNA$/ } $overlapped_feature->get_tagset_values("ncRNA_class" ) ) {
                print_verbose ( $o, "Skipping ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.") due to overlap with ".$overlapped_feature->primary_tag." (".$overlapped_feature->start."..".$overlapped_feature->end.")" );
                return 1;
            }
        }

        #########################INVERSE CONTAINMENT#########################
        foreach my $overlapped_feature ( get_overlapped_features ( $o, $check_feature, "overlaps", ["CDS"] ) ) {
            if ( $overlapped_feature->contains( $check_feature, "ignore" ) ) {
                # new CDS is inside a larger old CDS
                #                           <--------------- CDS1
                #--------------------------------------------------------------------------------------> CDS2
                # or
                #                                                                       ---------------> CDS1
                #--------------------------------------------------------------------------------------> CDS2
                print_verbose ( $o, "Removing previously annotated ".$overlapped_feature->primary_tag." (".$overlapped_feature->start."..".$overlapped_feature->end.") fully contained within ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.")" );
                delete_feature ( $o, $overlapped_feature);
            }
        }
    }
    # all good
    return 0;
}





######################################################## DATABASE #########################################################

sub get_nucleotide_sequence_of_feature {
    my ( $o, $feature) = @_;
    # get the sequence of the object (unfortunately it is not implemented in Bio:DB::*)
    my $sequence_object = $o->{"dbh_uuid"}->fetch_sequence( -seq_id=>$feature->seq_id, -start=>$feature->start, -end=>$feature->end, -bioseq=>1 );
    # set the id of the new sequence to the id of the feature (otherwise it has $seq_id:$start..$end format)
    $sequence_object->id( $feature->primary_id );
    # change the strand if needed
    if ( defined $feature->strand and $feature->strand == -1 ) {
        $sequence_object = $sequence_object->revcom;
    }
    return $sequence_object;
}

# wrapper for get_nucleotide_sequence_of_feature
sub get_protein_sequence_of_feature {
    my ( $o, $feature) = @_;
    my $sequence_object = get_nucleotide_sequence_of_feature( $o, $feature )->translate( -complete => 1, -codontable_id => 11 );
    return $sequence_object;
}

##################################################### IO FILE ######################################################
sub prepare_input {
  # prepares the input
    my ( $o, $s, $input_type ) = @_;
    prepare_sequences ( $o, $s ) if ( $input_type eq "sequences" );
    prepare_homology ( $o, $s ) if ( $input_type eq "homology" );
    prepare_annotation ( $o, $s ) if ( $input_type eq "annotation" );
    return;
}

sub prepare_sequences {
  # prepares an input file containing input sequences in FASTA format
  my ( $o ) = @_;
  my $output_filehandle = Bio::SeqIO->new(-file => ">".$o->{"job"}."/input_sequences", -format => "fasta");

  foreach my $seq_id ( sort keys %{$o->{"r"}} ) {
    # do not just dump $o->{"r"}; instead re-create to get rid of the description (breaks some tools, e.g. GeneMarkS)
    $output_filehandle->write_seq( Bio::Seq->new( -seq => $o->{"r"}->{$seq_id}->seq(),
                                                  -id => $seq_id ) );
  }
  exit_program ( $o, "Required file ".$o->{"job"}."/input_sequences does not exist.") if not (-e $o->{"job"}."/input_sequences" );
}

sub prepare_homology {
    my ( $o, $s ) = @_;
    # query file contains blast input if is has NOT been calculated previously
    my $query_filehandle = Bio::SeqIO->new(-file => ">".$o->{"job"}."/input_homology_".$o->{"cds-h"}, -format => "fasta" );
    foreach my $intergenic_and_cds_sequence ( get_intergenic_and_cds_sequences ( $o, $s ) ) {
        $intergenic_and_cds_sequence->id ( $intergenic_and_cds_sequence->id );
        if ( not exists $o->{"homology"}->{$intergenic_and_cds_sequence->id} ) {
            # set the "exists" tag
            $o->{"homology"}->{$intergenic_and_cds_sequence->id} = 1;
            # store the empty for future runs
            $query_filehandle->write_seq( $intergenic_and_cds_sequence );
        }
    }
}

sub prepare_annotation {
    my ( $o, $s ) = @_;
    my $query_filehandle = Bio::SeqIO->new(-file => ">".$o->{"job"}."/input_annotation", -format => "fasta" );
    foreach my $seq_id (sort keys %$s) {
        foreach my $feature ( $o->{"dbh_uuid"}->features( -seq_id => $seq_id, -type => ["CDS"] ) ) {
            my $sequence_object = get_protein_sequence_of_feature( $o, $feature );
            $sequence_object->id ( $sequence_object->id );
            $query_filehandle->write_seq( $sequence_object );
        }
    }
}

sub clean_up_description {
  # returns a description that follows the NCBI guidelines
  my ( $description ) = @_;
  # NCBI comments
  # [4] We prefer that the terms like 'homolog', 'paralog' or 'analog' not be
  # used as part of a protein name as this infers an evolutionary relationship
  # that has generally not been determined. Please use '-like protein' instead.
  # feature contains 'Homolog'
  # CDS     60-kDa SS-A/Ro ribonucleoprotein homolog        lcl|PBC10988:2544384-2546009    PBC10988_15970
  $description =~ s/\shomolog/-like/gi;
  $description =~ s/\sparalog/-like/gi;
  $description =~ s/\sanalog/-like/gi;
  # feature contains '@'
  # CDS     Acetaldehyde dehydrogenase @ Acetaldehyde dehydrogenase, ethanolamine utilization cluster       lcl|PBC10988:c161956-160472
  $description =~ s/\s\@.*//g;
  # [1] Some proteins are described as fragments or truncated.
  # 85 features contain 'Fragment'
  # CDS     MFS domain-containing protein (Fragment)        lcl|PBC10988:c2765-1416 PBC10988_0020
  $description =~ s/\s\(Fragment\)//gi;
  # feature contains 'Bacteroides'
  # /tmp/tmp.CBX0IaLU0p:CDS	Bacteroides aerotolerance operon BatA	lcl|sequence_1:c840585-839548	PBC10988_4360
  $description =~ s/^Bacteroides\s//gi;
  # feature contains 'disulphide', Replace with 'disulfide'
  $description =~ s/disulphide/disulfide/gi;
  # feature contains 'dyhydrogenase', Replace with 'dehydrogenase'
  $description =~ s/dyhydrogenase/dehydrogenase/gi;
  # feature contains 'methlytransferase', Replace with 'methyltransferase'
  $description =~ s/methlytransferase/methyltransferase/gi;
  # feature contains 'predicted', Replace with 'putative'
  $description =~ s/predicted/putative/gi;
  # feature contains 'probable', Replace with 'putative'
  $description =~ s/probable/putative/gi;
  # feature contains 'sulpho', Replace with 'sulfo'
  $description =~ s/sulpho/sulfo/gi;
  # feature contains 'SWIM zinc finger', Replace with 'SWIM zinc finger protein'
  $description =~ s/SWIM zinc finger/SWIM zinc finger protein/gi;
  # feature contains 'transposase and inactivated derivative', Replace with 'transposase'
  $description =~ s/transposase and inactivated derivative/Transposase/gi;
  # feature contains 'dimerisation', Replace with 'dimerization'
  $description =~ s/dimerisation/dimerization/gi;
  # feature contains 'sulphate', Replace with 'sulfate'
  $description =~ s/sulphate/sulfate/gi;
  # feature contains 'sulphide', Replace with 'sulfide'
  $description =~ s/sulphide/sulfide/gi;
  # feature contains 'utilisation', Replace with 'utilization'
  $description =~ s/utilisation/utilization/gi;
  # feature contains 'highly conserved'
  $description =~ s/^Highly conserved //gi;
  # feature starts with 'domain of unknown function', Replace with 'protein of unknown function'
  $description =~ s/^Domain of unknown function (DUF3291)$//;
  # Protein name ends with bracket and may contain organism name FEATURE: Prot: Bifunctional deaminase-reductase domain-containing protein [Kribbella flavida DSM] [gnl|PRJNA578297|GC106_64800:1-162] [gnl|PRJNA578297|GC106_64800: raw, aa len= 162]
  $description =~ s/ \[.+\]$//; # remove square brackets, but only from the end, avoid removing  3-oxoacyl-[acyl-carrier-protein] reductase or [Butirosin acyl-carrier protein]--L-glutamate ligase or 3-hydroxyacyl-[acyl-carrier-protein] dehydratase
  # feature contain 'faimly', Replace with 'family'.
  $description =~ s/faimly/family/gi;
  # feature contains 'doubtful'
  $description =~ s/, doubtful CDS$//gi;
  # feature contains 'doubtful'
  $description =~ s/ # A$//gi;
  # feature contains 'DNA for'
  $description =~ s/^DNA for //gi;
  # feature contains 'paralog'
  $description =~ s/paralog without usual motifs$//gi;

  return $description // "Hypothetical protein";
}

sub clean_up_uniprot {
  # returns a description that follows the NCBI guidelines
  my ( $description ) = @_;
  # WARNING: valid [SEQ_FEAT.ProteinNameEndsInBracket] Protein name ends with bracket and may contain organism name FEATURE: Prot: Nicotinate-nucleotide pyrophosphorylase [carboxylating] [lcl|sequence_1_108:1-297] [lcl|sequence_1_108: raw, aa len= 297]
  $description =~ s/\[/(/g;
  $description =~ s/\]/)/g;

  return $description // "Hypothetical protein";
}

# memoized sub parse_file
{
    my %filehandles;
    my $aragorn_sequence_name;
    my $glimmer_sequence_name;
    sub parse_file {
        my ( $o, $file, $action, $divider, $program ) = @_;
        # open file if it"s not open.
        if (not $filehandles{$file}) {
            open $filehandles{$file}, "<", $file;
        }
        # get-a-record
        if ( $action eq "record" ) {
            # consider double slash (or other) as a record division
            local $/ = $divider;
            while (my $line = readline ( $filehandles{$file})) {
                return $line;
            }
        }
        # get-a-line
        elsif ( $action eq "line" ) {
            while (my $line = readline ( $filehandles{$file})) {
                # chomp line first
                chomp $line;
                # skip empty lines
                next if ( $line =~ m/^$/);
                # skip comments
                next if ( $line =~ m/^#/);
                # skip header in tRNAscan-SE
                next if (( $program eq "trnascanse") && ( $line =~ m/^Sequence\s+tRNA|Name\s+tRNA|\-+\s+\-+/));
                # skip summary line in ARAGORN
                next if (( $program eq "aragorn") && ( $line =~ m/^\d+\s+gene/));
                # skip the line containg sequence name in ARAGORN, but store it for the next cycle
                next if (( $program eq "aragorn") && ( $line =~ m/^>(\S+)/) && ( $aragorn_sequence_name = $1));
                # skip the line containg sequence name in Glimmer, but store it for the next cycle
                next if (( $program eq "glimmer") && ( $line =~ m/^>(\S+)/) && ( $glimmer_sequence_name = $1));
                # skip the header line from PANNZER
                next if (( $program eq "pannzer") && ( $line =~ m/^qpid/));
                # warn on parse error
                print_log( $o, "Parse error in file $file, line $.") if not (my @line = split (/$divider/, $line));
                # make use of the previously stored $aragorn_sequence_name
                push @line, $aragorn_sequence_name if ( ( $aragorn_sequence_name ) && ( $program eq "aragorn" ) );
                # make use of the previously stored $glimmer_sequence_name
                push @line, $glimmer_sequence_name if ( ( $glimmer_sequence_name ) && ( $program eq "glimmer" ) );
                return \@line;
            }
        }
    }
}




sub get_command {
    my ( $o, $program ) = @_;
    if ( $program eq "rnammer" ) {
      return $o->{"cwd"}."/bin/rnammer/rnammer -m lsu,ssu,tsu -S bac -T /tmp ".$o->{"job"}."/input_sequences -gff ".$o->{"job"}."/rnammer";
    }
    elsif ( $program eq "trnascanse" ) {
      return $o->{"cwd"}."/bin/trnascanse/tRNAscan-SE -q -Q -B ".$o->{"job"}."/input_sequences -o ".$o->{"job"}."/trnascanse";
    }
    elsif ( $program eq "aragorn" ) {
      return $o->{"cwd"}."/bin/aragorn/aragorn -w -gc11 ".$o->{"job"}."/input_sequences -o ".$o->{"job"}."/aragorn";
    }
    elsif ( $program eq "infernal" ) {
      return $o->{"cwd"}."/bin/infernal/cmscan --tblout ".$o->{"job"}."/infernal --notextw --cut_tc --cpu 1 ".$o->{"cwd"}."/databases/rfam/rfam_bacteria.cm ".$o->{"job"}."/input_sequences";
    }
    elsif ( $program eq "glimmer" ) {
      return $o->{"cwd"}."/bin/glimmer/long-orfs --trans_table 11 --linear -n -t 1.15 ".$o->{"job"}."/input_sequences ".$o->{"job"}."/glimmer.longorfs; ".
             $o->{"cwd"}."/bin/glimmer/extract --nowrap -t ".$o->{"job"}."/input_sequences ".$o->{"job"}."/glimmer.longorfs | tee ".$o->{"job"}."/glimmer.train; ".
             $o->{"cwd"}."/bin/glimmer/build-icm -r ".$o->{"job"}."/glimmer.icm < ".$o->{"job"}."/glimmer.train; ".
             $o->{"cwd"}."/bin/glimmer/glimmer3 --trans_table 11 --linear --extend ".$o->{"job"}."/input_sequences ".$o->{"job"}."/glimmer.icm ".$o->{"job"}."/glimmer; ".
             "mv ".$o->{"job"}."/glimmer.predict ".$o->{"job"}."/glimmer";
    }
    elsif ( $program eq "prodigal" ) {
        # Error: Sequence must be 20000 characters (only 16084 read).
        #(Consider running with the '-p anon' option or finding more contigs from the same genome.)
        # on the other hand:
        # Error: Can't specify translation table with anon mode or training file.
        # therefore, either -p anon, or -g 11
        return $o->{"cwd"}."/bin/prodigal/prodigal -e 0 -q -f gff".( $o->{"input_size"} < 20000 ? " -p anon" : " -g 11 -i " ).$o->{"job"}."/input_sequences -o ".$o->{"job"}."/prodigal";
    }
    elsif ( $program eq "genemarks" ) {
        return $o->{"cwd"}."/bin/genemarks/gmsn.pl --format GFF --prok ".$o->{"job"}."/input_sequences --output ".$o->{"job"}."/genemarks";
    }
    elsif ( $program eq "blast" ) {
        return $o->{"cwd"}."/bin/blast/blastx -query_gencode 11 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen' -max_target_seqs 1 -db ".$o->{"cwd"}."/databases/uniprot/sprot_bacteria.fasta -query ".$o->{"job"}."/input_homology_".$o->{"cds-h"}." | tee -a ".$o->{"job"}."/output_homology_".$o->{"cds-h"};
    }
    elsif ( $program eq "diamond" ) {
        return $o->{"cwd"}."/bin/diamond/diamond blastx --query-gencode 11 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen --max-target-seqs 1 --query ".$o->{"job"}."/input_homology_".$o->{"cds-h"}." --db ".$o->{"cwd"}."/databases/uniprot/sprot_bacteria.fasta --out ".$o->{"job"}."/output_homology_".$o->{"cds-h"};
    }
    elsif ( $program eq "pannzer" ) {
        return "python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -H localhost -T localhost -d ".$o->{"cwd"}."/databases/pannzer -i ".$o->{"job"}."/input_annotation -o ,".$o->{"job"}."/output.desc,".$o->{"job"}."/output.go,".$o->{"job"}."/output.anno";
    }
    elsif ( $program eq "emapper" ) {
        return $o->{"cwd"}."/bin/emapper/emapper.py --override --cpu 1 --data_dir ".$o->{"cwd"}."/databases/emapper -m diamond -i ".$o->{"job"}."/input_annotation -o ".$o->{"job"}."/output"
    }
}

1;
