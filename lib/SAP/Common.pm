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
                get_line
                create_feature
                check_and_store_feature
                add_leading_zeros
                clean_up_description
            );

######################################################## PARAMETERS PART ########################################################
sub initialize_options {
    #setting default values
    my $o = {
        #modules
            #submodules
                #parameters
        #DNA feature prediction
        "dna" => 1,
            #CRISPR
            "dna-c" => 1,
                "dna-c-program" => "minced",
                "dna-c-score" => -inf,
            #tandem repeats
            "dna-t" => 0,
                "dna-t-program" => "trf",
                "dna-t-score" => -inf,
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
                "cds-i-score" => -inf,
            #homology
            "cds-h" => 1,
                "cds-h-program" => "diamond",
                "cds-h-score" => -inf,
            #annotation
            "cds-a" => 1,
                "cds-a-program" => "pannzer",
                "cds-a-score" => -inf,
        #additional parameters
        "program-name" => "ProSnap: Prokaryotic Sequence Annotation Pipeline",
        "program-version" => "20200518",
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
        "gcode" => "11",
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
        #DNA feature prediction
        "dna!"                                              => \$o->{"dna"},
            #CRISPR
            "dna-c!"                                        => \$o->{"dna-c"},
                "dna-c-program=s"                           => \$o->{"dna-c-program"},
                "dna-c-score=s"                             => \$o->{"dna-c-score"},
            #tandem repeats
            "dna-t!"                                        => \$o->{"dna-t"},
                "dna-t-program=s"                           => \$o->{"dna-t-program"},
                "dna-t-score=s"                             => \$o->{"dna-t-score"},
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
            #annotation
            "cds-a!"                                        => \$o->{"cds-a"},
                "cds-a-program=s"                           => \$o->{"cds-a-program"},
                "cds-a-score=s"                             => \$o->{"cds-a-score"},
        #meta data
        "taxid=i"                                           => \$o->{"taxid"},
        "gcode=i"                                           => \$o->{"gcode"},
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
    return $o;
}

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
        -gcode                  integer         -                                       NCBI genetic code table (e.g. 11)
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
        -cds-h-program          string          blast                                   blast, diamond
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
    # forward both stdout and stderr to respective log
    system ( $command . " 2>>". $o->{"job"}."/stdio 1>&2" );
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

    print_verbose ( $o, "Storing feature: ".$feature->primary_tag.":".$feature->start."-".$feature->end );

    # store feature on the record
    $o->{"r"}->{$feature->seq_id()}->add_SeqFeature( $feature );
    # mark unique feature location
    $o->{"feature_by_loc"}->{ $feature->seq_id().$feature->primary_tag.$feature->location()->to_FTstring() } = $feature;
    # store by feature unique id
    $o->{"feature_by_id"}->{ $feature->annotation()->{"unique_id"} } = $feature;

    return;
}

sub feature_rules {
  my ( $o, $feature ) = @_;
  if ( $feature->primary_tag =~ m/^CDS/ ) {

    # get the translation with complete = false
    my $translation = get_protein_sequence_of_feature( $o, $feature, 0 );
    # checking for proper stop codon
    if ( $translation !~ m/\*$/ ) {
      print_verbose ( $o, "Skipping CDS feature with no stop codon at the end: ".$feature->primary_tag.":".$feature->start."-".$feature->end );
      return 1;
    }
    # refresh the translation with complete = true
    $translation = get_protein_sequence_of_feature( $o, $feature, 1 );
    # checking for proper start codon
    if ( $translation !~ m/^M/ ) {
      print_verbose ( $o, "Skipping CDS feature with no start codon at the beginning: ".$feature->primary_tag.":".$feature->start."-".$feature->end );
      return 1;
    }

    # checking for internal stop codon
    if ( $translation =~ m/\*/ ) {
      print_verbose ( $o, "Skipping CDS feature with internal stop codons: ".$feature->primary_tag.":".$feature->start."-".$feature->end );
      return 1;
    }

  }
  # all good
  return 0;
}

sub get_adjoined_features {
  # returns all features that equal/contain/overlap with certain feature
  # accepts: $feature - a feature to check for equality/containment/overlap
  #          $rule - rule, according to which the check is being performed (equals, contains, overlaps)
  #          @types - a list of one or more feature types to include in the check
  # returns: a list of features of type @types that have a $rule relation (equal, contain, overlap) to the feature in question
  my ( $o, $feature, $rule, @types ) = @_;
  my @candidates;
  foreach my $type ( @types ) {
    foreach my $candidate ( $o->{"r"}->{ $feature->seq_id() }->get_SeqFeatures($type) ) {
      if ( $feature->$rule( $candidate ) ) {
        push @candidates, $candidate;
      }
    }
  }
  return @candidates;
}




sub create_seq_hash {
    my ( $o ) = @_;
    # to decide on sequence numbering below
    my $tmp_filehandle = Bio::SeqIO->new(-file => $o->{"i"}, -format => $o->{"input_format"} );
    while ( my $record = $tmp_filehandle->next_seq ) {
        $o->{"input_count"}++;
        $o->{"input_size"} = ( $o->{"input_size"} ? $o->{"input_size"} : 0 ) + $record->length;
    }
    $o->{"sequence_digit"} = length( $o->{"input_count"} );
    print_log( $o, "Input has ".$o->{"input_count"}." sequences, setting sequence numbering to ".$o->{"sequence_digit"}."-digit" );
    # the real part
    my $input_filehandle = Bio::SeqIO->new(-file => $o->{"i"}, -format => $o->{"input_format"} );

    while ( my $record = $input_filehandle->next_seq ) {

      # we have to change the name in the very beginning
      $record->id(add_leading_zeros(( keys %{$o->{"r"}} ) + 1, $o->{"sequence_digit"}));

      # put the input record for further reference
      $o->{"r"}->{$record->id} = $record;

      # set the sequence to circular if specified
      if ( $o->{"circular"}->{$record->id} ) {
        $record->is_circular( $o->{"circular"}->{$record->id} );
      }

      # transparent
      if ( $o->{"transparent"} ) {
        # store all features in case of transparent mode
        @{ $o->{"r"}->{$record->id}->{"transparent"} } = $record->get_all_SeqFeatures();
      }

      # remove all features, transparent or not
      $o->{"r"}->{$record->id}->remove_SeqFeatures();

    }

    $o->{"locus_digit"} = int (log( 100000 )/log(10) - 2);
    print_log( $o, "Input is ".$o->{"input_size"}." bp long, setting locus numbering to ".$o->{"locus_digit"}."-digit" );
}

sub overlap_rules {
    my ( $o, $check_feature ) = @_;

    # catch-all rule for unique location/primary tag combination
    if ( exists $o->{"feature_by_loc"}->{ $check_feature->seq_id().$check_feature->primary_tag.$check_feature->location()->to_FTstring() } ) {
      print_verbose ( $o, "Skipping exact location duplicate feature ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.")" );
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
      if ( get_adjoined_features ( $o, $check_feature, "overlaps", "tmRNA" ) ) {
        # a tRNA gene within the tmRNA gene is the same gene -> skip it
        print_verbose ( $o, "Skipping same gene annotated at the same place ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.")" );
        return 1;
      }
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

    elsif ( $check_feature->primary_tag eq "CDS" ) {
      # CDS features are most complex

      # CDS-CDS overlaps
      foreach my $feature ( get_adjoined_features ( $o, $check_feature, "overlaps", "CDS" ) ) {
        # derived from homology
        if ( $check_feature->annotation()->{"method"} eq "homology" ) {
          # doesn't matter whom contans whom -> if this is the same gene -> do not add a second one at the same place
          if ( ( $feature->contains( $check_feature ) or $check_feature->contains( $feature ) ) and $check_feature->annotation()->{"gene_id"} eq $feature->annotation()->{"gene_id"} ) {
            print_verbose ( $o, "Skipping same gene annotated at the same place ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.")" );
            return 1;
          }
        }
        # derived from ab initio
        if ( $check_feature->annotation()->{"method"} eq "ab initio" ) {
          # doesn't matter whom contans whom -> do not add a second one at the same place
          if ( $feature->contains( $check_feature ) or $check_feature->contains( $feature ) ) {
            print_verbose ( $o, "Skipping a CDS annotated at the same place ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.")" );
            return 1;
          }
        }
      }

      # CDS-RNA overlaps
      foreach my $feature ( get_adjoined_features ( $o, $check_feature, "overlaps", "rRNA", "tRNA", "tmRNA" ) ) {
        print_verbose ( $o, "Skipping CDS feature overlapping an RNA gene: ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.")" );
        return 1;
      }

      # CDS-repeat overlaps
      foreach my $feature ( get_adjoined_features ( $o, $check_feature, "overlaps", "repeat_region" ) ) {
        # if overlaps CRISPR, do not annotate
        if ( $feature->annotation()->{"type"} eq "CRISPR" ) {
          print_verbose ( $o, "Skipping CDS feature overlapping a CRISPR array: ".$check_feature->primary_tag." (".$check_feature->start."..".$check_feature->end.")" );
          return 1;
        }
      }

    }

    # all good
    return 0;
}





######################################################## DATABASE #########################################################

sub get_nucleotide_sequence_of_feature {
  # get the sequence of a feature
  my ( $o, $feature ) = @_;
  return $feature->spliced_seq()->seq();
}

sub get_protein_sequence_of_feature {
  # get the sequence of a feature
  my ( $o, $feature, $complete ) = @_;
  # lack of -complete passed to ->translate allows to expose terminal stop codons as asterisks (*)
  # -complete translation converts start codons into M
  my $sequence = $feature->spliced_seq()->translate( -offset => ( $feature->get_tagset_values( "codon_start" ) )[0] || 1,
                                                     -codontable_id => $o->{"gcode"},
                                                     -complete => $complete )->seq();

  # if this hit has a /transl_except
  if ( my $exception = $feature->annotation()->{"exception"} ) {
    if ( $exception->[0] eq $feature->annotation()->{"hit_id"} ) {
      if ( $exception->[3] eq "aa:Sec" ) {
        # modify to unusual codon tables
        substr( $sequence,$exception->[1]-1,1 ) = "U";
      }
      elsif ( $exception->[3] eq "aa:Pyl" ) {
        # modify to unusual codon tables
        substr( $sequence,$exception->[1]-1,1 ) = "O";
      }
    }
  }

  return $sequence;
}

##################################################### IO FILE ######################################################
sub prepare_input {
  # prepares the input
    my ( $o, $input_type ) = @_;
    prepare_sequences ( $o ) if ( $input_type eq "sequences" );
    prepare_annotation ( $o ) if ( $input_type eq "annotation" );
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

sub prepare_annotation {
  # prepares an input file containing protein sequences in FASTA format
  my ( $o, $s ) = @_;
  my $output_filehandle = Bio::SeqIO->new(-file => ">".$o->{"job"}."/input_annotation", -format => "fasta" );
  foreach my $seq_id ( sort keys %{$o->{"r"}} ) {
    foreach my $feature ( grep { ! $_->has_tag("product") } $o->{"r"}->{$seq_id}->get_SeqFeatures( "CDS" ) ) {
      my $sequence = get_protein_sequence_of_feature( $o, $feature, 1 );
      $output_filehandle->write_seq( Bio::Seq->new( -seq => $sequence,
                                                    -id => $feature->annotation()->{"unique_id"} ) );
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

sub parse_file {
  my ( $o, $file, $divider, $program ) = @_;
  # open file if it"s not open.
  open my $filehandle, "<", $file;
  my @lines;
  while ( my $line = readline ( $filehandle) ) {
    # chomp line first
    chomp $line;
    # skip empty lines
    next if ( $line =~ m/^$/);
    # skip comments
    next if ( $line =~ m/^#/);
    # skip header in tRNAscan-SE
    next if $program eq "trnascanse" and $line =~ m/^Sequence|^Name|^-/;
    # skip sequence lines in ARAGORN
    next if (( $program eq "aragorn") && ( $line !~ m/^>/));
    # skip the header line from PANNZER
    next if (( $program eq "pannzer") && ( $line =~ m/^qpid/));
    # warn on parse error
    print_log( $o, "Parse error in file $file, line $.") if not (my @line = split (/$divider/, $line));
    # push the line to the lines
    push @lines, \@line;
  }
  return map { parse_line( $o, $program, $_ ) } @lines;
}

sub parse_line {
  # returns line elements that are required depending on a type of a program
  my ( $o, $program, $l ) = @_;
  if ( $program eq "infernal" ) {
    #idx target name            accession query name           accession clan name mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc olp anyidx afrct1 afrct2 winidx wfrct1 wfrct2 description of target
    #--- ---------------------- --------- -------------------- --------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- --- ------ ------ ------ ------ ------ ------ ---------------------
    #[0] [1]                    [2]       [3]                  [4]       [5]       [6] [7]      [8]      [9]      [10]     [11]   [12]  [13] [14] [15]  [16]   [17]      [18][19][20]   [21]   [22]   [23]   [24]   [25]   [26]
    # get the RFAM information line
    my $i = get_line( $o, $o->{"cwd"}."/databases/rfam/family.txt", $l->[2], "\t" );
    # column names from http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.sql
    #rfam_acc rfam_id auto_wiki description author seed_source gathering_cutoff trusted_cutoff noise_cutoff comment previous_id cmbuild cmcalibrate cmsearch num_seed num_full num_genome_seq num_refseq type structure_source number_of_species number_3d_structures num_pseudonokts tax_seed ecmli_lambda ecmli_mu ecmli_cal_db ecmli_cal_hits maxl clen match_pair_node hmm_tau hmm_lambda created updated
    #[0]      [1]     [2]       [3]         [4]    [5]         [6]              [7]            [8]          [9]     [10]        [11]    [12]        [13]     [14]     [15]     [16]           [17]       [18] [19]             [20]              [21]                 [22]            [23]     [24]         [25]     [26]         [27]           [28] [29] [30]            [31]    [32]       [33]    [34]
    return { "seq_id" => $l->[3],
             "accession" => $l->[2],
             "hit_id" => $l->[1],
             "start" => $l->[11] eq "+" ? $l->[9] : $l->[10],
             "start_type" => "EXACT",
             # "start_type" => ( $l->[12] =~ /5'/ && $l->[11] eq "+" ) || ( $l->[12] =~ /3'/ && $l->[11] eq "-" ) ? "BEFORE" : "EXACT",
             # commented due to tbl2asn warnings
             # [SEQ_FEAT.PartialProblem] PartialLocation: Start does not include first/last residue of sequence FEATURE: misc_feature: /inference=profile:infernal:1.1.3:rfam:RF01766; cspA thermoregulator [lcl|sequence_1:c>2150785-<2150428] [lcl|sequence_1: raw, dna len= 5697240]
             "end" => $l->[11] eq "+" ? $l->[10] : $l->[9],
             "end_type" => "EXACT",
             # "end_type" => ( $l->[12] =~ /3'/ && $l->[11] eq "+" ) || ( $l->[12] =~ /5'/ && $l->[11] eq "-" ) ? "AFTER" : "EXACT",
             # commented due to tbl2asn warnings
             # [SEQ_FEAT.PartialProblem] PartialLocation: Start does not include first/last residue of sequence FEATURE: misc_feature: /inference=profile:infernal:1.1.3:rfam:RF01766; cspA thermoregulator [lcl|sequence_1:c>2150785-<2150428] [lcl|sequence_1: raw, dna len= 5697240]
             "strand" => $l->[11] eq "+" ? 1 : -1,
             "score" => $l->[16],
             "product" => ( join " ", @$l[26..$#$l] ),
             "type" => $i->[18],
             "clan" => $l->[5],
             "comment" => $i->[9],
           }
  }
  if ( $program eq "rnammer" ) {
    # seqname           source                      feature     start      end   score   +/-  frame  attribute
    # ---------------------------------------------------------------------------------------------------------
    # [0]               [1]                         [2]         [3]        [4]   [5]     [6]  [7]    [8]
    return { "seq_id" => $l->[0],
             "start" => $l->[3],
             "start_type" => "EXACT", # rnammer can only produce exact coordinates
             "end" => $l->[4],
             "end_type" => "EXACT", # rnammer can only produce exact coordinates
             "strand" => $l->[6] eq "+" ? 1 : -1,
             "score" => $l->[5],
             "type" => "rRNA",
             "product" => $l->[8] =~ s/s_rRNA/S subunit ribosomal rRNA/ri, # make the product name nicer by using the /r option
           }
  }
  if ( $program eq "trnascanse" ) {
    #Sequence                tRNA    Bounds  tRNA    Anti    Intron Bounds   Inf
    #Name            tRNA #  Begin   End     Type    Codon   Begin   End     Score   Note
    #--------        ------  -----   ------  ----    -----   -----   ----    ------  ------
    #1               1       229152  229228  Ile     GAT     0       0       75      pseudo
    #[0]             [1]     [2]     [3]     [4]     [5]     [6]     [7]     [8]     [9]
    return { "seq_id" => $l->[0] =~ s/\s+//ri, # make use of the /r option, sequence name has trailing spaces
             "start" => $l->[2] > $l->[3] ? $l->[3] : $l->[2],
             "start_type" => "EXACT", # trnascanse can only produce exact coordinates
             "end" => $l->[2] > $l->[3] ? $l->[2] : $l->[3],
             "end_type" => "EXACT", # trnascanse can only produce exact coordinates
             "strand" => $l->[2] > $l->[3] ? -1 : 1,
             "score" => $l->[8],
             "type" => "tRNA",
             "product" => $l->[4] eq "Undet" ? "tRNA-Xxx" : "tRNA-".$l->[4]."(".$l->[5].")", # SEQ_FEAT.MissingTrnaAA
                                                                                             # Explanation: The amino acid that the tRNA carries is not included.
                                                                                             # Suggestion: Include the amino acid as the product of the tRNA. If the amino acid of a tRNA is unknown, use tRNA-Xxx as the product.
                                                                                             # From: https://www.ncbi.nlm.nih.gov/genbank/genome_validation/
           }
  }
  if ( $program eq "aragorn" ) {
    #>1-60 tmRNA c[3604541,3604903]
    #>1-98 tRNA-Phe(gaa) c[5387551,5387626]
    #[0]   [1]           [2]
    return { "seq_id" => ( $l->[0] =~ m/^>(\d+)-\d+$/ )[0],
             "start" => ( $l->[2] =~ m/^c?\[(\d+),\d+\]$/ )[0],
             "start_type" => "EXACT", # aragorn can only produce exact coordinates
             "end" => ( $l->[2] =~ m/^c?\[\d+,(\d+)\]$/ )[0],
             "end_type" => "EXACT", # aragorn can only produce exact coordinates
             "strand" => $l->[2] =~ m/^c/ ? -1 : 1,
             "score" => 0,
             "product" => $l->[1],
           }
  }

  if ( $program eq "prodigal" ) {
    ##gff-version  3
    # Sequence Data: seqnum=1;seqlen=84908;seqhdr="1"
    # Model Data: version=Prodigal.v2.6.3;run_type=Metagenomic;model="20|Escherichia_coli_UMN026|B|50.7|11|1";gc_cont=50.70;transl_table=11;uses_sd=1
    #1       Prodigal_v2.6.3 CDS     3       704     76.7    +       0       ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.546;conf=100.00;score=76.68;cscore=73.46;sscore=3.22;rscore=0.00;uscore=0.00;tscore=3.22;
    #[0]     [1]             [2]     [3]     [4]     [5]     [6]     [7]     [8]

    return { "seq_id" => $l->[0],
             "start" => $l->[3],
             #"start_type" => $l->[8] =~ /partial=1\d/ ? "BEFORE" : "EXACT",
             # commented due to tbl2asn warnings
             # [SEQ_FEAT.PartialProblem] PartialLocation: Start does not include first/last residue of sequence FEATURE: misc_feature: /inference=profile:infernal:1.1.3:rfam:RF01766; cspA thermoregulator [lcl|sequence_1:c>2150785-<2150428] [lcl|sequence_1: raw, dna len= 5697240]
             "start_type" => "EXACT",
             "end" => $l->[4],
             #"end_type" => $l->[8] =~ /partial=\d1/ ? "AFTER" : "EXACT",
             # commented due to tbl2asn warnings
             # [SEQ_FEAT.PartialProblem] PartialLocation: Start does not include first/last residue of sequence FEATURE: misc_feature: /inference=profile:infernal:1.1.3:rfam:RF01766; cspA thermoregulator [lcl|sequence_1:c>2150785-<2150428] [lcl|sequence_1: raw, dna len= 5697240]
             "end_type" => "EXACT",
             "strand" => $l->[6] eq "+" ? 1 : -1,
             "score" => $l->[5],
             "method" => "ab initio",
             "type" => "CDS",
             "tags" => {
                          "transl_table" => $o->{"g_code"},
                          "codon_start" => 1,
              }
           }
  }

  if ( $program eq "genemarks" ) {
    ##gff-version 2
    ##source-version GeneMark.hmm_PROKARYOTIC 3.36
    ##date: Tue Feb 25 20:19:42 2020
    # Sequence file name: /annotate/public/jobs/67E87C78-5803-11EA-ABE0-A900B28B97A2/input_sequences
    # Model file name: GeneMark_hmm_combined.mod
    # RBS: true
    # Model information: GeneMarkS_gcode_11
    #1       GeneMark.hmm    CDS     3       98      1.290982        +       0       gene_id=1
    #[0]     [1]             [2]     [3]     [4]     [5]             [6]     [7]     [8]
    return { "seq_id" => $l->[0],
             "start" => $l->[3],
             "start_type" => "EXACT", # genemarks can only produce exact coordinates
             "end" => $l->[4],
             "end_type" => "EXACT", # genemarks can only produce exact coordinates
             "strand" => $l->[6] eq "+" ? 1 : -1,
             "score" => $l->[5],
             "method" => "ab initio",
             "type" => "CDS",
             "tags" => {
                          "transl_table" => $o->{"g_code"},
                          "codon_start" => 1,
              }
           }
  }

  if ( $program eq "glimmer" ) {
    #1      orf00001      343     2799  +1    10.26
    #[0]    [1]           [2]     [3]   [4]   [5]
    return { "seq_id" => $l->[0],
             "start" => $l->[2] > $l->[3] ? $l->[3] : $l->[2],
             "start_type" => "EXACT", # glimmer can only produce exact coordinates
             "end" => $l->[2] > $l->[3] ? $l->[2] : $l->[3],
             "end_type" => "EXACT", # glimmer can only produce exact coordinates
             "strand" => $l->[2] > $l->[3] ? -1 : 1,
             "score" => $l->[5],
             "method" => "ab initio",
             "type" => "CDS",
             "tags" => {
                          "transl_table" => $o->{"g_code"},
                          "codon_start" => 1,
              }
           }
  }

  if ( $program eq "diamond" or $program eq "blast" ) {
    # get the EXCEPTION information line
    #accession transl_except
    #--------- ----------------
    #[0]       [1]
    my $i = get_line( $o, $o->{"cwd"}."/databases/uniprot/exceptions.txt", $l->[1], "\t" );
    # print Dumper($l->[1]);
    #qseqid         sseqid                  qstart      qend        sstart   send       slen      pident    bitscore     stitle
    #1              sp|A7ZUK1|RPOB_ECO24    5146755     5150780     1        1342       1342      100.0     2638.2       sp|A7ZUK1|RPOB_ECO24 DNA-directed RNA polymerase subunit beta OS=Escherichia coli O139:H28 (strain E24377A / ETEC) OX=331111 GN=rpoB PE=1 SV=1
    #[0]            [1]                     [2]         [3]         [4]      [5]        [6]       [7]       [8]          [9]
    # always extend one less aa than $l->[4] to account for 1-based start location in bio
    my $extension_5 = ( $l->[4] - 1 ) * 3;
    # always extend one more aa than $l->[8] - $l->[7] to account for and include stop codon
    my $extension_3 = ( $l->[6] - $l->[5] + 1 ) * 3;
    # start depends on query orientation
    my $start = $l->[2] > $l->[3] ? $l->[3] : $l->[2];
    # so does start extension
    my $extension_start = $l->[2] > $l->[3] ? $extension_3 : $extension_5;
    # so does end
    my $end = $l->[2] > $l->[3] ? $l->[2] : $l->[3];
    # and its extension
    my $extension_end = $l->[2] > $l->[3] ? $extension_5 : $extension_3;
    # never extend start beyond sequence
    $start = ( $start - $extension_start ) < 1 ? $start : ( $start - $extension_start );
    # never extend end beyond sequence
    $end = ( $end + $extension_end ) > $o->{"r"}->{$l->[0]}->length ? $end : ( $end + $extension_end );
    return {
             "seq_id" => $l->[0],
             "hit_id" => $l->[1],
             "gene_id" => ( $l->[1] =~ m/^\S+\|(\S+)_\S+$/ )[0],
             "start" => $start,
             "start_type" => "EXACT", # diamond can only produce exact coordinates
             "end" => $end,
             "end_type" => "EXACT", # diamond can only produce exact coordinates
             "strand" => $l->[2] > $l->[3] ? -1 : 1,
             "score" => $l->[8],
             "product" => clean_up_uniprot( ( $l->[9] =~ m/^\S+ (.*) OS=.*$/ )[0] ),
             "gene_name" => ( $l->[9] =~ m/ GN=(\S+)/ )[0] // undef,
             "exception" => $i,
             "type" => "CDS",
             "method" => "homology",
             "tags" => {
                          "transl_table" => $o->{"g_code"},
                          "codon_start" => 1,
              }
           }
  }

  if ( $program eq "pannzer" ) {
    #qpid    cluster_GSZ             cluster_RM1sum  cluster_size    cluster_desccount     RM2              val_avg         jac_avg                 desc                                                            genename
    #337     1188.6501920470491      88.1362508589   100             8585                  1.16895577275    0.19527405526   2.23321196985e-05       Bifunctional aspartate kinase/homoserine dehydrogenase I        ThrA
    #[0]     [1]                     [2]             [3]             [4]                   [5]              [6]             [7]                     [8]                                                             [9]
    return { "seq_id" => $l->[0],
             "product" => $l->[8],
             "gene_name" => $l->[9],
             "score" => $l->[1],
             "type" => "annotation",
           }
  }

  if ( $program eq "emapper" ) {
    #query_name     seed_eggNOG_ortholog    seed_ortholog_evalue    seed_ortholog_score     predicted_gene_name     GO_terms                  KEGG_KOs          BiGG_reactions       Annotation_tax_scope         OGs                bestOG|evalue|score     COG_cat          eggNOG_annot
    #9397           749414.SBI_06869        1.6e-59                 234.6                   RMLC                    GO:0000271,GO:0003674     K01790,K13316     TDPDRE               NOG[107]                     COG1898@NOG        NA|NA|NA                M, H             DTDP-4-dehydrorhamnose 3,5-epimerase
    #[0]            [1]                     [2]                     [3]                     [4]                     [5]                       [6]               [7]                  [8]                          [9]                [10]                    [11]             [12]
    return { "seq_id" => $l->[0],
             "product" => $l->[12],
             "score" => $l->[3],
             "type" => "annotation",
           }
  }

  if ( $program eq "trf" ) {
    #start    end     period     num_copies     consenus_size       adjacent_identity_overall     adjacent_indels_overall       score       A%      C%      G%      T%      entropy           sequences
    #777012   777170  78         2.0            79                  88                            1                             239         41      20      31      6       1.77              AGAAGAAAGCGGCTGCTGAAAAGGCAGCAGCTGATA...
    #[0]      [1]     [2]        [3]            [4]                 [5]                           [6]                           [7]         [8]     [9]     [10]    [11]    [12]              [13+]
    return { "seq_id" => $l->[0],
              "start" => $l->[0],
              "start_type" => "EXACT", # trf can only produce exact coordinates
              "end" => $l->[1],
              "end_type" => "EXACT", # trf can only produce exact coordinates
              "strand" => 0,
              "type" => "tandem",
              "score" => $l->[7],
           }
  }

  if ( $program eq "pilercr" ) {
    #Pos            Repeat      %id           Spacer     Left flank    Repeat                           Spacer
    #3760999        29          100.0         32         ATTACCTGAT    .............................    GGTAGTACGCGCCTCCGGACGTTTTTATGTCG
    #[0]            [1]         [2]           [3]        [4]           [5]                              [6]
    #Array          Sequence    Position      Length  # Copies  Repeat  Spacer   +  Consensus
    #    1                 1     2877884         580        10      29      32   +  CGGTTTATCCCCGCTGGCGCGGGGAACTC
    #    [0]               [1]   [2]             [3]        [4]     [5]     [6] [7] [8]
    return {  "seq_id" => $l->[7],
              "primary_tag" => $l->[2] eq "repeat_region" ? "repeat_region" : "misc_feature",
              "start" => $l->[0],
              "start_type" => "EXACT", # trf can only produce exact coordinates
              "end" => $l->[1],
              "end_type" => "EXACT", # trf can only produce exact coordinates
              "strand" => 0,
              "type" => "CRISPR",
              "score" => $l->[7],
           }
  }

  if ( $program eq "minced" ) {
    #1       minced:0.4.2    repeat_region   3760816 3761332 9       .       .       ID=CRISPR1;rpt_type=direct;rpt_family=CRISPR;rpt_unit_seq=CGGTTTATCCCCGCTGGCGCGGGGAACAC
    #1       minced:0.4.2    repeat_unit     3760816 3760844 1       .       .       Parent=CRISPR1;ID=DR.CRISPR1.1
    #[0]     [1]             [2]             [3]     [4]     [5]     [6]     [7]     [8]
    return { "seq_id" => $l->[0],
             "start" => $l->[3],
             "primary_tag" => $l->[2] eq "repeat_region" ? "repeat_region" : "misc_feature",
             "start_type" => "EXACT", # minced can only produce exact coordinates
             "end" => $l->[4],
             "end_type" => "EXACT", # minced can only produce exact coordinates
             "strand" => 0,
             "score" => 0,
             "type" => "CRISPR",
             "tags" => $l->[2] eq "repeat_region" ?
                       {
                          "rpt_type" => "direct",
                          "rpt_family" => "CRISPR",
                          "rpt_unit_seq" => ( $l->[8] =~ m/rpt_unit_seq=(.+)/ )[0],
                        } : { "note" => "repeat unit" },
           }
  }

}

sub get_line {
  my ( $o, $file, $id, $divider ) = @_;
  # open file if it"s not open.
  open my $filehandle, "<", $file;
  while ( my $line = readline ( $filehandle) ) {
    # chomp line first
    chomp $line;
    my @line = split (/$divider/, $line);
    if ( $line[0] eq $id ) {
      return \@line;
    }
        }
    }




sub get_command {
    my ( $o, $program ) = @_;
    if ( $program eq "rnammer" ) {
      return $o->{"cwd"}."/bin/rnammer/rnammer -m lsu,ssu,tsu -S bac -T /tmp ".$o->{"job"}."/input_sequences -gff ".$o->{"job"}."/rnammer";
    }
    elsif ( $program eq "trnascanse" ) {
      return $o->{"cwd"}."/bin/trnascanse/tRNAscan-SE -q -B ".$o->{"job"}."/input_sequences -o ".$o->{"job"}."/trnascanse";
    }
    elsif ( $program eq "aragorn" ) {
      return $o->{"cwd"}."/bin/aragorn/aragorn -fon -gc".$o->{"gcode"}." ".$o->{"job"}."/input_sequences -o ".$o->{"job"}."/aragorn";
    }
    elsif ( $program eq "infernal" ) {
      return $o->{"cwd"}."/bin/infernal/cmscan --tblout ".$o->{"job"}."/infernal --fmt 2 --anytrunc --oskip --oclan --clanin ".$o->{"cwd"}."/databases/rfam/Rfam.clanin --notextw --cut_ga ".$o->{"cwd"}."/databases/rfam/Rfam.cm ".$o->{"job"}."/input_sequences";
    }
    elsif ( $program eq "glimmer1" ) {
      return $o->{"cwd"}."/bin/glimmer/long-orfs --trans_table ".$o->{"gcode"}." --linear -n -t 1.15 ".$o->{"job"}."/input_sequences ".$o->{"job"}."/glimmer.longorfs";
    }
    elsif ( $program eq "glimmer2" ) {
      return $o->{"cwd"}."/bin/glimmer/extract --nowrap -t ".$o->{"job"}."/input_sequences ".$o->{"job"}."/glimmer.longorfs | ".$o->{"cwd"}."/bin/glimmer/build-icm -r ".$o->{"job"}."/glimmer.icm";
    }
    elsif ( $program eq "glimmer3" ) {
      return $o->{"cwd"}."/bin/glimmer/glimmer3 --trans_table ".$o->{"gcode"}." --linear ".$o->{"job"}."/input_sequences ".$o->{"job"}."/glimmer.icm ".$o->{"job"}."/glimmer"
    }
    elsif ( $program eq "prodigal" ) {
        # Error: Sequence must be 20000 characters (only 16084 read).
        # (Consider running with the '-p anon' option or finding more contigs from the same genome.)
        # Warning:  ideally Prodigal should be given at least 100000 bases for training.
        # You may get better results with the -p meta option.
        return $o->{"cwd"}."/bin/prodigal/prodigal -q -m -f gff -g ".$o->{"gcode"}."".( $o->{"input_size"} < 100000 ? " -p meta" : "" )." -i ".$o->{"job"}."/input_sequences -o ".$o->{"job"}."/prodigal";
    }
    elsif ( $program eq "genemarks" ) {
        return $o->{"cwd"}."/bin/genemarks/gmsn.pl --format GFF --prok ".$o->{"job"}."/input_sequences --output ".$o->{"job"}."/genemarks";
    }
    elsif ( $program eq "blast" ) {
        return $o->{"cwd"}."/bin/blast/blastx -query_gencode ".$o->{"gcode"}." -outfmt '6 qseqid sseqid qstart qend sstart send slen pident bitscore stitle' -max_target_seqs 1000000 -db ".$o->{"cwd"}."/databases/uniprot/uniprot_sprot.fasta -query ".$o->{"job"}."/input_sequences -out ".$o->{"job"}."/blast";
    }
    elsif ( $program eq "diamond" ) {
        return $o->{"cwd"}."/bin/diamond/diamond blastx --masking 0 --min-score 30 --id 95 --range-culling -F 15 --subject-cover 95 --query-gencode ".$o->{"gcode"}." --outfmt 6 qseqid sseqid qstart qend sstart send slen pident bitscore stitle --max-target-seqs 0 --query ".$o->{"job"}."/input_sequences --db ".$o->{"cwd"}."/databases/uniprot/uniprot_sprot.fasta --out ".$o->{"job"}."/diamond";
    }
    elsif ( $program eq "pannzer" ) {
        return "python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -H localhost -T localhost -d ".$o->{"cwd"}."/databases/pannzer -i ".$o->{"job"}."/input_annotation -o ,".$o->{"job"}."/pannzer,,";
    }
    elsif ( $program eq "emapper" ) {
        return "PATH=\$PATH:".$o->{"cwd"}."/bin/diamond ".$o->{"cwd"}."/bin/emapper/emapper.py --override --cpu 1 --data_dir ".$o->{"cwd"}."/databases/emapper -m diamond -i ".$o->{"job"}."/input_annotation -o ".$o->{"job"}."/emapper";
    }
    elsif ( $program eq "trf" ) {
        return $o->{"cwd"}."/bin/trf/trf ".$o->{"job"}."/input_sequences 2 7 7 80 10 50 500 -h -ngs | tee ".$o->{"job"}."/trf";
    }
    elsif ( $program eq "pilercr" ) {
        return $o->{"cwd"}."/bin/pilercr/pilercr -in ".$o->{"job"}."/input_sequences -out ".$o->{"job"}."/pilercr -noinfo -quiet";
    }
    elsif ( $program eq "minced" ) {
        return $o->{"cwd"}."/bin/minced/minced -gffFull ".$o->{"job"}."/input_sequences ".$o->{"job"}."/minced.tmp ".$o->{"job"}."/minced";
    }

}

1;
