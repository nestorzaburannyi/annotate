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
    my ( $o ) = @_;
    print_log( $o, "Running tRNAscan-SE ".$o->{"trnascanse-version"}."..." );
    my ( $c ) = get_command ( $o, "trnascanse" );
    run_program ( $o, $c );
}

sub run_aragorn {
    my ( $o ) = @_;
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
    if ( $o->{"rna-tm-program"} eq "aragorn" ) {
      # skip non-tmRNA aragorn predictons
      next if not ( $l->{"product"} =~ m/^tmRNA/ );
      # product tag
      $feature->add_tag_value ( "product", $l->{"product"} );
      # inference tag
      $feature->add_tag_value ( "inference", "profile:".$o->{"rna-tm-program"}.":".$o->{$o->{"rna-tm-program"}."-version"} );
    }
    elsif ( $o->{"rna-tm-program"} eq "infernal" ) {
      # skip non-tmRNA infernal predictons
      next if not ( $l->{"hit_id"} =~ m/tmRNA$/ );
      # product tag
      $feature->add_tag_value ( "product", $l->{"product"} );
      # inference tag
      $feature->add_tag_value ( "inference", "profile:".$o->{"rna-tm-program"}.":".$o->{$o->{"rna-tm-program"}."-version"}.":rfam:".$l->{"accession"} );
    }

    # check and store tmRNA sequence feature
    check_and_store_feature ( $o, $feature );
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

    # program-specific block
    if ( $o->{"rna-nc-program"} eq "infernal" ) {
      # skip rRNA infernal predictons
      next if ( $l->{"type"} =~ m/^Gene; rRNA;$/ );
      # skip tRNA infernal predictons
      next if ( $l->{"type"} =~ m/^Gene; tRNA;$/ );
      # skip tmRNA infernal predictons
      next if ( $l->{"hit_id"} =~ m/tmRNA$/ );

      # inference tag
      $feature->add_tag_value ( "inference", "profile:".$o->{"rna-nc-program"}.":".$o->{$o->{"rna-nc-program"}."-version"}.":rfam:".$l->{"accession"} );

      # RNaseP clan
      if ( $l->{"clan"} eq "CL00002" ) {
        $feature->add_tag_value ( "ncRNA_class", $l->{"hit_id"} );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      # SRP clan
      elsif ( $l->{"clan"} eq "CL00003" ) {
        $feature->add_tag_value ( "ncRNA_class", $l->{"hit_id"} );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      # Annotate these as ncRNA with ncRNA_class="ribozyme"
      # even though it is "Gene; riboswitch;", annotate is as a ribozyme
      elsif ( $l->{"accession"} eq "RF00234" ) { # glmS glucosamine-6-phosphate activated ribozyme
          $feature->add_tag_value ( "ncRNA_class", "ribozyme" );
          $feature->add_tag_value ( "product", $l->{"product"} );
      }
      # Riboswitches should be annotated as 'regulatory' features when there is a known bound moiety and they
      # aren't defined as ribozymes. They must also include the mandatory attribute /regulatory_class="riboswitch",
      # and should include the /bound_moiety qualifier.
      elsif ( $l->{"accession"} eq "RF00050" ) { # FMN riboswitch (RFN element)
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "flavin mononucleotide" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF00059" ) { # TPP riboswitch (THI element)
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "thiamine/thiamin pyrophosphate" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"clan"} eq "CL00101" # Cobalamin riboswitch clan
           or $l->{"accession"} eq "RF00174" # Cobalamin riboswitch
           or $l->{"accession"} eq "RF01482" # AdoCbl riboswitch
           or $l->{"accession"} eq "RF01689" # AdoCbl variant RNA
            ) {
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "adenosylcobalamin" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF00504" ) { # glycine riboswitch / gcvT element
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "glycine" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF02912" ) { # AAC AAD 5' leader riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "aminoglycoside" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF00162" # SAM riboswitch (S box leader)
           or $l->{"accession"} eq "RF00521" # SAM riboswitch (alpha-proteobacteria)
           or $l->{"accession"} eq "RF00634" # S-adenosyl methionine (SAM) riboswitch
           or $l->{"accession"} eq "RF01725" # SAM-I/IV variant riboswitch
           or $l->{"accession"} eq "RF01767" # SMK box translational riboswitch
           or $l->{"accession"} eq "RF01826" # SAM-V riboswitch
           or $l->{"accession"} eq "RF02885" # SAM-VI riboswitch
          ) {
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "S-adenosylmethionine" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF01727" ) { # SAM/SAH riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "S-adenosylmethionine and/or S-adenosylhomocysteine" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF01057" ) { # S-adenosyl-L-homocysteine riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "S-adenosylhomocysteine" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF00167" ) { # Purine riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "guanine and/or adenine" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF00442" ) { # Guanidine-I riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "guanidine" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }

      elsif ( $l->{"accession"} eq "RF00168" ) { # Lysine riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "lysine" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF00522" # PreQ1 riboswitch
           or $l->{"accession"} eq "RF01054" # preQ1-II (pre queuosine) riboswitch
           or $l->{"accession"} eq "RF02680" # PreQ1-III riboswitch
            ) {
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "pre-queuosine1" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF01831" ) { # THF riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "tetrahydrofolate" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF01055" ) { # Moco (molybdenum cofactor) riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "molybdenum cofactor" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF01786" ) { # Cyclic di-GMP-II riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "cyclic di-GMP" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF01056" ) { # Magnesium Sensor
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "magnesium" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF01734" ) { # Fluoride riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "fluoride" );
        $feature->add_tag_value ( "note", $l->{"product"} );

      }
      elsif ( $l->{"accession"} eq "RF02683" ) { # NiCo riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "nickel/cobalt" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"accession"} eq "RF01739" ) { # Glutamine riboswitch
        $feature->primary_tag( "regulatory" );
        $feature->add_tag_value ( "regulatory_class", "riboswitch" );
        $feature->add_tag_value ( "bound_moiety", "glutamine" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      # Riboswitches with unknown functions or bound moieties are annotated as misc_features with the name in a note.
      # for now, includes:
      # RF00080 # yybP-ykoY manganese riboswitch
      # RF00379 # ydaO/yuaA leader
      # RF01510 # M. florum riboswitch
      # RF01750 # ZMP/ZTP riboswitch
      # RF03057 # nhaA-I RNA
      # RF03058 # sul1 RNA
      # RF03071 # DUF1646 RNA
      # RF03072 # raiA RNA
      elsif ( $l->{"type"} eq "Gene; riboswitch" ) {
          $feature->primary_tag( "misc_feature" );
          $feature->add_tag_value ( "note", $l->{"product"} );
      }
      # Annotate these as antisense ncRNA with a ncRNA_class="antisense"
      elsif ( $l->{"type"} eq "Gene; antisense;" ) {
        $feature->add_tag_value ( "ncRNA_class", "antisense" );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; antitoxin" ) {
        $feature->add_tag_value ( "ncRNA_class", "antitoxin" );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; CRISPR" ) {
        $feature->add_tag_value ( "ncRNA_class", "CRISPR" );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; lncRNA" ) {
        $feature->add_tag_value ( "ncRNA_class", "lncRNA" );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      # Annotate microRNA as ncRNA with ncRNA_class="miRNA"
      elsif ( $l->{"type"} eq "Gene; miRNA" ) {
        $feature->add_tag_value ( "ncRNA_class", "miRNA" );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; ribozyme" ) {
        $feature->add_tag_value ( "ncRNA_class", "ribozyme" );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; snRNA" ) {
        $feature->add_tag_value ( "ncRNA_class", "snRNA" );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; sRNA" ) {
        $feature->add_tag_value ( "ncRNA_class", "sRNA" );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      # Catalytic introns should be annotated as ncRNAs with ncRNA_class="autocatalytically_spliced_intron"
      elsif ( $l->{"type"} eq "Intron;" ) {
        $feature->add_tag_value ( "ncRNA_class", "autocatalytically_spliced_intron" );
        $feature->add_tag_value ( "product", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; frameshift_element" ) {
        $feature->primary_tag( "misc_feature" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; IRES;" ) {
        $feature->primary_tag( "misc_feature" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; leader;" ) {
        $feature->primary_tag( "misc_feature" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
      elsif ( $l->{"type"} eq "Gene; thermoregulator;" ) {
          $feature->primary_tag( "misc_feature" );
          $feature->add_tag_value ( "note", $l->{"product"} );
      }
      else {
        $feature->primary_tag("misc_feature" );
        $feature->add_tag_value ( "note", $l->{"product"} );
      }
    }

    # check and store ncRNA/misc_feature/misc_binding sequence feature
    check_and_store_feature ( $o, $feature );
  }
}

1;
