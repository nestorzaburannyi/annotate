package SAP::Update;

use strict;
use warnings;
use SAP::Common;
use Net::FTP;
use LWP::Simple;
use SWISS::Entry;
use Archive::Extract;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(update);

sub update {
    my ( $o ) = @_;
    print_log( $o, "Initializing / updating databases..." );
    mkdir( $o->{"cwd"}."/databases" );
    sqlite ( $o );
    taxonomy ( $o );
    rfam ( $o );
    sprot ( $o );
    pannzer ( $o );
    emapper ( $o );
    print_log( $o, "Updating databases complete." );
    exit 0;
}

sub sqlite {
    # empty database, just in case we decide to store something, someday in it
    my ( $o ) = @_;
    print_log( $o, "Initializing SQL database..." );
    mkdir( $o->{"cwd"}."/databases/sqlite" );
    my $database = DBI->connect("dbi:SQLite:".$o->{"cwd"}."/databases/sqlite/database");
}

sub taxonomy {
    my ( $o ) = @_;
    mkdir ( $o->{"cwd"}."/databases/taxonomy" );
    print_log( $o, "Downloading NCBI Taxonomy database..." );
    download_file ( $o, "ftp", "ftp.ncbi.nlm.nih.gov", "pub/taxonomy/taxdump.tar.gz", $o->{"cwd"}."/databases/taxonomy/taxdump.tar.gz" );
    decompress_file ( $o, $o->{"cwd"}."/databases/taxonomy/taxdump.tar.gz", $o->{"cwd"}."/databases/taxonomy/" );
    $o->{"dbh_taxonomy"} = Bio::DB::Taxonomy->new( -force => 1, -source => "flatfile", -directory => $o->{"cwd"}."/databases/taxonomy", -nodesfile => $o->{"cwd"}."/databases/taxonomy/nodes.dmp", -namesfile => $o->{"cwd"}."/databases/taxonomy/names.dmp" );
    return;
}

sub rfam {
    my ( $o ) = @_;
    mkdir ( $o->{"cwd"}."/databases/rfam" );
    print_log( $o, "Downloading RFAM taxonomy information..." );
    download_and_uncompress_file ( $o, "ftp", "ftp.ebi.ac.uk", "pub/databases/Rfam/CURRENT/database_files/family_ncbi.txt.gz", $o->{"cwd"}."/databases/rfam/family_ncbi.txt" );
    print_log( $o, "Parsing RFAM taxonomy information..." );
    my %rfam_records;
    while (my $line = parse_file( $o, $o->{"cwd"}."/databases/rfam/family_ncbi.txt", "line", "\t", "" )) {
        # if all the possible taxonomies are present, no need to iterate the class again
        next if ( exists $rfam_records{$line->[2]} and exists $rfam_records{$line->[2]}{"taxonomy"} and exists $rfam_records{$line->[2]}{"bacteria"} and exists $rfam_records{$line->[2]}{"archaea"} );
        if ( not exists $rfam_records{$line->[2]}{"bacteria"} and taxon_id_belongs_to ( $o, $line->[0], 2 ) ) {
            $rfam_records{$line->[2]}{"taxonomy"} = 1;
            $rfam_records{$line->[2]}{"bacteria"} = 1;
        }
        if ( not exists $rfam_records{$line->[2]}{"archaea"} and taxon_id_belongs_to ( $o, $line->[0], 2158 ) ) {
            $rfam_records{$line->[2]}{"taxonomy"} = 1;
            $rfam_records{$line->[2]}{"archaea"} = 1;
        }
    }
    print_log( $o, "Downloading RFAM family information..." );
    download_and_uncompress_file ( $o, "ftp", "ftp.ebi.ac.uk", "pub/databases/Rfam/CURRENT/database_files/family.txt.gz", $o->{"cwd"}."/databases/rfam/family.txt" );
    print_log( $o, "Parsing RFAM family information..." );
    while ( my $line = parse_file( $o, $o->{"cwd"}."/databases/rfam/family.txt", "line", "\t", "" ) ) {
        next if not ( exists $rfam_records{$line->[0]} and exists $rfam_records{$line->[0]}{"taxonomy"} );
        $rfam_records{$line->[0]}{"name"} = $line->[1];
        $rfam_records{$line->[0]}{"product"} = $line->[3];
        $rfam_records{$line->[0]}{"type"} = $line->[18];
    }
    print_log( $o, "Downloading RFAM database..." );
    download_and_uncompress_file ( $o, "ftp", "ftp.ebi.ac.uk", "pub/databases/Rfam/CURRENT/Rfam.cm.gz", $o->{"cwd"}."/databases/rfam/Rfam.cm" );
    print_log( $o, "Parsing RFAM database..." );
    open my $bacteria_filehandle, ">", $o->{"cwd"}."/databases/rfam/rfam_bacteria.cm" or die "Could not open ".$o->{"cwd"}."/databases/rfam/rfam_bacteria.cm - $!";
    open my $archaea_filehandle, ">", $o->{"cwd"}."/databases/rfam/rfam_archaea.cm" or die "Could not open ".$o->{"cwd"}."/databases/rfam/rfam_archaea.cm - $!";
    while ( my $line = parse_file( $o, $o->{"cwd"}."/databases/rfam/Rfam.cm", "record", "\n//\n", "" ) ) {
        my ( $spacer, $accession ) = ( $line =~ m/ACC\s(\s+)(\S+)\n/ ) or next;
        next if not ( exists $rfam_records{$accession} and exists $rfam_records{$accession}{"taxonomy"} );
        $line =~ s/ACC\s+$accession\n/ACC $spacer$accession\nDESC$spacer$rfam_records{$accession}{"type"} $rfam_records{$accession}{"product"}\n/;
        print {$bacteria_filehandle} $line if ( exists $rfam_records{$accession}{"bacteria"} );
        print {$archaea_filehandle} $line if ( exists $rfam_records{$accession}{"archaea"} );
        }
    close( $bacteria_filehandle );
    close( $archaea_filehandle );
    print_log( $o, "Building RFAM databases..." );
    system( $o->{"cwd"}."/bin/infernal/cmpress -F ".$o->{"cwd"}."/databases/rfam/rfam_bacteria.cm" );
    system( $o->{"cwd"}."/bin/infernal/cmpress -F ".$o->{"cwd"}."/databases/rfam/rfam_archaea.cm" );
    return;
}

sub sprot {
    my ( $o ) = @_;
    mkdir ( $o->{"cwd"}."/databases/uniprot" );
    print_log( $o, "Downloading Uniprot/Swissprot..." );
    download_and_uncompress_file ( $o, "ftp", "ftp.uniprot.org", "pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz", $o->{"cwd"}."/databases/uniprot/uniprot_sprot.dat" );
    print_log( $o, "Parsing Uniprot/Swissprot..." );
    # consider double slash as a record division
    local $/ = "\n//\n";
    # opening input file
    open my $input_filehandle, "<", $o->{"cwd"}."/databases/uniprot/uniprot_sprot.dat" or die "Could not open ".$o->{"cwd"}."/databases/uniprot/uniprot_sprot.dat - $!";
    open my $bacteria_filehandle, ">", $o->{"cwd"}."/databases/uniprot/sprot_bacteria.fasta" or die "Could not open ".$o->{"cwd"}."/databases/uniprot/sprot_bacteria.fasta - $!";
    open my $archaea_filehandle, ">", $o->{"cwd"}."/databases/uniprot/sprot_archaea.fasta" or die "Could not open ".$o->{"cwd"}."/databases/uniprot/sprot_archaea.fasta - $!";
    while ( my $line = <$input_filehandle> ) {
        my $entry = SWISS::Entry->fromText( $line );
        # do not consider partial proteins as reliable, unless requested
        $o->{"update-sprot-complete-only"} = 1;
        next if ( ( $entry->isFragment or $entry->DEs->hasFragment ) and ( $o->{"update-sprot-complete-only"} ) );
        # get the protein evidence
        my $pe = substr( $entry->PE->{"text"}, 0, 1 ) or die "Could not get protein evidence tag from:\n $line\n\nThis should not have happened!";
        # do not consider proteins with no protein evidence as reliable, unless requested
        $o->{"update-sprot-protein-evidence"} = 1;
        next if ( $o->{"update-sprot-protein-evidence"} > $pe );
        # get the taxid
        my $taxon_id = ( $entry->OXs->NCBI_TaxID->elements )[0]->text;
        my $taxon;
        # limit db to certain taxons
        if ( taxon_id_belongs_to ( $o, $taxon_id, 2 ) ) {
            $taxon = "bacteria";
        }
        elsif ( taxon_id_belongs_to ( $o, $taxon_id, 2157 ) ) {
            $taxon = "archaea";
        }
        else {
            next;
        }
        # get the ID
        my $protein_id = $entry->ID or die "Could not get protein ID from:\n $line\n\nThis should not have happened!";
        # get the sequence
        my $sequence = $entry->SQ or die "Could not get protein sequence from:\n $line\n\nThis should not have happened!";
        # set the non-standard tags
        my $non_std = $line =~ m/ID   IF3_/ ? "_UNUSUAL_START_ATT" : "";
        $non_std = $line =~ m/ID   PCNB_/ ? "$non_std"."_UNUSUAL_START_ATT" : $non_std;
        $non_std = $line =~ m/NON_STD\s+\d+\s+\d+\s+Selenocysteine/ ? "$non_std"."_UNUSUAL_NONSTOP_TGA" : $non_std;
        $non_std = $line =~ m/NON_STD\s+\d+\s+\d+\s+Pyrrolysine/ ? "$non_std"."_UNUSUAL_NONSTOP_TAG" : $non_std;
        # write the sequence
        if ( $taxon eq "bacteria" ) {
            print { $bacteria_filehandle } ">$protein_id$non_std\n$sequence\n";
        }
        elsif ( $taxon eq "archaea" ) {
            print { $archaea_filehandle } ">$protein_id$non_std\n$sequence\n";
        }
    }
    close ( $input_filehandle );
    close ( $bacteria_filehandle );
    close ( $archaea_filehandle );
    # build sprot databases
    system( $o->{"cwd"}."/bin/blast/makeblastdb -dbtype prot -in ".$o->{"cwd"}."/databases/uniprot/sprot_bacteria.fasta" );
    system( $o->{"cwd"}."/bin/blast/makeblastdb -dbtype prot -in ".$o->{"cwd"}."/databases/uniprot/sprot_archaea.fasta" );
    system( $o->{"cwd"}."/bin/diamond/diamond makedb --in ".$o->{"cwd"}."/databases/uniprot/sprot_bacteria.fasta --db ".$o->{"cwd"}."/databases/uniprot/sprot_bacteria.fasta" );
    system( $o->{"cwd"}."/bin/diamond/diamond makedb --in ".$o->{"cwd"}."/databases/uniprot/sprot_archaea.fasta --db ".$o->{"cwd"}."/databases/uniprot/sprot_archaea.fasta" );
    return;
}

sub pannzer {
    my ( $o ) = @_;
    print_log( $o, "Updating Pannzer..." );
    # update instructions taken from http://ekhidna2.biocenter.helsinki.fi/sanspanz/
    mkdir ( $o->{"cwd"}."/databases/pannzer" );
    print_log( $o, "Downloading GO..." );
    system( "wget 'https://www.uniprot.org:443/taxonomy/?query=*&compress=yes&format=tab' -O - | gzip -d > ".$o->{"cwd"}."/databases/pannzer/taxonomy-all.tab" );
    # EC and KEGG mappings to GO
    print_log( $o, "EC and KEGG mappings to GO..." );
    system( "wget http://geneontology.org/external2go/ec2go -O - | perl ".$o->{"cwd"}."/bin/pannzer/uniprot/external2go.pl > ".$o->{"cwd"}."/databases/pannzer/ec2go.tab" );
    system( "wget http://geneontology.org/external2go/kegg2go -O - | perl ".$o->{"cwd"}."/bin/pannzer/uniprot/external2go.pl > ".$o->{"cwd"}."/databases/pannzer/kegg2go.tab" );
    # get, format and index uniprot
    print_log( $o, "Downloading uniprot..." );
    system( "wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O ".$o->{"cwd"}."/databases/pannzer/uniprot_sprot.fasta.gz" );
    system( "wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz -O ".$o->{"cwd"}."/databases/pannzer/uniprot_trembl.fasta.gz" );
    print_log( $o, "Indexing uniprot..." );
    system( "perl ".$o->{"cwd"}."/bin/pannzer/saisformatdb.pl ".$o->{"cwd"}."/databases/pannzer/uniprot ".$o->{"cwd"}."/databases/pannzer/uniprot_sprot.fasta.gz ".$o->{"cwd"}."/databases/pannzer/uniprot_trembl.fasta.gz" );
    # update counts: uniprot.phr is indexed by SANSparallel and contains the description lines of the SANSparallel sequence database
    print_log( $o, "Getting counts..." );
    system( "perl -pe 's/^\\S+\\s+//' ".$o->{"cwd"}."/databases/pannzer/uniprot.phr| perl -pe 's/ \\w+\\=.*//' | python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -m Cleandesc -f tab -c desc -b desc | cut -f 2 > ".$o->{"cwd"}."/databases/pannzer/x" );
    system( "sort ".$o->{"cwd"}."/databases/pannzer/x | uniq -c > ".$o->{"cwd"}."/databases/pannzer/uniprot.desc.uc.counts" );
    system( "perl -pe 's/ +/\n/g' ".$o->{"cwd"}."/databases/pannzer/x | sort | uniq -c > ".$o->{"cwd"}."/databases/pannzer/uniprot.word.uc.counts" );
    system( "sort ".$o->{"cwd"}."/databases/pannzer/x | uniq -c > ".$o->{"cwd"}."/databases/pannzer/uniprot.desc.uc.counts" );
    system( "sort ".$o->{"cwd"}."/databases/pannzer/x | uniq -c > ".$o->{"cwd"}."/databases/pannzer/uniprot.desc.uc.counts" );
    system( "perl ".$o->{"cwd"}."/bin/pannzer/uniprot/sumcol.pl 1 < ".$o->{"cwd"}."/databases/pannzer/uniprot.desc.uc.counts > ".$o->{"cwd"}."/databases/pannzer/nprot" );
    system( "perl ".$o->{"cwd"}."/bin/pannzer/uniprot/sumcol.pl 1 < ".$o->{"cwd"}."/databases/pannzer/uniprot.word.uc.counts > ".$o->{"cwd"}."/databases/pannzer/nwordtotal" );
    # GO dictionaries
    system( "wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz -O - | gzip -d > ".$o->{"cwd"}."/databases/pannzer/goa_uniprot_all.gaf" );
    system( "wget http://geneontology.org/ontology/go-basic.obo -O ".$o->{"cwd"}."/databases/pannzer/go-basic.obo" );
    # GO dictionaries using gorelatives
    # GO hierarchy and information content of GO terms
    system( "python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -m obo -i ".$o->{"cwd"}."/databases/pannzer/go-basic.obo -o '".$o->{"cwd"}."/obo.tab,'" );
    system( "cut -f 2,5,7 ".$o->{"cwd"}."/databases/pannzer/goa_uniprot_all.gaf | python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -m gaf2propagated_with_evidencecode -f tab -c 'qpid goid evidence_code' -o ',".$o->{"cwd"}."/databases/pannzer/godict.txt' 2> ".$o->{"cwd"}."/databases/pannzer/err" );
    # workaround because of hard-coded obo.tab path in gaf2propagated_with_evidencecode
    system( "mv ".$o->{"cwd"}."/obo.tab"." ".$o->{"cwd"}."/databases/pannzer/obo.tab" );
    system( "cut -f 2 ".$o->{"cwd"}."/databases/pannzer/godict.txt | sort -T ".$o->{"cwd"}."/databases/pannzer | uniq -c | perl -pe 's/^\\s+//' | perl -pe 's/ /\\t/' > ".$o->{"cwd"}."/databases/pannzer/godict_counts" );
    system( "python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -m BayesIC -i ".$o->{"cwd"}."/databases/pannzer/godict_counts -f tab -c 'qpid propagated' --eval_OBOTAB ".$o->{"cwd"}."/databases/pannzer/obo.tab -o ',".$o->{"cwd"}."/databases/pannzer/obo_with_ic.tab'" );
    # GO propagated parent list using gorelatives
    system( "grep 'id: GO:' ".$o->{"cwd"}."/databases/pannzer/go-basic.obo | perl -pe 's/\\S*id: GO://' > ".$o->{"cwd"}."/databases/pannzer/go.list" );
    system( $o->{"cwd"}."/bin/pannzer/uniprot/gorelatives -b ".$o->{"cwd"}."/databases/pannzer/go-basic.obo -q ".$o->{"cwd"}."/databases/pannzer/go.list -r isa,partof,altid -d parents -l > ".$o->{"cwd"}."/databases/pannzer/go_data" );
    system( "perl ".$o->{"cwd"}."/bin/pannzer/uniprot/generate_godict.pl ".$o->{"cwd"}."/databases/pannzer/go_data ".$o->{"cwd"}."/databases/pannzer/obo_with_ic.tab ".$o->{"cwd"}."/databases/pannzer/ec2go.tab ".$o->{"cwd"}."/databases/pannzer/kegg2go.tab > ".$o->{"cwd"}."/databases/pannzer/mergeGO.out" );
    # all done, start the servers
    system( "nohup python ".$o->{"cwd"}."/bin/pannzer/DictServer.py -H localhost -d ".$o->{"cwd"}."/databases/pannzer&" );
    system( "nohup mpirun --allow-run-as-root -np \$(nproc --all) -output-filename ".$o->{"cwd"}."/bin/pannzer/log ".$o->{"cwd"}."/bin/pannzer/server ".$o->{"cwd"}."/databases/pannzer/uniprot 54321 uniprot&" );
}

sub emapper {
  my ( $o ) = @_;
  print_log( $o, "Updating emapper..." );
  mkdir( $o->{"cwd"}."/databases/emapper" );
  system( $o->{"cwd"}."/bin/emapper/download_eggnog_data.py -y -f --data_dir ".$o->{"cwd"}."/databases/emapper/ none" );
}

sub download_file {
    my ( $o, $protocol, $host, $path, $target) = @_;
    if ( $protocol eq "ftp" ) {
        my $ftp = Net::FTP->new( $host, Passive => 1 ) or exit_program ( $o, "Cannot connect to $protocol://$host: $@" );
        $ftp->login or exit_program ( $o, "Cannot login, ".( $ftp->message));
        $ftp->binary or exit_program ( $o, "Cannot set binary mode, ".( $ftp->message));
        $ftp->get( $path, $target) or exit_program ( $o, "Cannot get file $protocol://$host/$path, ".( $ftp->message ) );
        $ftp->quit;
    }
    elsif ( $protocol eq "http" ) {
        my $ua = LWP::UserAgent->new;
        my $response = $ua->get( "$protocol://$host/$path" );
        if ($response->is_success) {
            open my $output_filehandle, ">", "$target" or die "Could not open $target - $!";
            print {$output_filehandle} $response->decoded_content;
            close $output_filehandle;
        }
        else {
            exit_program ( $o, "Cannot get file $protocol://$host/$path, ".(  $response->status_line ) );
        }
    }
    return;
}

sub decompress_file {
    my ( $o, $source, $target) = @_;
    my $ae = Archive::Extract->new (archive => $source);
    $ae->extract (to => $target) or exit_program ( $o, "Decompression failed: $ae->error" );
    return;
}

sub download_and_uncompress_file {
  # a wrapper to combine download_file() and decompress_file() in one call
  my ( $o, $protocol, $host, $file, $target) = @_;
  download_file ( $o, $protocol, $host, $file, $target.".gz" );
  decompress_file ( $o, $target.".gz", $target);
  return;
}

1;
