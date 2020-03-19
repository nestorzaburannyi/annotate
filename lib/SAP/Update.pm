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
  my ( $o, $cwd ) = @_;
  # set the option to cwd
  $o->{"cwd"} = $cwd;
  # update job has a special job uuid
  $o->{"job_uuid"} = "update";
  # the the base for output folder for the log files
  $o->{"job"} = $o->{"cwd"};
  print_log( $o, "Initializing / updating databases..." );
  mkdir( $o->{"cwd"}."/databases" );
  taxonomy ( $o );
  rfam ( $o );
  sprot ( $o );
  pannzer ( $o );
  emapper ( $o );
  print_log( $o, "Updating databases complete." );
  exit 0;
}

sub taxonomy {
  my ( $o ) = @_;
  mkdir ( $o->{"cwd"}."/databases/taxonomy" );
  print_log( $o, "Downloading NCBI Taxonomy database..." );
  download_file ( $o, "ftp", "ftp.ncbi.nlm.nih.gov", "pub/taxonomy/taxdump.tar.gz", $o->{"cwd"}."/databases/taxonomy/taxdump.tar.gz" );
  decompress_file ( $o, $o->{"cwd"}."/databases/taxonomy/taxdump.tar.gz", $o->{"cwd"}."/databases/taxonomy/" );
  Bio::DB::Taxonomy->new( -force => 1, -source => "flatfile", -directory => $o->{"cwd"}."/databases/taxonomy", -nodesfile => $o->{"cwd"}."/databases/taxonomy/nodes.dmp", -namesfile => $o->{"cwd"}."/databases/taxonomy/names.dmp" );
  return;
}

sub rfam {
  my ( $o ) = @_;
  mkdir ( $o->{"cwd"}."/databases/rfam" );
  print_log( $o, "Downloading RFAM family information..." );
  download_and_uncompress_file ( $o, "ftp", "ftp.ebi.ac.uk", "pub/databases/Rfam/CURRENT/database_files/family.txt.gz", $o->{"cwd"}."/databases/rfam/family.txt" );
  print_log( $o, "Downloading RFAM database..." );
  download_and_uncompress_file ( $o, "ftp", "ftp.ebi.ac.uk", "pub/databases/Rfam/CURRENT/Rfam.cm.gz", $o->{"cwd"}."/databases/rfam/Rfam.cm" );
  print_log( $o, "Building RFAM databases..." );
  run_program( $o, $o->{"cwd"}."/bin/infernal/cmpress -F ".$o->{"cwd"}."/databases/rfam/Rfam.cm" );
  print_log( $o, "Downloading RFAM clan information..." );
  download_and_uncompress_file ( $o, "ftp", "ftp.ebi.ac.uk", "pub/databases/Rfam/CURRENT/Rfam.clanin", $o->{"cwd"}."/databases/rfam/Rfam.clanin" );
  return;
}

sub sprot {
  my ( $o ) = @_;
  mkdir ( $o->{"cwd"}."/databases/uniprot" );
  print_log( $o, "Downloading Uniprot/Swissprot..." );
  download_and_uncompress_file ( $o, "ftp", "ftp.uniprot.org", "pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz", $o->{"cwd"}."/databases/uniprot/uniprot_sprot.dat" );
  download_and_uncompress_file ( $o, "ftp", "ftp.uniprot.org", "pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz", $o->{"cwd"}."/databases/uniprot/uniprot_sprot.fasta" );
  print_log( $o, "Parsing Uniprot/Swissprot..." );
  # consider double slash as a record division
  local $/ = "\n//\n";
  # opening input file
  open my $input_filehandle, "<", $o->{"cwd"}."/databases/uniprot/uniprot_sprot.dat" or die "Could not open ".$o->{"cwd"}."/databases/uniprot/uniprot_sprot.dat - $!";
  open my $output_filehandle, ">", $o->{"cwd"}."/databases/uniprot/exceptions.txt" or die "Could not open ".$o->{"cwd"}."/databases/uniprot/exceptions.txt - $!";
  while ( my $record = <$input_filehandle> ) {
    my $entry = SWISS::Entry->fromText( $record );
    # get the id
    my $id = $entry->ID or die "Could not get ID from:\n $record\n\nThis should not have happened!";
    # get the accession
    my $ac = $entry->AC or die "Could not get AC from:\n $record\n\nThis should not have happened!";
    # find the non-standard tags
    while ( $record =~ m/NON_STD\s+(\d+)\s+(\d+)\s+(Se)leno(c)ysteine/g ) {
      # id from to type
      print $output_filehandle "sp|$ac|$id\t$1\t$2\taa:$3$4\n";
    }
    while ( $record =~ m/NON_STD\s+(\d+)\s+(\d+)\s+(Py)rro(l)ysine/g ) {
      # id from to type
      print $output_filehandle "sp|$ac|$id\t$1\t$2\taa:$3$4\n";
    }
  }
  # build sprot databases
  run_program( $o, $o->{"cwd"}."/bin/blast/makeblastdb -dbtype prot -in ".$o->{"cwd"}."/databases/uniprot/uniprot_sprot.fasta" );
  run_program( $o, $o->{"cwd"}."/bin/diamond/diamond makedb --in ".$o->{"cwd"}."/databases/uniprot/uniprot_sprot.fasta --db ".$o->{"cwd"}."/databases/uniprot/uniprot_sprot.fasta" );
  return;
}

sub pannzer {
    my ( $o ) = @_;
    print_log( $o, "Updating Pannzer..." );
    # update instructions taken from http://ekhidna2.biocenter.helsinki.fi/sanspanz/
    mkdir ( $o->{"cwd"}."/databases/pannzer" );
    print_log( $o, "Downloading GO..." );
    run_program( $o, "wget 'https://www.uniprot.org:443/taxonomy/?query=*&compress=yes&format=tab' -O - | gzip -d > ".$o->{"cwd"}."/databases/pannzer/taxonomy-all.tab" );
    # EC and KEGG mappings to GO
    print_log( $o, "EC and KEGG mappings to GO..." );
    run_program( $o, "wget http://geneontology.org/external2go/ec2go -O - | perl ".$o->{"cwd"}."/bin/pannzer/uniprot/external2go.pl > ".$o->{"cwd"}."/databases/pannzer/ec2go.tab" );
    run_program( $o, "wget http://geneontology.org/external2go/kegg2go -O - | perl ".$o->{"cwd"}."/bin/pannzer/uniprot/external2go.pl > ".$o->{"cwd"}."/databases/pannzer/kegg2go.tab" );
    # get, format and index uniprot
    print_log( $o, "Downloading uniprot..." );
    run_program( $o, "wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O ".$o->{"cwd"}."/databases/pannzer/uniprot_sprot.fasta.gz" );
    run_program( $o, "wget ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz -O ".$o->{"cwd"}."/databases/pannzer/uniprot_trembl.fasta.gz" );
    print_log( $o, "Indexing uniprot..." );
    run_program( $o, "perl ".$o->{"cwd"}."/bin/pannzer/saisformatdb.pl ".$o->{"cwd"}."/databases/pannzer/uniprot ".$o->{"cwd"}."/databases/pannzer/uniprot_sprot.fasta.gz ".$o->{"cwd"}."/databases/pannzer/uniprot_trembl.fasta.gz" );
    # update counts: uniprot.phr is indexed by SANSparallel and contains the description lines of the SANSparallel sequence database
    print_log( $o, "Getting counts..." );
    run_program( $o, "perl -pe 's/^\\S+\\s+//' ".$o->{"cwd"}."/databases/pannzer/uniprot.phr| perl -pe 's/ \\w+\\=.*//' | python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -m Cleandesc -f tab -c desc -b desc | cut -f 2 > ".$o->{"cwd"}."/databases/pannzer/x" );
    run_program( $o, "sort ".$o->{"cwd"}."/databases/pannzer/x | uniq -c > ".$o->{"cwd"}."/databases/pannzer/uniprot.desc.uc.counts" );
    run_program( $o, "perl -pe 's/ +/\n/g' ".$o->{"cwd"}."/databases/pannzer/x | sort | uniq -c > ".$o->{"cwd"}."/databases/pannzer/uniprot.word.uc.counts" );
    run_program( $o, "sort ".$o->{"cwd"}."/databases/pannzer/x | uniq -c > ".$o->{"cwd"}."/databases/pannzer/uniprot.desc.uc.counts" );
    run_program( $o, "sort ".$o->{"cwd"}."/databases/pannzer/x | uniq -c > ".$o->{"cwd"}."/databases/pannzer/uniprot.desc.uc.counts" );
    run_program( $o, "perl ".$o->{"cwd"}."/bin/pannzer/uniprot/sumcol.pl 1 < ".$o->{"cwd"}."/databases/pannzer/uniprot.desc.uc.counts > ".$o->{"cwd"}."/databases/pannzer/nprot" );
    run_program( $o, "perl ".$o->{"cwd"}."/bin/pannzer/uniprot/sumcol.pl 1 < ".$o->{"cwd"}."/databases/pannzer/uniprot.word.uc.counts > ".$o->{"cwd"}."/databases/pannzer/nwordtotal" );
    # GO dictionaries
    run_program( $o, "wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz -O - | gzip -d > ".$o->{"cwd"}."/databases/pannzer/goa_uniprot_all.gaf" );
    run_program( $o, "wget http://geneontology.org/ontology/go-basic.obo -O ".$o->{"cwd"}."/databases/pannzer/go-basic.obo" );
    # GO dictionaries using gorelatives
    # GO hierarchy and information content of GO terms
    run_program( $o, "python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -m obo -i ".$o->{"cwd"}."/databases/pannzer/go-basic.obo -o '".$o->{"cwd"}."/obo.tab,'" );
    run_program( $o, "cut -f 2,5,7 ".$o->{"cwd"}."/databases/pannzer/goa_uniprot_all.gaf | python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -m gaf2propagated_with_evidencecode -f tab -c 'qpid goid evidence_code' -o ',".$o->{"cwd"}."/databases/pannzer/godict.txt' 2> ".$o->{"cwd"}."/databases/pannzer/err" );
    # workaround because of hard-coded obo.tab path in gaf2propagated_with_evidencecode
    run_program( $o, "mv ".$o->{"cwd"}."/obo.tab"." ".$o->{"cwd"}."/databases/pannzer/obo.tab" );
    run_program( $o, "cut -f 2 ".$o->{"cwd"}."/databases/pannzer/godict.txt | sort -T ".$o->{"cwd"}."/databases/pannzer | uniq -c | perl -pe 's/^\\s+//' | perl -pe 's/ /\\t/' > ".$o->{"cwd"}."/databases/pannzer/godict_counts" );
    run_program( $o, "python ".$o->{"cwd"}."/bin/pannzer/runsanspanz.py -m BayesIC -i ".$o->{"cwd"}."/databases/pannzer/godict_counts -f tab -c 'qpid propagated' --eval_OBOTAB ".$o->{"cwd"}."/databases/pannzer/obo.tab -o ',".$o->{"cwd"}."/databases/pannzer/obo_with_ic.tab'" );
    # GO propagated parent list using gorelatives
    run_program( $o, "grep 'id: GO:' ".$o->{"cwd"}."/databases/pannzer/go-basic.obo | perl -pe 's/\\S*id: GO://' > ".$o->{"cwd"}."/databases/pannzer/go.list" );
    run_program( $o, $o->{"cwd"}."/bin/pannzer/uniprot/gorelatives -b ".$o->{"cwd"}."/databases/pannzer/go-basic.obo -q ".$o->{"cwd"}."/databases/pannzer/go.list -r isa,partof,altid -d parents -l > ".$o->{"cwd"}."/databases/pannzer/go_data" );
    run_program( $o, "perl ".$o->{"cwd"}."/bin/pannzer/uniprot/generate_godict.pl ".$o->{"cwd"}."/databases/pannzer/go_data ".$o->{"cwd"}."/databases/pannzer/obo_with_ic.tab ".$o->{"cwd"}."/databases/pannzer/ec2go.tab ".$o->{"cwd"}."/databases/pannzer/kegg2go.tab > ".$o->{"cwd"}."/databases/pannzer/mergeGO.out" );
    # all done, start the servers
    run_program( $o, "nohup python ".$o->{"cwd"}."/bin/pannzer/DictServer.py -H localhost -d ".$o->{"cwd"}."/databases/pannzer&" );
    run_program( $o, "nohup mpirun --allow-run-as-root -np \$(nproc --all) -output-filename ".$o->{"cwd"}."/bin/pannzer/log ".$o->{"cwd"}."/bin/pannzer/server ".$o->{"cwd"}."/databases/pannzer/uniprot 54321 uniprot&" );
}

sub emapper {
  my ( $o ) = @_;
  print_log( $o, "Updating emapper..." );
  mkdir( $o->{"cwd"}."/databases/emapper" );
  run_program( $o, $o->{"cwd"}."/bin/emapper/download_eggnog_data.py -y -f --data_dir ".$o->{"cwd"}."/databases/emapper/ none" );
}

sub download_file {
  # downloads a file, either with FTP or HTTP
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
  # decompresses a gzip file
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
