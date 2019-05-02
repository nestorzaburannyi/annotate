package SWISS::CCcofactor;

use vars qw($AUTOLOAD @ISA %fields);

use Carp;
use strict;
use SWISS::TextFunc;
use SWISS::ListBase;


BEGIN {
  @ISA = ( 'SWISS::ListBase' );
  
  %fields = (
          form    => undef,  
	      note    => undef
	      #list
	    );
}


sub new {
  my $ref = shift;
  my $class = ref($ref) || $ref;
  my $self = new SWISS::ListBase;
  $self->rebless($class);
  return $self;
}


sub fromText {
    my $class = shift;
    my $textRef = shift;
    my $self = new SWISS::CCcofactor;
    my $text = $$textRef;
    $self->initialize();
  
    $text =~ s/ *-!- COFACTOR: +//;
    $text =~ s/; {2,}/; /g;
    $text =~ s/, {2,}/, /g;
    
    if ( $text !~ /Name=/ && $text !~ /Note=/ ) { # old format!
    	  my $note = $text;
    	( my $note_ev = $1 ) if $note =~ s/($SWISS::TextFunc::evidencePattern)\.?//;
    	$self->{ note } = [ $note.".", $note_ev ] if $note;
    }
    else { # new format
        ( my $form    = $1 )           if $text =~ s/^([^:]+): (?=N)//;
        ( my $note    = $1 ) =~ s/;$// if $text =~ s/ ?Note=(.+)$//;
          my $note_ev = $1             if $note =~ s/($SWISS::TextFunc::evidencePattern)//;
        $self->{ form }  = $form if $form;
        $self->{ note }  = [ $note, $note_ev ] if $note;
        foreach my $name ( split / +(?=Name=)/, $text ) {
            my $ev = $name =~ s/ +Evidence=($SWISS::TextFunc::evidencePattern);// ? " ".$1 : undef;
            $self->add( [ $name, $ev ] );
        }	
    }

    $self->{_dirty} = 0;
    return $self;
}


sub topic {
	return "COFACTOR";
}


sub toString {
    my ( $self ) = @_;
    
    my $form = $self->{ form } ? " $self->{ form }:" : "";
    my $text = "CC   -!- COFACTOR:$form\n";
    
    foreach my $name_ev ( @{ $self->{ list } } ) {
        my ( $name, $ev ) = @$name_ev;
        my   $line = $name;
        if ( $ev ) { $ev=~s/^ +//; $line .= " Evidence=".$ev.";" }
        $text .= SWISS::TextFunc->wrapOn( 'CC       ',"CC         ", $SWISS::TextFunc::lineLength, $line, "(?<=;) ", "(?<=,) ", $SWISS::TextFunc::textWrapPattern1 );
    }
    if ( $self->{ note } ) {
    	my ( $note, $note_ev ) = @{ $self->{ note } };
        $note_ev ||= "";
        $note      = "Note=" . $note . $note_ev . ";";
        $text     .= SWISS::TextFunc->wrapOn( 'CC       ',"CC       ", $SWISS::TextFunc::lineLength, $note );
    }
    
    return $text; 
}


sub comment {
    my ( $self, $with_ev ) = @_;
    
    my $text = "";
    foreach my $name_ev ( @{ $self->{ list } } ) {
        my ( $name, $ev ) = @$name_ev;
        $text .= ( $text ? " " : "" ) . $name;
        if ( $ev && $with_ev ) { $ev=~s/^ +//; $text .= " Evidence=".$ev.";" }
    }
    if ( $self->{ note } ) {
    	my ( $note, $note_ev ) = @{ $self->{ note } };
    	$note_ev = "" unless $with_ev && $note_ev;
        $note      = "Note=" . $note . $note_ev . ";";
        $text .= ( $text ? " " : "" ) . $note;
    }
    
    return $text;
}


#sub sort {
#    my $self = shift;
#    $self->{ list } = sort { lc( $a->[0] ) cmp lc( $b->[0] ) } @{ $self->{ list } };
#}


1;

__END__

=head1 Name

SWISS::CCcofactor

=head1 Description

B<SWISS::CCcofactor> represents a comment on the topic 'COFACTOR'
within a Swiss-Prot or TrEMBL entry as specified in the user manual
http://www.expasy.org/sprot/userman.html .  Comments on other topics are stored
in other types of objects, such as SWISS::CC (see SWISS::CCs for more information).

Collectively, comments of all types are stored within a SWISS::CCs container
object.

=head1 Inherits from

SWISS::ListBase.pm

=head1 Attributes

=over

=item topic

The topic of this comment ('COFACTOR').

=item comment

The "text" version of this comment (without evidences and new lines).

=item note

The note and evidence (Note= in new format or full description in old format)
reference to an array of [ $note, $note_ev ] (strings)

=item elements

An array of [name_str, evidence_tags_str], if any.

=back
=head1 Methods

=head2 Standard methods

=over

=item new

=item fromText

=item toString

Returns a string representation of this comment.

=back
