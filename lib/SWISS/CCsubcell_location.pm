package SWISS::CCsubcell_location;

use vars qw($AUTOLOAD @ISA @_properties %fields);

use Carp;
use strict;
use SWISS::TextFunc;
use SWISS::BaseClass;

BEGIN {
    @ISA = ('SWISS::BaseClass');

    %fields = (
        locations => undef,
        note      => undef
    );
        
}


sub new {
  my $ref = shift;
  my $class = ref($ref) || $ref;
  my $self = new SWISS::BaseClass;
  $self->rebless($class);
  return $self;
}


sub fromText {
    my $class   = shift;
    my $textRef = shift;
    my $self    = new SWISS::CCsubcell_location;

    my $text = $$textRef;
    $self->initialize();

    $self->{ locations } = undef;
    $self->{ note }      = undef;
    
    $text =~ s/^ *-!- SUBCELLULAR LOCATION: +//;
    # e.g.: isoform1:Cell membrane; Lipid-anchor, GPI-anchor. Nucleus membrane {ECO:0000269|PubMed:15282802}.
    #     = form....:component....; topology................. a 2nd location...
    
    my ( $core, $note ) = split /\.? Note=/, $text;
    foreach my $lstr ( split /\. /, $core ) {
    	my $form = $lstr =~ s/^(.+?)\: // ? $1 : "";
    	my @locs = split /; /, $lstr;
    	push( @{ $self->{ locations } }, 
    	   {
    	       'form'        => $form, 
    	       'component'   => _parse_txt_ev( $locs[0] ),
    	       'topology'    => $locs[1] ? _parse_txt_ev( $locs[1] ) : undef,
    	       'orientation' => $locs[2] ? _parse_txt_ev( $locs[2] ) : undef
    	   }
    	);
    }
    $self->{note} = _parse_txt_ev( $note ) if $note;
    
    return $self;
}


sub _parse_txt_ev {
    my $txt = shift;
    
    my ( $evidence ) = $txt =~ /($SWISS::TextFunc::evidencePattern)/m;
    if ( $evidence ) {
        my $quotedEvidence = quotemeta $evidence;
        $txt =~ s/$quotedEvidence//m;
    }
    $txt =~ s/\.$//;
    return [ $txt, $evidence ];
}


sub topic {
    return "SUBCELLULAR LOCATION";
}


sub toString {
    my $self = shift;
    
    my $text = "CC   -!- SUBCELLULAR LOCATION: " . $self->comment( "true" );
    
    $text = SWISS::TextFunc->wrapOn( '', "CC       ", $SWISS::TextFunc::lineLength, $text);
    
    return $text;
}


sub comment {
    my ( $self, $with_ev ) = @_;
    
    my $text = "";
    
    foreach my $location ( @{ $self->{locations} } ) {
    	$text .= " " if $text;
    	$text .= $location->{ form }.": " if $location->{ form };
    	my $component    = $location->{ component }->[0];
    	my $comp_ev      = $location->{ component }->[1] || "";
    	$text .= $component;
    	$text .= $comp_ev if $with_ev;
    	if ( $location->{ topology } ) {
            my $topology = $location->{ topology }->[0];
            my $topo_ev  = $location->{ topology }->[1] || "";
            $text .= "; ".$topology;
            $text .= $topo_ev if $with_ev;		
    	}
    	if ( $location->{ orientation } ) {
            my $orientation = $location->{ orientation }->[0];
            my $orien_ev    = $location->{ orientation }->[1] || "";
            $text .= "; ".$orientation;
            $text .= $orien_ev if $with_ev;   
        }
        $text .= ".";
    }
    
    if ( $self->{ note } ) {
    	$text .= " " if $text;
        my $note      = $self->{ note }->[0];
        my $note_ev   = $self->{ note }->[1];
        my $is_new_ev = defined $note_ev ? $note_ev =~ /ECO:/ : "";
        $text .= "Note=".$note;
        if ( $note_ev ) {
        	my $extra = $is_new_ev ? "." : "";
        	$text .= $extra.$note_ev
        }
        $text .= ".";  
    }
    
    return $text;
}


sub locations {
	my ( $self, $value ) = @_;
	
    if ( defined $value ) {
        $self->{ locations } = $value;
    }
    	
	return $self->{ locations };
}


sub note {
    my ( $self, $value, $ev ) = @_;
    
    if ( defined $value ) {
        my $new_ev = $ev ? $ev : $self->{ note }->[0];
        $self->{ note } = [ $value, $new_ev ];
    }
    
    return $self->{'note'};
}


1;

__END__

=head1 Name

SWISS::CCsubcell_location.pm

=head1 Description

B<SWISS::CCdisease> represents a comment on the topic 'SUBCELLULAR LOCATION'
within a Swiss-Prot or TrEMBL entry as specified in the user manual
http://www.expasy.org/sprot/userman.html .  Comments on other topics are stored
in other types of objects, such as SWISS::CC (see SWISS::CCs for more information).

Collectively, comments of all types are stored within a SWISS::CCs container
object.

=head1 Inherits from

SWISS::BaseClass.pm

=head1 Attributes

=over

=item topic

The topic of this comment ('SUBCELLULAR LOCATION').

=back
=head1 Methods

=item locations

The locations; reference to an array [ 
                {   'form'        => $in_form,
                    'component'   => [ $component,   $component_ev ],
                    'topology'    => [ $topology,    $topology_ev ],
                    'orientation' => [ $orientation, $orientation_ev ]
                }, ... ]

=item locations( $new_locations )

Set locations to $new_locations, that should be an array: [ 
                {   'form'        => $in_form,
                    'component'   => [ $component,   $component_ev ],
                    'topology'    => [ $topology,    $topology_ev ],
                    'orientation' => [ $orientation, $orientation_ev ]
                }, ... ]

=item note

The note and evidence of the disease (Note= in new CC disease format or full description in old format)
reference to an array of [ $note, $note_ev ] (strings)

=item note( [ $new_note, $new_note_ev ] )

Set note to [ $new_note, $new_note_ev ]

=item comment

The "text" version of this comment.

=head2 Standard methods

=over

=item new

=item fromText

=item toString

Returns a string representation of this comment.

=back
