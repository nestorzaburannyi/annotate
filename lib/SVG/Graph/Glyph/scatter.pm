package SVG::Graph::Glyph::scatter;

use base SVG::Graph::Glyph;
use strict;

=head1 NAME

SVG::Graph::Glyph::scatter - SVG::Graph::Glyph::scatter object

=head1 METHODS

=head2 draw

 Title   : draw
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub draw{
   my ($self,@args) = @_;

   my $id = 'n'.sprintf("%07d",int(rand(9999999)));
   my $group = $self->svg->group(id=>"scatter$id");

   foreach my $datum($self->group->data){
	 my $cx = (($datum->x - $self->group->xmin) * $self->xscale) + $self->xoffset;
	 my $cy = (($self->xsize - ($datum->y - $self->group->ymin) * $self->yscale)) + $self->yoffset;

	 $group->circle(cx=>$cx,cy=>$cy,r=>3,style=>{$self->_style});
   }
}

1;
