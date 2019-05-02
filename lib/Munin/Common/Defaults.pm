use warnings;
use strict;

# If you change the class path take a look in get_defaults too, please!
package Munin::Common::Defaults;

# $Id$

use English qw(-no_match_vars);
use File::Basename qw(dirname);

# This file's package variables are changed during the build process.

# This variable makes only sense in development environment
my $COMPONENT_ROOT = dirname(__FILE__) . '/../../..';


our $DROPDOWNLIMIT     = 1;

our $MUNIN_PREFIX     = q{/usr};
our $MUNIN_CONFDIR    = q{/etc/munin};
our $MUNIN_BINDIR     = q{/usr/bin};
our $MUNIN_SBINDIR    = q{/usr/sbin};
our $MUNIN_DOCDIR     = q{/usr/doc};
our $MUNIN_LIBDIR     = q{/usr/share/munin};
our $MUNIN_HTMLDIR    = q{/var/cache/munin/www};
our $MUNIN_CGIDIR     = q{/usr/lib/munin/cgi};
our $MUNIN_CGITMPDIR     = '';
our $MUNIN_DBDIR      = q{/var/lib/munin};
our $MUNIN_PLUGSTATE  = q{/var/lib/munin-node/plugin-state};
our $MUNIN_SPOOLDIR   = q{/var/lib/munin-async};
our $MUNIN_MANDIR     = q{/usr/share/man};
our $MUNIN_LOGDIR     = q{/var/log/munin};
our $MUNIN_STATEDIR   = q{/var/run/munin};
our $MUNIN_USER       = q{munin};
our $MUNIN_GROUP      = q{munin};
our $MUNIN_PLUGINUSER = q{nobody};
our $MUNIN_VERSION    = q{2.0.33-1};
our $MUNIN_PERL       = q{/usr/bin/perl};
our $MUNIN_PERLLIB    = q{/usr/share/perl5};
our $MUNIN_GOODSH     = q{/bin/sh};
our $MUNIN_BASH       = q{/bin/bash};
our $MUNIN_PYTHON     = q{/usr/bin/env python};
our $MUNIN_RUBY       = q{/usr/bin/env ruby};
our $MUNIN_OSTYPE     = q{linux};
our $MUNIN_HOSTNAME   = q{localhost.localdomain};
our $MUNIN_MKTEMP     = q{mktemp -p /tmp/ MKTEMP     = };
our $MUNIN_HASSETR    = q{1};


sub get_defaults {
    my ($class) = @_;
    
    ## no critic

    no strict 'refs';
    my $defaults = {};
    for my $g (keys %{"Munin::Common::Defaults::"}) {
        next unless $g =~ /MUNIN_/;
        $defaults->{$g} = ${*$g{'SCALAR'}};
    }

    ## use critic

    return $defaults;
}


sub export_to_environment {
    my ($class) = @_;

    my %defaults = %{$class->get_defaults()};
    while (my ($k, $v) = each %defaults) {
        $ENV{$k} = $v;
    }

    return
}


1;


__END__


=head1 NAME

Munin::Common::Defaults - Default values defined by installation
scripts.


=head1 PACKAGE VARIABLES

See L<http://munin-monitoring.org/wiki/MuninInstallProcedure> for
more information on the variables provided by this package.


=head1 METHODS

=over

=item B<get_defaults>

  \%defaults = $class->get_defaults()

Returns all the package variables as key value pairs in a hash.

=item B<export_to_environment>

  $class = $class->export_to_environment()

Export all the package variables to the environment.

=back

