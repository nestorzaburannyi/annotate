package #
Date::Manip::TZ::eukiro00;
# Copyright (c) 2008-2016 Sullivan Beck.  All rights reserved.
# This program is free software; you can redistribute it and/or modify it
# under the same terms as Perl itself.

# This file was automatically generated.  Any changes to this file will
# be lost the next time 'tzdata' is run.
#    Generated on: Thu Dec  1 11:59:53 EST 2016
#    Data version: tzdata2016j
#    Code version: tzcode2016j

# This module contains data from the zoneinfo time zone database.  The original
# data was obtained from the URL:
#    ftp://ftp.iana.org/tz

use strict;
use warnings;
require 5.010000;

our (%Dates,%LastRule);
END {
   undef %Dates;
   undef %LastRule;
}

our ($VERSION);
$VERSION='6.57';
END { undef $VERSION; }

%Dates         = (
   1    =>
     [
        [ [1,1,2,0,0,0],[1,1,2,3,18,48],'+03:18:48',[3,18,48],
          'LMT',0,[1919,6,30,23,59,59],[1919,7,1,3,18,47],
          '0001010200:00:00','0001010203:18:48','1919063023:59:59','1919070103:18:47' ],
     ],
   1919 =>
     [
        [ [1919,7,1,0,0,0],[1919,7,1,3,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1930,6,20,20,59,59],[1930,6,20,23,59,59],
          '1919070100:00:00','1919070103:00:00','1930062020:59:59','1930062023:59:59' ],
     ],
   1930 =>
     [
        [ [1930,6,20,21,0,0],[1930,6,21,1,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1981,3,31,19,59,59],[1981,3,31,23,59,59],
          '1930062021:00:00','1930062101:00:00','1981033119:59:59','1981033123:59:59' ],
     ],
   1981 =>
     [
        [ [1981,3,31,20,0,0],[1981,4,1,1,0,0],'+05:00:00',[5,0,0],
          '+05',1,[1981,9,30,18,59,59],[1981,9,30,23,59,59],
          '1981033120:00:00','1981040101:00:00','1981093018:59:59','1981093023:59:59' ],
        [ [1981,9,30,19,0,0],[1981,9,30,23,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1982,3,31,19,59,59],[1982,3,31,23,59,59],
          '1981093019:00:00','1981093023:00:00','1982033119:59:59','1982033123:59:59' ],
     ],
   1982 =>
     [
        [ [1982,3,31,20,0,0],[1982,4,1,1,0,0],'+05:00:00',[5,0,0],
          '+05',1,[1982,9,30,18,59,59],[1982,9,30,23,59,59],
          '1982033120:00:00','1982040101:00:00','1982093018:59:59','1982093023:59:59' ],
        [ [1982,9,30,19,0,0],[1982,9,30,23,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1983,3,31,19,59,59],[1983,3,31,23,59,59],
          '1982093019:00:00','1982093023:00:00','1983033119:59:59','1983033123:59:59' ],
     ],
   1983 =>
     [
        [ [1983,3,31,20,0,0],[1983,4,1,1,0,0],'+05:00:00',[5,0,0],
          '+05',1,[1983,9,30,18,59,59],[1983,9,30,23,59,59],
          '1983033120:00:00','1983040101:00:00','1983093018:59:59','1983093023:59:59' ],
        [ [1983,9,30,19,0,0],[1983,9,30,23,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1984,3,31,19,59,59],[1984,3,31,23,59,59],
          '1983093019:00:00','1983093023:00:00','1984033119:59:59','1984033123:59:59' ],
     ],
   1984 =>
     [
        [ [1984,3,31,20,0,0],[1984,4,1,1,0,0],'+05:00:00',[5,0,0],
          '+05',1,[1984,9,29,21,59,59],[1984,9,30,2,59,59],
          '1984033120:00:00','1984040101:00:00','1984092921:59:59','1984093002:59:59' ],
        [ [1984,9,29,22,0,0],[1984,9,30,2,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1985,3,30,21,59,59],[1985,3,31,1,59,59],
          '1984092922:00:00','1984093002:00:00','1985033021:59:59','1985033101:59:59' ],
     ],
   1985 =>
     [
        [ [1985,3,30,22,0,0],[1985,3,31,3,0,0],'+05:00:00',[5,0,0],
          '+05',1,[1985,9,28,21,59,59],[1985,9,29,2,59,59],
          '1985033022:00:00','1985033103:00:00','1985092821:59:59','1985092902:59:59' ],
        [ [1985,9,28,22,0,0],[1985,9,29,2,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1986,3,29,21,59,59],[1986,3,30,1,59,59],
          '1985092822:00:00','1985092902:00:00','1986032921:59:59','1986033001:59:59' ],
     ],
   1986 =>
     [
        [ [1986,3,29,22,0,0],[1986,3,30,3,0,0],'+05:00:00',[5,0,0],
          '+05',1,[1986,9,27,21,59,59],[1986,9,28,2,59,59],
          '1986032922:00:00','1986033003:00:00','1986092721:59:59','1986092802:59:59' ],
        [ [1986,9,27,22,0,0],[1986,9,28,2,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1987,3,28,21,59,59],[1987,3,29,1,59,59],
          '1986092722:00:00','1986092802:00:00','1987032821:59:59','1987032901:59:59' ],
     ],
   1987 =>
     [
        [ [1987,3,28,22,0,0],[1987,3,29,3,0,0],'+05:00:00',[5,0,0],
          '+05',1,[1987,9,26,21,59,59],[1987,9,27,2,59,59],
          '1987032822:00:00','1987032903:00:00','1987092621:59:59','1987092702:59:59' ],
        [ [1987,9,26,22,0,0],[1987,9,27,2,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1988,3,26,21,59,59],[1988,3,27,1,59,59],
          '1987092622:00:00','1987092702:00:00','1988032621:59:59','1988032701:59:59' ],
     ],
   1988 =>
     [
        [ [1988,3,26,22,0,0],[1988,3,27,3,0,0],'+05:00:00',[5,0,0],
          '+05',1,[1988,9,24,21,59,59],[1988,9,25,2,59,59],
          '1988032622:00:00','1988032703:00:00','1988092421:59:59','1988092502:59:59' ],
        [ [1988,9,24,22,0,0],[1988,9,25,2,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1989,3,25,21,59,59],[1989,3,26,1,59,59],
          '1988092422:00:00','1988092502:00:00','1989032521:59:59','1989032601:59:59' ],
     ],
   1989 =>
     [
        [ [1989,3,25,22,0,0],[1989,3,26,2,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1989,9,23,22,59,59],[1989,9,24,2,59,59],
          '1989032522:00:00','1989032602:00:00','1989092322:59:59','1989092402:59:59' ],
        [ [1989,9,23,23,0,0],[1989,9,24,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1990,3,24,22,59,59],[1990,3,25,1,59,59],
          '1989092323:00:00','1989092402:00:00','1990032422:59:59','1990032501:59:59' ],
     ],
   1990 =>
     [
        [ [1990,3,24,23,0,0],[1990,3,25,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1990,9,29,22,59,59],[1990,9,30,2,59,59],
          '1990032423:00:00','1990032503:00:00','1990092922:59:59','1990093002:59:59' ],
        [ [1990,9,29,23,0,0],[1990,9,30,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1991,3,30,22,59,59],[1991,3,31,1,59,59],
          '1990092923:00:00','1990093002:00:00','1991033022:59:59','1991033101:59:59' ],
     ],
   1991 =>
     [
        [ [1991,3,30,23,0,0],[1991,3,31,3,0,0],'+04:00:00',[4,0,0],
          '+04',0,[1992,3,28,21,59,59],[1992,3,29,1,59,59],
          '1991033023:00:00','1991033103:00:00','1992032821:59:59','1992032901:59:59' ],
     ],
   1992 =>
     [
        [ [1992,3,28,22,0,0],[1992,3,29,2,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1992,9,26,22,59,59],[1992,9,27,2,59,59],
          '1992032822:00:00','1992032902:00:00','1992092622:59:59','1992092702:59:59' ],
        [ [1992,9,26,23,0,0],[1992,9,27,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1993,3,27,22,59,59],[1993,3,28,1,59,59],
          '1992092623:00:00','1992092702:00:00','1993032722:59:59','1993032801:59:59' ],
     ],
   1993 =>
     [
        [ [1993,3,27,23,0,0],[1993,3,28,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1993,9,25,22,59,59],[1993,9,26,2,59,59],
          '1993032723:00:00','1993032803:00:00','1993092522:59:59','1993092602:59:59' ],
        [ [1993,9,25,23,0,0],[1993,9,26,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1994,3,26,22,59,59],[1994,3,27,1,59,59],
          '1993092523:00:00','1993092602:00:00','1994032622:59:59','1994032701:59:59' ],
     ],
   1994 =>
     [
        [ [1994,3,26,23,0,0],[1994,3,27,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1994,9,24,22,59,59],[1994,9,25,2,59,59],
          '1994032623:00:00','1994032703:00:00','1994092422:59:59','1994092502:59:59' ],
        [ [1994,9,24,23,0,0],[1994,9,25,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1995,3,25,22,59,59],[1995,3,26,1,59,59],
          '1994092423:00:00','1994092502:00:00','1995032522:59:59','1995032601:59:59' ],
     ],
   1995 =>
     [
        [ [1995,3,25,23,0,0],[1995,3,26,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1995,9,23,22,59,59],[1995,9,24,2,59,59],
          '1995032523:00:00','1995032603:00:00','1995092322:59:59','1995092402:59:59' ],
        [ [1995,9,23,23,0,0],[1995,9,24,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1996,3,30,22,59,59],[1996,3,31,1,59,59],
          '1995092323:00:00','1995092402:00:00','1996033022:59:59','1996033101:59:59' ],
     ],
   1996 =>
     [
        [ [1996,3,30,23,0,0],[1996,3,31,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1996,10,26,22,59,59],[1996,10,27,2,59,59],
          '1996033023:00:00','1996033103:00:00','1996102622:59:59','1996102702:59:59' ],
        [ [1996,10,26,23,0,0],[1996,10,27,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1997,3,29,22,59,59],[1997,3,30,1,59,59],
          '1996102623:00:00','1996102702:00:00','1997032922:59:59','1997033001:59:59' ],
     ],
   1997 =>
     [
        [ [1997,3,29,23,0,0],[1997,3,30,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1997,10,25,22,59,59],[1997,10,26,2,59,59],
          '1997032923:00:00','1997033003:00:00','1997102522:59:59','1997102602:59:59' ],
        [ [1997,10,25,23,0,0],[1997,10,26,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1998,3,28,22,59,59],[1998,3,29,1,59,59],
          '1997102523:00:00','1997102602:00:00','1998032822:59:59','1998032901:59:59' ],
     ],
   1998 =>
     [
        [ [1998,3,28,23,0,0],[1998,3,29,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1998,10,24,22,59,59],[1998,10,25,2,59,59],
          '1998032823:00:00','1998032903:00:00','1998102422:59:59','1998102502:59:59' ],
        [ [1998,10,24,23,0,0],[1998,10,25,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[1999,3,27,22,59,59],[1999,3,28,1,59,59],
          '1998102423:00:00','1998102502:00:00','1999032722:59:59','1999032801:59:59' ],
     ],
   1999 =>
     [
        [ [1999,3,27,23,0,0],[1999,3,28,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[1999,10,30,22,59,59],[1999,10,31,2,59,59],
          '1999032723:00:00','1999032803:00:00','1999103022:59:59','1999103102:59:59' ],
        [ [1999,10,30,23,0,0],[1999,10,31,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2000,3,25,22,59,59],[2000,3,26,1,59,59],
          '1999103023:00:00','1999103102:00:00','2000032522:59:59','2000032601:59:59' ],
     ],
   2000 =>
     [
        [ [2000,3,25,23,0,0],[2000,3,26,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2000,10,28,22,59,59],[2000,10,29,2,59,59],
          '2000032523:00:00','2000032603:00:00','2000102822:59:59','2000102902:59:59' ],
        [ [2000,10,28,23,0,0],[2000,10,29,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2001,3,24,22,59,59],[2001,3,25,1,59,59],
          '2000102823:00:00','2000102902:00:00','2001032422:59:59','2001032501:59:59' ],
     ],
   2001 =>
     [
        [ [2001,3,24,23,0,0],[2001,3,25,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2001,10,27,22,59,59],[2001,10,28,2,59,59],
          '2001032423:00:00','2001032503:00:00','2001102722:59:59','2001102802:59:59' ],
        [ [2001,10,27,23,0,0],[2001,10,28,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2002,3,30,22,59,59],[2002,3,31,1,59,59],
          '2001102723:00:00','2001102802:00:00','2002033022:59:59','2002033101:59:59' ],
     ],
   2002 =>
     [
        [ [2002,3,30,23,0,0],[2002,3,31,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2002,10,26,22,59,59],[2002,10,27,2,59,59],
          '2002033023:00:00','2002033103:00:00','2002102622:59:59','2002102702:59:59' ],
        [ [2002,10,26,23,0,0],[2002,10,27,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2003,3,29,22,59,59],[2003,3,30,1,59,59],
          '2002102623:00:00','2002102702:00:00','2003032922:59:59','2003033001:59:59' ],
     ],
   2003 =>
     [
        [ [2003,3,29,23,0,0],[2003,3,30,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2003,10,25,22,59,59],[2003,10,26,2,59,59],
          '2003032923:00:00','2003033003:00:00','2003102522:59:59','2003102602:59:59' ],
        [ [2003,10,25,23,0,0],[2003,10,26,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2004,3,27,22,59,59],[2004,3,28,1,59,59],
          '2003102523:00:00','2003102602:00:00','2004032722:59:59','2004032801:59:59' ],
     ],
   2004 =>
     [
        [ [2004,3,27,23,0,0],[2004,3,28,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2004,10,30,22,59,59],[2004,10,31,2,59,59],
          '2004032723:00:00','2004032803:00:00','2004103022:59:59','2004103102:59:59' ],
        [ [2004,10,30,23,0,0],[2004,10,31,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2005,3,26,22,59,59],[2005,3,27,1,59,59],
          '2004103023:00:00','2004103102:00:00','2005032622:59:59','2005032701:59:59' ],
     ],
   2005 =>
     [
        [ [2005,3,26,23,0,0],[2005,3,27,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2005,10,29,22,59,59],[2005,10,30,2,59,59],
          '2005032623:00:00','2005032703:00:00','2005102922:59:59','2005103002:59:59' ],
        [ [2005,10,29,23,0,0],[2005,10,30,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2006,3,25,22,59,59],[2006,3,26,1,59,59],
          '2005102923:00:00','2005103002:00:00','2006032522:59:59','2006032601:59:59' ],
     ],
   2006 =>
     [
        [ [2006,3,25,23,0,0],[2006,3,26,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2006,10,28,22,59,59],[2006,10,29,2,59,59],
          '2006032523:00:00','2006032603:00:00','2006102822:59:59','2006102902:59:59' ],
        [ [2006,10,28,23,0,0],[2006,10,29,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2007,3,24,22,59,59],[2007,3,25,1,59,59],
          '2006102823:00:00','2006102902:00:00','2007032422:59:59','2007032501:59:59' ],
     ],
   2007 =>
     [
        [ [2007,3,24,23,0,0],[2007,3,25,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2007,10,27,22,59,59],[2007,10,28,2,59,59],
          '2007032423:00:00','2007032503:00:00','2007102722:59:59','2007102802:59:59' ],
        [ [2007,10,27,23,0,0],[2007,10,28,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2008,3,29,22,59,59],[2008,3,30,1,59,59],
          '2007102723:00:00','2007102802:00:00','2008032922:59:59','2008033001:59:59' ],
     ],
   2008 =>
     [
        [ [2008,3,29,23,0,0],[2008,3,30,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2008,10,25,22,59,59],[2008,10,26,2,59,59],
          '2008032923:00:00','2008033003:00:00','2008102522:59:59','2008102602:59:59' ],
        [ [2008,10,25,23,0,0],[2008,10,26,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2009,3,28,22,59,59],[2009,3,29,1,59,59],
          '2008102523:00:00','2008102602:00:00','2009032822:59:59','2009032901:59:59' ],
     ],
   2009 =>
     [
        [ [2009,3,28,23,0,0],[2009,3,29,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2009,10,24,22,59,59],[2009,10,25,2,59,59],
          '2009032823:00:00','2009032903:00:00','2009102422:59:59','2009102502:59:59' ],
        [ [2009,10,24,23,0,0],[2009,10,25,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2010,3,27,22,59,59],[2010,3,28,1,59,59],
          '2009102423:00:00','2009102502:00:00','2010032722:59:59','2010032801:59:59' ],
     ],
   2010 =>
     [
        [ [2010,3,27,23,0,0],[2010,3,28,3,0,0],'+04:00:00',[4,0,0],
          '+04',1,[2010,10,30,22,59,59],[2010,10,31,2,59,59],
          '2010032723:00:00','2010032803:00:00','2010103022:59:59','2010103102:59:59' ],
        [ [2010,10,30,23,0,0],[2010,10,31,2,0,0],'+03:00:00',[3,0,0],
          '+03',0,[2011,3,26,22,59,59],[2011,3,27,1,59,59],
          '2010103023:00:00','2010103102:00:00','2011032622:59:59','2011032701:59:59' ],
     ],
   2011 =>
     [
        [ [2011,3,26,23,0,0],[2011,3,27,3,0,0],'+04:00:00',[4,0,0],
          '+04',0,[2014,10,25,21,59,59],[2014,10,26,1,59,59],
          '2011032623:00:00','2011032703:00:00','2014102521:59:59','2014102601:59:59' ],
     ],
   2014 =>
     [
        [ [2014,10,25,22,0,0],[2014,10,26,1,0,0],'+03:00:00',[3,0,0],
          '+03',0,[9999,12,31,0,0,0],[9999,12,31,3,0,0],
          '2014102522:00:00','2014102601:00:00','9999123100:00:00','9999123103:00:00' ],
     ],
);

%LastRule      = (
);

1;
