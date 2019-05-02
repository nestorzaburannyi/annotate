c------------------------------------------------------------------------------
c Copyright (c) 2012, liisa holm
c All rights reserved.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice, this
c    list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyright notice,
c    this list of conditions and the following disclaimer in the documentation
c    and/or other materials provided with the distribution.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
c ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
c WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
c DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
c ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
c (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
c LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
c ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
c (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c------------------------------------------------------------------------------

c       Ussage: sais8 uniprot. 
c       in: uniprot.pin, uniprot.psq
c       out: uniprot.SA [8 bytes]
c
c       gcc -c sais8.c
c       gfortran -O3 sais8.f sais8.o -o sais8
c
        program sais8
        implicit none
        integer, parameter :: long = selected_int_kind(range(1)*2)
        character(len=356) sname,fn_SA,fn_pin,fn_psq
        integer*8 spos,naa,p,i,offset,k
        integer ok,irec,l,numarg
        character, dimension(:), allocatable :: T
        integer*8, dimension(:), allocatable :: SA
        integer, parameter :: bufsize=2000000000,bufsizeout=bufsize/8

        numarg=iargc()
        if(numarg.lt.1) stop 'USAGE: sais8 database-name'
        call getarg(1,sname)
        fn_psq=trim(sname)//'.psq'
        fn_pin=trim(sname)//'.pin'
        fn_SA=trim(sname)//'.SA'
        ! get naa
        open(10,file=fn_pin,status='old')
10      read(10,*,end=19) p,l
        goto 10
19      close(10)
        naa=p+l
        write(*,*) '# naa =',naa
        
        allocate(T(naa),stat=ok)
        if (ok.ne.0) stop "Memory allocation failed"
        allocate(SA(naa),stat=ok)
        if (ok.ne.0) stop "Memory allocation failed"
        call openfile(10,fn_psq,1*bufsize)
        offset=0
        irec=0
        do while(offset.lt.naa)
                irec=irec+1
                k=min(bufsize,naa-offset)
                write(*,*) '# reading',irec,offset,naa,bufsize,k
                read(10,rec=irec) (T(i),i=offset+1,offset+k)
                offset=offset+k
        end do
        close(10)  
        write(*,*) '#T done',naa,offset
        call sais_from_fortran(T,SA,naa)
        write(*,*) '#SA done',naa
        ! output SA
        call openfile(10,fn_SA,8*bufsizeout)
        offset=0
        irec=0
        do while(offset.lt.naa)
                irec=irec+1
                k=min(bufsizeout,naa-offset)
                write(*,*) '# writing',irec,offset,naa,bufsizeout,k
                write(10,rec=irec) (SA(i),i=offset+1,offset+k)

                ! output 10 checkpoints
                do i=offset+1,offset+k,k/10
                  spos=SA(i)
                  write(*,*) '#SA',spos,(T(p),p=spos+1,min(naa,spos+20))
                end do

                offset=offset+k
        end do
        close(10)
        write(*,*) '# output done ',trim(fn_SA),naa,offset

        end program sais8
c
c======================================================================
c
        subroutine openfile(unit10,fn,recordlength)
c
c       open binary file for direct access
c
        implicit none
        integer unit10,recordlength
        character*256 fn

        open(unit=unit10,
     $          file=fn,
     $          form='unformatted',
     $          access='direct',
     $          recl=recordlength)

        return
        end
c
c====================================================================
c

