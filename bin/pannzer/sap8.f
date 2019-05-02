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
c       in: uniprot.ISA [8 bytes], uniprot.pin
c       out: uniprot.SAP [4 bytes] uniprot.SRES [2 bytes]
c
c       USAGE: kb7_31 uniprot.pin uniprot.ISA uniprot.SAP uniprot.sres
c
c       algo: SAP[ISA[seq]]=prot(seq)
c
c       gfortran -O3 sap8.f -o sap8
c
        program kb7_31
        implicit none
        INTEGER, PARAMETER :: long = SELECTED_INT_KIND(18)
        INTEGER, PARAMETER :: giga=1000000000 ! 1 G
        INTEGER(long), PARAMETER :: MEMORY=500_long*giga ! 500 G
        integer(long) :: seq,naa,offset,ptr,k,i,chunksize
        integer bufsize,ibuf,irec,seqlen,jrec,nprot
        integer, dimension(:), allocatable :: sap,sres
        parameter(bufsize=250000000) ! 250 M = 2 Gbytes
        character*256 sname,fn_pin,fn_ISA, fn_SAP,fn_sres
        integer(long) :: isabuffer(bufsize),lex_from,lex_to
        integer iprot,ires,ok,outbuf,outbuffer(bufsize),divisor
        integer*2 outbuffer_sres(bufsize)

        ! mandatory arguments
        if(iargc().lt.1) then
          write(*,*) 'sap8 generates protein labels for suffix array'
          write(*,*) ' '
          write(*,*) 'usage: sap8 database-name '
          stop
        end if
        call getarg(1,sname)
        fn_pin=trim(sname)//'.pin'
        fn_ISA=trim(sname)//'.ISA'
        fn_SAP=trim(sname)//'.SAP'
        fn_sres=trim(sname)//'.SRES'
        write(*,*) '# inputs: ',fn_pin,fn_ISA
        write(*,*) '# outputs: ',fn_SAP,fn_sres
        naa=0 

        ! open input ISA
        call openfile(10,fn_ISA,8*bufsize)

        ! open output SAP
        call openfile(20,fn_SAP,4*bufsize)
        ! open output sres
        call openfile(40,fn_sres,2*bufsize)

        ! open input pin
        open(30,file=fn_pin,status='old')

        ! read uniprot.pin, fill SAP[ISA[seq]]=prot(seq)
        ! get naa
        nprot=0
10      read(30,*,end=19) ptr,seqlen 
        nprot=nprot+1
        goto 10
19      continue
        naa=ptr+seqlen+1
        close(30)

        ! split database to fit MEMORY
        ! sap(chunksize)+isabuffer(bufsize)+outbuffer(bufsize).lt.MEMORY
        chunksize=naa
        divisor=1
        do while(chunksize*8_long+bufsize*14_long.gt.MEMORY) 
                divisor=divisor+1
                chunksize=1+naa/divisor
        end do
        write(*,*) '# chunksize,naa,chunks',chunksize,naa,divisor
        if(chunksize*divisor.lt.naa) stop 'chunksize bug'
                
        ! dynamic memory allocation
        allocate(sap(chunksize),stat=ok)
        if(ok /= 0) stop 'sap memory allocation failed'
        allocate(sres(chunksize),stat=ok)
        if(ok /= 0) stop 'sres memory allocation failed'

        ! scan chunks
        lex_from=0
        lex_to=chunksize
        irec=0 ! output
        outbuf=0
	offset=0
        do while(lex_from.lt.naa)
          write(*,*) '#processing chunk',lex_from,lex_to,naa
          ! rewind
          open(30,file=fn_pin,status='old')
          jrec=0 ! input
          ibuf=bufsize
          iprot=0
          do iprot=1,nprot
20          read(30,*,end=29) ptr,seqlen
            seq=ptr
            !write(*,*) '#processing',iprot,ptr,seqlen
            do ires=1,seqlen+1 ! aaseq plus newline
                seq=seq+1
                ibuf=ibuf+1
                if(ibuf.gt.bufsize) then
                        jrec=jrec+1
                        k=min(bufsize,naa-offset)
	  write(*,*) '#to read',jrec,offset,k,offset+k,naa,iprot
                        read(10,rec=jrec,err=29) (isabuffer(i),i=1,k)
	                offset=offset+k
                        ibuf=1
                end if
                i=isabuffer(ibuf)-lex_from
                if(i.gt.0.and.i.le.chunksize) then
                        sap(i)=iprot
                        sres(i)=ires
          !if(iprot.eq.1) write(*,*) '#sap',iprot,ires,seq,i,i+lex_from
                end if
            end do
            if(mod(iprot,100000).eq.0) 
     $          write(*,*) iprot,naa,ptr,seqlen,seq,float(seq)/naa
          end do
29        continue
          close(30)

          ! buffered output SAP; flush full records
          do i=1,chunksize
                if(outbuf.ge.bufsize) then ! flush
                  irec=irec+1
                  write(20,rec=irec) (outbuffer(k),k=1,bufsize)
                  write(40,rec=irec) (outbuffer_sres(k),k=1,bufsize)
                  outbuf=0
                end if
                outbuf=outbuf+1
                outbuffer(outbuf)=sap(i)
                outbuffer_sres(outbuf)=sres(i) 
          end do
          lex_from=lex_from+chunksize
          lex_to=lex_to+chunksize
        end do

        deallocate(sap)

        ! output last record
        if(outbuf.gt.0) then
                irec=irec+1
                write(20,rec=irec) (outbuffer(k),k=1,outbuf)
                write(40,rec=irec) (outbuffer_sres(k),k=1,outbuf)
        end if

        close(10)
        close(20)
        close(40)

        end program kb7_31

c
c====================================================================
c
        subroutine openfile(unit10,fn,recordlength)
c
c       open binary file for direct access
c
        implicit none
        integer unit10,recordlength
        character*(*) fn

        open(unit=unit10,
     $          file=fn,
     $          form='unformatted',
     $          access='direct',
     $          recl=recordlength)

        return
        end
c
c======================================================================
c

