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
c       in: unirot.pin uniprot.SA [8 bytes] 
c       out: uniprot.ISA [8 bytes]
c
c       USAGE: kb7_30 uniprot
c
c       algo: ISA[SA[lex]]=lex
c
c       gfortran -O3 isa8.f -o isa8
c
        program kb7_30
        implicit none
        INTEGER, PARAMETER :: long = SELECTED_INT_KIND(18)
        INTEGER, PARAMETER :: unit10=10
        INTEGER, PARAMETER :: giga=1000000000 ! 1 G
        INTEGER(long), PARAMETER :: MEMORY=500_long*giga ! 500 G 
        integer(long) :: lex,naa,offset,k,i,seq,ptr,seqlen
        integer(long), dimension(:), allocatable :: isa
        integer bufsize,ibuf,irec,srec,ok,jrec
        parameter(bufsize=50000000) ! 250 M * 8 bytes = 2 Gbytes
        character*256 sname,fn_SA,fn_ISA,fn_pin
        integer(long) :: sabuffer(bufsize),getsize,chunksize
        integer(long) :: seq_from,seq_to,outbuffer(bufsize)
        integer outbuf,divisor

        ! mandatory arguments
        if(iargc().lt.1) then
                write(*,*) 'isa8 generates inverse suffix array'
                write(*,*) ' '
                write(*,*) 'usage: isa8 database-name '
                stop
        end if
        call getarg(1,sname)
        fn_pin=trim(sname)//'.pin'
        fn_SA=trim(sname)//'.SA'
        fn_ISA=trim(sname)//'.ISA'

        ! get naa from uniprot.pin
        open(unit10,file=fn_pin,status='old')
10      read(unit10,*,end=19) ptr,seqlen
        goto10
19      continue
        naa=ptr+seqlen+1
        close(unit10)

        ! open input SA
        call openfile(10,fn_SA,8*bufsize)

        ! open output ISA
        call openfile(20,fn_ISA,8*bufsize)

        ! split database to fit MEMORY
        ! isa(chunksize)+outbuffer(bufsize)+sabuffer(bufsize) .lt. MEMORY
        chunksize=naa
        divisor=1
        do while(chunksize*8+bufsize*8*2.gt.MEMORY) 
                divisor=divisor+1
                chunksize=1+naa/divisor
        end do
        write(*,*) '# chunksize,naa,chunks',chunksize,naa,divisor,MEMORY
        if(chunksize*divisor.lt.naa) stop 'chunksize bug'
 
        ! dynamic memory allocation
        allocate(isa(chunksize), stat=ok)
        if(ok /= 0) stop 'isa memory allocation failed'

        ! scan chunks
        seq_from=0
        seq_to=chunksize
        irec=0 ! output
        outbuf=0
        do while(seq_from.lt.naa) 
          write(*,*) '#processing chunk',seq_from,seq_to,naa
          write(*,*) '# lex / naa',lex,naa

          ! scan SA(lex), fill ISA[SA[lex]]=lex
          ! assign isa[seq_from:seq_to]
          ibuf=bufsize
          jrec=0 ! input
          lex=0
          do while(lex.lt.naa)
                lex=lex+1
                ibuf=ibuf+1
                if(ibuf.gt.bufsize) then
                        jrec=jrec+1
                        k=min(bufsize,naa-lex)
                        write(*,*) '# to read',jrec,k
                        read(10,rec=jrec) (sabuffer(i),i=1,k)
                        ibuf=1
                end if
                if(mod(lex,100000000).eq.0) then
                        write(*,*) '#processing',lex,float(lex)/naa
                end if
                ! chunk coordinates
                seq=1+sabuffer(ibuf)-seq_from ! C to F numbering
                if(seq.gt.0.and.seq.le.chunksize) then
                        isa(seq)=lex
                end if
          end do

          ! buffered output ISA: flush full records
          do i=1,chunksize
                if(outbuf.ge.bufsize) then ! flush
                        irec=irec+1
                        write(20,rec=irec) (outbuffer(k),k=1,bufsize)
                        outbuf=0
                end if
                outbuf=outbuf+1
                outbuffer(outbuf)=isa(i)
          end do
          seq_from=seq_from+chunksize
          seq_to=seq_to+chunksize
        end do

        ! free memory
        deallocate(isa)

        ! output last record
        if(outbuf.gt.0) then
                irec=irec+1
                write(20,rec=irec) (outbuffer(k),k=1,outbuf)
        end if

        close(10)
        close(20)

        end program kb7_30

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
c
c====================================================================
c
        function getsize(fn)
        implicit none
        character*(*) fn
        integer*8 getsize
        integer l,i
        character*256 cmd

        ! use system call to get size of file in bytes
        l=len_trim(fn)
        cmd(1:6)='ls -l '
        cmd(7:l+6)=fn
        cmd(l+7:l+51)=
     $          " | perl -pe 's/\ +/\t/g' | cut -f 5 > tmp_sans"
        do i=l+52,256
                cmd(i:i)=' '
        end do

        write(*,*) '# call system ',cmd
        call system(cmd(1:len_trim(cmd)))
        open(unit=10,file='tmp_sans',status='old')
        read(10,*) getsize
        close(10)

        return
        end
c
c====================================================================
c
