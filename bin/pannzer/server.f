c------------------------------------------------------------------------------
c Copyright (c) 2015, Liisa Holm
c All rights reserved.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions aremet:
c
c 1. Redistributions of source code must retain the above copyrightnotice, this
c    list of conditions and the following disclaimer.
c 2. Redistributions in binary form must reproduce the above copyrightnotice,
c    this list of conditions and the following disclaimer in thedocumentation
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
c
c
c=======================================================================
c
        module maxheap
        implicit none

        private
        public popheap,pushheap,popheap8,pushheap8,testheap

        ! purpose: hash (key,value) pairs, return key with maximum value
        !
        ! extrinsic variables:
        ! integer n = number of keys
        ! integer d(m) = indices to values; dimension m>=n
        ! integer|real v(m) = values; dimension m>=n
        !
        ! public methods:
        ! (a) integer keys, integer values
        ! popheap(n,d,v,maxkey,maxvalue)
        !       - removes key with maximum value from heap
        !       - result: maxkey, maxvalue
        ! pushheap(key,val,n,d,v)
        !       - add (key,val) to heap
        ! (b) integer keys, real values
        ! popheapr(), pushheapr()
        ! (c) integer*8 keys, integer values
        ! popheap8(), pushheap8()

        contains
c
c-----------------------------------------------------------------------
c
        subroutine popheap8(n,d,v,maxitem)
        implicit none
        integer n,maxitem,v(maxitem),i,j
        integer*8 d(maxitem)

        ! remove key with maximum value from heap
        call swap8(1,n,n,d,maxitem)
        n=n-1
        ! repair heap
        i=1
        do while(i.lt.n)
                j=i*2
                if(j.gt.n) exit
                if(j+1.le.n.and.v(d(j+1)).gt.v(d(j))) j=j+1
                if(v(d(i)).ge.v(d(j))) exit ! heap ok
                call swap8(j,i,n,d,maxitem)
                i=j
        end do

        end subroutine popheap8
c
c-----------------------------------------------------------------------
c
        subroutine pushheap8(val,n,d,v,maxitem)
        implicit none
        integer val,maxitem,n,v(maxitem),p,i
        integer*8 d(maxitem)
        n=n+1
        d(n)=n
        v(n)=val
        ! repair heap
        i=n
        p=i/2 ! integer division
        do while(p.ge.1.and.v(d(p)).lt.v(d(i)))
                call swap8(i,p,n,d,maxitem)
                i=p
                p=i/2
        end do

        end subroutine pushheap8

c
c-----------------------------------------------------------------------
c
        subroutine swap8(i,j,n,d,maxitem)
        implicit none
        integer i,j,n,maxitem
        integer*8 d(maxitem),x

        x=d(i)
        d(i)=d(j)
        d(j)=x

        end subroutine swap8
c
c-----------------------------------------------------------------------
c
        subroutine popheap(n,d,v,maxitem)
        implicit none
        integer n,maxitem,d(maxitem),v(maxitem),maxvalue,i,j

        ! remove key with maximum value from heap
        call swap(1,n,n,d)
        n=n-1
        ! repair heap
        i=1
        do while(i.lt.n)
                j=i*2
                if(j.gt.n) exit
                if(j+1.le.n.and.v(d(j+1)).gt.v(d(j))) j=j+1
                if(v(d(i)).ge.v(d(j))) exit ! heap ok
                call swap(j,i,n,d) 
                i=j
        end do

        end subroutine popheap
c
c-----------------------------------------------------------------------
c
        subroutine pushheap(val,n,d,v,maxitem)
        implicit none
        integer val,n,maxitem,d(maxitem),v(maxitem),p,i

        n=n+1
        d(n)=n
        v(n)=val
        if(n.eq.1) return
        ! repair heap
        i=n
        p=i/2 ! integer division
        do while(v(d(p)).lt.v(d(i)))
                call swap(i,p,n,d)
                i=p
                p=i/2
                if(p.lt.1) exit ! loop
        end do

        end subroutine pushheap
c
c-----------------------------------------------------------------------
c
        subroutine swap(i,j,n,d)
        implicit none
        integer i,j,n,d(n),x

        x=d(i)
        d(i)=d(j)
        d(j)=x

        end subroutine swap
c
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c
        subroutine testheap()
        implicit none
        integer, parameter :: m=100
        integer input(m),output(m),n,d(m),v(m),i,j,val

        n=0
        do i=1,m
                val=int(m*rand())
                input(i)=val
                call pushheap(val,n,d,v,m)
                !write(*,*) '#push',i,d(i),n,(v(d(j)),j=1,n)
        end do
        do i=1,m
                output(i)=v(d(1))
                call popheap(n,d,v,m)
                !write(*,*) '#pop',i,d(i),n,(v(d(j)),j=1,n)
                !write(*,*) '#n,rank,value',n,i,val
        end do
        write(*,*) '#input',input
        write(*,*) '#output',output

        end subroutine testheap

        end module maxheap
c
c=======================================================================
c
!        program test
!        use maxheap
!        call testheap()
!        end program test
c
c====================================================================
c
        module altschul
        ! evalue (BLOSUM62/gapopen -10/gapelon -2)
        !real, parameter ::
        !lambda=0.316,lnkappa=log(0.133),entropy=0.403
        real, parameter :: lambda=0.291,lnkappa=log(0.075),entropy=0.403
        real, parameter :: MAXEVALUE=1.0
        integer hsplen, neff
        integer*8 meff,neffmeff

        private

        public get_bitscore, get_evalue, get_bitscorecutoff
        public prepare_altschul, getblosum_inline
        public mapASCII_NUC, getblosum_inline_NUC

        contains

c
c====================================================================
c
        subroutine prepare_altschul(qseqlen,naa_tot,nprot_tot)
        implicit none
        integer qseqlen,nprot_tot
        integer*8 naa_tot

        hsplen=lnkappa+log(qseqlen+0.1)+log(naa_tot/1000.0)+
     $          log(1000.0/entropy)
        neff=max(1,qseqlen-hsplen)
        meff=naa_tot-hsplen*nprot_tot
        neffmeff=neff*meff

        end subroutine prepare_altschul
c
c====================================================================
c
        real function get_bitscorecutoff(evaluecutoff) result(bitscore)
        implicit none
        real evaluecutoff

        bitscore=log(float(neff))+log(float(meff))-log(evaluecutoff)

        end function get_bitscorecutoff
c
c====================================================================
c
        real function get_bitscore(blosumscore) result(bitscore)
        implicit none
        integer blosumscore

        bitscore=lambda*blosumscore-lnkappa ! in nats

        end function get_bitscore
c
c====================================================================
c
        double precision function get_evalue(bitscore) result(evalue)
        implicit none
        real bitscore

        evalue=neffmeff*exp(-bitscore)

        end function get_evalue
c
c====================================================================
c
        subroutine getblosum_inline(blosum,alfmap)
        implicit none
        integer blosum(0:20,0:20),alfmap(128)
        character x,alfabet(0:20)
        integer i,j,k,l,s(23)
        character(len=62) b(22)
        data b /
     $  "# Entries for the BLOSUM62 matrix at a scale of ln(2)/2.0.",
     $  "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V",
     $  "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0",
     $  "R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3",
     $  "N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3",
     $  "D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3",
     $  "C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1",
     $  "Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2",
     $  "E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2",
     $  "G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3",
     $  "H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3",
     $  "I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3",
     $  "L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1",
     $  "K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2",
     $  "M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1",
     $  "F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1",
     $  "P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2",
     $  "S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2",
     $  "T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0",
     $  "W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3",
     $  "Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1",
     $  "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4"/

        do i=0,20
                do j=0,20
                        blosum(i,j)=-4
                end do
        end do
        blosum(0,0)=1
        alfabet(0)='X'
        read(b(2),*) (alfabet(i),i=1,20)
        !A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
        do i=1,20
                k=alfmap(iachar(alfabet(i)))
                read(b(i+2),*) x,(s(j),j=1,20)
                do j=1,20
                        l=alfmap(iachar(alfabet(j)))
                        blosum(k,l)=s(j)
                end do
        end do

        return
        end subroutine getblosum_inline
c
c====================================================================
c
        subroutine mapASCII_NUC(alfmap)
c
c       alphabet size is 21 (range 0:20), \n is 0
c       valid characters ATGCSWRYKMBVHDN
c       all other characters mapped to 0
c
        implicit none
        integer i,alfmap(128)

        write(*,*) '# This is mapASCII_NUC'
        do i=1,128
                alfmap(i)=0
        end do
        alfmap(65)=1 ! A = 1
        alfmap(84)=2 ! T = 2
        alfmap(71)=3 ! G = 3
        alfmap(67)=4 ! C = 4
        alfmap(83)=5 ! S = 5
        alfmap(87)=6 ! W = 6
        alfmap(82)=7 ! R = 7
        alfmap(89)=8 ! Y = 8
        alfmap(75)=9 ! K = 9
        alfmap(77)=10 ! M = 10
        alfmap(66)=11 ! B = 11
        alfmap(86)=12 ! V = 12
        alfmap(72)=13 ! H = 13
        alfmap(68)=14 ! D = 14
        alfmap(79)=15 ! N = 15
        write(*,*) alfmap

        return
        end subroutine mapASCII_NUC
c
c====================================================================
c
        subroutine getblosum_inline_NUC(blosum,alfmap)
        implicit none
        integer blosum(0:20,0:20),alfmap(128)
        character x,alfabet(0:20)
        integer i,j,k,l,s(23)
        character(len=62) b(17)
        data b /
     $  "# This matrix was created by Todd Lowe   12/10/92",
     $  "    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N",
     $  "A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2",
     $  "T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2",
     $  "G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2",
     $  "C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2",
     $  "S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1",
     $  "W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1",
     $  "R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1",
     $  "Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1",
     $  "K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1",
     $  "M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1",
     $  "B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1",
     $  "V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1",
     $  "H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1",
     $  "D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1",
     $  "N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1"/

        write(*,*) '# This is getblosum_inline_NUC'
        do i=0,20
                do j=0,20
                        blosum(i,j)=-4
                end do
        end do
        blosum(0,0)=0
        alfabet(0)='X'
        read(b(2),*) (alfabet(i),i=1,15)
        !A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
        do i=1,15
                k=alfmap(iachar(alfabet(i)))
                read(b(i+2),*) x,(s(j),j=1,15)
                do j=1,15
                        l=alfmap(iachar(alfabet(j)))
                        blosum(k,l)=s(j)
                end do
        end do
        write(*,*) '# done'

        return
        end subroutine getblosum_inline_NUC
c
c======================================================================
c
        end module altschul
c
c====================================================================
c
c
c====================================================================
c
        subroutine blosum2pssm(qtxt,qseqlen,blosum,PSSM,alfmap)
        implicit none
        integer qseqlen,blosum(0:20,0:20),PSSM(0:20,qseqlen)
        integer alfmap(128)
        integer*1 qtxt(qseqlen)
        integer i,j,k

        do i=1,qseqlen
                k=alfmap(qtxt(i))
                do j=0,20
                        PSSM(j,i)=blosum(j,k)
                end do
                !write(*,*) '#2PSSM ',char(qtxt(i)),i,(PSSM(j,i),j=0,20)
        end do

        end subroutine blosum2pssm
c
c====================================================================
c
        subroutine mapASCII(alfmap)
c
c       alphabet size is 21 (range 0:20), \n is 0
c       valid characters ACDEFGHIKLMNPQRSTVWY
c       all other characters mapped to 0
c
        implicit none
        integer i,alfmap(128)

        do i=1,128
                alfmap(i)=0
        end do
        alfmap(65)=1 ! A = 1
        do i=67,73 ! CDEFGHI = 2-8
                alfmap(i)=i-65
        end do
        do i=75,78 ! KLMN = 9-12
                alfmap(i)=i-66
        end do
        do i=80,84 ! PQRST = 13-17
                alfmap(i)=i-67
        end do
        alfmap(86)=18 ! V
        alfmap(87)=19 ! W
        alfmap(89)=20 ! Y

        return
        end subroutine mapASCII
c
c====================================================================
c
        subroutine getblosum(blosum,alfmap,unit10,infile)
        implicit none
        integer blosum(0:20,0:20),alfmap(128),unit10
        character*(*) infile
        character x,alfabet(0:20)
        integer i,j,k,l,s(23)

        do i=0,20
                do j=0,20
                        blosum(i,j)=-4
                end do
        end do
        blosum(0,0)=1
        open(unit10,file=infile,status='old')
        read(unit10,*) x ! skip header
        alfabet(0)='X'
        read(unit10,*) (alfabet(i),i=1,23)
        !A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
        do i=1,23
                k=alfmap(iachar(alfabet(i)))
                read(10,*) x,(s(j),j=1,23)
                do j=1,23
                        l=alfmap(iachar(alfabet(j)))
                        blosum(k,l)=s(j)
                end do
        end do
        close(unit10)

        return
        end
c
c======================================================================
c
c====================================================================
c
        subroutine real_Qsort(i,j,n,d,v)
c
c       non-recursive Quicksort
c       sorts integer array v which has n elements from smallest to largest
c       d contains indices to v
c
        implicit none
        integer i,j,n,d(n)
        real v(n)
        integer stacksize
        parameter(stacksize=100)
c
        integer p,top,bottom,stack(stacksize), nstack, realpartition

        stack(1)=j
        stack(2)=i
        nstack=2
        do while(nstack.ne.0)
                top=stack(nstack)
                bottom=stack(nstack-1)
                nstack=nstack-2
                do while(top.lt.bottom)
                        p=realpartition(top, bottom, n, d, v)
                        if((p-top).gt.(bottom-p)) then
                                stack(nstack+1)=p-1
                                stack(nstack+2)=top
                                top=p+1
                                nstack=nstack+2
                                if(nstack.gt.stacksize)
     $  stop ' FATAL ERROR: Stack overflow in Qsort '
                        else
                                stack(nstack+1)=bottom
                                stack(nstack+2)=p+1
                                bottom=p-1
                                nstack=nstack+2
                                if(nstack.gt.stacksize)
     $  stop ' FATAL ERROR: Stack overflow in Qsort '
                        end if
                end do
        end do

        return
        end

c ======================================================================
        function realpartition(i,j,n,d,v)
        implicit none
        integer i,j,n,realpartition, d(n)
        real v(n)
c
        integer upper, lower, saved

        upper=i
        lower=j
        saved=d(i)

        do while (upper.ne.lower)
                do while((upper.lt.lower).and.(v(saved).le.v(d(lower))))
                        lower=lower-1
                end do
                if(upper.ne.lower) d(upper)=d(lower)
                do while((upper.lt.lower).and.(v(saved).ge.v(d(upper))))
                        upper=upper+1
                end do
                if(upper.ne.lower) d(lower)=d(upper)
        end do
        d(upper)=saved
        realpartition=upper

        return
        end
c
c====================================================================
c

c====================================================================
c
        subroutine j_index_Qsort(i,j,n,d,v)
c
c       non-recursive Quicksort
c       sorts integer array v which has n elements from smallest to largest
c       d contains indices to v
c
        implicit none
        integer i,j,n,d(n)
        integer v(n)
        integer stacksize
        parameter(stacksize=100)
c
        integer p,top,bottom,stack(stacksize), nstack, partition

        stack(1)=j
        stack(2)=i
        nstack=2
        do while(nstack.ne.0)
                top=stack(nstack)
                bottom=stack(nstack-1)
                nstack=nstack-2
                do while(top.lt.bottom)
                        p=partition(top, bottom, n, d, v)
                        if((p-top).gt.(bottom-p)) then
                                stack(nstack+1)=p-1
                                stack(nstack+2)=top
                                top=p+1
                                nstack=nstack+2
                                if(nstack.gt.stacksize)
     $  stop ' FATAL ERROR: Stack overflow in Qsort '
                        else
                                stack(nstack+1)=bottom
                                stack(nstack+2)=p+1
                                bottom=p-1
                                nstack=nstack+2
                                if(nstack.gt.stacksize)
     $  stop ' FATAL ERROR: Stack overflow in Qsort '
                        end if
                end do
        end do

        return
        end

c ======================================================================
        function partition(i,j,n,d,v)
        implicit none
        integer i,j,n,partition, d(n)
        integer v(n)
c
        integer upper, lower, saved

        upper=i
        lower=j
        saved=d(i)

        do while (upper.ne.lower)
                do while((upper.lt.lower).and.(v(saved).le.v(d(lower))))
                        lower=lower-1
                end do
                if(upper.ne.lower) d(upper)=d(lower)
                do while((upper.lt.lower).and.(v(saved).ge.v(d(upper))))
                        upper=upper+1
                end do
                if(upper.ne.lower) d(lower)=d(upper)
        end do
        d(upper)=saved
        partition=upper

        return
        end
c
c====================================================================
c

c
c======================================================================
c
        module sanslocaldata
        implicit none
        integer, parameter :: MAXPROC=128, MAXRES=64000, CYCMOD=32000
        integer*1,dimension(:),allocatable :: TXT
        integer*1, dimension(:), allocatable :: PHR
        integer, dimension(:), allocatable :: SAP
        integer*8, dimension(:), allocatable :: ptr, ptr_phr
        integer*2, dimension(:), allocatable :: SRES, CYCLIC
        integer*8 naa,naa_tot, nphr, nphr_tot
        integer nprot,nprot_tot,MPI_COMM_SANS

        private
        public load_data,deload_data,scatter_txt,scatter_protlen
        public TXT,PHR,SAP,ptr,ptr_phr,SRES,CYCLIC,naa,nprot,nphr
        public CYCMOD,MPI_COMM_SANS,MAXPROC
        public naa_tot, nprot_tot

        contains
c
c======================================================================
c
        subroutine load_data(sname,rank,nproc)
        implicit none
        integer nproc,rank
        character(len=256) sname
        integer ok

        if(nproc.gt.MAXPROC) stop 'increase MAXPROC'
        ! space allocation: naa, nprot are LOCAL 
        call get_naa_nprot(sname,rank,nproc)
        allocate(TXT(naa+1),stat=ok) ! 1 byte
        allocate(PHR(nphr+1),stat=ok) ! 1 byte
        allocate(SAP(naa+1),stat=ok) ! 4 bytes
        allocate(SRES(naa+1),stat=ok) ! 2 bytes
        allocate(CYCLIC(naa+1),stat=ok) ! 2 bytes
        allocate(ptr(nprot+1),stat=ok) ! 8 bytes
        allocate(ptr_phr(nprot+1),stat=ok) ! 8 bytes

        end subroutine load_data
c
c====================================================================
c
        subroutine deload_data()
        implicit none

        deallocate(TXT)
        deallocate(PHR)
        deallocate(SAP)
        deallocate(SRES)
        deallocate(CYCLIC)
        deallocate(ptr)
        deallocate(ptr_phr)

        end subroutine deload_data
c
c====================================================================
c
        subroutine scatter_txt(sname,rank,nproc)
        implicit none
        include 'mpif.h'
        integer i,j,k,l,m,n,r,ok,irec,ierr,rank,nproc
        integer*1,dimension(:),allocatable::map
        integer*1, dimension(:), allocatable:: tmptxt,sendbuf,recvbuf
        integer*8 p,offset,local_offset
        character(len=256) sname,fn_pin,fn_psq
        integer, parameter :: bufsize=2000000000 ! 2 G
        integer counts(MAXPROC),displs(MAXPROC),wher(MAXPROC)
        integer*8, dimension(:), allocatable :: ptr_tot,ptr_phr_tot,
     $          ptr_buf,ptr_phr_buf

        !write(*,*)'#txt',rank,naa,naa_tot,nprot,nprot_tot,trim(sname)
        ! read ptr(), fill map
        allocate(ptr(nprot+1),stat=ok)
        allocate(ptr_phr(nprot+1),stat=ok)
        if(rank.eq.0) then
          allocate(ptr_tot(nprot_tot+nproc),stat=ok)
          allocate(ptr_phr_tot(nprot_tot+nproc),stat=ok)
          allocate(ptr_buf(nprot_tot+nproc),stat=ok)
          allocate(ptr_phr_buf(nprot_tot+nproc),stat=ok)
          fn_pin=trim(sname)//'.pin'
          open(10,file=fn_pin,status='old')
          do i=1,nprot_tot
                read(10,*,end=19) ptr_tot(i),l,ptr_phr_tot(i),m
          end do
19        close(10)
          ptr_tot(nprot_tot+1)=ptr_tot(nprot_tot)+l
          ptr_phr_tot(nprot_tot+1)=ptr_phr_tot(nprot_tot)+m
          write(*,*) '#ptr_tot',(ptr_tot(i),i=1,10)
          write(*,*) '#ptr_phr_tot',(ptr_phr_tot(i),i=1,10)
          ! reorder ptr_buf, ptr_phr_buf
          ptr_buf=0
          ptr_phr_buf=0
          do r=1,nproc
                j=(r-1)*nprot
                do i=r,nprot_tot,nproc
                  j=j+1
                  if(i.eq.r) then
                    ptr_buf(j)=0
                    ptr_phr_buf(j)=0
                  else
                    ptr_buf(j)=ptr_buf(j-1)+
     $                  ptr_tot(i-nproc+1)-ptr_tot(i-nproc)
                    ptr_phr_buf(j)=ptr_phr_buf(j-1)+
     $                  ptr_phr_tot(i-nproc+1)-ptr_phr_tot(i-nproc)
                  end if
                  !write(*,*) '#i',r,i,j,k,ptr_buf(j),ptr_phr_buf(j)
                end do
          end do
          ! allocate map
          allocate(map(max(naa_tot,nphr_tot)),stat=ok)
          do i=1,nprot_tot
                m=1+mod(i-1,nproc)
                ! map txt
                do p=ptr_tot(i)+1,ptr_tot(i+1)
                        if(p.lt.1) cycle
                        map(p)=m
                end do
          end do
          ! read full text
          allocate(tmptxt(max(naa_tot,nphr_tot)),stat=ok)
          irec=0
          offset=0
          fn_psq=trim(sname)//'.psq'
          call openfile(10,fn_psq,1*bufsize)
          do while(offset.lt.naa_tot)
                irec=irec+1
                k=min(bufsize,naa_tot-offset)
                write(*,*) '#read psq',irec,offset,k,naa_tot
                read(10,rec=irec) (tmptxt(offset+i),i=1,k)
                offset=offset+k
          end do
          close(10)
          ! allocate sendbuf
          allocate(sendbuf(bufsize),stat=ok)
        end if
        ! distribute pointers
        call MPI_SCATTER(ptr_buf,nprot,MPI_INTEGER8,ptr,nprot,
     $          MPI_INTEGER8,0,MPI_COMM_SANS,ierr)
        call MPI_SCATTER(ptr_phr_buf,nprot,MPI_INTEGER8,ptr_phr,nprot,
     $          MPI_INTEGER8,0,MPI_COMM_SANS,ierr)
        if(ptr(nprot).lt.1) ptr(nprot)=ptr(nprot-1)
        ptr(nprot+1)=ptr(nprot)
        ! check
        write(*,*) '#ptr',rank,(ptr(i),i=1,10),
     $          (ptr(i),i=nprot-10,nprot+1)
        write(*,*) '#ptr_phr',rank,(ptr_phr(i),i=1,10),
     $          (ptr_phr(i),i=nprot-10,nprot+1)

        ! reorder sendbuf, scatterv buffer
        allocate(recvbuf(bufsize),stat=ok)
        offset=0
        local_offset=0
        do while(offset.lt.naa_tot)
          !write(*,*) '#offset',rank,offset,naa_tot,local_offset,naa
          if(rank.eq.0) then
                ! get counts,displs
                counts=0
                k=min(bufsize,naa_tot-offset)
                do p=offset+1,offset+k
                        r=map(p)
                        if(r.lt.1) cycle
                        counts(r)=counts(r)+1
                end do
                !write(*,*) '#counts',offset,(counts(r),r=1,nproc)
                displs(1)=0
                do r=2,nproc
                        displs(r)=displs(r-1)+counts(r-1)
                end do
                do r=1,nproc
                        wher(r)=displs(r)
                end do
                ! copy tmptxt to sendbuf
                do p=offset+1,offset+k
                        r=map(p)
                        if(r.lt.1) cycle
                        wher(r)=wher(r)+1
                        sendbuf(wher(r))=tmptxt(p)
                end do
                !write(*,*) '#wher',(wher(r),r=1,nproc)
          end if
          call MPI_SCATTER(counts,1,MPI_INTEGER,n,1,MPI_INTEGER,
     $          0,MPI_COMM_SANS,ierr)
          call MPI_SCATTERV(sendbuf,counts,displs,MPI_INTEGER1,
     $          recvbuf,n,MPI_INTEGER1,0,MPI_COMM_SANS,ierr)
          write(*,*) '#scv',rank,n,local_offset,naa
          ! copy received to TXT
          do i=1,n
                TXT(local_offset+i)=recvbuf(i)
          end do
          local_offset=local_offset+n
          offset=offset+bufsize
        end do
        ! check
        write(*,*) '#TXT',rank,naa,(char(TXT(p)),p=1,100)
        write(*,*) '#...',rank,naa,(char(TXT(p)),p=naa-100,naa)
        if(rank.eq.0) then
          ! map phr
          do i=1,nprot_tot
                m=1+mod(i-1,nproc)
                do p=ptr_phr_tot(i)+1,ptr_phr_tot(i+1)
                        if(p.lt.1) cycle
                        map(p)=m
                end do
          end do
          ! read full phr to tmptxt
          irec=0
          offset=0
          fn_psq=trim(sname)//'.phr'
          call openfile(10,fn_psq,1*bufsize)
          do while(offset.lt.nphr_tot)
                irec=irec+1
                k=min(bufsize,nphr_tot-offset)
                write(*,*) '#read phr',irec,offset,k,nphr_tot
                read(10,rec=irec) (tmptxt(offset+i),i=1,k)
                offset=offset+k
          end do
        end if
        ! reorder sendbuf, scatterv buffer
        offset=0
        local_offset=0
        do while(offset.lt.nphr_tot)
          !write(*,*) '#offset',rank,offset,nphr_tot,local_offset,nphr
          if(rank.eq.0) then
                ! get counts,displs
                counts=0
                k=min(bufsize,nphr_tot-offset)
                do p=offset+1,offset+k
                        r=map(p)
                        if(r.lt.1) cycle
                        counts(r)=counts(r)+1
                end do
                !write(*,*) '#counts',offset,(counts(r),r=1,nproc)
                displs(1)=0
                do r=2,nproc
                        displs(r)=displs(r-1)+counts(r-1)
                end do
                do r=1,nproc
                        wher(r)=displs(r)
                end do
                ! copy tmptxt to sendbuf
                do p=offset+1,offset+k
                        r=map(p)
                        if(r.lt.1) cycle
                        wher(r)=wher(r)+1
                        sendbuf(wher(r))=tmptxt(p)
                end do
                !write(*,*) '#wher',(wher(r),r=1,nproc)
          end if
          call MPI_BARRIER(MPI_COMM_SANS,ierr)
          call MPI_SCATTER(counts,1,MPI_INTEGER,n,1,MPI_INTEGER,
     $          0,MPI_COMM_SANS,ierr)
          call MPI_SCATTERV(sendbuf,counts,displs,MPI_INTEGER1,
     $          recvbuf,n,MPI_INTEGER1,0,MPI_COMM_SANS,ierr)
          write(*,*) '#scv',rank,n,local_offset,nphr
          ! copy received to PHR
          do i=1,n
                PHR(local_offset+i)=recvbuf(i)
          end do
          local_offset=local_offset+n
          offset=offset+bufsize
        end do
        ! check
!       write(*,*) '#PHR',rank,naa,(char(PHR(p)),p=1,100)
!        write(*,*) '#...',rank,naa,(char(PHR(p)),p=nphr-100,nphr)

        ! free space
        if(rank.eq.0) then
                deallocate(tmptxt)
                deallocate(map)
                deallocate(sendbuf)
                deallocate(ptr_tot)
                deallocate(ptr_phr_tot)
        end if
        deallocate(recvbuf)

        end subroutine scatter_txt
c
c====================================================================
c
        subroutine scatter_protlen(sname,rank,nproc)
        implicit none
        include 'mpif.h'
        integer nproc,rank
        character(len=256) fn_pin,sname,fn_psq,fn_SAP,fn_SRES
        integer*8 spos,p,q,offset,offset_tot,lex_local,lex
        integer i,j,k,l,n,r,ok,ierr,displs(MAXPROC),irec,nsend(MAXPROC)
        integer s,counts(MAXPROC),sprot,wher(MAXPROC)
        integer, parameter :: bufsize=500000000 ! 500 M
        integer*2, dimension(:), allocatable :: sresbuf,cycbuf2,sresbuf2
        integer*2, dimension(:), allocatable :: recvbuf2,recvbuf22
        integer, dimension(:), allocatable :: sapbuf,sapbuf2,recvbuf4

        ! SAP(mod sprot),SRES(mod sprot),CYCLIC=mod(lex,64K)
        if(rank.eq.0) then
                irec=0
                fn_SAP=trim(sname)//'.SAP'
                fn_SRES=trim(sname)//'.SRES'
                call openfile(10,fn_SAP,4*bufsize)
                call openfile(20,fn_SRES,2*bufsize)
                allocate(sapbuf(bufsize),stat=ok)
                allocate(sresbuf(bufsize),stat=ok)
                allocate(sapbuf2(bufsize),stat=ok)
                allocate(sresbuf2(bufsize),stat=ok)
                allocate(cycbuf2(bufsize),stat=ok)
        end if
        allocate(recvbuf4(bufsize),stat=ok)
        allocate(recvbuf2(bufsize),stat=ok)
        allocate(recvbuf22(bufsize),stat=ok)
        lex=0
        lex_local=0
        do while(lex.lt.naa_tot)
                write(*,*) '#lex',rank,lex_local,naa,lex,naa_tot
                if(rank.eq.0) then
                        irec=irec+1
                        k=min(bufsize,naa_tot-lex)
                        write(*,*) '#read sap',irec,lex,k,naa_tot
                        read(10,rec=irec) (sapbuf(i),i=1,k)
                        read(20,rec=irec) (sresbuf(i),i=1,k)
                        ! reorder buffers for scatterv
                        ! counts(r)
                        counts=0
                        do i=1,k
                                sprot=sapbuf(i)
                                if(sprot.lt.1) cycle
                                 r=1+mod(sprot-1,nproc)
                               counts(r)=counts(r)+1
                        end do
                        s=0
                        do r=1,nproc
                                s=s+counts(r)
                        end do
                        write(*,*) '#counts sum',s
                        displs(1)=0
                        do r=2,nproc
                                displs(r)=displs(r-1)+counts(r-1)
                        end do
                        ! fill sapbuf2,sresbuf2,cycbuf2
                        do r=1,nproc
                                wher(r)=displs(r)
                        end do
                        sapbuf2=0
                        sresbuf2=0
                        cycbuf2=0
                        do i=1,k
                                sprot=sapbuf(i)
                                if(sprot.lt.1) cycle
                                r=1+mod(sprot-1,nproc)
                                wher(r)=wher(r)+1
                                ! sapbuf2=local_sprot
                                sapbuf2(wher(r))=1+(sprot-1)/nproc
                                sresbuf2(wher(r))=sresbuf(i)
                                cycbuf2(wher(r))=mod(lex+i,CYCMOD)
                        end do
                        ! counts, displs, wher are INTEGER !
                        write(*,*) '#counts',(counts(i),i=1,nproc)
                        write(*,*) '#displs',(displs(i),i=1,nproc)
                        write(*,*) '#wher',(wher(i),i=1,nproc)
                end if
                call MPI_BARRIER(MPI_COMM_SANS,ierr)
                call MPI_SCATTER(counts,1,MPI_INTEGER,n,1,MPI_INTEGER,
     $                  0,MPI_COMM_SANS,ierr)
                !write(*,*) '#scv',rank,n,lex_local,naa,lex,naa_tot
                call MPI_SCATTERV(sapbuf2,counts,displs,MPI_INTEGER,
     $                  recvbuf4,n,MPI_INTEGER,0,MPI_COMM_SANS,ierr)
                call MPI_SCATTERV(sresbuf2,counts,displs,MPI_INTEGER2,
     $                  recvbuf2,n,MPI_INTEGER2,0,MPI_COMM_SANS,ierr)
                call MPI_SCATTERV(cycbuf2,counts,displs,MPI_INTEGER2,
     $                  recvbuf22,n,MPI_INTEGER2,0,MPI_COMM_SANS,ierr)
                ! append recvbuf to SAP,SRES,CYCLIC
                do i=1,n
                        lex_local=lex_local+1
                        if(lex_local.gt.naa) then
                                write(*,*) '#bug',rank,lex_local,naa,
     $  i,recvbuf4(i),recvbuf2(i),recvbuf22(i),n
                                cycle
                        end if
                        SAP(lex_local)=recvbuf4(i)
                        SRES(lex_local)=recvbuf2(i)
                        CYCLIC(lex_local)=recvbuf22(i)
                end do
                write(*,*) '#n',rank,lex_local,n,naa     
                lex=lex+bufsize
        end do
        ! check
        write(*,*) '#SAP',rank,(SAP(p),p=1,10),(SAP(p),p=naa-10,naa)
        write(*,*) '#SRES',rank,(SRES(p),p=1,10),(SRES(p),p=naa-10,naa)
        write(*,*) '#C',rank,(CYCLIC(p),p=1,10),(CYCLIC(p),p=naa-10,naa)

        if(rank.eq.0) then
                close(10)
                close(20)
                deallocate(sapbuf)
                write(*,*) '#sapbuf deallocated'
                deallocate(sresbuf)
                write(*,*) '#sresbuf deallocated'
                deallocate(sapbuf2)
                write(*,*) '#sapbuf2 deallocated'
                deallocate(sresbuf2)
                write(*,*) '#sresbuf2 deallocated'
                deallocate(cycbuf2)
                write(*,*) '#cycbuf2 deallocated'
        end if
        deallocate(recvbuf4)
        write(*,*) '#recvbuf4 deallocated'
        deallocate(recvbuf2)
        write(*,*) '#recvbuf2 deallocated'
        deallocate(recvbuf22)
        write(*,*) '#recvbuf22 deallocated'

        end subroutine scatter_protlen
c
c======================================================================
c
        subroutine get_naa_nprot(sname,rank,nproc)
        implicit none
        include 'mpif.h'
        character(len=256) sname,fn_pin
        integer rank,nproc,ierr,i,l,m
        integer*8 p,q,naa_list(MAXPROC),nphr_list(MAXPROC)

        if(rank.eq.0) then
                fn_pin=trim(sname)//'.pin'
                open(10,file=fn_pin,status='old')
                nprot_tot=0
                nphr_tot=0
                naa_tot=0
                naa_list=0
                nphr_list=0
10              read(10,*,end=19) p,l,q,m
                nprot_tot=nprot_tot+1
                naa_tot=naa_tot+l+1 ! add newline
                nphr_tot=nphr_tot+m+1 ! add newline
                i=mod(nprot_tot,nproc)
                if(i.eq.0) i=nproc
                naa_list(i)=naa_list(i)+l+1
                nphr_list(i)=nphr_list(i)+m+1
                goto 10
19              close(10)
                write(*,*) '# total nprot, naa,nphr',nprot_tot,naa_tot,
     $                  (naa_list(i),i=1,nproc),(nphr_list(i),i=1,nproc)
                nprot=1+nprot_tot/nproc
        end if
        call MPI_BCAST(nprot,1,MPI_INTEGER,0,MPI_COMM_SANS,ierr)
        call MPI_SCATTER(naa_list,1,MPI_INTEGER8,naa,1,MPI_INTEGER8,0,
     $          MPI_COMM_SANS,ierr)
        call MPI_SCATTER(nphr_list,1,MPI_INTEGER8,nphr,1,MPI_INTEGER8,0,
     $          MPI_COMM_SANS,ierr)
        call MPI_BCAST(naa_tot,1,MPI_INTEGER8,0,MPI_COMM_SANS,ierr)
        call MPI_BCAST(nprot_tot,1,MPI_INTEGER,0,MPI_COMM_SANS,ierr)
        call MPI_BCAST(nphr_tot,1,MPI_INTEGER8,0,MPI_COMM_SANS,ierr)
        write(*,*) '# local nprot,naa',rank,nprot,nprot_tot,naa,naa_tot,
     $          nphr,nphr_tot

        end subroutine get_naa_nprot
c
c======================================================================
c
c
c======================================================================
c
        end module sanslocaldata
c
c======================================================================
c
c
c======================================================================
c
        module myhash
        implicit none

        integer, parameter :: maxitem=10000000
        integer ncollision, nkey, shortlist(maxitem),nshort

        contains

        subroutine linkhash(sprot,diag,val,linklist_key,linklist_next,
     $          linklist_value,linklist_first,M,nprot)
        implicit none
        integer sprot,diag,val,M,linklist_key(M),linklist_next(M)
        integer linklist_value(M),nprot,linklist_first(nprot),ix

        if(nkey.ge.maxitem) return
        ! key=sprot; value=diag=sres-qres
        ix=linklist_first(sprot)
        !write(*,*) '#linkhash',sprot,diag,ix,nkey,nshort
        if(ix.eq.0) then ! first key
                nkey=nkey+1
                linklist_first(sprot)=nkey
                linklist_next(nkey)=0
                linklist_value(nkey)=val
                linklist_key(nkey)=diag
                nshort=nshort+1
                shortlist(nshort)=sprot
                return
        end if
        ! new key added to start of chain
        nkey=nkey+1
        linklist_next(nkey)=linklist_first(sprot)
        linklist_first(sprot)=nkey
        linklist_value(nkey)=val
        linklist_key(nkey)=diag
        return

        end subroutine linkhash

        end module myhash
c
c======================================================================
c
c
c====================================================================
c
        module keyspace
        use maxheap, only : pushheap, popheap ! maxheap8.f
        implicit none
        ! keys on heap
        integer base(7)
        data base /1,21,441,9261,194481,4084101,85766121/
        integer,parameter:: MAXKEY4=194481 ! 21^4
        integer,parameter :: MAXKEYSCORE=128
        integer :: keymap(4,MAXKEY4) ! [score,qres,qprot,mutpos]
        ! scoring
        integer,parameter :: maxres=64000
        integer*1 atxt(maxres)
        ! heap
        integer heap_n,heap_input(MAXKEY4),heap_v(MAXKEY4)
        integer heap_d(MAXKEY4)

        private
        public get_keymap, keymap, MAXKEY4, get_ix

        contains
c
c----------------------------------------------------------------------
c
        subroutine get_ix(atxt,qseqlen,ix)
        implicit none
        integer qseqlen,ix(qseqlen)
        integer*1 atxt(qseqlen)
        integer qres

        if(qseqlen.le.2) return
        ix(1)=base(4)*atxt(1)+base(3)*atxt(2)+base(2)*atxt(3)+atxt(4)
        do qres=2,qseqlen-3
                ix(qres)=mod(ix(qres-1),base(4))*base(2)+atxt(qres+3)
        end do
        do qres=qseqlen-2,qseqlen
                ix(qres)=0
        end do

        end subroutine get_ix
c
c----------------------------------------------------------------------
c
        ! call with rank=0, nproc=1 for serial version
        subroutine seedkeys(qprot,kmer,qtxt,qseqlen,PSSM,alfmap)
        implicit none
        integer qprot,kmer,qseqlen,alfmap(128)
        integer qres,score,i,ix(maxres)
        logical ok(maxres)
        integer*1 a,qtxt(qseqlen)
        integer PSSM(0:20,maxres),mutpos

        ! alfmap of qseq, exclude X-words
        do qres=1,qseqlen
                ok(qres)=.true.
                atxt(qres)=alfmap(QTXT(qres)) ! map 128 ASCII characters to 0-20
                if(atxt(qres).lt.1) then
                        do i=max(1,qres-kmer+1),min(qseqlen,qres+kmer-1)
                                ok(qres)=.false.
                        end do
                end if
        end do
        call get_ix(atxt,qseqlen,ix)
        !write(*,*) '#xwords',(ok(i),i=1,qseqlen)
        mutpos=0
        do qres=1,qseqlen-kmer+1
                if(.not.ok(qres)) cycle
                score=0
                do i=1,kmer
                        a=atxt(qres+i-1)
                        score=score+PSSM(a,qres+i-1) ! self
!       write(*,*) '#score',qprot,qres,i,a,PSSM(a,qres+i-1),score
                end do
!                write(*,*) '#qix',qtxt(qres),atxt(qres),score,
!     $                  PSSM(atxt(qres),qres),qres,ix(qres)
                ! push to heap
                call savekey(ix(qres),score,qres,qprot,mutpos)
        end do

        end subroutine seedkeys
c
c----------------------------------------------------------------------
c
        subroutine get_keymap(qprot,scorecutoff,kmer,alfmap,
     $          MPI_COMM_SANS,PSSM,qtxt,qseqlen)
        implicit none
        include 'mpif.h'
        integer offset,x,i,j,k,keyix,aaix,ix,slice,kmer,alfmap(128)
        integer scorecutoff,ierr,MPI_COMM_SANS
        integer PSSM(0:20,maxres)
        integer ixout,scoreout,mutpos,q,qpos,qprot,score,qseqlen
        integer*1 qtxt(qseqlen)

        !write(*,*) '#get_keymap',qprot,scorecutoff,kmer,qseqlen
        ! unmutated query words
        heap_n=0
        call seedkeys(qprot,kmer,qtxt,qseqlen,PSSM,alfmap)
        ! mutated words scoring above scorecutoff
        do while(heap_n.gt.0)
                j=heap_d(1)
                ix=heap_input(j)
                score=heap_v(j) ! keyscore
                !write(*,*) '#pop',qprot,ix,score,heap_n
                call popheap(heap_n,heap_d,heap_v,MAXKEY4)
                ! key data saved in keymap on pushheap()
                ! mutate words
                call pointmutate(ix,scorecutoff,kmer,PSSM,qseqlen,qprot)
        end do
        ! keymap is result

        end subroutine get_keymap
c
c----------------------------------------------------------------------
c
        subroutine pointmutate(ix,scorecutoff,kmer,PSSM,qseqlen,qprot)
        implicit none
        integer ix,ixout,i,j,k,mutpos,aaix,qpos,qprot,score,qseqlen
        integer scoreout,scorecutoff,slice,kmer,PSSM(0:20,qseqlen)

        qprot=keymap(3,ix)
        qpos=keymap(2,ix)
        score=keymap(1,ix)
!        write(*,*) '#pointmutate',ix,scorecutoff,
!     $
!     qprot,qpos,score,keys_mutpos(ix),(atxt(qpos+i),i=1,kmer)
        do mutpos=keymap(4,ix)+1,kmer
                j=qpos+mutpos-1 ! mutated qpos
                k=atxt(j) ! original aaix
                if(k.le.0) cycle
                i=base(kmer-mutpos+1) ! mutpos multiplier
                do aaix=1,20
                      if(aaix.eq.k) cycle ! no self
                      ixout=ix+i*(aaix-k)
                      scoreout=score-PSSM(k,j)+PSSM(aaix,j)
!         write(*,*) '#mut',mutpos,k,aaix,ix,ixout,score,scoreout,
!     $          PSSM(k,qpos+mutpos-1),PSSM(aaix,qpos+mutpos-1)
                      if(scoreout.lt.scorecutoff) cycle
                      if(keymap(3,ixout).eq.qprot.and.
     $                  keymap(1,ixout).ge.scoreout) cycle
                      ! save key if better
                      if(keymap(3,ixout).ne.qprot.or.
     $                  keymap(1,ixout).lt.scoreout) then
                        call savekey(ixout,scoreout,qpos,qprot,mutpos)
                      end if
                end do
        end do

        end subroutine pointmutate
c
c----------------------------------------------------------------------
c
        subroutine savekey(ix,score,qpos,qprot,mutpos)
        implicit none
        integer ix,score,qpos,qprot,mutpos

        ! check dimensions
        if(heap_n.ge.MAXKEY4) return ! overflow exit
        ! push to heap
        heap_input(heap_n+1)=ix
        keymap(1,ix)=score
        keymap(2,ix)=qpos
        keymap(3,ix)=qprot
        keymap(4,ix)=mutpos
        call pushheap(score,heap_n,heap_d,heap_v,MAXKEY4)
        !write(*,*) '# keymap',ix,score,qpos,qprot,mutpos

        end subroutine savekey
c
c----------------------------------------------------------------------
c
        end module keyspace
c
c======================================================================
c


c
c====================================================================
c
        module sans
        use altschul
        use maxheap, only : pushheap, popheap ! maxheap8.f
        use myhash, only : linkhash,nkey,nshort,shortlist
        use sanslocaldata, only: TXT,naa,nprot,SAP,SRES,CYCLIC,ptr,
     $          PHR,CYCMOD,ptr_phr,MPI_COMM_SANS,MAXPROC,
     $          naa_tot,nprot_tot
        use keyspace, only : keymap, get_keymap, MAXKEY4, get_ix
        implicit none
        integer qseqlen,naccept,nreject
        integer, parameter :: MAXVOTE=10000,MAXRES=64000,MAXH=1000
        integer, parameter :: tubewidth=10, MAXLINE=MAXH*MAXRES
        integer*1 qtxt(MAXRES)
        integer, parameter :: MAXDIAG=30000000 ! 30 M
        integer diaglist(3,MAXDIAG) ! (sprot,summa,diagonal=qres-sres)
        integer heap_v(MAXDIAG),heap_n,offset,heap_vote(MAXDIAG)
        integer heap_input(MAXDIAG),heap_diag(MAXDIAG),heap_d(MAXDIAG)
        integer, parameter :: LCP_DEPTH=127

        ! HACK: milestones independent of H
        integer, parameter :: HM=100

        ! HACK: SAP bounds for bsearch
        integer*8 SAP_left, SAP_rite

        ! hash diagonals per sprot
        integer, dimension(:),allocatable :: linklist_first
        integer linklist_key(MAXDIAG),linklist_next(MAXDIAG)
        integer linklist_value(MAXDIAG)
        !>> altschul module
        integer PSSM(0:20,MAXRES),blosum(0:20,0:20),alfmap(128)
        logical :: lfasta=.true.

        public 

        contains
c
c====================================================================
c
        subroutine make_qtxt_QSA(qseq,QSA,lreverse)
        implicit none
        character(len=MAXRES) rseq,qseq
        logical lreverse
        integer i,j,QSA(MAXRES)

        ! prepare reverse/forward qtxt
        if(lreverse) then
                j=0
                do i=qseqlen,1,-1
                        j=j+1
                        rseq(j:j)=qseq(i:i)
                end do
                call sais_from_fortran(rseq,QSA,qseqlen)
                do i=1,qseqlen
                        qtxt(i)=iachar(rseq(i:i))
                end do
        else
                call sais_from_fortran(qseq,QSA,qseqlen)
                do i=1,qseqlen
                        qtxt(i)=iachar(qseq(i:i))
                end do
        end if

        end subroutine make_qtxt_QSA
c
c====================================================================
c
        subroutine local_votes(QSA,vote_hist,naccept,
     $          rank,HX,W,MIN_SUMMA)
        implicit none
        include 'mpif.h'
        logical lreverse
        integer rank,ierr,qres,qprot,HX,W,MIN_SUMMA
        integer vote_hist(MAXVOTE),naccept
        integer i,j,sprot, s,ok,summa,QSA(MAXRES)
        integer*8 p,lex,diag,m_prev,mlist(MAXRES)
        integer*1 a,b,c,d,e,f
        integer*8 spos,m,m_left(MAXRES),m_rite(MAXRES)
        integer*8 llist(MAXRES),rlist(MAXRES)

        ! collect votes
        nkey=0 ! index to hash shortlist
        nshort=0
        m_prev=SAP_left
        do j=1,qseqlen ! test qpos in lexical order
                mlist(j)=0 ! case masked
                qres=QSA(j)
                if(qres.gt.qseqlen-10) cycle
                if(qres.lt.1) cycle
                ! skip composition bias
                a=qtxt(qres+1)
                if(a.eq.88) cycle ! X
                c=qtxt(qres+3)
                e=qtxt(qres+5)
                if(a.eq.c.and.c.eq.e) cycle
                ! binary search for insert position in SA
                !m=bsearch(qres,m_prev,naa-1,rank) ! uses qtxt
                m=bsearch(qres,m_prev,SAP_rite,rank) ! uses qtxt
                m_prev=m ! to be safe
                ! test min-max brackets
                mlist(j)=m
                if(m.lt.1) m=1
                if(m.ge.naa) m=naa-1
                llist(j)=CYCLIC(m)
                rlist(j)=CYCLIC(m+1)
                if(llist(j).gt.rlist(j)) rlist(j)=rlist(j)+CYCMOD
        end do
        call MPI_ALLREDUCE(llist,m_left,qseqlen,MPI_INTEGER8,
     $          MPI_MAX,MPI_COMM_SANS,ierr)
        call MPI_ALLREDUCE(rlist,m_rite,qseqlen,MPI_INTEGER8,
     $          MPI_MIN,MPI_COMM_SANS,ierr)
        do j=1,qseqlen
                m=mlist(j)
                if(m.lt.1) cycle
                if(m.gt.naa) cycle
                qres=QSA(j)
                if(qres.lt.1) cycle
                if(qres.gt.qseqlen) cycle
                lex=m+1
                if(lex.lt.1) cycle
                if(lex.gt.naa) cycle
                if(.false..and.qres.lt.qseqlen-20)
     $          write(*,*) '#m ',rank,(char(qtxt(qres+i)),i=0,19),
     $   qprot,qres,lex,naa,llist(j),rlist(j),m_left(j),m_rite(j),
     $   m_rite(j)-m_left(j),' ',
     $   (char(TXT(p)),p=ptr(SAP(m))+SRES(m),ptr(SAP(m))+SRES(m)+19)
                ! search left
                do while(CYCLIC(lex).gt.m_left(j)-HX)
                        lex=lex-1
                        if(lex.lt.1) exit ! loop
                        sprot=SAP(lex)
                        if(sprot.lt.1) cycle
                        if(sprot.gt.nprot) cycle
                        s=SRES(lex)
                        call linkhash(sprot,
     $         qres-s,1,linklist_key,linklist_next,
     $         linklist_value,linklist_first,MAXDIAG,nprot)
                        if(CYCLIC(lex).gt.CYCLIC(lex+1))
     $                                  m_left(j)=m_left(j)+CYCMOD
                end do
                ! search rite
                lex=mlist(j)
                if(lex.lt.1) cycle
                if(lex.gt.naa) cycle
                do while(CYCLIC(lex).lt.m_rite(j)+HX)
                        lex=lex+1
                        if(lex.ge.naa) exit ! loop
                        sprot=SAP(lex)
                        if(sprot.lt.1) cycle
                        if(sprot.gt.nprot) cycle
                        s=SRES(lex)
                        call linkhash(sprot,
     $         qres-s,1,linklist_key,linklist_next,
     $         linklist_value,linklist_first,MAXDIAG,nprot)
                        if(CYCLIC(lex).lt.CYCLIC(lex-1))
     $                                  m_rite(j)=m_rite(j)-CYCMOD
                end do
        end do
        ! return diaglist(3,:), vote_hist(MAXVOTE)
        call diag_sum_linklist(MIN_SUMMA,W,vote_hist,naccept,rank)
        ! clear hash table
        do i=1,nshort
                linklist_first(shortlist(i))=0
        end do
 
        end subroutine local_votes
c
c====================================================================
c
        subroutine do_all(socket,dbversion,programtag,nucl)
        implicit none
        include 'mpif.h'
        integer rank,ierr,ok,nproc,socket
        double precision starttime,t,t0,dt
        real fdr_cutoff
        character(len=MAXRES) qstring
        integer request,stat(MPI_STATUS_SIZE)
        logical flag,lterminate,nucl
        integer, parameter :: secs=0, usecs=500000, TIMEOUT = 60 ! snooze 0.5 s
        character(len=256) blosumfile,dbversion,programtag
        integer i,j,uid

        ! allocate hash tables
        allocate(linklist_first(nprot),stat=ok)
        call MPI_COMM_SIZE(MPI_COMM_SANS,nproc,ierr)
        call MPI_COMM_RANK(MPI_COMM_SANS,rank,ierr)
        ! alignment scoring
        if(nucl) then
                call mapASCII_NUC(alfmap)
                call getblosum_inline_NUC(blosum,alfmap)
        else
                call mapASCII(alfmap)
                call getblosum_inline(blosum,alfmap)
        end if
        ! hack
        SAP_left=1
        SAP_rite=naa
        do while(SAP(SAP_left).le.0.and.SAP_left.lt.naa) 
                SAP_left=SAP_left+1
        end do
        do while(SAP(SAP_rite).le.0.and.SAP_rite.gt.1)
                SAP_rite=SAP_rite-1
        end do
        write(*,*) '#SAP_left, SAP_rite',SAP_left,SAP_rite,naa
        ! input loop
        starttime=MPI_WTIME()
        t0=starttime
        uid=0
        linklist_first=0 ! clear hash
        keymap=0
        lterminate=.false.
        do while(.not.lterminate)
                ! all farm nodes snooze until receive query from socket
                call MPI_IRECV(qseqlen,1,MPI_INTEGER,MPI_ANY_SOURCE,
     $                  socket,MPI_COMM_WORLD,request,ierr)
                !write(*,*) '#received query',qseqlen,ierr
                do while(.not.lterminate)
                  call MPI_TEST(request,flag,stat,ierr)
                  if(flag) then ! data was receivedn
                        !write(*,*) '#receiving query',qprot,qseqlen
                        ! master receives qseq from socket;
                        ! bcast(qseq,master,mpi_comm_sans)
                        if(rank.eq.0) call MPI_RECV(qstring,qseqlen,
     $  MPI_CHARACTER,socket,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
                        call MPI_BCAST(qstring,qseqlen,MPI_CHARACTER,
     $                          0,MPI_COMM_SANS,ierr)
                        !write(*,*) '#query:',rank,qstring(1:qseqlen)
                        ! test for termination: qstring .eq. 'SHUTDOWN!'
                        if(qstring(1:qseqlen).eq.'SHUTDOWN!') then
                                lterminate=.true.
                                exit ! inner loop
                        end if
                        ! send <sbjctlist> to socket
                        uid=uid+1
                        call do_one(socket,qstring,rank,nproc,uid,
     $                          dbversion,programtag)
                        ! reset timeout
                        t0=MPI_WTIME()
                        exit ! loop
                  else ! test for timeout
                        dt=MPI_WTIME()-t0
                        if(dt.gt.TIMEOUT) then
                            !write(*,*) '#snoozing',qprot,dt
                            call microtimer(secs,usecs)
                        end if
                  end if
                end do
        end do
        ! master must send reply to socket on shutdown!
        if(rank.eq.0) then
          call MPI_SEND(qseqlen,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
          call MPI_SEND(qstring,qseqlen,MPI_CHARACTER,0,
     $                  0,MPI_COMM_WORLD,ierr)
        end if
        t=MPI_WTIME()
        write(*,*) '# ',rank,uid,' querys in ',t-starttime,' seconds'

        end subroutine do_all
c
c--------------------------------------------------------------------
c
        subroutine verifast(qprot,rank,nproc,qseq,qhdr,H,HX,W,MIN_SUMMA,
     $          votelist_size,protocol,dbversion,programtag)
        implicit none
        include 'mpif.h'
        integer rank,ierr,qprot,H,HX,W,MIN_SUMMA,nproc,votelist_size
        integer vote_cutoff,vote_hist(MAXVOTE),hist(MAXVOTE)
        integer i,summa,lmsg
        character(len=MAXRES) qseq,qstring,qhdr
        character(len=MAXRES) a
        character line(MAXLINE)
        character resultlines(MAXLINE)
        double precision evalue,t0,t
        integer QSA(MAXRES),protocol,x
        character(len=256) dbversion,programtag

        t0=MPI_WTIME()
        ! set EERIE vote cutoff
        call voting(qseq,.true.,H,HX,W,MIN_SUMMA,rank,vote_cutoff)
        ! add "safety margin"
        if(protocol.le.-99) then
                x=sqrt(float(vote_cutoff))
        else
                x=-protocol ! protcol is negative
        end if
        vote_cutoff=vote_cutoff+x 
        if(vote_cutoff.gt.MIN_SUMMA) MIN_SUMMA=vote_cutoff
        !write(*,*) '#verifast-reverse',rank,MIN_SUMMA,naccept

        ! increase vote_cutoff if too many hits
        call voting(qseq,.false.,votelist_size,HX,W,MIN_SUMMA,rank,
     $          vote_cutoff)
        !write(*,*) '#verifast-forward',rank,vote_cutoff,naccept

        ! direct output without alignment
        call open_queryblock(qprot,rank,vote_cutoff,qhdr,qseq,
     $          line,offset,dbversion,programtag,H,999.9)
        evalue=1.0
        do i=1,naccept
                summa=diaglist(2,i)
                if(summa.lt.vote_cutoff) cycle ! unsorted
                call write_sbjct(.true.,lfasta,
     $                  diaglist(1,i),summa,0,0,0,0.0,evalue,
     $                  diaglist(3,i),line,offset,0,0,0,0)
        end do
        ! collect all results on master and send to socket
        t=MPI_WTIME()-t0
        call close_queryblock(qprot,rank,nproc,offset,line,t)
        return

        end subroutine verifast
c
c--------------------------------------------------------------------
c
        subroutine voting(qseq,lreverse,H,HX,W,MIN_SUMMA,rank,
     $          vote_cutoff)
        implicit none
        include 'mpif.h'
        integer H,rank,vote_cutoff,HX,W,MIN_SUMMA
        logical lreverse
        character(len=MAXRES) qseq
        integer QSA(MAXRES),vote_hist(MAXVOTE),hist(MAXVOTE),htot,ierr

        call make_qtxt_QSA(qseq,QSA,lreverse) ! qtxt,qseqlen global
        call local_votes(QSA,vote_hist,naccept,rank,HX,W,MIN_SUMMA)
        call MPI_REDUCE(vote_hist, hist, MAXVOTE,MPI_INTEGER,
     $                  MPI_SUM, 0, MPI_COMM_SANS,IERR)
        if(rank.eq.0) call get_votecutoff(H,hist,vote_cutoff,MIN_SUMMA,
     $                  htot,.false.)
        call MPI_BCAST(vote_cutoff,1,MPI_INTEGER,0,MPI_COMM_SANS,ierr)

        end subroutine voting
c
c--------------------------------------------------------------------
c
        subroutine next_sbjct_heap(sprot,vote,tuplesum,diagonal)
        implicit none
        integer sprot,vote,tuplesum,diagonal
        integer j

        ! pop best
        j=heap_d(1)
        if(j.le.0) then
                write(*,*) '#WARNING: bad entry in heap',j,heap_n
                sprot=0
                tuplesum=0
                diagonal=0
                vote=0
                return
        end if
        sprot=heap_input(j)
        tuplesum=heap_v(j) ! tuplesum
        diagonal=heap_diag(j)
        vote=heap_vote(j) ! SANS-vote
        call popheap(heap_n,heap_d,heap_v,MAXDIAG)

        end subroutine next_sbjct_heap
c
c--------------------------------------------------------------------
c        
        subroutine do_one(socket,qstring,rank,nproc,uid, dbversion,
     $          programtag)
        implicit none
        include 'mpif.h'
        character(len=256) dbversion,programtag
        integer rank,ierr,qres,qprot,H,HX,W,MIN_SUMMA,nproc,socket,uid
        integer MINKEYSCORE
        integer vote_cutoff,vote_hist(MAXVOTE),hist(MAXVOTE)
        double precision starttime,t,searchtime
        real fdr(MAXVOTE),st,sf,sp,evalue_cutoff
        integer i,j,k,n,sprot, s,ok,summa,rhist(MAXVOTE),lmsg
        integer*8 p,lex,diag,m_prev,mlist(MAXRES)
        character(len=MAXRES) qseq,qstring,qhdr
        character(len=MAXRES) a
        character line(MAXLINE)
        character resultlines(MAXLINE)
        integer displs(MAXPROC),offsets(MAXPROC)
        integer diagonal,bestscore,bestnide,bestlali
        integer htot,rtot,ntot,newtot,vote, tuplesum
        real bitscore
        double precision evalue
        integer n1,QSA(MAXRES)
        integer niter,ntest,R,protocol,votelist_size
        integer score_cutoff,score_cutoff_prev,nvote
        integer, parameter :: maxiter=10
        integer milestones(maxiter),nnew,totnnew,n0

        ! extract runtime parameters from query string
        ! MINKEYSCORE, H, R, HX, W, MIN_SUMMA, tubewidth,
        ! EVALUE_CUTOFF
        read(qstring(1:qseqlen),*,err=99,end=99) qprot,H,HX,W,MIN_SUMMA,
     $          MINKEYSCORE,evalue_cutoff,R,votelist_size,protocol,
     $          qseq,qhdr
        lfasta=(qprot.ge.0) ! hack: no sseq output if qprot is negative
        ! check qseq  is [A-Za-z]
        qseqlen=len_trim(qseq)
        if(qseqlen.lt.6) goto 99 ! too short is illformed 
        do i=1,qseqlen
                k=iachar(qseq(i:i))
                if(k.lt.65.or.k.gt.122.or.(k.gt.90.and.k.lt.97)) then
                        qseq(i:i)='X'
                end if
        end do
        ! verifast procedure is special case: alpha < 0
        if(protocol.lt.0) then
              call verifast(qprot,rank,nproc,qseq,qhdr,H,HX,W,MIN_SUMMA,
     $                  votelist_size,protocol,dbversion,programtag)
              return
        end if
        t=MPI_WTIME()
        ! slow, verislow use tupleheap
!        if(rank.eq.0) write(*,*) '# do_one',rank,qprot,qseqlen,
!     $          R,votelist_size,protocol,H,HX,MINKEYSCORE,tubewidth,
!     $          evalue_cutoff
        ! vote cutoff; naccept entries saved in diaglist(3,:)
        call voting(qseq,.false.,votelist_size,HX,W,MIN_SUMMA,
     $          rank,vote_cutoff)
!        if(rank.eq.0) write(*,*)
!     $          '#voting time',MPI_WTIME()-t,vote_cutoff,naccept
        ! fast,slow,verislow use alignment; PSSM used by tuplesum
        call prepare_altschul(qseqlen,naa_tot,nprot_tot)
        call blosum2pssm(qtxt,qseqlen,blosum,PSSM,alfmap) 
         ! call keymap sequentially; do calculation on all nodes
         if(protocol.ge.1) then
!                t=MPI_WTIME()
                call get_keymap(uid,MINKEYSCORE,4,alfmap,
     $                  MPI_COMM_SANS,PSSM,qtxt,qseqlen)
!                if(rank.eq.0) write(*,*) '#keymap time',
!     $                  MPI_WTIME()-t,naccept
         end if
         ! diaglist(3,1:naccept) saved in local_votes
         ! filter1=vote; degenerate tail
!         t=MPI_WTIME()
         heap_n=0
         vote_hist=0
         do i=1,naccept
                vote=diaglist(2,i)
                if(vote.lt.vote_cutoff) cycle ! unsorted diaglist
                sprot=diaglist(1,i)
                diagonal=diaglist(3,i)
                select case(protocol)
                        case(0)
                                tuplesum=vote
                        case(1) ! vote diagonals are bad!
                                call get_bestdiag(uid,sprot,tubewidth,
     $                                  tuplesum,diagonal)
                        case(2) ! tuple diagonals are much better!
                                call get_bestdiag0(uid,sprot,tubewidth,
     $                                  tuplesum,diagonal)
                end select
!         write(*,*)'#best',sprot,summa,tuplesum,diaglist(3,i),diagonal
                if(tuplesum.lt.1) cycle 
                heap_input(heap_n+1)=sprot
                heap_diag(heap_n+1)=diagonal
                heap_vote(heap_n+1)=vote ! SANS-vote
                call pushheap(tuplesum,heap_n,heap_d,heap_v,MAXDIAG)
!                write(*,*) '#push',uid,sprot,vote,tuplesum,heap_n
                ! count tuplesum histogram for milestones
                j=min(MAXVOTE,tuplesum)
                vote_hist(j)=vote_hist(j)+1 ! SANS-vote or tuplesum
                if(heap_n.ge.MAXDIAG) exit ! loop
         end do
        call get_milestones(rank,vote_hist,max(H,HM),maxiter,milestones)
!         if(rank.eq.0) write(*,*) '#keyspace time',MPI_WTIME()-t,heap_n

         ! output result
         call open_queryblock(qprot,rank,vote_cutoff,qhdr,qseq,
     $          line,offset,dbversion,programtag,H,evalue_cutoff)
         ! test for termination: total H accepts or R rejects
         naccept=0
         nreject=0
         ! milestones so each rank goes to equal depth each iteration
         score_cutoff=get_bitscorecutoff(evalue_cutoff)
         score_cutoff_prev=0
         niter=0
         ntest=0
         n0=heap_n
         vote_hist=0 ! test convergence of Hth alignment score
         htot=0
         rtot=0
         ntot=0
         newtot=0
         do while(.true.)
           nnew=0
           niter=niter+1
           if(niter.gt.maxiter) exit ! loop
           write(*,*) '#iter',niter,milestones(niter),naccept,nreject,
     $          heap_n,score_cutoff_prev
           do while(heap_n.gt.0)
                ntest=ntest+1
                summa=-1
                call align_next_sbjct(evalue_cutoff,vote_hist,
     $                  line,offset,summa,nnew,score_cutoff_prev)
                if(summa.lt.milestones(niter)) exit ! sbjct loop
           end do
           ! test convergence
!           write(*,*) '#test',heap_n,H,R,naccept,nreject,nnew,
!     $ ntest,score_cutoff_prev,htot
           call get_totals(naccept,nreject,heap_n,nnew,
     $          htot,rtot,ntot,newtot)
           if(htot.ge.H.and.protocol.eq.0) exit ! terminate: H accepts
           if(rtot.ge.R) exit ! terminate: too many rejects
           if(ntot.le.0) exit ! terminate: no more hits in heap
           if(newtot.le.0) exit ! terminate: no change in top-H
!           if(terminate(protocol,H,R,heap_n,nnew)) exit ! loop
           call MPI_REDUCE(vote_hist, hist, MAXVOTE,MPI_INTEGER,
     $          MPI_SUM, 0, MPI_COMM_SANS,IERR)
           if(rank.eq.0)  call get_votecutoff(H,hist,score_cutoff,
     $          MIN_SUMMA,htot,.false.)
           call MPI_BCAST(score_cutoff,1,MPI_INTEGER,0,
     $          MPI_COMM_SANS,ierr)
           score_cutoff_prev=score_cutoff
         end do
         ! collect all results on master and send to socket
         searchtime=MPI_WTIME()-t
         call close_queryblock(qprot,rank,nproc,offset,line,searchtime)
         if(rank.eq.0) write(*,*) '#time',searchtime,qprot,niter,
     $          htot,rtot,ntot,newtot
         return

         ! error exit; master sends error message to socket
99       call illformed(qprot,rank,qstring)
         return 

         end subroutine do_one
c
c--------------------------------------------------------------------
c
         subroutine align_next_sbjct(evalue_cutoff,
     $          vote_hist,line,offset,summa,nnew,score_cutoff_prev)
         implicit none
         integer sprot,vote,summa,diagonal,vote_hist(MAXVOTE),j
         integer bestscore,bestnide,bestlali,offset
         integer score_cutoff_prev,nnew
         character line(MAXLINE)
         real bitscore,evalue_cutoff
         double precision evalue
         integer qfrom,qto,sfrom,sto

!         write(*,*) '#this is align_next_sbjct'
         ! pop best
         call next_sbjct_heap(sprot,vote,summa,diagonal)
!         write(*,*) '#next_sbjct',sprot,vote,summa,diagonal,heap_n
         ! test bitscore
         call tube_dp(sprot,diagonal,tubewidth,bestscore,
     $                  bestnide,bestlali,qfrom,qto,sfrom,sto)
!         write(*,*) '#bestscore',bestscore
         if(bestscore.gt.score_cutoff_prev) nnew=nnew+1
         if(bestscore.lt.1) return
         j=min(MAXVOTE,bestscore)
         vote_hist(j)=vote_hist(j)+1
         bitscore=get_bitscore(bestscore)
         evalue=get_evalue(bitscore)
!         write(*,*) '#pop',sprot,diagonal,summa,
!     $                bitscore,evalue,naccept,nreject,heap_n
         ! compare rejects to score_cutoff
         if(evalue.gt.evalue_cutoff) then
                nreject=nreject+1
         else
                naccept=naccept+1
                call write_sbjct(.false.,lfasta,sprot,summa,bestscore,
     $                  bestnide,bestlali,bitscore/log(2.0),evalue,
     $                  diagonal,line,offset,qfrom,qto,sfrom,sto)
         end if

         end subroutine align_next_sbjct
c
c--------------------------------------------------------------------
c
        subroutine send_result(qprot,lmsg,resultlines)
        implicit none
        include 'mpif.h'
        integer lmsg,ierr,qprot
        character resultlines(lmsg)

!        write(*,*) '# Result on server:',qprot,lmsg,char(10),
!     $          (resultlines(i),i=1,lmsg),'//'
         !write(*,*) '# length of result is',qprot,lmsg
         if(lmsg>MAXLINE) lmsg=MAXLINE
         ! socket is MPI_COMM_WORLD rank 0
         call MPI_SEND(lmsg,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
         !write(*,*) '#sent result header',lmsg,ierr
         call MPI_SEND(resultlines,lmsg,MPI_CHARACTER,0,
     $                  0,MPI_COMM_WORLD,ierr)
         !write(*,*) '# sent result data',lmsg,ierr

        end subroutine send_result
c
c--------------------------------------------------------------------
c
        subroutine get_totals(naccept,nreject,heap_n,nnew,
     $          htot,rtot,ntot,newtot)
        implicit none
        include 'mpif.h'
        integer naccept,nreject,heap_n,nnew,htot,rtot,ntot,newtot
        integer datin(4),datout(4),ierr

        datin(1)=naccept
        datin(2)=nreject
        datin(3)=heap_n
        datin(4)=nnew
        call MPI_ALLREDUCE(datin,datout,4, 
     $            MPI_INTEGER,MPI_SUM,MPI_COMM_SANS,ierr)
        htot=datout(1)
        rtot=datout(2)
        ntot=datout(3)
        newtot=datout(4)
        end subroutine get_totals
c
c--------------------------------------------------------------------
c
        logical function terminate(protocol,H,R,n,nnew) result(lterm)
        implicit none
        include 'mpif.h'
        integer datin(4),datout(4),ierr,n,H,R,protocol,nnew

        lterm=.true.
        datin(1)=naccept
        datin(2)=nreject
        datin(3)=n ! heap_n
        datin(4)=nnew
        call MPI_ALLREDUCE(datin,datout,4, 
     $            MPI_INTEGER,MPI_SUM,MPI_COMM_SANS,ierr)
        !write(*,*) '#term',datout,H,R
        if(protocol.lt.0.and.datout(1).ge.H) return ! terminate: global naccept.ge.H
        if(datout(2).ge.R) return ! terminate: global nreject.ge.R
        if(datout(3).le.0) return ! terminate: no more hits in heap on any rank
        if(datout(4).eq.0) return ! terminate: no change in top-H list on iteration
        lterm=.false.

        end function terminate
c
c--------------------------------------------------------------------
c
        subroutine illformed(qprot,rank,qstring)
        implicit none
        integer qprot,lmsg,rank
        character(len=MAXRES) qstring
        character resultlines(MAXLINE)
        character(len=MAXRES) a

        ! only master sends result to socket
        if(rank.ne.0) return
        write(a,*) '<QUERY nid=',qprot,'>',char(10),
     $  '#ill-formed query: ',qstring(1:qseqlen),char(10),
     $  '</QUERY>',char(10)
        lmsg=0
        call writeresult(a,len_trim(a)+1,resultlines,lmsg)
        call send_result(qprot,lmsg,resultlines)

        end subroutine illformed
c
c--------------------------------------------------------------------
c      sandwiched between open_queryblock() and close_queryblock()
c
        subroutine write_sbjct(lshort,lfasta,sprot,vote,summa,nide,lali,
     $         bitscore,evalue,diagonal,line,offset,qfrom,qto,sfrom,sto)
        implicit none
        integer sprot,summa,nide,lali,diagonal,offset,sseqlen,vote
        real bitscore,pide
        double precision evalue
        character line(MAXLINE)
        character(len=MAXRES) a
        integer*8 p
        logical lshort,lfasta
        integer qfrom,qto,sfrom,sto

        sseqlen=ptr(sprot+1)-ptr(sprot)-1
        if(lshort) then ! verifast: sprot,vote,diagonal
                write(a,
     $  '(a12,i8,1x,2(a6,i8),1x,2a1)')
     $          '<SBJCT VOTE=',vote,
     $          ' DIAG=',diagonal,
     $          ' LSEQ=',sseqlen,
     $          '>',char(10)
        else
                pide=float(nide)/max(1,lali)
                write(a,
     $  '(a12,i8,a6,i8,a6,f6.2,a6,i6,a6,f8.1,a8,e16.8,2(a6,i8),'//
     $          '2(a7,i6,a5,i6),1x,2a1)')
     $          '<SBJCT VOTE=',vote,
     $          ' TUPS=',summa,
     $          ' PIDE=',pide,
     $          ' LALI=',lali,
     $          ' BITS=',bitscore,
     $          ' EVALUE=',evalue,
     $          ' DIAG=',diagonal,
     $          ' LSEQ=',sseqlen,
     $          ' QFROM=',qfrom,
     $          ' QTO=',qto,
     $          ' SFROM=',sfrom,
     $          ' STO=',sto,
     $          '>',char(10)
        end if
        call writeresult(a,len_trim(a)+1,line,offset)
        line(offset)='>' ! fasta header tag
        do p=ptr_phr(sprot)+1,ptr_phr(sprot+1)
                offset=offset+1
                line(offset)=char(PHR(p))
        end do
        if(lfasta) then
                do p=ptr(sprot)+1,ptr(sprot+1)
                        offset=offset+1
                        line(offset)=char(TXT(p))
                end do
        else ! empty sequence line
                offset=offset+1
                line(offset)=char(10)
        end if
        write(a,'(a8,a1)') '</SBJCT>',char(10)
        call writeresult(a,9,line,offset)

         end subroutine write_sbjct
c
c--------------------------------------------------------------------
c       query header must be written before sbjcts
c 
        subroutine open_queryblock(qprot,rank,vote_cutoff,
     $          qhdr,qseq,line,offset,dbversion,programtag,
     $          H,evalue_cutoff)
        implicit none
        integer H
        real evalue_cutoff
        integer qprot,vote_cutoff,offset,cdcutoff,rank
        character line(MAXLINE)
        character(len=MAXRES) qseq,qhdr
        character(len=256) dbversion,a,b,c,programtag

        offset=0
        if(rank.ne.0) return 
        write(a,'(a11,i6,a13,i12,a13,i12,a6,i8,1x,a1)') '<QUERY nid=',
     $          qprot,' vote_cutoff=',vote_cutoff,' LSEQ=',qseqlen,'>'
        call writeresult(a,len_trim(a),line,offset)
        call writeresult(char(10),1,line,offset)
        write(b,*) naa_tot
        write(c,*) nprot_tot
        a='<DATABASE= '//trim(dbversion)//' letters= '//trim(b)
     $          //' sequences= '//trim(c)//' >'//char(10)
        call writeresult(a,len_trim(a),line,offset)
        a=trim(programtag)//char(10)
        call writeresult(a,len_trim(a),line,offset)
        write(a,*) '<PARAM H=',H,' EVALUE_CUTOFF=',evalue_cutoff,'>'
     $          ,char(10)
        call writeresult(a,len_trim(a),line,offset)
        write(a,*) '<timeInfo start="',timestamp(),'"/>',char(10)
        call writeresult(a,len_trim(a),line,offset)
        ! echo query-fasta
        ! Jalview, Muscle accept repeated protein-identifier
        call writeresult(qhdr,len_trim(qhdr),line,offset)
        call writeresult(char(10),1,line,offset)
        call writeresult(qseq,len_trim(qseq),line,offset)
        call writeresult(char(10),1,line,offset)

        end subroutine open_queryblock
c
c--------------------------------------------------------------------
c
        subroutine close_queryblock(qprot,rank,nproc,offset,line,t)
        implicit none
        include 'mpif.h'
        integer offset,lmsg,nproc,rank,qprot
        double precision t
        character(len=MAXRES) a
        character line(MAXLINE)
        character resultlines(MAXLINE)
        integer displs(MAXPROC),offsets(MAXPROC)
        integer i,ierr

        call MPI_GATHER(offset,1,MPI_INTEGER,offsets,1,
     $          MPI_INTEGER,0, MPI_COMM_SANS, ierr)
        if(rank.eq.0) then
                displs(1)=0
                do i=2,nproc
                        displs(i)=displs(i-1)+offsets(i-1)
                       if(displs(i)+offsets(i).gt.MAXLINE)
     $  offsets(i)=max(0,MAXLINE-displs(i))
                end do
        end if
        ! in case offsets were modified
        call MPI_SCATTER(offsets,1,MPI_INTEGER,offset,1,MPI_INTEGER,
     $          0,MPI_COMM_SANS,ierr)
        call MPI_GATHERV(line,offset,MPI_CHARACTER,resultlines,
     $          offsets,displs,MPI_CHARACTER,0, MPI_COMM_SANS, ierr)
        lmsg=displs(nproc)+offsets(nproc)
        ! close query block
        if(rank.eq.0) then
                ! output search time
                write(a,*) '<timeInfo search= ',t,' >',char(10)
                call writeresult(a,len_trim(a),resultlines,lmsg)
                write(a,*) '<timeInfo end="',timestamp(),'"/>',char(10)
                call writeresult(a,len_trim(a),resultlines,lmsg)
                ! message terminator
                write(a,'(a8,a1)') '</QUERY>',char(10)
                call writeresult(a,9,resultlines,lmsg)
               ! only master sends result to socket
                call send_result(qprot,lmsg,resultlines)
        end if

        end subroutine close_queryblock
c
c--------------------------------------------------------------------
c 
        function timestamp() result(a)
        implicit none
        integer*4 today(3), now(3)
        character(len=30) a

        call idate(today)   ! today(1)=day, (2)=month, (3)=year
        call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
        write ( a, 1000 )  today(2), today(1), today(3), now
 1000 format ('Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
        
        end function timestamp
c
c--------------------------------------------------------------------
c       sum tuples in votediagonal window
c
       subroutine get_bestdiag(uid,sprot,tw,summa,diagonal)
       implicit none
       integer sprot,diagonal,sres,d,six,sseqlen,ix(maxres),uid
       integer*8 spos
       integer*1 s(maxres)
       integer summa,tw,mind,maxd

       sres=0
       sseqlen=ptr(sprot+1)-ptr(sprot)
       spos=ptr(sprot)
       do sres=1,sseqlen
                spos=spos+1
                s(sres)=alfmap(TXT(spos))
       end do
       call get_ix(s,sseqlen,ix)

       mind=diagonal-tw
       maxd=diagonal+tw
       summa=0
       do sres=1,sseqlen-3
                ! diagonal=qres-sres
                six=ix(sres)
                if(six.lt.1) cycle
                if(keymap(3,six).ne.uid) cycle
                d=keymap(2,six)-sres
                if(d.lt.mind) cycle
                if(d.gt.maxd) cycle
                summa=summa+keymap(1,six) ! 4-tuple matches in diag-band
       end do

       end subroutine get_bestdiag
c
c--------------------------------------------------------------------
c 
       subroutine get_bestdiag0(uid,sprot,tw,bests,bestdiag)
       implicit none
       integer sprot,diagonal,sres,d,six,i,j,sseqlen,ix(maxres),uid
       integer*8 spos
       integer*1 s(maxres)
       integer summa(-maxres:maxres),tw,bests,bestdiag,x,mind,maxd

       sres=0
       sseqlen=ptr(sprot+1)-ptr(sprot)
       spos=ptr(sprot)
       do sres=1,sseqlen
                spos=spos+1
                s(sres)=alfmap(TXT(spos))
       end do
       call get_ix(s,sseqlen,ix)

       do d=-sseqlen-tw,max(qseqlen,sseqlen)+tw
                summa(d)=0
       end do
       mind=sseqlen
       maxd=-sseqlen
       do sres=1,sseqlen-3
                ! diagonal=qres-sres
                six=ix(sres)
                if(six.lt.1) cycle
                if(keymap(3,six).ne.uid) cycle
                d=keymap(2,six)-sres
                summa(d)=summa(d)+keymap(1,six) ! 4-tuple matches in diag-band
!                write(*,*) '#d',sprot,sres,s(sres),six,d,summa(d),
!     $                  (keymap(j,six),j=1,2)
                if(d.lt.mind) mind=d
                if(d.gt.maxd) maxd=d
       end do
       ! find best window
       x=0
       do d=mind,mind+tw-1
                x=x+summa(d)
       end do
       bests=0
       bestdiag=0
       do d=mind,maxd
                x=x-summa(d-tw)+summa(d+tw)
                if(x.gt.bests) then
                        bests=x
                        bestdiag=d
                end if
       end do
!       write(*,*) '#bestdiag result',uid,sprot,mind,maxd,bests,bestdiag

       end subroutine get_bestdiag0
c
c--------------------------------------------------------------------
c 
        subroutine get_votecutoff(H,hist,vote_cutoff,MIN_SUMMA,HTOT,
     $          verbose)
        integer hist(MAXVOTE),vote_cutoff,MIN_SUMMA,H
        logical verbose
        integer HTOT,i,j

        j=0
        !write(*,*) '#get_votecutoff',H,MIN_SUMMA
        do i=MAXVOTE,MIN_SUMMA,-1
                if(verbose) write(*,*) '#votes',i,hist(i),j
                if(hist(i).le.0) cycle
                j=j+hist(i)
                if(j.ge.H) exit
        end do
        vote_cutoff=i
        HTOT=j

        end subroutine get_votecutoff
c
c--------------------------------------------------------------------
c 
        subroutine get_milestones(rank,vote_hist,H,maxiter,milestones)
        implicit none
        include 'mpif.h'
        integer vote_hist(MAXVOTE),hist(MAXVOTE),maxiter,rank
        integer milestones(maxiter),i,j,H,trgt,s,ierr,v0

        !write(*,*) '#get_milestones',rank,H,maxiter
        call MPI_REDUCE(vote_hist,hist,MAXVOTE, MPI_INTEGER,MPI_SUM,
     $          0,MPI_COMM_SANS,IERR)
        milestones=0
        if(rank.eq.0) then
                trgt=H
                v0=MAXVOTE
                do i=1,maxiter
                        s=0
                        !write(*,*) '#iter',i,v0,trgt
                        do j=v0,1,-1
                                s=s+hist(j)
                                if(s.lt.trgt) cycle
                                milestones(i)=j
                                v0=j-1
                                exit ! loop
                        end do
                        trgt=trgt+H
                 end do
                !write(*,*) '#iter done',v0,trgt
        end if
        call MPI_BCAST(milestones,maxiter,MPI_INTEGER,0,MPI_COMM_SANS,
     $          IERR)

        end subroutine get_milestones
c
c--------------------------------------------------------------------
c 
        subroutine writeresult(string,l,resultlines,lmsg)
        implicit none
        integer l,lmsg,i
        character(len=l) string
        character resultlines(MAXLINE)

        do i=1,min(l,MAXLINE-lmsg)
                lmsg=lmsg+1
                resultlines(lmsg)=string(i:i)
        end do

        end subroutine writeresult
c
c--------------------------------------------------------------------
c
        subroutine diag_sum_linklist(MIN_SUMMA,W,vote_hist,naccept,rank)
        implicit none
        integer i,summa,ix,sprot,j,k,n,d(MAXDIAG),v(MAXDIAG),beg,bestsum
        integer vsum(MAXRES),MIN_SUMMA,W,s,bests,vote_hist(MAXVOTE)
        integer naccept,rank
        integer*8 diag,bestdiag

        heap_n=0
        naccept=0
        vote_hist=0
        !write(*,*) '#diag_sum_linklist',rank,nshort
        do j=1,nshort
                sprot=shortlist(j)
                ix=linklist_first(sprot)
                n=0
                ! collect list of diagonals (dots)
                do while(ix.gt.0)
                        n=n+1
                        d(n)=n
                        v(n)=linklist_key(ix) ! diagonal
                        ix=linklist_next(ix)
                end do
                ! sort diagonals
                if(n.gt.1) call j_index_Qsort(1,n,n,d,v)
                !write(*,*) '#sorted',rank,n,(v(d(i)),i=1,n)
                ! rolling sum over diagonals; remember best
                i=1
                beg=1
                bestsum=0
                bestdiag=0 ! highest scoring diagonal
                bests=0
                do while(i.lt.n.and.heap_n.lt.MAXDIAG)
                        ix=d(i)
                        if(ix.lt.1) exit ! loop
                        i=i+1
                        do while(v(ix)-v(d(beg)).gt.W.and.beg.lt.n)
                                beg=beg+1
                        end do
                        k=d(i)
                        do while(k.gt.0.and.v(ix).eq.v(k).and.i.lt.n)
                                i=i+1 ! same diagonal
                                k=d(i)
                        end do
                        summa=i-beg+1 ! diagonal weight == 1
                        if(summa.gt.bestsum) then
                                bestsum=summa
                                bestdiag=(v(d(beg))+v(d(i)))/2
!      write(*,*) '#best',rank,summa,beg,i,v(d(beg)),v(d(i)),bestdiag,
!     $  naccept
                        end if
                end do
                if(bestsum.lt.MIN_SUMMA) cycle
                ! save best only per sprot
                naccept=naccept+1
                diaglist(1,naccept)=sprot
                diaglist(2,naccept)=bestsum
                diaglist(3,naccept)=bestdiag ! not used
                if(bestsum.gt.MAXVOTE) bestsum=MAXVOTE
                vote_hist(bestsum)=vote_hist(bestsum)+1
                if(naccept.ge.MAXDIAG) exit ! loop
        end do
        
        end subroutine diag_sum_linklist
c
c--------------------------------------------------------------------
c
        recursive function bsearch(qres,left,rite,rank) result(m)
        implicit none
        integer*8 left,rite,m,i
        integer qres,qseqlen,l,rank
        integer*8 spos_m,spos_l,spos_r
        integer*1 q,s

        m=0
        if(left.lt.1) return
        if(rite.gt.naa) return 
        m=(left+rite)/2 
        spos_m=ptr(SAP(m))+SRES(m)-1
        spos_l=ptr(SAP(left))+SRES(left)-1
        spos_r=ptr(SAP(rite))+SRES(rite)-1
        l=lcp(qres,spos_m)+1
        q=qtxt(qres+l)
        s=TXT(spos_m+l)
        if(rite.le.left) then
                if(q.lt.s) m=m-1 ! always return txt(m)<=qtxt
                return
                write(*,*) '#bsearch result',rank,qres,left,m,
     $          rite,spos_l,spos_m,spos_r,
     $          l,' ',(char(qtxt(i)),i=qres+1,qres+10),' ',
     $          (char(max(32,TXT(i))),i=spos_l+1,spos_l+10),' ',
     $          (char(max(32,TXT(i))),i=spos_m+1,spos_m+10),' ',
     $          (char(max(32,TXT(i))),i=spos_r+1,spos_r+10),
     $          (q.lt.s),(q.gt.s)
                 return
        end if
        if(.false.) then
                write(*,*) '#bsearch',rank,qres,left,m,
     $          rite,spos_l,spos_m,spos_r,
     $          SAP(m),SRES(m),
     $          SAP(left),SRES(left),SAP(rite),SRES(rite),
     $          l,' ',(char(qtxt(i)),i=qres+1,qres+10),' ',
     $          (char(max(32,TXT(i))),i=spos_l+1,spos_l+10),' ',
     $          (char(max(32,TXT(i))),i=spos_m+1,spos_m+10),' ',
     $          (char(max(32,TXT(i))),i=spos_r+1,spos_r+10),
     $          (q.lt.s),(q.gt.s)
        end if

        if(q.lt.s) then ! query is in left half
                m=bsearch(qres,left,m-1,rank)
        else if(q.gt.s) then ! query is in right half
                m=bsearch(qres,m+1,rite,rank)
        else ! equal
                return
        end if

        end function bsearch
c
c--------------------------------------------------------------------
c
        integer function lcp(qres,spos)
        implicit none
        integer qres
        integer*8 spos

        lcp=0
        do while(qtxt(qres+lcp+1).eq.TXT(spos+lcp+1))
                lcp=lcp+1
                if(lcp.gt.LCP_DEPTH) exit
                if(qres+lcp.ge.qseqlen) exit
                if(spos+lcp.ge.naa) exit
                if(TXT(spos+lcp).lt.32) exit
        end do

        end function lcp
c
c
c--------------------------------------------------------------------
c       diagonal=qres-sres
c
        subroutine tube_dp(sprot,diagonal,tubewidth,bestscore,
     $          bestnide,bestlali,bestqfrom,bestqto,bestsfrom,beststo)
        implicit none
        integer, parameter :: maxtubewidth=1000
        integer qprot,sprot,bestscore,diagonal,tubewidth
        integer bestqfrom,bestqto,bestsfrom,beststo
        integer bestnide,bestlali
        integer insert(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer delete(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer match(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer nide_match(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer lali_match(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer nide_insert(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer lali_insert(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer nide_delete(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer lali_delete(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer qfrom_match(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer qfrom_insert(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer qfrom_delete(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer sfrom_match(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer sfrom_insert(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer sfrom_delete(0:1,-maxtubewidth-1:maxtubewidth+1)
        integer x,qpos,testdiag,ide
        integer i,j,curr,prev,s(MAXRES),dat(3)
        integer*8 s0,diag,spos,sleft,srite
        integer, parameter :: GAPOPEN=-10,GAPELON=-2,infinite=99999999
        integer*1 stxt(maxres)

        ! check dimensions
        if(tubewidth.gt.maxtubewidth) stop 'increase maxtubewidth'

        s0=ptr(sprot)
        !write(*,*) '#tube',sprot,diagonal,ptr(sprot),ptr(sprot+1)
        ! get best trace
        qpos=1
        sleft=0
        srite=ptr(sprot+1)-ptr(sprot)-1
        ! initialize ARRAYS
         delete=0
         match=0
         insert=0
         lali_match=0
         lali_insert=0
         lali_delete=0
         nide_match=0
         nide_insert=0
         nide_delete=0
         qfrom_match=0
         qfrom_insert=0
         qfrom_delete=0
         sfrom_match=0
         sfrom_insert=0
         sfrom_delete=0
         i=0
         do spos=ptr(sprot)+1,ptr(sprot+1)-1
                i=i+1
                if(i.gt.MAXRES) exit
                ! longest sequences 74488 aa and 74074 aa
                ! in uniprot.Oct2018
                j=TXT(spos)
                stxt(i)=j
                s(i)=alfmap(j)       
         end do

        ! iterate over qpos; i=mod(qpos,2)=index to match/insert/delete
        bestscore=0
        bestnide=0
        bestlali=0
        bestqfrom=1
        bestqto=0
        bestsfrom=1
        beststo=0
        prev=1
        curr=0
        do qpos=1,QSEQLEN
                ! swap prev, curr
                x=prev
                prev=curr
                curr=x
                ! jump between diagonals [-tubewidth:tubewidth]
                do testdiag=-tubewidth,tubewidth
                        spos=qpos-diagonal+testdiag
                        if(spos.le.sleft) cycle
                        if(spos.ge.srite) cycle
                        delete(curr,testdiag)=0 
                        match(curr,testdiag)=0
                        insert(curr,testdiag)=0
                        ! update match/insert/delete score @ (qpos,spos)
                        ide=0
                        if(qtxt(qpos).eq.stxt(spos)) ide=1
                        dat(1)=match(prev,testdiag)
                        dat(2)=delete(prev,testdiag)
                        dat(3)=insert(prev,testdiag)
                        j=max3(dat)
                        match(curr,testdiag)=PSSM(s(spos),qpos)+dat(j)
                        select case(j)
                                case(1) ! M->M
                                        nide_match(curr,testdiag)=
     $                                  nide_match(prev,testdiag)+ide
                                        lali_match(curr,testdiag)=
     $                                  lali_match(prev,testdiag)+1
                                        qfrom_match(curr,testdiag)=
     $                                  qfrom_match(prev,testdiag)
                                        sfrom_match(curr,testdiag)=
     $                                  sfrom_match(prev,testdiag)
                                case(2) ! D->M
                                        nide_match(curr,testdiag)=
     $                                  nide_delete(prev,testdiag)+ide
                                        lali_match(curr,testdiag)=
     $                                  lali_delete(prev,testdiag)+1
                                        qfrom_match(curr,testdiag)=
     $                                  qfrom_delete(prev,testdiag)
                                        sfrom_match(curr,testdiag)=
     $                                  sfrom_delete(prev,testdiag)
                                case(3) ! I->M
                                        nide_match(curr,testdiag)=
     $                                  nide_insert(prev,testdiag)+ide
                                        lali_match(curr,testdiag)=
     $                                  lali_insert(prev,testdiag)+1
                                        qfrom_match(curr,testdiag)=
     $                                  qfrom_insert(prev,testdiag)
                                        sfrom_match(curr,testdiag)=
     $                                  sfrom_insert(prev,testdiag)
                        end select
                        ! case we started a new trace at this match
                        if(qfrom_match(curr,testdiag).eq.0) then
                                qfrom_match(curr,testdiag)=qpos
                                sfrom_match(curr,testdiag)=spos
                        end if
                        dat(1)=match(curr,testdiag-1)+GAPOPEN !+GAPELON
                        dat(2)=insert(curr,testdiag-1)+GAPELON
                        if(dat(1).ge.dat(2)) then
                                j=1
                        else
                                j=2
                        end if
                        insert(curr,testdiag)=dat(j)
                        select case(j)
                                case(1) ! M->I
                                        nide_insert(curr,testdiag)=
     $                                  nide_match(curr,testdiag-1)
                                        lali_insert(curr,testdiag)=
     $                                  lali_match(curr,testdiag-1)
                                        qfrom_insert(curr,testdiag)=
     $                                  qfrom_match(curr,testdiag-1)
                                        sfrom_insert(curr,testdiag)=
     $                                  sfrom_match(curr,testdiag-1)
                                case(2) ! I->I
                                        nide_insert(curr,testdiag)=
     $                                  nide_insert(curr,testdiag-1)
                                        lali_insert(curr,testdiag)=
     $                                  lali_insert(curr,testdiag-1)
                                        qfrom_insert(curr,testdiag)=
     $                                  qfrom_insert(curr,testdiag-1)
                                        sfrom_insert(curr,testdiag)=
     $                                  sfrom_insert(curr,testdiag-1)
                        end select
                        dat(1)=match(prev,testdiag+1)+GAPOPEN !+GAPELON
                        dat(2)=delete(prev,testdiag+1)+GAPELON
                        if(dat(1).ge.dat(2)) then
                                j=1
                        else
                                j=2
                        end if
                        delete(curr,testdiag)=dat(j)
                        select case(j)
                                case(1) ! M->D
                                        nide_delete(curr,testdiag)=
     $                                  nide_match(prev,testdiag+1)
                                        lali_delete(curr,testdiag)=
     $                                  lali_match(prev,testdiag+1)
                                        qfrom_delete(curr,testdiag)=
     $                                  qfrom_match(prev,testdiag+1)
                                        sfrom_delete(curr,testdiag)=
     $                                  sfrom_match(prev,testdiag+1)
                                case(2) ! D->D
                                        nide_delete(curr,testdiag)=
     $                                  nide_delete(prev,testdiag+1)
                                        lali_delete(curr,testdiag)=
     $                                  lali_delete(prev,testdiag+1)
                                        qfrom_delete(curr,testdiag)=
     $                                  qfrom_delete(prev,testdiag+1)
                                        sfrom_delete(curr,testdiag)=
     $                                  sfrom_delete(prev,testdiag+1)
                        end select
                        ! update bestscore
                        if(match(curr,testdiag).gt.bestscore) then
                                bestscore=match(curr,testdiag)
                                bestnide=nide_match(curr,testdiag)
                                bestlali=lali_match(curr,testdiag)
                                bestqfrom=qfrom_match(curr,testdiag)
                                bestsfrom=sfrom_match(curr,testdiag)
                                bestqto=qpos
                                beststo=spos
                        end if
                        ! smith-waterman
                        if(match(curr,testdiag).lt.0) then
                                match(curr,testdiag)=0
                                lali_match(curr,testdiag)=0
                                nide_match(curr,testdiag)=0
                                qfrom_match(curr,testdiag)=qpos
                                sfrom_match(curr,testdiag)=spos
                        end if
                end do
!                write(*,*)'#',qpos,char(qtxt(qpos)),' ',
!     $                  (char(TXT(s0+spos)),
!     $                  PSSM(alfmap(TXT(s0+spos)),qpos),
!     $                  spos=max(sleft+1,qpos-diagonal-tubewidth),
!     $                  min(srite-1,qpos-diagonal+tubewidth)),bestscore
        end do
       ! write(*,*) '#alignment score',bestscore

        end subroutine tube_dp
c
c----------------------------------------------------------------------
c
        integer function max3(dat) result(k)
        implicit none
        integer dat(3)

        if(dat(1).ge.dat(2).and.dat(1).ge.dat(3)) then
                k=1
        else if(dat(2).ge.dat(3)) then
                k=2
        else
                k=3
        end if

        end function max3
c
c----------------------------------------------------------------------
c

        end module sans
c
c======================================================================
c


c>> >2 G chunks bug (integer*8 indices)
c>> ? bug in scatterv: n>naa (x_8.1.6, x_16.1.06, x_32.1.06, x_64.1.38)
c
c       module add openmpi-x86_64
c       mpif90 -O3 maxheap8.f sans_module.f x.f Cinxia/sais.o -o x
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
c======================================================================
c
        subroutine test_txt_phr()
        use sanslocaldata
        implicit none
        integer*8 i,j,l

         do i=1,10
          l=ptr(i+1)-ptr(i)-1
          write(*,*) '#TXT',i,l,ptr(i),(char(TXT(ptr(i)+j)),j=1,l)
          l=ptr_phr(i+1)-ptr_phr(i)-1
        write(*,*) '#PHR',i,l,ptr_phr(i),(char(PHR(ptr_phr(i)+j)),j=1,l)
         end do
         do i=nprot-10,nprot-1
          l=ptr(i+1)-ptr(i)-1
          write(*,*) '#TXT',i,l,ptr(i),(char(TXT(ptr(i)+j)),j=1,l)
          l=ptr_phr(i+1)-ptr_phr(i)-1
        write(*,*) '#PHR',i,l,ptr_phr(i),(char(PHR(ptr_phr(i)+j)),j=1,l)
         end do

        end subroutine test_txt_phr     
c
c======================================================================
c
        subroutine test_SAP(lex)
        use sanslocaldata
        implicit none
        integer*8 lex,spos,p
        integer sprot,s

        sprot=SAP(lex)
        s=SRES(lex)
        spos=ptr(sprot)+s
        if(spos.gt.0.and.spos.lt.naa-20)
     $  write(*,*) '#SAP',lex,sprot,s,spos,(char(TXT(p)),p=spos,spos+20)

        end subroutine test_SAP
c
c======================================================================
c
        program sais_mpi
        use sanslocaldata
        use sans
        implicit none
        include 'mpif.h'
        integer rank,nproc,ierr,col,socket,master,SERVER_PORT
        integer*8 lex
        character(len=256) sname,dbversion,programtag
        logical nucl
        integer i

        nucl=.false.
        programtag=
     $  '<program name=SANSparallel version=2.1 citation=PMID:25855811>'
        if(iargc().lt.3) stop 
     $    'USAGE: server database-prefix port-no database-version [NUC]'
        do i=2,iargc()
                if(i.eq.2) then
                        call getarg(2,sname)
                        read(sname,*) SERVER_PORT 
                else if(i.eq.3) then
                        call getarg(3,sname)
                        read(sname,*) dbversion
                else if(i.eq.4) then
                        call getarg(4,sname)
                        nucl=(sname.eq.'NUC')
                        write(*,*) '# Mucleotide data ',nucl
                end if
        end do
        ! start mpi
        call mpi_init(ierr)
        call mpi_comm_size(mpi_comm_world,nproc,ierr)
        call mpi_comm_rank(mpi_comm_world,rank,ierr)

        ! split communicators 
        if(rank.gt.0) then
                col=1
        else
                col=2
        end if
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,col,rank,MPI_COMM_SANS,ierr)
        ! socket sends query/transmitss results
        socket=0
        master=1
        call mpi_comm_size(mpi_comm_world,nproc,ierr)
        call mpi_comm_rank(mpi_comm_world,rank,ierr)
        if(rank.eq.socket) then
                ! match MPI_BARRIER of mpi_comm-sans group
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                write(*,*) '#call socket',socket,nproc
                call socket_loop(SERVER_PORT)
                write(*,*) '#socket closed',socket
        else ! mpi_comm_sans group

                call mpi_comm_size(mpi_comm_sans,nproc,ierr)
                call mpi_comm_rank(mpi_comm_sans,rank,ierr)
                write(*,*) '#I am sans node',rank
                ! load data in MPI_COMM_SANS!
                if(rank.eq.0) call getarg(1,sname)
                call load_data(sname,rank,nproc)
                ! generate local ptr()
                call scatter_txt(sname,rank,nproc)
                ! test local
                call test_txt_phr()
                call scatter_protlen(sname,rank,nproc)
                if(.true.) then
                        do lex=1000000,naa-1,1000000
                                call test_SAP(lex)
                                call test_SAP(lex+1)
                        end do
                end if
                write(*,*) '#data loaded',rank
                ! match MPI_BARRIER of socket
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                call do_all(socket,dbversion,programtag,nucl)
                write(*,*) '#exit do_all loop'
                call deload_data()
                write(*,*) '#data deloaded'
        end if
        call MPI_FINALIZE(ierr)

        end program sais_mpi
c
c======================================================================
c

