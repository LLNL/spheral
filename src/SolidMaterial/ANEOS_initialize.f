C-------------------------------------------------------------------------------
C Spheral's initialization methods for the ANEOS EOS package.
C 
C Note I'm trying to reverse engineer ANEOS in all of this, so take it with a 
C caveat!
C
C Arguments:
C    filename : name of the ANEOS input file, normally something like ANEOS.INPUT
C    num      : number of materials 
C    izetl    : array of num EOS numbers

C-------------------------------------------------------------------------------
      subroutine ANEOS_initialize(in_filename, out_filename, num, izetl)

      implicit none

C ANEOS common blocks
      integer KLST, KINP
      COMMON /FILEOS/ KLST, KINP

C Declare the arguments.
      character*256 in_filename, out_filename, lineinp
      integer num, izetl(21)

C Open the files.
      open(10, file=in_filename, status='old')
      open(12, file=out_filename)
C      open(12, file=out_filename, status='new')
      kinp = 10
      klst = 12

C Call the main ANEOS initialization method.
      call ANEOS2(1, num, 0, izetl)

      end

C-------------------------------------------------------------------------------
C call_ANEOS
C
C Provides a wrapper around the ANEOS method used to compute the EOS response.
C-------------------------------------------------------------------------------
      subroutine call_ANEOS(matnum, T, rho, P, e, s, cv, dpdt, dpdr)

      implicit none

      integer matnum
      real*8 T, rho, P, e, s, cv, dpdt, dpdr

c$$$      integer kpa
c$$$      real*8 fkro, cs
c$$$      call ANEOS(T, rho, P, e, s, cv, dpdt, dpdr, fkro, cs, kpa, matnum)
c$$$
c$$$      print *, " --> ", T, rho, P, e, s, cv, dpdt, dpdr, fkro, cs, kpa

C     the following contains sqrt(t(i)) for vector aneos entry
C     ipsqts points to current value
      integer IPSQTS,MATBUF
      PARAMETER (MATBUF=64)
      real*8 sqts(MATBUF)
      COMMON /ANESQT/ SQTS,IPSQTS

      ipsqts = 1
      sqts(ipsqts) = dsqrt(T)
      call ANEOS1(T, rho, P, e, s, cv, dpdt, dpdr, matnum)

      end

C-------------------------------------------------------------------------------
C get_ANEOS_atomicWeight
C
C Use the stored ANEOS data to compute the atomic weight of the given material.
C This has to be called after the initialization method above.
C
C I've cribbed and reimplemented the stuff in ANEOS2 for DIN(29) since it doesn't
C seem to be otherwise accessible.
C-------------------------------------------------------------------------------
c$$$      subroutine get_ANEOS_atomicWeight(matnum)
c$$$
c$$$      implicit none
c$$$
c$$$C Crap straight from ANEOS
c$$$      PARAMETER (MAXMAT=10)
c$$$      PARAMETER(NINPUT=48) !NINPUT must be a multiple of 8!
c$$$      COMMON /ANES/  ACK(99*MAXMAT),ZZS(30*MAXMAT),COT(30*MAXMAT)
c$$$     1 ,FNI(30*MAXMAT),RCT(MAXMAT+1),TCT(MAXMAT+1),RSOL(100*MAXMAT)
c$$$     2 ,RVAP(100*MAXMAT),TTWO(100*MAXMAT),SAVER(92),BOLTS,EIP(4370)
c$$$     3 ,LOCSV(MAXMAT+1),LOCKP(MAXMAT+1),LOCKPL(MAXMAT+1),NOTRAD
c$$$
c$$$      result = 0.0
c$$$      do 100 i = iz, izi
c$$$         IKK = ZZS(I)
c$$$         IKJ = IKK+(IKK*(IKK+1))/2
c$$$         TEIP = EIP(IKJ-IKK)           !Atomic weight of species I
c$$$         result = result + COT(I)*TEIP !mean atomic weight
c$$$ 100  continue
