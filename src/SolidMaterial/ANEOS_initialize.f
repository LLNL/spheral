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
C call_ANEOS1
C
C Provides a wrapper around the ANEOS1 method to insert stuff in the common 
C block.
C-------------------------------------------------------------------------------
      subroutine call_ANEOS1(T, rho, P, E, S, CV, DPDT, DPDR, zbar, L)

C Common block tomfoolery.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MATBUF=64)
      COMMON /ANESQT/ SQTS(MATBUF),IPSQTS
      COMMON /ANEEL/ TEVX,RHOX,ABARX,ZBARM,T32X,FNX
     &  ,PE,EE,SE,CVE,DPTE,DPRE
     &  ,NMATSX,IIZX

      ipsqts = 1
      sqts(ipsqts) = dsqrt(T)
      call ANEOS1(T, rho, P, E, S, CV, DPDT, DPDR, L)

C We have to dig the atomic weight (ZBARM) out of the common block.  Hope this is right!
      zbar = ZBARM;

      end
