      logical function test_big(x,i_1,i_2,i_3,big)
c     
      implicit none
c     
      integer i,i_1,i_2,i_3
      real x(i_2),big
c     
      test_big=.false.
      do i=i_1,i_2,i_3
       test_big=abs(x(i)) .gt. big
       if(test_big) return
      end do
c     
      return
      end
