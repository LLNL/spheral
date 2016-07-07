      logical function check_high(point,array,minimum,
     >  dimension)
c
      implicit none
c
      integer dimension
      integer point,array(dimension),minimum
c
      check_high=.false.
      if(point .gt. 0) check_high=array(point) .ge. minimum
      return
      end

