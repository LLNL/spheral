      subroutine wrapper(pointer,wrap)
c
      implicit none
c
      integer pointer,wrap,min_p
c
      min_p=min(abs(pointer),abs(pointer-wrap),abs(pointer+wrap))
c
      if(min_p .eq. abs(pointer-wrap)) then
       pointer=pointer-wrap
      else if(min_p .eq. abs(pointer+wrap)) then
       pointer=pointer+wrap
      end if
c
      return
      end
