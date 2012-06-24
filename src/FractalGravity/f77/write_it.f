      subroutine write_int(file,array,length)
c
      implicit none
c
      integer file,length
      integer array(length)
c
      write(file)length,array
c
c      print*,'integer ',length
      return
      end
c
      subroutine write_real(file,array,length)
c
      implicit none
c
      integer file,length
      real array(length)
c
      write(file)length,array
c
c      print*,'real ',length
      return
      end
c
      subroutine write_log(file,array,length)
c
      implicit none
c
      integer file,length
      logical array(length)
c
      write(file)length,array
c
c      print*,'log ',length
      return
      end
c
