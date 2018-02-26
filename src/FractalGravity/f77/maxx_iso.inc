      integer iso_maxx,iso,green_length
      parameter (iso=0)
      parameter (iso_maxx=iso*
     > ((2*grid_length_max)**2*(2*grid_length_max+2)-1)+1)
      parameter (green_length=
     > ((grid_length_max+1)**3-1)*iso+1)
