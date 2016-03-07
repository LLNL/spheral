! this file contains routines to invert the helmholtz eos
!
! routine invert_helm_pt is used when the pressure and temperature are given
! routine invert_helm_pt_quiet is as above, but supresses all error messages
! routine invert_helm_pd is used when the pressure and density are given
! routine invert_helm_et is used when the energy and temperature are given
! routine invert_helm_ed is used when the energy and density are given
! routine invert_helm_st is used when the entropy and temperature are given
! routine invert_helm_st_quiet is as above, but supresses all error messages





      subroutine invert_helm_pt
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)


! given the pressure, temperature, and composition
! find everything else

! it is assumed that ptot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))
      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,400

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpd_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_pt'
      write(6,*) 'pipeline element',j
      write(6,01) 'pwant  =',eoswrk03(j),' temp =',temp_row(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)
      stop 'could not find a density in routine invert_helm_pt'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      return
      end







      subroutine invert_helm_pt_quiet
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)


! given the pressure, temperature, and composition
! find everything else

! it is assumed that ptot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.

! this version is quiet on all errors


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpd_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
!      write(6,*)
!      write(6,*) 'newton-raphson failed in routine invert_helm_pt'
!      write(6,*) 'pipeline element',j
!      write(6,01) 'pres  =',eoswrk03(j)
! 01   format(1x,5(a,1pe16.8))
!      write(6,01) 'error =',eoswrk01(j),
!     1            '  eostol=',eostol,'  fpmin =',fpmin
!      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
!      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
!      write(6,*)
!      stop 'could not find a density in routine invert_helm_pt'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      return
      end







      subroutine invert_helm_pd
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)


! given the pressure, density, and composition
! find everything else

! it is assumed that ptot_row(j), den_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input temp_row(j) conatins a guess for the temperature,
! on output temp_row(j) contains the converged temperature.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision tmpold,tmp,f,df,tmpnew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = temp_row(j)
      end do


! do the first newton loop with all elements in the pipe
      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpt_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       tmp    = temp_row(j)
       tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
       eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
       temp_row(j)  = min(1.0d14,max(tmpnew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpt_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        tmp    = temp_row(j)
        tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
        eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
        temp_row(j)  = min(1.0d13,max(tmpnew,1.0d3))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_pd'
      write(6,*) 'pipeline element',j
      write(6,01) 'pwant  =',eoswrk03(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'tmp   =',temp_row(j),'  tmpold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)
      stop 'could not find a temperature in routine invert_helm_pd'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      return
      end






      subroutine invert_helm_pd_quiet
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)


! given the pressure, density, and composition
! find everything else

! it is assumed that ptot_row(j), den_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input temp_row(j) conatins a guess for the temperature,
! on output temp_row(j) contains the converged temperature.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision tmpold,tmp,f,df,tmpnew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = temp_row(j)
      end do


! do the first newton loop with all elements in the pipe
      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpt_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       tmp    = temp_row(j)
       tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
       eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
       temp_row(j)  = min(1.0d14,max(tmpnew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpt_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        tmp    = temp_row(j)
        tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
        eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
        temp_row(j)  = min(1.0d13,max(tmpnew,1.0d3))

! end of netwon loop
       end do


! we did not converge if we land here
!      write(6,*)
!      write(6,*) 'newton-raphson failed in routine invert_helm_pd'
!      write(6,*) 'pipeline element',j
!      write(6,01) 'pwant  =',eoswrk03(j)
! 01   format(1x,5(a,1pe16.8))
!      write(6,01) 'error =',eoswrk01(j),
!     1            '  eostol=',eostol,'  fpmin =',fpmin
!      write(6,01) 'tmp   =',temp_row(j),'  tmpold=',eoswrk04(j)
!      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
!      write(6,*)
!      stop 'could not find a temperature in routine invert_helm_pd'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      return
      end






      subroutine invert_helm_et
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)


! given the specific internal energy, temperature, and composition,
! find everything else

! it is assumed that etot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = etot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      do j = jlo_eos, jhi_eos

       f     = etot_row(j)/eoswrk03(j) - 1.0d0
       df    = ded_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

        f     = etot_row(j)/eoswrk03(j) - 1.0d0
        df    = ded_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_et'
      write(6,*) 'pipeline element',j
      write(6,01) 'ewant  =',eoswrk03(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)


      stop 'could not find a density in routine invert_helm_et'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)

      return
      end





      subroutine invert_helm_ed
!, temp_row, den_row, etot_row,
!     1     abar_row, zbar_row)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)
      integer npart

! given the specific internal energy density, density, and composition
! find everything else

! it is assumed that etot_row(j), den_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input temp_row(j) conatins a guess for the temperature,
! on output temp_row(j) contains the converged temperature.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision tmpold,tmp,f,df,tmpnew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)

! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = etot_row(j)
       eoswrk04(j) = temp_row(j)
      end do


! do the first newton loop with all elements in the pipe

      call helmeos

      do j = jlo_eos, jhi_eos

       f     = etot_row(j)/eoswrk03(j) - 1.0d0
       df    = det_row(j)/eoswrk03(j)

       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       tmp    = temp_row(j)
       tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
       eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new temperature, keep it within the table limits
! and make sure we have at least 1 iteration in Newton loop
       temp_row(j)  = min(1.0d13,max(tmpnew,1.000001*small_temp))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20


        if (temp_row(j) .le. small_temp ) then
           temp_row(j)=small_temp
           goto 20
        end if

        jlo_eos = j
        jhi_eos = j

        call helmeos

        f     = etot_row(j)/eoswrk03(j) - 1.0d0
        df    = det_row(j)/eoswrk03(j)

        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        tmp    = temp_row(j)
        tmpnew = min(max(0.5d0*tmp,tmp - eoswrk02(j)),2.0d0*tmp)

! compute the error
        eoswrk01(j)  = abs((tmpnew - tmp)/tmp)

! store the new density, keep it within the table limits
        temp_row(j)  = min(1.0d13,max(tmpnew,small_temp))

! end of newton loop
       end do

! we did not converge if we land here
!      write(6,*)
!      write(6,*) 'newton-raphson failed in routine invert_helm_ed'
!      write(6,*) 'pipeline element',j
!      write(6,01) 'ewant  =',eoswrk03(j)
! 01   format(1x,5(a,1pe16.8))
!      write(6,01) 'error =',eoswrk01(j),
!     1            '  eostol=',eostol,'  fpmin =',fpmin
!      write(6,01) 'tmp   =',temp_row(j),'  tmpold=',eoswrk04(j),
!     1     'den=', den_row(j)
!      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
!      write(6,*)
!
!
!      call pretty_eos_out('helm:  ')
!      temp_row(j)=-1.

!      stop 'could not find a temperature in routine invert_helm_ed'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos

      return
      end





      subroutine invert_helm_st
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)


! given the entropy, temperature, and composition
! find everything else

! it is assumed that stot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.

! this version is quiet on all errors


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = stot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
      call helmeos

      do j = jlo_eos, jhi_eos

       f     = stot_row(j)/eoswrk03(j) - 1.0d0
       df    = dsd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos

        f     = stot_row(j)/eoswrk03(j) - 1.0d0
        df    = dsd_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
      write(6,*)
      write(6,*) 'newton-raphson failed in routine invert_helm_st'
      write(6,*) 'pipeline element',j
      write(6,01) 'entr  =',eoswrk03(j)
 01   format(1x,5(a,1pe16.8))
      write(6,01) 'error =',eoswrk01(j), &
                  '  eostol=',eostol,'  fpmin =',fpmin
      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
      write(6,*)
      stop 'could not find a density in routine invert_helm_st'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos

      return
      end



      subroutine invert_helm_st_quiet
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'

! given the entropy, temperature, and composition
! find everything else

! it is assumed that stot_row(j), temp_row(j), abar_row(j),
! zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have
! been set before calling this routine.

! on input den_row(j) conatins a guess for the density,
! on output den_row(j) contains the converged density.

! To get the greatest speed advantage, the eos should be fed a
! large pipe of data to work on.

! this version is quiet on all errors


! local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8, &
                        fpmin  = 1.0d-14)


! initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = stot_row(j)
       eoswrk04(j) = den_row(j)
      end do


! do the first newton loop with all elements in the pipe
      call helmeos

      do j = jlo_eos, jhi_eos

       f     = stot_row(j)/eoswrk03(j) - 1.0d0
       df    = dsd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

! limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

      enddo



! now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. &
            abs(eoswrk02(j)) .le. fpmin) goto 20

        jlo_eos = j
        jhi_eos = j

        call helmeos

        f     = stot_row(j)/eoswrk03(j) - 1.0d0
        df    = dsd_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

! limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

! compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

! store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

! end of netwon loop
       end do


! we did not converge if we land here
!      write(6,*)
!      write(6,*) 'newton-raphson failed in routine invert_helm_pt'
!      write(6,*) 'pipeline element',j
!      write(6,01) 'pres  =',eoswrk03(j)
! 01   format(1x,5(a,1pe16.8))
!      write(6,01) 'error =',eoswrk01(j),
!     1            '  eostol=',eostol,'  fpmin =',fpmin
!      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
!      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
!      write(6,*)
!      stop 'could not find a density in routine invert_helm_pt'



! land here if newton loop converged, back for another pipe element
 20    continue
      end do



! call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos

      return
      end






      subroutine wrapper_invert_helm_ed(npart, in_den_row, &
       in_etot_row, in_abar_row, in_zbar_row, in_temp_row, in_ptot_row, &
       in_small_temp, in_vsound_row, in_gam1_row, in_stot_row)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      double precision &
                in_den_row(nrowmax),  in_etot_row(nrowmax), &
                in_abar_row(nrowmax), in_zbar_row(nrowmax), &
                in_temp_row(nrowmax), in_ptot_row(nrowmax), & 
		in_vsound_row(nrowmax), in_gam1_row(nrowmax), &
                in_stot_row(nrowmax)
      double precision in_small_temp
      integer npart, i, j

!    Set the pipeline limits from 1 to npart
      jlo_eos=1
      jhi_eos=npart
      if (jhi_eos .gt. nrowmax) then
         write(6,*) 'Too many particles supplied, adjust nrowmax!!!'
      end if

!    Set the minimum temperature to avoid numerical trouble
      small_temp=max(1d3,in_small_temp)
!      small_temp=in_small_temp

!    Copy the local variables into the common block
      do i=1,npart
         den_row(i)=in_den_row(i)
         etot_row(i)=in_etot_row(i)
         abar_row(i)=in_abar_row(i)
         zbar_row(i)=in_zbar_row(i)
         temp_row(i)=in_temp_row(i)
         gam1_row(i)=in_gam1_row(i)
         stot_row(i)=in_stot_row(i)
      enddo

      call invert_helm_ed

!    Copy back into the temperature and pressure arrays
!      call pretty_eos_out('helm:  ')
      do i=1,npart
         in_temp_row(i)=temp_row(i)
         in_ptot_row(i)=ptot_row(i)
		 in_vsound_row(i)=cs_row(i)
         in_gam1_row(i)=gam1_row(i)
         in_stot_row(i)=stot_row(i)
      enddo

      end


