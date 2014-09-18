!!$      program teos
!!$      implicit none
!!$      save
!!$      include 'vector_eos.dek'
!!$      include 'table_data.dek'
!!$
!!$!
!!$! tests the helmholtz eos routine
!!$!
!!$! ionmax  = number of isotopes in the network
!!$! xmass   = mass fractions
!!$! ymass   = molar fractions
!!$! aion    = number of nucleons
!!$! zion    = number of protons
!!$
!!$      integer          ionmax, itermax, i, j, k
!!$      parameter        (ionmax=3)
!!$      double precision xmass(ionmax),ymass(ionmax), &
!!$                       aion(ionmax),zion(ionmax),temp,den,abar,zbar
!!$      double precision Tn, Tnplusone, fn, dfdxn, etarget, errtol
!!$!      double precision temp_row(nrowmax), den_row(nrowmax),
!!$!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)
!!$
!!$! set the mass fractions, z's and a's of the composition
!!$! hydrogen
!!$      xmass(1) = 0.75d0
!!$      aion(1)  = 1.0d0
!!$      zion(1)  = 1.0d0
!!$
!!$! helium
!!$      xmass(2) = 0.23d0
!!$      aion(2)  = 4.0d0
!!$      zion(2)  = 2.0d0
!!$
!!$! carbon 12
!!$      xmass(3) = 0.02d0
!!$      aion(3)  = 12.0d0
!!$      zion(3)  = 6.0d0
!!$
!!$
!!$! get abar, zbar and a few other composition variables
!!$      call azbar(xmass,aion,zion,ionmax, &
!!$                 ymass,abar,zbar)
!!$
!!$      abar=13.7143
!!$      zbar=6.85714
!!$!   i     radius            mass              density           energy            temp
!!$! 16   .9599699824E+04  2.4848507176E-11  3.4999995804E+07  2.9321170671E+17  1.0000000000E+07  1.37
!!$
!!$
!!$! set the input vector
!!$      temp_row(1) = 1.e7
!!$      den_row(1)  = 3.56598997E+06
!!$      abar_row(1) = abar
!!$      zbar_row(1) = zbar
!!$
!!$! here the pipeline is only 1 element long
!!$      jlo_eos = 1
!!$      jhi_eos = 1
!!$
!!$
!!$      call init_helm_table
!!$
!!$
!!$! call the eos
!!$      call helmeos
!!$!(temp_row, den_row, etot_row, abar_row, zbar_row)
!!$
!!$
!!$! write out the results
!!$      call pretty_eos_out('helm:  ')
!!$!, temp_row, den_row, etot_row,
!!$!     1     abar_row, zbar_row)
!!$
!!$
!!$!    NOW USE THIS KNOWN TEST TO SEE WHETHER NEWTON RAPHSON CAN GET THE
!!$!    INPUT T BACK. As far as I understand, the specific energy is
!!$!    stored in the variable "ener" inside helmeos, which then gets
!!$!    stored in "etot_row(j)".
!!$
!!$!    The Newton Raphson scheme is a very simple and fast converging
!!$!    scheme, based on the fact that starting from a a guess xn, the
!!$!    value xnplusone=xn-f(xn)/dfdx(xn) should be closer to the root of
!!$!    f(x).
!!$
!!$!    In this case, I want to find T, such that u(T)=u,
!!$!    i.e. f(T)=u-u(T). Thus, for this scheme, I will also need
!!$!    dfdT=-dudT, which is also provided by the code, and stored in
!!$!    det_row(j).
!!$
!!$!    Initial guess for T. Since we know the input, let's just use
!!$!    something close
!!$!      do k=1,10000
!!$
!!$
!!$
!!$!      j=1
!!$!      errtol=1e-6
!!$!      Tnplusone=temp_row(j)*1.25
!!$!      etarget=etot_row(j)
!!$!
!!$!      itermax=100
!!$!      do i=1,itermax
!!$!    New iteration, replace n+1 by n
!!$!         Tn=Tnplusone
!!$!         temp_row(j)=Tn
!!$!         call helmeos
!!$!            call pretty_eos_out('helm:  ')
!!$!, temp_row, den_row, etot_row,
!!$!     1           abar_row, zbar_row)
!!$
!!$
!!$!(temp_row, den_row, etot_row, abar_row, zbar_row)
!!$!    Reset f and df/dx
!!$!         fn=etarget-etot_row(j)
!!$!         dfdxn=-det_row(j)
!!$!    Need an if statement here to prevent dividing by 0
!!$!         Tnplusone=Tn-fn/dfdxn
!!$!         write(6,*), Tnplusone, tnplusone-tn, errtol*tn
!!$!         if (abs(Tnplusone-Tn) .lt. errtol*Tn) then
!!$!
!!$!            exit
!!$!         end if
!!$!      enddo
!!$!      enddo
!!$
!!$
!!$
!!$      stop 'normal termination'
!!$      end






      subroutine azbar(xmass,aion,zion,ionmax, &
                       ymass,abar,zbar)
      implicit none
      save

! this routine calculates composition variables for an eos routine

! input:
! mass fractions     = xmass(1:ionmax)
! number of nucleons = aion(1:ionmax)
! charge of nucleus  = zion(1:ionmax)
! number of isotopes = ionmax

! output:
! molar abundances        = ymass(1:ionmax),
! mean number of nucleons = abar
! mean nucleon charge     = zbar


! declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax), &
                       ymass(ionmax),abar,zbar,zbarxx,ytot1

      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
      end







      subroutine pretty_eos_out
!(whose, temp_row, den_row, etot_row,
!     1     abar_row, zbar_row)
      implicit none
      save
      include 'vector_eos.dek'
!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)
!
! writes a pretty output for the eos tester
!
! declare
      integer     j
      character*7 whose

! popular formats
01    format(1x,t2,a,t11,'total',t24,'ion',t34,'e- & e+', &
             t46,'radiation',t58,'coulomb')
02    format(1x,t2,a,1p6e12.4)
03    format(1x,t2,a6,1pe12.4,t22,a6,1pe12.4, &
               t42,a6,1pe12.4,t62,a6,1pe12.4)



      do j=jlo_eos,jhi_eos


! the input
      write(6,03) 'temp =',temp_row(1),'den  =',den_row(1), &
                  'abar =',abar_row(1),'zbar =',zbar_row(1)
      write(6,*) ' '


! and the output
! first the totals from each of the components
      write(6,01)  whose
      write(6,02) 'pres =', &
                  ptot_row(j),pion_row(j),pele_row(j), &
                  prad_row(j),pcou_row(j)
      write(6,02) 'ener =', &
                  etot_row(j),eion_row(j),eele_row(j), &
                  erad_row(j),ecou_row(j)
      write(6,02) 'entr =', &
                  stot_row(j),sion_row(j),sele_row(j), &
                  srad_row(j),scou_row(j)

! derivatives of the totals with respect to the input variables
      write(6,*)  ' '
      write(6,03) 'dp/dd=',dpd_row(j),'dp/dt=',dpt_row(j), &
                  'dp/da=',dpa_row(j),'dp/dz=',dpz_row(j)
      write(6,03) 'de/dd=',ded_row(j),'de/dt=',det_row(j), &
                  'de/da=',dea_row(j),'de/dz=',dez_row(j)
      write(6,03) 'ds/dd=',dsd_row(j),'ds/dt=',dst_row(j), &
                  'ds/da=',dsa_row(j),'ds/dz=',dsz_row(j)


! derivatives of the electron-positron compoenets with
! respect to the input variables
      write(6,*) ' '
      write(6,03) 'dpepd=',dpepd_row(j),'dpept=',dpept_row(j), &
                  'dpepa=',dpepa_row(j),'dpepz=',dpepz_row(j)
      write(6,03) 'deepd=',deepd_row(j),'deept=',deept_row(j), &
                  'deepa=',deepa_row(j),'deepz=',deepz_row(j)
      write(6,03) 'dsepd=',dsepd_row(j),'dsept=',dsept_row(j), &
                  'dsepa=',dsepa_row(j),'dsepz=',dsepz_row(j)


! the thermodynamic consistency relations, these should all be
! at the floating poiint limit of zero
      write(6,*) ' '
      write(6,03) 'maxw1=',dse_row(j),'maxw2=',dpe_row(j), &
                  'maxw3=',dsp_row(j)


! number density of electrons, poistrons, matter electrons, and ions
      write(6,03) 'xne  =',xne_row(j),'xnp  =',xnp_row(j), &
                  'xnem =',xnem_row(j),'xni  =',xni_row(j)


! derivatibves of the electron number density with
! respect to the input variables
      write(6,03) 'dxned=',dxned_row(j),'dxnet=',dxnet_row(j), &
                  'dxnea=',dxnea_row(j),'dxnez=',dxnez_row(j)


! electron chemical potential, positron chemical potential
! and derivatives of electron chemical potential with respect
! to the input variables
      write(6,03) 'eta  =',etaele_row(j),'etap =',etapos_row(j)
      write(6,03) 'detad=',detad_row(j),'detat=',detat_row(j), &
                  'detaa=',detaa_row(j),'detaz=',detaz_row(j)


! specific heats, and ratio of electostatic to thermal energy
      write(6,03) 'cp   =',cp_row(j),'cv   =',cv_row(j), &
                  'plasg=',plasg_row(j)

! the 3 gammas and the sound speed
      write(6,03) 'gam1 =',gam1_row(j),'gam2 =',gam2_row(j), &
                  'gam3 =',gam3_row(j),'csond=',cs_row(j)
      write(6,*) ' '

      enddo
      return
      end


      subroutine get_helm_table(of, ofd, oft, ofdd, oftt, ofdt, ofddt, &
       ofdtt, ofddtt, odpdf, odpdfd, odpdft, odpdfdt, oef, oefd, oeft, &
       oefdt, oxf, oxfd, oxft, oxfdt)

      implicit none
      save
      include 'vector_eos.dek'
      include 'table_data.dek'

! for the helmholtz free energy tables
      double precision of(imax,jmax),ofd(imax,jmax), &
                       oft(imax,jmax),ofdd(imax,jmax),oftt(imax,jmax), &
                       ofdt(imax,jmax),ofddt(imax,jmax), &
                       ofdtt(imax,jmax),ofddtt(imax,jmax)

! for the pressure derivative with density ables
      double precision odpdf(imax,jmax),odpdfd(imax,jmax), &
                       odpdft(imax,jmax),odpdfdd(imax,jmax), &
                       odpdftt(imax,jmax),odpdfdt(imax,jmax)

! for chemical potential tables
      double precision oef(imax,jmax),oefd(imax,jmax), &
                       oeft(imax,jmax),oefdd(imax,jmax), &
                       oeftt(imax,jmax),oefdt(imax,jmax)

! for the number density tables
      double precision oxf(imax,jmax),oxfd(imax,jmax), &
                       oxft(imax,jmax),oxfdd(imax,jmax), &
                       oxftt(imax,jmax),oxfdt(imax,jmax)
      double precision outf
      integer i,j

!    read the helmholtz free energy table
      do j=1,jmax
         do i=1,imax
            of(i,j)      = f(i,j)
            ofd(i,j)     = fd(i,j)
            oft(i,j)     = ft(i,j)
            ofdd(i,j)    = fdd(i,j)
            oftt(i,j)    = ftt(i,j)
            ofdt(i,j)    = fdt(i,j)
            ofddt(i,j)   = fddt(i,j)
            ofdtt(i,j)   = fdtt(i,j)
            ofddtt(i,j)  = fddtt(i,j)
            odpdf(i,j)   = dpdf(i,j)
            odpdfd(i,j)  = dpdfd(i,j)
            odpdft(i,j)  = dpdft(i,j)
            odpdfdt(i,j) = dpdfdt(i,j)
            oef(i,j)     = ef(i,j)
            oefd(i,j)    = efd(i,j)
            oeft(i,j)    = eft(i,j)
            oefdt(i,j)   = efdt(i,j)
            oxf(i,j)     = xf(i,j)
            oxfd(i,j)    = xfd(i,j)
            oxft(i,j)    = xft(i,j)
            oxfdt(i,j)   = xfdt(i,j)
         enddo
      enddo

      end


      subroutine set_helm_table(of, ofd, oft, ofdd, oftt, ofdt, ofddt, &
       ofdtt, ofddtt, odpdf, odpdfd, odpdft, odpdfdt, oef, oefd, oeft, &
       oefdt, oxf, oxfd, oxft, oxfdt)

      implicit none
      save
      include 'vector_eos.dek'
      include 'table_data.dek'

! for the helmholtz free energy tables
      double precision of(imax,jmax),ofd(imax,jmax), &
                       oft(imax,jmax),ofdd(imax,jmax),oftt(imax,jmax), &
                       ofdt(imax,jmax),ofddt(imax,jmax), &
                       ofdtt(imax,jmax),ofddtt(imax,jmax)

! for the pressure derivative with density ables
      double precision odpdf(imax,jmax),odpdfd(imax,jmax), &
                       odpdft(imax,jmax),odpdfdd(imax,jmax), &
                       odpdftt(imax,jmax),odpdfdt(imax,jmax)

! for chemical potential tables
      double precision oef(imax,jmax),oefd(imax,jmax), &
                       oeft(imax,jmax),oefdd(imax,jmax), &
                       oeftt(imax,jmax),oefdt(imax,jmax)

! for the number density tables
      double precision oxf(imax,jmax),oxfd(imax,jmax), &
                       oxft(imax,jmax),oxfdd(imax,jmax), &
                       oxftt(imax,jmax),oxfdt(imax,jmax)
      double precision outf
      integer i,j

!    read the helmholtz free energy table
      do j=1,jmax
         do i=1,imax
            f(i,j)     = of(i,j)
            fd(i,j)    = ofd(i,j)
            ft(i,j)    = oft(i,j)
            fdd(i,j)   = ofdd(i,j)
            ftt(i,j)   = oftt(i,j)
            fdt(i,j)   = ofdt(i,j)
            fddt(i,j)  = ofddt(i,j)
            fdtt(i,j)  = ofdtt(i,j)
            fddtt(i,j) = ofddtt(i,j)
            dpdf(i,j)  = odpdf(i,j)
            dpdfd(i,j) = odpdfd(i,j)
            dpdft(i,j) = odpdft(i,j)
            dpdfdt(i,j)= odpdfdt(i,j)
            ef(i,j)    = oef(i,j)
            efd(i,j)   = oefd(i,j)
            eft(i,j)   = oeft(i,j)
            efdt(i,j)  = oefdt(i,j)
            xf(i,j)    = oxf(i,j)
            xfd(i,j)   = oxfd(i,j)
            xft(i,j)   = oxft(i,j)
            xfdt(i,j)  = oxfdt(i,j)
         enddo
      enddo

      end


      subroutine init_helm_table
      implicit none
      save
      include 'vector_eos.dek'
      include 'table_data.dek'
      integer i, j


!    open the table
      open(unit=2,file='helm_table.dat',status='old')

!    read the helmholtz free energy table
      do j=1,jmax
         do i=1,imax
            read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j), &
                 fddt(i,j),fdtt(i,j),fddtt(i,j)
         enddo
      enddo

!    read the pressure derivative with density table
      do j=1,jmax
         do i=1,imax
            read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
         enddo
      enddo

!    read the electron chemical potential table
      do j=1,jmax
         do i=1,imax
            read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
         enddo
      enddo

!    read the number density table
      do j=1,jmax
         do i=1,imax
            read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
         enddo
      enddo

      close(unit=2)

      end





      subroutine helmeos
!(temp_row, den_row, etot_row, abar_row,
!     1     zbar_row)
      implicit none
      save
      include 'vector_eos.dek'
      include 'table_data.dek'

!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)

! given a temperature temp [K], density den [g/cm**3], and a composition
! characterized by abar and zbar, this routine returns most of the other
! thermodynamic quantities. of prime interest is the pressure [erg/cm**3],
! specific thermal energy [erg/gr], the entropy [erg/g/K], along with
! their derivatives with respect to temperature, density, abar, and zbar.
! other quantites such the normalized chemical potential eta (plus its
! derivatives), number density of electrons and positron pair (along
! with their derivatives), adiabatic indices, specific heats, and
! relativistically correct sound speed are also returned.
!
! this routine assumes planckian photons, an ideal gas of ions,
! and an electron-positron gas with an arbitrary degree of relativity
! and degeneracy. interpolation in a table of the helmholtz free energy
! is used to return the electron-positron thermodynamic quantities.
! all other derivatives are analytic.
!
! references: cox & giuli chapter 24 ; timmes & swesty apj 1999


! declare
      double precision pi,amu,kerg,clight,avo,qe,h,ssol,asol
      parameter       (pi      = 3.1415926535897932384d0, &
                        amu    = 1.6605402d-24, &
                        kerg   = 1.380658d-16, &
                        clight = 2.99792458d10, &
                        avo    = 6.0221367d23, &
                        qe     = 4.8032068d-10, &
                        h      = 6.6260755d-27, &
                        ssol   = 5.67051d-5, &
                        asol   = 4.0d0 * ssol / clight)

      integer          i,j
      double precision x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida, &
                       dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt, &
                       dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt, &
                       deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt, &
                       dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion, &
                       sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd, &
                       dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp, &
                       gam1,gam2,gam3,chit,chid,nabad,sound,etaele, &
                       detadt,detadd,xnefer,dxnedt,dxnedd,s, &
                       temp,den,abar,zbar,ytot1,ye, &
                       sioncon,forth,forpi,kergavo,ikavo,asoli3,light2

      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h), &
                        forth   = 4.0d0/3.0d0, &
                        forpi   = 4.0d0 * pi, &
                        kergavo = kerg * avo, &
                        ikavo   = 1.0d0/kergavo, &
                        asoli3  = asol/3.0d0, &
                        light2  = clight * clight)

! for the abar derivatives
      double precision dpradda,deradda,dsradda, &
                       dpionda,deionda,dsionda, &
                       dpepda,deepda,dsepda, &
                       dpresda,denerda,dentrda, &
                       detada,dxneda


! for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz, &
                       dpiondz,deiondz,dsiondz, &
                       dpepdz,deepdz,dsepdz, &
                       dpresdz,denerdz,dentrdz, &
                       detadz,dxnedz

! for the interpolations
      integer          iat,jat
      double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi, &
                       tsav,dsav,free,df_d,df_t,df_dd,df_tt,df_dt
      double precision dth,dt2,dti,dt2i,dd,dd2,ddi,dd2i,xt,xd,mxt,mxd, &
                       si0t,si1t,si2t,si0mt,si1mt,si2mt, &
                       si0d,si1d,si2d,si0md,si1md,si2md, &
                       dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, &
                       dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, &
                       ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, &
                       ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md, &
                       z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, &
                       dpsi2,ddpsi2,din,h5,fi(36), &
                       xpsi0,xdpsi0,xpsi1,xdpsi1,h3, &
                       w0t,w1t,w2t,w0mt,w1mt,w2mt, &
                       w0d,w1d,w2d,w0md,w1md,w2md

! for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax), &
                       dti_sav(jmax),dt2i_sav(jmax), &
                       dd_sav(imax),dd2_sav(imax), &
                       ddi_sav(imax),dd2i_sav(imax)


! for the coulomb corrections
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd, &
                       plasg,plasgdd,plasgdt,plasgda,plasgdz, &
                       a1,b1,c1,d1,e1,a2,b2,c2, &
                       ecoul,decouldd,decouldt,decoulda,decouldz, &
                       pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                       scoul,dscouldd,dscouldt,dscoulda,dscouldz, &
                       tmelt,tfermi,rhocond,z2,x1,x2,third,esqu
      parameter        (a1    = -0.898004d0, &
                        b1    =  0.96786d0, &
                        c1    =  0.220703d0, &
                        d1    = -0.86097d0, &
                        e1    =  2.5269d0, &
                        a2    =  0.29561d0, &
                        b2    =  1.9885d0, &
                        c2    =  0.288675d0, &
                        third = 1.0d0/3.0d0, &
                        esqu  = qe * qe)


! for initialization
      integer          ifirst
      data             ifirst/0/


! quintic hermite polynomial statement functions
! psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)

! psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)

! psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)

! biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= &
             fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt &
           + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t &
           + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt &
           + fi(13) *w1d*w0t   + fi(14) *w1md*w0t &
           + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt &
           + fi(17) *w2d*w0t   + fi(18) *w2md*w0t &
           + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt &
           + fi(21) *w1d*w1t   + fi(22) *w1md*w1t &
           + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt &
           + fi(25) *w2d*w1t   + fi(26) *w2md*w1t &
           + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt &
           + fi(29) *w1d*w2t   + fi(30) *w1md*w2t &
           + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt &
           + fi(33) *w2d*w2t   + fi(34) *w2md*w2t &
           + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



! cubic hermite polynomial statement functions
! psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)

! psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


! bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = &
             fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t &
           + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt &
           + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t &
           + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt &
           + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t &
           + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt &
           + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t &
           + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



! popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))


! do this stuff once
      if (ifirst .eq. 0) then
       ifirst = 1

! open the table
       call init_helm_table

! read the helmholtz free energy table
       tlo   = 3.0d0
       thi   = 13.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -12.0d0
       dhi   = 15.0d0
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
        enddo
       enddo

! construct the temperature and density deltas and their inverses
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
       enddo


      end if



! start of vectorization loop, normal execution starts here
      eosfail = .false.
      do j=jlo_eos,jhi_eos

!       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
!       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = ytot1 * zbar

! initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


! radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad*tempi
       dpradda = 0.0d0
       dpraddz = 0.0d0

       erad    = 3.0d0 * prad*deni
       deraddd = -erad*deni
       deraddt = 3.0d0 * dpraddt*deni
       deradda = 0.0d0
       deraddz = 0.0d0


       srad    = (prad*deni + erad)*tempi
       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
       dsradda = 0.0d0
       dsraddz = 0.0d0

! ion section:
       xni     = avo * ytot1 * den
       dxnidd  = avo * ytot1
       dxnida  = -xni * ytot1

       pion    = xni * kt
       dpiondd = dxnidd * kt
       dpiondt = xni * kerg
       dpionda = dxnida * kt
       dpiondz = 0.0d0

       eion    = 1.5d0 * pion*deni
       deiondd = (1.5d0 * dpiondd - eion)*deni
       deiondt = 1.5d0 * dpiondt*deni
       deionda = 1.5d0 * dpionda*deni
       deiondz = 0.0d0


       x       = abar*abar*sqrt(abar) * deni/avo
       s       = sioncon * temp
       z       = x * s * sqrt(s)
       y       = log(z)
       sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
       dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
                  - kergavo * deni * ytot1
       dsiondt = (dpiondt*deni + deiondt)*tempi - &
                 (pion*deni + eion) * tempi*tempi &
                 + 1.5d0 * kergavo * tempi*ytot1
       x       = avo*kerg/abar
       dsionda = (dpionda*deni + deionda)*tempi &
                 + kergavo*ytot1*ytot1* (2.5d0 - y)
       dsiondz = 0.0d0


! electron-positron section:
! assume complete ionization
       xnem    = xni * zbar

! enter the table with ye*den
       din = ye*den

! bomb proof the input
       if (temp .gt. t(jmax)) then
        write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
        write(6,*) 'temp too hot, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (temp .lt. t(1)) then
        !write(6,01) 'temp=',temp,' t(1)=',t(1)
        !write(6,*) 'temp too cold, off grid'
        !write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .gt. d(imax)) then
        write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
        write(6,*) 'ye*den too big, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .lt. d(1)) then
        write(6,01) 'ye*den=',din,' d(1)=',d(1)
        write(6,*) 'ye*den too small, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if

! hash locate this temperature and density
       jat = int((log10(temp) - tlo)*tstpi) + 1
       jat = max(1,min(jat,jmax-1))
       iat = int((log10(din) - dlo)*dstpi) + 1
       iat = max(1,min(iat,imax-1))


! access the table locations only once
       fi(1)  = f(iat,jat)
       fi(2)  = f(iat+1,jat)
       fi(3)  = f(iat,jat+1)
       fi(4)  = f(iat+1,jat+1)
       fi(5)  = ft(iat,jat)
       fi(6)  = ft(iat+1,jat)
       fi(7)  = ft(iat,jat+1)
       fi(8)  = ft(iat+1,jat+1)
       fi(9)  = ftt(iat,jat)
       fi(10) = ftt(iat+1,jat)
       fi(11) = ftt(iat,jat+1)
       fi(12) = ftt(iat+1,jat+1)
       fi(13) = fd(iat,jat)
       fi(14) = fd(iat+1,jat)
       fi(15) = fd(iat,jat+1)
       fi(16) = fd(iat+1,jat+1)
       fi(17) = fdd(iat,jat)
       fi(18) = fdd(iat+1,jat)
       fi(19) = fdd(iat,jat+1)
       fi(20) = fdd(iat+1,jat+1)
       fi(21) = fdt(iat,jat)
       fi(22) = fdt(iat+1,jat)
       fi(23) = fdt(iat,jat+1)
       fi(24) = fdt(iat+1,jat+1)
       fi(25) = fddt(iat,jat)
       fi(26) = fddt(iat+1,jat)
       fi(27) = fddt(iat,jat+1)
       fi(28) = fddt(iat+1,jat+1)
       fi(29) = fdtt(iat,jat)
       fi(30) = fdtt(iat+1,jat)
       fi(31) = fdtt(iat,jat+1)
       fi(32) = fdtt(iat+1,jat+1)
       fi(33) = fddtt(iat,jat)
       fi(34) = fddtt(iat+1,jat)
       fi(35) = fddtt(iat,jat+1)
       fi(36) = fddtt(iat+1,jat+1)


! various differences
       xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
       xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
       mxt = 1.0d0 - xt
       mxd = 1.0d0 - xd

! the six density and six temperature basis functions
       si0t =   psi0(xt)
       si1t =   psi1(xt)*dt_sav(jat)
       si2t =   psi2(xt)*dt2_sav(jat)

       si0mt =  psi0(mxt)
       si1mt = -psi1(mxt)*dt_sav(jat)
       si2mt =  psi2(mxt)*dt2_sav(jat)

       si0d =   psi0(xd)
       si1d =   psi1(xd)*dd_sav(iat)
       si2d =   psi2(xd)*dd2_sav(iat)

       si0md =  psi0(mxd)
       si1md = -psi1(mxd)*dd_sav(iat)
       si2md =  psi2(mxd)*dd2_sav(iat)

! derivatives of the weight functions
       dsi0t =   dpsi0(xt)*dti_sav(jat)
       dsi1t =   dpsi1(xt)
       dsi2t =   dpsi2(xt)*dt_sav(jat)

       dsi0mt = -dpsi0(mxt)*dti_sav(jat)
       dsi1mt =  dpsi1(mxt)
       dsi2mt = -dpsi2(mxt)*dt_sav(jat)

       dsi0d =   dpsi0(xd)*ddi_sav(iat)
       dsi1d =   dpsi1(xd)
       dsi2d =   dpsi2(xd)*dd_sav(iat)

       dsi0md = -dpsi0(mxd)*ddi_sav(iat)
       dsi1md =  dpsi1(mxd)
       dsi2md = -dpsi2(mxd)*dd_sav(iat)

! second derivatives of the weight functions
       ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
       ddsi1t =   ddpsi1(xt)*dti_sav(jat)
       ddsi2t =   ddpsi2(xt)

       ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
       ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
       ddsi2mt =  ddpsi2(mxt)

!       ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
!       ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
!       ddsi2d =   ddpsi2(xd)

!       ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
!       ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
!       ddsi2md =  ddpsi2(mxd)


! the free energy
       free  = h5(iat,jat, &
               si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
               si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density
       df_d  = h5(iat,jat, &
               si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, &
               dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

! derivative with respect to temperature
       df_t = h5(iat,jat, &
               dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
               si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to density**2
!       df_dd = h5(iat,jat,
!     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
!     2         ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

! derivative with respect to temperature**2
       df_tt = h5(iat,jat, &
             ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, &
               si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

! derivative with respect to temperature and density
       df_dt = h5(iat,jat, &
               dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, &
               dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



! now get the pressure derivative with density, chemical potential, and
! electron positron number densities
! get the interpolation weight functions
       si0t   =  xpsi0(xt)
       si1t   =  xpsi1(xt)*dt_sav(jat)

       si0mt  =  xpsi0(mxt)
       si1mt  =  -xpsi1(mxt)*dt_sav(jat)

       si0d   =  xpsi0(xd)
       si1d   =  xpsi1(xd)*dd_sav(iat)

       si0md  =  xpsi0(mxd)
       si1md  =  -xpsi1(mxd)*dd_sav(iat)


! derivatives of weight functions
       dsi0t  = xdpsi0(xt)*dti_sav(jat)
       dsi1t  = xdpsi1(xt)

       dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
       dsi1mt = xdpsi1(mxt)

       dsi0d  = xdpsi0(xd)*ddi_sav(iat)
       dsi1d  = xdpsi1(xd)

       dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
       dsi1md = xdpsi1(mxd)


! look in the pressure derivative only once
       fi(1)  = dpdf(iat,jat)
       fi(2)  = dpdf(iat+1,jat)
       fi(3)  = dpdf(iat,jat+1)
       fi(4)  = dpdf(iat+1,jat+1)
       fi(5)  = dpdft(iat,jat)
       fi(6)  = dpdft(iat+1,jat)
       fi(7)  = dpdft(iat,jat+1)
       fi(8)  = dpdft(iat+1,jat+1)
       fi(9)  = dpdfd(iat,jat)
       fi(10) = dpdfd(iat+1,jat)
       fi(11) = dpdfd(iat,jat+1)
       fi(12) = dpdfd(iat+1,jat+1)
       fi(13) = dpdfdt(iat,jat)
       fi(14) = dpdfdt(iat+1,jat)
       fi(15) = dpdfdt(iat,jat+1)
       fi(16) = dpdfdt(iat+1,jat+1)

! pressure derivative with density
       dpepdd  = h3(iat,jat, &
                       si0t,   si1t,   si0mt,   si1mt, &
                       si0d,   si1d,   si0md,   si1md)
       dpepdd  = max(ye * dpepdd,0.0d0)



! look in the electron chemical potential table only once
       fi(1)  = ef(iat,jat)
       fi(2)  = ef(iat+1,jat)
       fi(3)  = ef(iat,jat+1)
       fi(4)  = ef(iat+1,jat+1)
       fi(5)  = eft(iat,jat)
       fi(6)  = eft(iat+1,jat)
       fi(7)  = eft(iat,jat+1)
       fi(8)  = eft(iat+1,jat+1)
       fi(9)  = efd(iat,jat)
       fi(10) = efd(iat+1,jat)
       fi(11) = efd(iat,jat+1)
       fi(12) = efd(iat+1,jat+1)
       fi(13) = efdt(iat,jat)
       fi(14) = efdt(iat+1,jat)
       fi(15) = efdt(iat,jat+1)
       fi(16) = efdt(iat+1,jat+1)


! electron chemical potential etaele
       etaele  = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to density
       x       = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
       detadd  = ye * x

! derivative with respect to temperature
       detadt  = h3(iat,jat, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to abar and zbar
      detada = -x * din * ytot1
      detadz =  x * den * ytot1



! look in the number density table only once
       fi(1)  = xf(iat,jat)
       fi(2)  = xf(iat+1,jat)
       fi(3)  = xf(iat,jat+1)
       fi(4)  = xf(iat+1,jat+1)
       fi(5)  = xft(iat,jat)
       fi(6)  = xft(iat+1,jat)
       fi(7)  = xft(iat,jat+1)
       fi(8)  = xft(iat+1,jat+1)
       fi(9)  = xfd(iat,jat)
       fi(10) = xfd(iat+1,jat)
       fi(11) = xfd(iat,jat+1)
       fi(12) = xfd(iat+1,jat+1)
       fi(13) = xfdt(iat,jat)
       fi(14) = xfdt(iat+1,jat)
       fi(15) = xfdt(iat,jat+1)
       fi(16) = xfdt(iat+1,jat+1)

! electron + positron number densities
      xnefer   = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to density
      x        = h3(iat,jat, &
                     si0t,   si1t,   si0mt,   si1mt, &
                    dsi0d,  dsi1d,  dsi0md,  dsi1md)
      x = max(x,0.0d0)
      dxnedd   = ye * x

! derivative with respect to temperature
      dxnedt   = h3(iat,jat, &
                    dsi0t,  dsi1t,  dsi0mt,  dsi1mt, &
                     si0d,   si1d,   si0md,   si1md)

! derivative with respect to abar and zbar
      dxneda = -x * din * ytot1
      dxnedz =  x  * den * ytot1



! the desired electron-positron thermodynamic quantities

! dpepdd at high temperatures and low densities is below the
! floating point limit of the subtraction of two large terms.
! since dpresdd doesn't enter the maxwell relations at all, use the
! bicubic interpolation done above instead of this one
       x       = din * din
       pele    = x * df_d
       dpepdt  = x * df_dt
!       dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
       s       = dpepdd/ye - 2.0d0 * din * df_d
       dpepda  = -ytot1 * (2.0d0 * pele + s * din)
       dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


       x       = ye * ye
       sele    = -df_t * ye
       dsepdt  = -df_tt * ye
       dsepdd  = -df_dt * x
       dsepda  = ytot1 * (ye * df_dt * din - sele)
       dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


       eele    = ye*free + temp * sele
       deepdt  = temp * dsepdt
       deepdd  = x * df_d + temp * dsepdd
       deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
       deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




! coulomb section:
! initialize


        pcoul    = 0.0d0
        dpcouldd = 0.0d0
        dpcouldt = 0.0d0
        dpcoulda = 0.0d0
        dpcouldz = 0.0d0
        ecoul    = 0.0d0
        decouldd = 0.0d0
        decouldt = 0.0d0
        decoulda = 0.0d0
        decouldz = 0.0d0
        scoul    = 0.0d0
        dscouldd = 0.0d0
        dscouldt = 0.0d0
        dscoulda = 0.0d0
        dscouldz = 0.0d0


! uniform background corrections only
! from yakovlev & shalybkov 1989
! lami is the average ion seperation
! plasg is the plasma coupling parameter
        z        = forth * pi
        s        = z * xni
        dsdd     = z * dxnidd
        dsda     = z * dxnida

        lami     = 1.0d0/s**third
        inv_lami = 1.0d0/lami
        z        = -third * lami
        lamidd   = z * dsdd/s
        lamida   = z * dsda/s

        plasg    = zbar*zbar*esqu*ktinv*inv_lami
        z        = -plasg * inv_lami
        plasgdd  = z * lamidd
        plasgda  = z * lamida
        plasgdt  = -plasg*ktinv * kerg
        plasgdz  = 2.0d0 * plasg/zbar


! yakovlev & shalybkov 1989 equations 82, 85, 86, 87
         if (plasg .ge. 1.0) then
          x        = plasg**(0.25d0)
          y        = avo * ytot1 * kerg
          ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
          pcoul    = third * den * ecoul
          scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
                    + d1 * (log(plasg) - 1.0d0) - e1)

          y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
          decouldd = y * plasgdd
          decouldt = y * plasgdt + ecoul/temp
          decoulda = y * plasgda - ecoul/abar
          decouldz = y * plasgdz

          y        = third * den
          dpcouldd = third * ecoul + y*decouldd
          dpcouldt = y * decouldt
          dpcoulda = y * decoulda
          dpcouldz = y * decouldz


          y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
          dscouldd = y * plasgdd
          dscouldt = y * plasgdt
          dscoulda = y * plasgda - scoul/abar
          dscouldz = y * plasgdz


! yakovlev & shalybkov 1989 equations 102, 103, 104
         else if (plasg .lt. 1.0) then
          x        = plasg*sqrt(plasg)
          y        = plasg**b2
          z        = c2 * x - third * a2 * y
          pcoul    = -pion * z
          ecoul    = 3.0d0 * pcoul/den
          scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

          s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
          dpcouldd = -dpiondd*z - pion*s*plasgdd
          dpcouldt = -dpiondt*z - pion*s*plasgdt
          dpcoulda = -dpionda*z - pion*s*plasgda
          dpcouldz = -dpiondz*z - pion*s*plasgdz

          s        = 3.0d0/den
          decouldd = s * dpcouldd - ecoul/den
          decouldt = s * dpcouldt
          decoulda = s * dpcoulda
          decouldz = s * dpcouldz

          s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
          dscouldd = s * plasgdd
          dscouldt = s * plasgdt
          dscoulda = s * plasgda - scoul/abar
          dscouldz = s * plasgdz
         end if



! bomb proof
        x   = prad + pion + pele + pcoul
        y   = erad + eion + eele + ecoul
        z   = srad + sion + sele + scoul
!        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then

! .    SPECIAL SETTING FOR COLD WHITE DWARF:
!       DON'T TURN COULOMB CORRECTIONS OFF!
        if (x .le. 0.0 .or. y .le. 0.0) then

!         write(6,*)
!         write(6,*) 'coulomb corrections are causing a negative pressure'
!         write(6,*) 'setting all coulomb corrections to zero'
!         write(6,*)

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if





! sum all the components
       pres    = prad + pion + pele + pcoul
       ener    = erad + eion + eele + ecoul
       entr    = srad + sion + sele + scoul

       dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd
       dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
       dpresda = dpradda + dpionda + dpepda + dpcoulda
       dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz

       denerdd = deraddd + deiondd + deepdd + decouldd
       denerdt = deraddt + deiondt + deepdt + decouldt
       denerda = deradda + deionda + deepda + decoulda
       denerdz = deraddz + deiondz + deepdz + decouldz

       dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
       dentrdt = dsraddt + dsiondt + dsepdt + dscouldt
       dentrda = dsradda + dsionda + dsepda + dscoulda
       dentrdz = dsraddz + dsiondz + dsepdz + dscouldz


! the temperature and density exponents (c&g 9.81 9.82)
! the specific heat at constant volume (c&g 9.92)
! the third adiabatic exponent (c&g 9.93)
! the first adiabatic exponent (c&g 9.97)
! the second adiabatic exponent (c&g 9.105)
! the specific heat at constant pressure (c&g 9.98)
! and relativistic formula for the sound speed (c&g 14.29)
       zz    = pres*deni
       zzi   = den/pres
       chit  = temp/pres * dpresdt
       chid  = dpresdd*zzi
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + light2)*zzi
       sound = clight * sqrt(gam1/z)


! maxwell relations; each is zero if the consistency is perfect
       x   = den * den
       dse = temp*dentrdt/denerdt - 1.0d0
       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
       dsp = -dentrdd*x/dpresdt - 1.0d0


! store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd
        dpa_row(j)    = dpresda
        dpz_row(j)    = dpresdz

        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd
        dea_row(j)    = denerda
        dez_row(j)    = denerdz

        stot_row(j)   = entr
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd
        dsa_row(j)    = dentrda
        dsz_row(j)    = dentrdz

        prad_row(j)   = prad
        erad_row(j)   = erad
        srad_row(j)   = srad

        pion_row(j)   = pion
        eion_row(j)   = eion
        sion_row(j)   = sion
        xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda
        dpepz_row(j)  = dpepdz

        eele_row(j)   = eele
        epos_row(j)   = 0.0d0
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda
        deepz_row(j)  = deepdz

        sele_row(j)   = sele
        spos_row(j)   = 0.0d0
        dsept_row(j)  = dsepdt
        dsepd_row(j)  = dsepdd
        dsepa_row(j)  = dsepda
        dsepz_row(j)  = dsepdz

        xnem_row(j)   = xnem
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxnedt
        dxned_row(j)  = dxnedd
        dxnea_row(j)  = dxneda
        dxnez_row(j)  = dxnedz
        xnp_row(j)    = 0.0d0

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        etapos_row(j) = 0.0d0

        pcou_row(j)   = pcoul
        ecou_row(j)   = ecoul
        scou_row(j)   = scoul
        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        cs_row(j)     = sound

! end of vectorization loop
      enddo


!      call pretty_eos_out('inside helmeos:  ')
!, temp_row, den_row, etot_row,
!     1     abar_row, zbar_row)
      return
      end




      subroutine wrapper_helmeos(npart, in_den_row, &
       in_etot_row, in_abar_row, in_zbar_row, in_temp_row, &
       in_ptot_row)
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'
      double precision &
                in_den_row(nrowmax),in_etot_row(nrowmax), &
                in_abar_row(nrowmax),in_zbar_row(nrowmax), &
                in_temp_row(nrowmax), in_ptot_row(nrowmax)
!      double precision temp_row(nrowmax), den_row(nrowmax),
!     1     etot_row(nrowmax), abar_row(nrowmax), zbar_row(nrowmax)
      integer npart, i, j

!    Set the pipeline limits from 1 to npart
      jlo_eos=1
      jhi_eos=npart
      if (jhi_eos .gt. nrowmax) then
         write(6,*) 'Too many particles supplied, adjust nrowmax!!!'
      end if

!    Copy the local variables into the common block
      do i=1,npart
         den_row(i)=in_den_row(i)
         etot_row(i)=in_etot_row(i)
         abar_row(i)=in_abar_row(i)
         zbar_row(i)=in_zbar_row(i)
         temp_row(i)=in_temp_row(i)
      enddo

      call helmeos
!(temp_row, den_row, etot_row, abar_row, zbar_row)
!    Copy back into the energy array output

!      call pretty_eos_out('helm:  ')
!, temp_row, den_row, etot_row,
!     1     abar_row, zbar_row)

      do i=1,npart
         in_etot_row(i)=etot_row(i)
         in_ptot_row(i)=ptot_row(i)
      enddo

      end
