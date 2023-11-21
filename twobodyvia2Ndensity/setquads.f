c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     now only contains Setquad12()
c            -- Setquad3() was only for spectator integration, which is absent in density 
c     twoSmax/twoMz dependence: none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Setquad12(th12,Nth12,phi12,Nphi12,
     &     Nordth12,Nthbins12,Nordphi12,Nphibins12,
     &     AngularType12,angweight12,Nanggrid12,verbosity)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set up quadratures for the radial, theta, & phi integrations
c      
c     These routines set up the quadratures for the radial, theta, & phi
c     integrations of the (12) systems.
c     
c     Written by D.P.-11/97
c     hgrie 20 June 2014: modified for use of LebedevLaikov or Gaussian
c     integration for theta & phi separately,
c     for solid angle integral in (12) system 
c     combined wth*wphi*sin(th) (weight of angles theta and phi) into
c     one array angweight12(,) -- so the sum of all weights is 4\pi.
c     
c     hgrie May 2017: split single routine
c      setquad(th3,wth3,Nth3,th12,Nth12,phi12,Nphi12,
c     &     Nordth3,Nthbins3,Nordth12,Nthbins12,Nordphi12,Nphibins12,
c     &     AngularType12,angweight12,Nanggrid12,verbosity)
c      into separate routines for (12)
c
c     setquad12(th12,Nth12,phi12,Nphi12,
c     &     Nordth12,Nthbins12,Nordphi12,Nphibins12,
c     &     AngularType12,angweight12,Nanggrid12,verbosity)
c
c     and spectator (3)
c      
c      setquad3(th3,wth3,Nth3,Nordth3,Nthbins3,verbosity)
c      
c     since spectator (3) integrations not needed in 2N density approach. 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This routine sets up the quadratures for the radial, theta, & phi
c     integrations, as well as those for the numerical zero to one int.
c     used to calculate the invariant functions A1 through A6.
c     
c     Written by D.P.-11/97
c     hgrie 20 June 2014: modified for use of LebedevLaikov or Gaussian
c     integration for theta & phi separately,
c     for solid angle integral in (12) system 
c     combined wth*wphi*sin(th) (weight of angles theta and phi) into
c     one array angweight12(,) -- so the sum of all weights is 4\pi.
c     
c     
      implicit none
      include '../common-densities/params.def'
      include '../common-densities/constants.def'
      include '../common-densities/calctype.def'
c     
c     *******************************************************************
c     
c     OUTPUT VARIABLES:
c     
c     pq,wp-radial quadratures and weights, set up on [0,infty]
c     Np-total number of radial quadratures
c     thq,wth-theta quadratures and weights, set up on [0,PI]
c     Nth-total number of theta quadratures
c     phi12,wphi-phi quadratures and weights, set up on [0,2 PI]
c     Nphi12-total number of phi quadratures
c     
      real*8 th12(Nangmax),wth12(Nangmax)
      real*8 phi12(Nangmax),wphi(Nangmax)
      integer Nth12,Nphi12
c     
c     *******************************************************************
c     
c     INPUT VARIABLES:
c     
c     Nordp,Npbins-number of quadratures/bin and number of bins for
c     radial integration
c     Radialtype-determines nature of mapping from [0,infty] to [0,1]
c     Nordth,Nthbins-number of quadratures/bin and number of bins 
c     for theta integration
c     Nordphi,Nphibins-number of quadratures/bin and number of bins 
c     for phi integration
c     densitytype-used to determine whether we need to cut off the p
c     quadratures
c     hgrie 20 June 2014:
c     AngularType12-determines nature of angular integration in (12) system
c     Nanggrid12: number of points on solid angle grid of Lebedev-Laikov grid
c     angweight12: combined weight of angular integrations, _including_ sin^2(theta)
c     so sum of all weights is 4*Pi
c     angwgth: local variable passed from LebedevLaikov routine 
c     
      integer Nordth12,Nthbins12,Nordphi12,Nphibins12

      integer AngularType12,Nanggrid12
      integer verbosity
      real*8 angweight12(Nangmax,Nangmax)
      real*8 angwgth(Nangmax)

      integer ith,iphi
      
c     
c     *******************************************************************
c     
c     hgrie 20 June 2014: set up  two-dimensional angular integration mesh
c     CalculateIntegralI2 ((12)-integration) needs theta & phi mesh
c     Fill array with ZEROES
      angweight12=0.0E0
      if (AngularType12.eq.1) then ! Gauss-Legendre integration in angles
         call AnglePtsWts(Nordth12,Nthbins12,Nangmax,0.d0,PI,th12,wth12,Nth12,verbosity)
         call AnglePtsWts(Nordphi12,Nphibins12,Nangmax,0.d0,2.0d0*PI,phi12,wphi,Nphi12,verbosity)
         do ith=1,Nth12
            do iphi=1,Nphi12
               angweight12(ith,iphi) = wth12(ith) * wphi(iphi) * dsin(th12(ith))
            end do
         end do
      else if (AngularType12.eq.2) then ! LebedevLaikov integration in angles
         call LebedevLaikovWts(th12,phi12,angwgth,Nanggrid12,verbosity)
         do ith=1,Nanggrid12
            angweight12(ith,ith) = angwgth(ith)
         end do
c     be brutal and just re-define Nth12 to be Nanggrid12
         Nth12 = Nanggrid12
      end if   
c     
      if (verbosity.eq.1) then  ! sum should be 4\pi, for total solid angle
         write(*,*) "Sum of all angular weights (expect 4Pi) = ",SUM(angweight12)/(4*Pi),"*4Pi"
      endif
      if (verbosity.eq.1000) continue
c     end hgrie mod       
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
