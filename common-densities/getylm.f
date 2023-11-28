cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for One/Twobody Contributions to Few-Nucleon Processes Calculated Via 1N/2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              getsphericalharmonics : Calculates the spherical harmonics Ylm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, identical to file of same name in common-densities/ of Compton density code v2.0 hgrie Oct 2022
c           New documentation -- kept only documentation of changes in Compton if relevant/enlightening for this code. 
c           No back-compatibility 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getsphericalharmonics(Yl,l,theta,phi)
c
c**********************************************************************
c
c     Calculates the spherical harmonics Yl for m from -l to l for a giv. l and a given theta and phi
c    
c**********************************************************************
c
      implicit none
c
c**********************************************************************
c
      include 'params.def'
      include 'constants.def'
c
c**********************************************************************
c   OUTPUT VARIABLES:
c
c      complex*16 Yl(-Lmax:Lmax)
      complex*16,intent(out) :: Yl(-5:5)
c
c Yl(ml)-spherical harmonics of angular momentum l at angle (theta,phi)
c
c----------------------------------------------------------------------
c   INPUT VARIABLES:
c
      integer,intent(in) :: l
      real*8,intent(in)  :: theta,phi
c
c     theta,phi-input solid angle
c     l-angular momentum
c
c----------------------------------------------------------------------
c   LOCAL VARIABLES:
c
      integer ml
c
c     ml-projection of L on z-axis, loop variable to zero Yl array
c
c**********************************************************************
c
      Yl=c0
      if (l.eq.0) then
         Yl(0)=1.d0/dsqrt(4.d0*Pi)
      else if (l.eq.1) then
         Yl(1)=-1.d0/2.d0*dsqrt(3.d0/(2.d0*Pi))*dsin(theta)*exp(ci*phi)
         Yl(0)=1.d0/2.d0*dsqrt(3.d0/Pi)*dcos(theta)
      else if (l.eq.2) then
         Yl(2)=1.d0/4.d0*dsqrt(15.d0/(2.d0*Pi))*dsin(theta)**2*
     &        exp(2.d0*ci*phi)
         Yl(1)=-1.d0/2.d0*dsqrt(15.d0/(2.d0*Pi))*dsin(theta)*
     &                                         dcos(theta)*exp(ci*phi)
        Yl(0)=dsqrt(5.d0/Pi)*(2.d0*dcos(theta)**2 - dsin(theta)**2)/4.d0
      else if (l.eq.3) then
         Yl(3)=-1.d0/8.d0*dsqrt(35.d0/Pi)*dsin(theta)**3*
     &        exp(3.d0*ci*phi)
         Yl(2)=1.d0/4.d0*dsqrt(105.d0/(2.d0*Pi))*dcos(theta)*
     &        dsin(theta)**2*exp(2.d0*ci*phi)
         Yl(1)=-1.d0/8.d0*dsqrt(21.d0/Pi)*exp(ci*phi)*
     &        (4.d0*dcos(theta)**2 - dsin(theta)**2)*dsin(theta)
         Yl(0)=1.d0/4.d0*dsqrt(7.d0/Pi)*dcos(theta)*
     &        (2.d0*dcos(theta)**2-3.d0*dsin(theta)**2)
      else if (l.eq.4) then
         Yl(4)=3.d0/16.d0*dsqrt(35.d0/(2.d0*Pi))*exp(4.d0*ci*phi)*
     &        dsin(theta)**4
         Yl(3)=-3.d0/8.d0*dsqrt(35.d0/Pi)*dsin(theta)**3*dcos(theta)*
     &        exp(3.d0*ci*phi)
        Yl(2)=3.d0/8.d0*dsqrt(5.d0/(2.d0*Pi))*(7.d0*dcos(theta)**2-1.d0)
     &        *dsin(theta)**2*exp(2.d0*ci*phi)
         Yl(1)=-3.d0/8.d0*dsqrt(5.d0/Pi)*exp(ci*phi)*
     &        (7.d0*dcos(theta)**2 - 3.d0)*dsin(theta)*dcos(theta)
         Yl(0)=3.d0/(16.d0*dsqrt(Pi))*(3.d0 - 30.d0*dcos(theta)**2 
     &        + 35.d0*dcos(theta)**4)
      else if (l.eq.5) then
        Yl(5)=-3.d0/32.d0*exp(5.d0*ci*phi)*
     &        dsqrt(77.d0/Pi)*dsin(theta)**5
        Yl(4)=3.d0/16.d0*exp(4.d0*ci*phi)*
     &       dsqrt(385.d0/(2.d0*Pi))*dcos(theta)*dsin(theta)**4
        Yl(3)=-1.d0/32.d0*exp(3.d0*ci*phi)*
     &      dsqrt(385.d0/Pi)*(9.d0*dcos(theta)**2 - 1.d0)*dsin(theta)**3
        Yl(2)=1.d0/8.d0*exp(2.d0*ci*phi)*
     &       dsqrt(1155.d0/(2.d0*Pi))*(3.d0*dcos(theta)**2 - 1.d0)*
     &                                  dsin(theta)**2*dcos(theta)
        Yl(1)=-1.d0/16.d0*exp(ci*phi)*
     &       dsqrt(165.d0/(2.d0*Pi))*(21.d0*dcos(theta)**4 - 
     &                   14.d0*dcos(theta)**2 + 1.d0)*dsin(theta)
        Yl(0)=1.d0/16.d0*dsqrt(11.d0/Pi)*(63.d0*dcos(theta)**5
     &       - 70.d0*dcos(theta)**3 + 15.d0*dcos(theta))
c      else if (l.eq.6) then
c        Yl(2)=1.d0/64.d0*exp(2.d0*ci*phi)*
c     &     dsqrt(1365.d0/Pi)*(33.d0*dcos(theta)**4 - 18.d0*dcos(theta)**2
c     &                                + 1.d0)*dsin(theta)**2
c        Yl(1)=-1.d0/16.d0*exp(ci*phi)*
c     &       dsqrt(273.d0/(2.d0*Pi))*dcos(theta)*(33.d0*dcos(theta)**4 - 
c     &                   30.d0*dcos(theta)**2 + 5.d0)*dsin(theta)
c        Yl(0)=1.d0/32.d0*dsqrt(13.d0/Pi)*(231.d0*dcos(theta)**6
c     &       - 315.d0*dcos(theta)**4 + 105.d0*dcos(theta)**2 - 5.d0)
c      else if (l.eq.7) then
c        Yl(2)=3.d0/64.d0*exp(2.d0*ci*phi)*dcos(theta)*
c     &     dsqrt(35.d0/Pi)*(143.d0*dcos(theta)**4 - 110.d0*dcos(theta)**2
c     &                                + 15.d0)*dsin(theta)**2
c        Yl(1)=-1.d0/64.d0*exp(ci*phi)*dsin(theta)*
c     &       dsqrt(105.d0/(2.d0*Pi))*(429.d0*dcos(theta)**6
c     &            - 495.d0*dcos(theta)**4 + 135.d0*dcos(theta)**2 - 5.d0)
c        Yl(0)=1.0/32.0*dsqrt(15.d0/Pi)*(429.0*dcos(theta)**7 - 
c     &        693.0*dcos(theta)**5 + 315.0*dcos(theta)**3 - 
c     &        35.0*dcos(theta))
      end if
c----------------------------------------------------------------------
c 
c   Use relation of complex conjugate to calculate Yl for -l to -1
c
c----------------------------------------------------------------------
      do ml=1,l
         Yl(-ml)=(-1.d0)**ml*(Real(Yl(ml)) - ci*Imag(Yl(ml)))
      end do
      return
      end
c
c**********************************************************************
      
