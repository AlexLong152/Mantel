cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for One/Twobody Contributions to Few-Nucleon Processes Calculated Via 1N/2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              calcphotonmomenta : set up photon momenta and relative angle theta
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calcphotonmomenta(k,kth,kphi,t,kp,kpth,kpphi,omega,
     &     Q,Qth,Qphi,kgamma,thetacm,verbosity)
c     
c**********************************************************************
c     
c     Sets photon momenta: angle thetacm between the two vectors,
c     with the initial momentum aligned along the z-axis.
c     
c**********************************************************************
c     
      implicit none
      include 'constants.def'
c     
c----------------------------------------------------------------------
c     
c     Output variables:
c     
c     
      real*8,intent(out) :: k,kth,kphi,kp,kpth,kpphi,t,omega
      real*8,intent(out) :: Q,Qth,Qphi
c     
c----------------------------------------------------------------------
c     
c     Input variables:
c     
      real*8,intent(in)  :: kgamma,thetacm
      integer,intent(in) :: verbosity
c     
c----------------------------------------------------------------------
c     
c     Local variables:
c      
      real*8 kx,ky,kz,kpx,kpy,kpz
      real*8 Qx,Qy,Qz
c     
c**********************************************************************
c     
      t=-2.d0*kgamma**2*(1 - dcos(thetacm))
      k=kgamma
      kth=0.d0
      kphi=0.d0
      kx=k*dsin(kth)
      ky=0.d0
      kz=k*dcos(kth)
      kp=kgamma
      kpth=thetacm
      kpphi=0.d0
      kpx=kp*dsin(kpth)
      kpy=0.d0
      kpz=kp*dcos(kpth)
      omega=kgamma
      Qx=kpx-kx
      Qy=kpy-ky
      Qz=kpz-kz
      Q=dsqrt(Qx**2 + Qy**2 + Qz**2)
      Qth=dacos(Qz/Q)
      Qphi=0.d0                 !since both kphi  and kpphi are zero
      if (verbosity.eq.1000) continue
      return
      end
      
