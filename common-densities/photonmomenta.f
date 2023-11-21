c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018: used to be part of 3HeCompt/common
c     now part of common-densities, backward compatibility deliberately broken
c     no changes yet
c
c**********************************************************************
c
      subroutine calcphotonmomenta(k,kth,kphi,t,kp,kpth,kpphi,omega,
     &     Q,Qth,Qphi,kgamma,thetacm,verbosity)
c     
c**********************************************************************
c     
c     Sets photon momenta for Chi PT calculation: in other words
c     it just sets up an angle thetacm between the two vectors,
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
      
