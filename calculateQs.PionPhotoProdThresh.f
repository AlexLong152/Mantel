cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of KERNEL code for Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              CalculateQs      : compute vectors needed by kernel/diagrams
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO DO:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     v1.0 Nov 2023: New, loosely based on subroutine which was part of common-densities/calcmomenta.f of Compton density code v2.0 hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     COMMENTS:
c     specific to process considered
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine CalculateQs(qx,qy,qz,q12x,q12y,q12z,qpx,qpy,qpz,
     &     qp12x,qp12y,qp12z,qppx,qppy,qppz,qpp12x,qpp12y,qpp12z,
     &     qpppx,qpppy,qpppz,qppp12x,qppp12y,qppp12z,
     &     qsq,qpsq,qppsq,qpppsq,q12sq,qp12sq,qpp12sq,qppp12sq,px,py,pz,
     &     ppx,ppy,ppz,
     &     k,thetacm,verbosity)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:
c
      real*8,intent(in)  :: px,py,pz,ppx,ppy,ppz,k,thetacm
      integer,intent(in) :: verbosity
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:
c
      real*8,intent(out) :: qx,qy,qz,qpx,qpy,qpz,qppx,qppy,qppz,qsq,qpsq,qppsq
      real*8,intent(out) :: q12x,q12y,q12z,qp12x,qp12y,qp12z,qpp12x,qpp12y,qpp12z
      real*8,intent(out) :: q12sq,qp12sq,qpp12sq
      real*8,intent(out) :: qpppx,qpppy,qpppz, qpppsq
      real*8,intent(out) :: qppp12x,qppp12y,qppp12z, qppp12sq
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES: none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      qx=px - ppx + k*dsin(thetacm)/2.d0
      qy=py - ppy
      qz=pz - ppz + k*(1.d0 + dcos(thetacm))/2.d0
      qpppx=px + ppx - k*dsin(thetacm)/2.d0
      qpppy=py + ppy
      qpppz=pz + ppz - k*(1.d0 + 3.d0*dcos(thetacm))/6.d0
      qppp12x=-px - ppx - k*dsin(thetacm)/2.d0
      qppp12y=-py - ppy
      qppp12z=-pz - ppz - k*(1.d0 + 3.d0*dcos(thetacm))/6.d0
      q12x=-px + ppx + k*dsin(thetacm)/2.d0
      q12y=-py + ppy
      q12z=-pz + ppz + k*(1.d0 + dcos(thetacm))/2.d0
      qpx=px - ppx + k*dsin(thetacm)/2.d0
      qpy=py - ppy
      qpz=pz - ppz + k*(dcos(thetacm)- 1.d0)/2.d0
      qppx=px - ppx - k*dsin(thetacm)/2.d0
      qppy=py - ppy 
      qppz=pz - ppz + k*(1.d0 - dcos(thetacm))/2.d0
      qp12x=-px + ppx + k*dsin(thetacm)/2.d0
      qp12y=-py + ppy
      qp12z=-pz + ppz + k*(dcos(thetacm) - 1.d0)/2.d0
      qpp12x=-px + ppx - k*dsin(thetacm)/2.d0
      qpp12y=-py + ppy
      qpp12z=-pz + ppz + k*(1.d0 - dcos(thetacm))/2.d0
      qsq=qx**2 + qy**2 + qz**2
      q12sq=q12x**2 + q12y**2 + q12z**2
      qpsq=qpx**2 + qpy**2 + qpz**2
      qppsq=qppx**2 + qppy**2 + qppz**2
      qp12sq=qp12x**2 + qp12y**2 + qp12z**2
      qpp12sq=qpp12x**2 + qpp12y**2 + qpp12z**2
      qpppsq=qpppx**2 + qpppy**2 + qpppz**2
      qppp12sq=qppp12x**2 + qppp12y**2 + qppp12z**2


      
      if (verbosity.eq.1000) continue
      return
      end
