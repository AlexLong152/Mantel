c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018: used to be part of 3HeCompt/common
c     now part of common-densities, backward compatibility deliberately broken
c     no changes yet
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Calculateiap12p(p12p,p12pth,p12pphi,p12,
     &          thq,phiq,Q,Qth,Qphi,verbosity)
c**********************************************************************
c
c     Routine which calculates p12p=p - Q/2
c
c**********************************************************************
c
      implicit none
      include 'params.def'
      include 'constants.def'
c
c
c**********************************************************************
c
c  Input variables:
c
      real*8 p12,thq,phiq,Q,Qth,Qphi
      integer verbosity
c
c----------------------------------------------------------------------
c
c  Output variables:
c
      real*8 p12p,p12pth,p12pphi
c
c----------------------------------------------------------------------
c
c  Local variables;
c
      real*8 p12x,p12y,p12z,p12px,p12py,p12pz,Qx,Qy,Qz
      real*8 eps
      parameter (eps=1.0d-8)
c
c**********************************************************************
c
      p12x=p12*dcos(phiq)*dsin(thq)
      p12y=p12*dsin(phiq)*dsin(thq)
      p12z=p12*dcos(thq)
      Qx=Q*dcos(Qphi)*dsin(Qth)
      Qy=Q*dsin(Qphi)*dsin(Qth)
      Qz=Q*dcos(Qth)
      p12px=p12x - 0.5d0*(Qx)
      p12py=p12y - 0.5d0*(Qy)
      p12pz=p12z - 0.5d0*(Qz)
      p12p=dsqrt(p12px**2 + p12py**2 + p12pz**2)
      p12pth=dacos(p12pz/p12p)
      if ((phiq .eq. 0.d0) .and. (Qphi .eq. 0.d0)) then
         p12pphi=0.d0
      else
         p12pphi=dacos(p12px/(p12p*dsin(p12pth)))
      end if 
      if (abs((p12p*dsin(p12pphi)*dsin(p12pth)-p12py)/p12p).gt.eps) then
        p12pphi=-p12pphi
      end if
      if (abs((p12p*dsin(p12pphi)*dsin(p12pth)-p12py)/p12p).gt.eps) then
         write (*,*) 'ERROR: unable to get correct angles'
         write (*,*) p12p*dsin(p12pphi)*dsin(p12pth),p12py
      end if
         if (verbosity.eq.1000) continue
      return
      end
c
c**********************************************************************
      subroutine Calculatep3p(p3p,p3pth,p3pphi,p3,
     &          thq,phiq,Q,Qth,Qphi,verbosity)
c**********************************************************************
c
c     Routine which calculates p12p=p + Q/3
c
c**********************************************************************
c
c
      implicit none
c
c**********************************************************************
c
c  Input variables:
c
      real*8 p3,thq,phiq,Q,Qth,Qphi
      integer verbosity
c
c----------------------------------------------------------------------
c
c  Output variables:
c
      real*8 p3p,p3pth,p3pphi
c
c----------------------------------------------------------------------
c
c  Local variables;
c
      real*8 p3x,p3y,p3z,p3px,p3py,p3pz,Qx,Qy,Qz
      real*8 eps
      parameter (eps=0.00000001)
c
c**********************************************************************
c
      p3x=p3*dcos(phiq)*dsin(thq)
      p3y=p3*dsin(phiq)*dsin(thq)
      p3z=p3*dcos(thq)
      Qx=Q*dcos(Qphi)*dsin(Qth)
      Qy=Q*dsin(Qphi)*dsin(Qth)
      Qz=Q*dcos(Qth)
      p3px=p3x + 1.0d0/3.0d0*(Qx)
      p3py=p3y + 1.0d0/3.0d0*(Qy)
      p3pz=p3z + 1.0d0/3.0d0*(Qz)
      p3p=dsqrt(p3px**2 + p3py**2 + p3pz**2)
      p3pth=dacos(p3pz/p3p)
      if ((phiq.eq.0.d0).and.(Qphi.eq.0.d0)) then
         p3pphi=0.d0
      else
         p3pphi=dacos(p3px/(p3p*dsin(p3pth)))
      end if 
      if (abs((p3p*dsin(p3pphi)*dsin(p3pth)-p3py)/p3p).gt.eps) then
        p3pphi=-p3pphi
      end if
      if (abs((p3p*dsin(p3pphi)*dsin(p3pth)-p3py)/p3p).gt.eps) then
         write (*,*) 'ERROR: unable to get correct angles'
         write (*,*) p3p*dsin(p3pphi)*dsin(p3pth),p3py
      end if
         if (verbosity.eq.1000) continue
      return
      end
      

c**********************************************************************
      subroutine Calculatep3phat(p3pth,p3pphi,p3p,p3,
     &          thq,phiq,Q,Qth,Qphi,verbosity)
c**********************************************************************
c
c     Routine which calculates p12p=p + Q/3
c
c**********************************************************************
c
c
      implicit none
      include 'constants.def'
c
c**********************************************************************
c
c  Input variables:
c
      real*8 p3,thq,phiq,Q,Qth,Qphi,p3p
      integer verbosity
c
c----------------------------------------------------------------------
c
c  Output variables:
c
      real*8 p3pth,p3pphi
c
c----------------------------------------------------------------------
c
c  Local variables;
c
      real*8 p3x,p3y,p3z,p3px,p3py,p3pz,Qx,Qy,Qz
      real*8 eps
      parameter (eps=0.001)
c
c**********************************************************************
c
      p3x=p3*dcos(phiq)*dsin(thq)
      p3y=p3*dsin(phiq)*dsin(thq)
      p3z=p3*dcos(thq)
      Qx=Q*dcos(Qphi)*dsin(Qth)/3.d0
      Qy=Q*dsin(Qphi)*dsin(Qth)/3.d0
      Qz=Q*dcos(Qth)/3.d0
      p3px=p3x+Qx
      p3py=p3y+Qy
      p3pz=p3z+Qz
c
       if(p3p.NE.0) THEN
         p3pth=acos(p3pz/p3p)
       ELSE
         p3pth=0.0
       end if
       if(p3px.GT.0.0) THEN
         p3pphi=atan(p3py/p3px)
       ELSE if(p3px.LT.0.0) THEN
         if(p3py.GT.0.0) THEN
           p3pphi=acos(-1.0)+atan(p3py/p3px)
         ELSE
           p3pphi=-acos(-1.0)+atan(p3py/p3px)
         end if
       ELSE                     ! p13px.EQ.0
         if(p3py.GT.0.0) THEN
           p3pphi=acos(-1.0)/2.0
         ELSE
           p3pphi=-acos(-1.0)/2.0
         end if           
       end if
c
         if (verbosity.eq.1000) continue
      return
      end
      
c**********************************************************************
      subroutine Calculatep12phat(p12pth,p12pphi,p12p,p12,
     &          thq,phiq,Q,Qth,Qphi,verbosity)
c**********************************************************************
c
c     Routine which calculates p12p=p -Q/2
c
c**********************************************************************
c
c
      implicit none
      include 'constants.def'
c
c**********************************************************************
c
c  Input variables:
c
      real*8 p12,thq,phiq,Q,Qth,Qphi,p12p
      integer verbosity
c
c----------------------------------------------------------------------
c
c  Output variables:
c
      real*8 p12pth,p12pphi
c
c----------------------------------------------------------------------
c
c  Local variables;
c
      real*8 p12x,p12y,p12z,p12px,p12py,p12pz,Qx,Qy,Qz
      real*8 eps
      parameter (eps=0.001)
c
c**********************************************************************
c
      p12x=p12*dcos(phiq)*dsin(thq)
      p12y=p12*dsin(phiq)*dsin(thq)
      p12z=p12*dcos(thq)
      Qx=Q*dcos(Qphi)*dsin(Qth)/2.d0
      Qy=Q*dsin(Qphi)*dsin(Qth)/2.d0
      Qz=Q*dcos(Qth)/2.d0
      p12px=p12x-Qx
      p12py=p12y-Qy
      p12pz=p12z-Qz
c      
       if(p12p.NE.0) THEN
         p12pth=acos(p12pz/p12p)
       ELSE
         p12pth=0.0
       end if
       if(p12px.GT.0.0) THEN
         p12pphi=atan(p12py/p12px)
       ELSE if(p12px.LT.0.0) THEN
         if(p12py.GT.0.0) THEN
           p12pphi=acos(-1.0)+atan(p12py/p12px)
         ELSE
           p12pphi=-acos(-1.0)+atan(p12py/p12px)
         end if
       ELSE                     ! p12px.EQ.0
         if(p12py.GT.0.0) THEN
           p12pphi=acos(-1.0)/2.0
         ELSE
           p12pphi=-acos(-1.0)/2.0
         end if           
       end if
c
         if (verbosity.eq.1000) continue
      return
      end
c**********************************************************************
      subroutine Calculatep3(p3,p3th,p3phi,p3mag,
     &          thq,phiq,k,kth,kphi,verbosity)
c**********************************************************************
c
c     Routine which calculates p12p=p + Q/3
c
c**********************************************************************
c
c
      implicit none
c
c**********************************************************************
c
c  Input variables:
c
      real*8 p3mag,thq,phiq,k,kth,kphi
      integer verbosity
c
c----------------------------------------------------------------------
c
c  Output variables:
c
      real*8 p3,p3th,p3phi
c
c----------------------------------------------------------------------
c
c  Local variables;
c
      real*8 p3x,p3y,p3z,p3xold,p3yold,p3zold,kx,ky,kz
      real*8 eps
      parameter (eps=0.00000001)
c
c**********************************************************************
c
      p3xold=p3mag*dcos(phiq)*dsin(thq)
      p3yold=p3mag*dsin(phiq)*dsin(thq)
      p3zold=p3mag*dcos(thq)
      kx=k*dcos(kphi)*dsin(kth)
      ky=k*dsin(kphi)*dsin(kth)
      kz=k*dcos(kth)
      p3x=p3xold + 2.0d0/3.0d0*(kx)
      p3y=p3yold + 2.0d0/3.0d0*(ky)
      p3z=p3zold + 2.0d0/3.0d0*(kz)
      p3=dsqrt(p3x**2 + p3y**2 + p3z**2)
      p3th=dacos(p3z/p3)
      if ((phiq.eq.0.d0).and.(kphi.eq.0.d0)) then
         p3phi=0.d0
      else
         p3phi=dacos(p3x/(p3*dsin(p3th)))
      end if 
      if (abs((p3*dsin(p3phi)*dsin(p3th)-p3y)/p3).gt.eps) then
        p3phi=-p3phi
      end if
      if (abs((p3*dsin(p3phi)*dsin(p3th)-p3y)/p3).gt.eps) then
         write (*,*) 'ERROR: unable to get correct angles'
         write (*,*) p3*dsin(p3phi)*dsin(p3th),p3y
      end if
      if (verbosity.eq.1000) continue
      return
      end
c**********************************************************************
c
      subroutine calculatepvector(qx,qy,qz,q,thq,phiq,verbosity)
c
c  Calcualte Cartesian components given spherical polar co-ordinates
c     
c**********************************************************************
c
      implicit none
c
c**********************************************************************
c
c  Input variables:
c
      real*8 qx,qy,qz
      integer verbosity
c
c**********************************************************************
c
c  Output variables:
c
      real*8 q,thq,phiq
c
c**********************************************************************
c
      qx=q*dsin(thq)*dcos(phiq)
      qy=q*dsin(thq)*dsin(phiq)
      qz=q*dcos(thq)
      if (verbosity.eq.1000) continue
      return
      end
c
c**********************************************************************
c
      subroutine calculateqs(qx,qy,qz,q12x,q12y,q12z,qpx,qpy,qpz,
     &     qp12x,qp12y,qp12z,qppx,qppy,qppz,qpp12x,qpp12y,qpp12z,
     &     qpppx,qpppy,qpppz,qppp12x,qppp12y,qppp12z,
     &     qsq,qpsq,qppsq,qpppsq,q12sq,qp12sq,qpp12sq,qppp12sq,px,py,pz,
     &     ppx,ppy,ppz,
     &     k,thetacm,verbosity)
c
      implicit none
c
c**********************************************************************
c
c  Input variables:
c
      real*8 px,py,pz,ppx,ppy,ppz,k,thetacm
      integer verbosity
c
c**********************************************************************
c
c  Output variables:
c
      real*8 qx,qy,qz,qpx,qpy,qpz,qppx,qppy,qppz,qsq,qpsq,qppsq
      real*8 q12x,q12y,q12z,qp12x,qp12y,qp12z,qpp12x,qpp12y,qpp12z
      real*8 q12sq,qp12sq,qpp12sq
      real*8 qpppx,qpppy,qpppz, qpppsq
      real*8 qppp12x,qppp12y,qppp12z, qppp12sq
c
c**********************************************************************      
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
