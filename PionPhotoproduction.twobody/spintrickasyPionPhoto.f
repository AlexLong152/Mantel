c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie 17 Nov 2023: split the following subroutines into new file spinstructures.f and renamed two for more intuitive names:

c         singlesigma
c         Calchold => doublesigma
c      
c     This way, spintricks*f only contains individual diagram
c     contributions and not these routines which are generally relevant for spin structures.
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie 17 Nov 2023: split the following subroutines into new file spinstructures.f and renamed two for more intuitive names:

c         singlesigma
c         Calchold => doublesigma
c      
c     This way, spintricks*f only contains individual diagram
c     contributions and not these routines which are generally relevant for spin structures.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     no changes yet
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contains: CalcCompton2BAxxasy()
c               CalcCompton2BAxyasy()
c               CalcCompton2BAyxasy()
c               CalcCompton2BAyyasy()
c               CalcCompton2BBasy()
c               CalcCompton2BCxasy
c               CalcCompton2BCyasy
c               CalcCompton2BDxasy
c               CalcCompton2BDyasy
c               CalcCompton2BEasy
c               CalcCompton2BFgasy
c               CalcCompton2BFg2asy
c               CalcCompton2Bhiasy
c               CalcCompton2Bhi2asy
c               CalcCompton2Bjmasy
c               CalcCompton2Bnoasy
c               CalcCompton2Bfgfg2niasy
c               CalcCompton2Bhiniasy
c               CalcCompton2Bhini2asy
c               CalcCompton2Baa2hihi2niasy
c               CalcCompton2Bd2jmniasy
c               CalcCompton2Bdnoniasy
c               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine CalcCompton2BAasy(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &     factor,
     &     Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram C for x->x and y->x.
c     
c     Note: 1=+,2=-
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16  Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bpx(0:1,-1:1,0:1,-1:1),Comp2Bpy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor,factor12,polnfacx,polnfacy,polnfac12x
      real*8 polnfac12y,k
      real*8 qpppx,qpppy,qppx,qppy,qppz
      real*8 qppp12x,qppp12y,qpp12x,qpp12y,qpp12z
      integer Ms,Msp,Sp,S
      integer verbosity
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     εx: mapped to xx
c      call singlesigmaasy(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)   
      call singlesigmaasy(hold, 1.d0, 0.d0, 0.d0, factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εy: mapped to xy
c      call singlesigmaasy(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)   
      call singlesigmaasy(hold, 0.d0, 1.d0, 0.d0, factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εz: mapped to yx
c      call singlesigmaasy(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)   
      call singlesigmaasy(hold, 0.d0, 0.d0, 1.d0, factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      
      if (verbosity.eq.1000) continue
      return
      end
c====================================================================

c====================================================================
c     
      subroutine CalcCompton2BBasy(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &     factor,
     &     Ax,Ay,Az,Bx,By,Bz, ! A.σ, B.ε
     &     Sp,S,verbosity)    
c     
c********************************************************************
c     
c     Calculates diagram B
c     
c********************************************************************
c     
      implicit none
c     
c********************************************************************
c     
      include '../common-densities/constants.def'
c     
c********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16 Comp2Bxx(0:1,-1:1,0:1,-1:1),Comp2Byx(0:1,-1:1,0:1,-1:1)
      complex*16  Comp2Bxy(0:1,-1:1,0:1,-1:1),Comp2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Comp2Bpx(0:1,-1:1,0:1,-1:1),Comp2Bpy(0:1,-1:1,0:1,-1:1)
      complex*16  hold(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     
c     INPUT VARIABLES:
c     
      real*8 thetacm,factor
      real*8 Ax,Ay,Az,Bx,By,Bz
      integer Ms,Msp,Sp,S
      integer verbosity
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     singlesigmaasy(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
c     calculates 2*S.A, where S=(sigma1-sigma2)/2
c      call singlesigmaasy(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)   
      call singlesigmaasy(hold, Ax, Ay, Az, factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
c     εx: mapped to xx   
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)*Bx
c     εy: mapped to xy   
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)*By
c     εz: mapped to yx   
            Comp2Byx(Sp,Msp,S,Ms)=Comp2Byx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)*Bz
         end do
      end do  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      
      if (verbosity.eq.1000) continue
      return
      end
