c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie 17 Nov 2023: split the following subroutines into new file spinstructures.f and renamed two for more intuitive names:

c         singlesigmaasy => singlesigmaasy
c         Calcholdasy    => doublesigmaasy
c      
c     This way, spintricks*f only contains individual diagram
c     contributions and not these routines which are generally relevant for spin structures.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     no changes yet
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contains:
C               CalcKernel2BAasy
C               CalcKernel2BBasy
c               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine CalcKernel2BAasy(Kernel2B,
     &     factor,
     &     Sp,S,extQnumlimit,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram A
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
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16,intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
c     
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: factor
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     INTERNAL VARIABLES:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     εx:
      call singlesigmaasy(hold,1.d0,0.d0,0.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(1,Sp,Msp,S,Ms) = Kernel2B(1,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εy: 
      call singlesigmaasy(hold,0.d0,1.d0,0.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(2,Sp,Msp,S,Ms) = Kernel2B(2,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εz:  
      call singlesigmaasy(hold,0.d0,0.d0,1.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(3,Sp,Msp,S,Ms) = Kernel2B(3,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      
      if (verbosity.eq.1000) continue
      return
      end
c====================================================================
c====================================================================
c     
      subroutine CalcKernel2BBasy(Kernel2B,
     &     factor,
     &     Ax,Ay,Az,Bx,By,Bz, ! A.σ, B.ε
     &     Sp,S,extQnumlimit,verbosity)    
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
c     INPUT/OUTPUT VARIABLE:
c     
      complex*16,intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
      
c********************************************************************
c     INPUT VARIABLES:
c     
      real*8,intent(in)  :: factor
      real*8,intent(in)  :: Ax,Ay,Az,Bx,By,Bz
      integer,intent(in) :: Sp,S
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: verbosity
c     
c********************************************************************
c     INTERNAL VARIABLES:
c      
      complex*16 hold(0:1,-1:1,0:1,-1:1)
      integer Msp,Ms
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call singlesigmaasy(hold,Ax,Ay,Az,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
c     εx:
            Kernel2B(1,Sp,Msp,S,Ms) = Kernel2B(1,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)*Bx
c     εy:
            Kernel2B(2,Sp,Msp,S,Ms) = Kernel2B(2,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)*By
c     εz:
            Kernel2B(3,Sp,Msp,S,Ms) = Kernel2B(3,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)*Bz
         end do
      end do  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      
      if (verbosity.eq.1000) continue
      return
      end
