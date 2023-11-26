c     hgrie Oct 2022: v2.0 fewbody-Compton
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie 17 Nov 2023: split the following subroutines into new file spinstructures.f and renamed two for more intuitive names:

c         singlesigma => singlesigmasym
c         Calchold    => doublesigmasym
c      
c     This way, spintricks*f only contains individual diagram
c     contributions and not these routines which are generally relevant for spin structures.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie 25 Sep 2023: correction in subroutine singlesigma()
c     which computes (σ1+σ2).A and enters in Compton at e²δ⁴ (i.e. affects no publication!)
c     Alex Long found a missing sign in line 2420, which should read
c          hold(1,-1,1,-1)= - hold(1,1,1,1) ! minus sign added!!!
c     Rationale and confirmation by Daniel documented in 
c          documentation/corrections/spintrick-singlesigma.signerror-corrected.20230925.pdf
c     which are scans of manuscript Compton Densities Approach pp. 65a-f.
c     These notes show computation of the MEs of (σ1±σ2).A and claim to show that these are
c     now correctly implemented.
c     These notes therefore also clarify that there is NO mistake in spintrickasy.f's
c     corresponding subroutine singlesigmaasy(),
c     and also checked that the implementation in the deuteron's twobodypolnamp.f,
c     subroutine CalcSdotA(), is correct.
c     Fortunately, we have not yet actually RUN anything of the twobody code at e²δ⁴, so
c     the coding mistake does NOT affect published or produced results.
c     Remember that the e²δ⁴ twobody parts were coded by Arman in analogy to the Beame/... 2003
c     paper for the deuteron,e xtended to inlcude the spin-asymmetric piece (σ1-σ2).A .
c     However, we have as of Sep 2023 not yet CHECKED that/if Arman's code is a correct implementation.
c     So this is the first mistake we find in that part. 
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     no changes yet
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contains:
c               CalcKernel2BAsym
c               CalcKernel2BBsym
c              
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DRP Feb 2017: check of all factors and extensive commenting. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Aug-Oct 2016/hgrie Feb 2017: Arman added OQ4 diagrams
c====================================================================
c     
      subroutine CalcKernel2BAsym(Kernel2B,
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
      call singlesigmasym(hold,1.d0,0.d0,0.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(1,Sp,Msp,S,Ms) = Kernel2B(1,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εy:
      call singlesigmasym(hold,0.d0,1.d0,0.d0,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Kernel2B(2,Sp,Msp,S,Ms) = Kernel2B(2,Sp,Msp,S,Ms) + factor*hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εz:
      call singlesigmasym(hold,0.d0,0.d0,1.d0,Sp,S,verbosity)
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
      subroutine CalcKernel2BBsym(Kernel2B,
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
      complex*16, intent(inout) :: Kernel2B(1:extQnumlimit,0:1,-1:1,0:1,-1:1)
c      complex*16 Kernel2Bpx(0:1,-1:1,0:1,-1:1),Kernel2Bpy(0:1,-1:1,0:1,-1:1)
c     
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
      call singlesigmasym(hold,Ax,Ay,Az,Sp,S,verbosity)
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
