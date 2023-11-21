c     hgrie Oct 2022: v2.0 fewbody-Compton
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie 17 Nov 2023: split the following subroutines into new file spinstructures.f and renamed two for more intuitive names:

c         singlesigma
c         Calchold => doublesigma
c      
c     This way, spintricks*f only contains individual diagram
c     contributions and not these routines which are generally relevant for spin structures.
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c              
c              CalcCompton2BAxx
c              CalcCompton2BAxy
c              CalcCompton2BAyx
c              CalcCompton2BAyy
c              CalcCompton2BB
c              CalcCompton2BCx
c              CalcCompton2BCy
c              CalcCompton2BDx
c              CalcCompton2BDy
c              CalcCompton2BE
c              CalcCompton2Bfg
c              CalcCompton2Bfg2
c              CalcCompton2Bhi
c              CalcCompton2Bhi2
c              CalcCompton2Bjm
c              CalcCompton2Bno
c              CalcCompton2Bfgfg2ni
c              CalcCompton2Bhini
c              CalcCompton2Bhini2
c              CalcCompton2Baa2hihi2ni
c              CalcCompton2Bd2jmni
c              CalcCompton2Bdnoni
c              
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DRP Feb 2017: check of all factors and extensive commenting. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Aug-Oct 2016/hgrie Feb 2017: Arman added OQ4 diagrams
c====================================================================
c     
      subroutine CalcCompton2BA(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
     &     factor,
     &     Sp,S,verbosity)
c     
c********************************************************************
c     
c     Calculates diagram A
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
c     singlesigma(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
c     calculates 2*S.A, where S=(sigma1+sigma2)/2
c      call singlesigma(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)   
      call singlesigma(hold, 1.d0, 0.d0, 0.d0, factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxx(Sp,Msp,S,Ms)=Comp2Bxx(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εy: mapped to xy
c      call singlesigma(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)   
      call singlesigma(hold, 0.d0, 1.d0, 0.d0, factor,Sp,S,verbosity)
      do Msp=-Sp,Sp
         do Ms=-S,S
            Comp2Bxy(Sp,Msp,S,Ms)=Comp2Bxy(Sp,Msp,S,Ms)+hold(Sp,Msp,S,Ms)
         end do
      end do  
c     εz: mapped to yx
c      call singlesigma(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)   
      call singlesigma(hold, 0.d0, 0.d0, 1.d0, factor,Sp,S,verbosity)
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
      subroutine CalcCompton2BB(Comp2Bxx,Comp2Byx,Comp2Bxy,Comp2Byy,
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
c     singlesigma(hold,Ax,Ay,Az,factor,Sp,S,verbosity)
c     calculates 2*S.A, where S=(sigma1+sigma2)/2
c      call singlesigma(hold,-qppx,-qppy,-qppz,factor,Sp,S,verbosity)   
      call singlesigma(hold, Ax, Ay, Az, factor,Sp,S,verbosity)
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
