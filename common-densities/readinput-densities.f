c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018, based on common/readinput.f
c
c     subroutines to read input files for density calculations:
c
c     ReadinputCommon: input for parameters which are identical for onebody and twobody
c     ReadinputOnebody: input for parameters for onebody
c     ReadinputTwobody: input for parameters for twobody
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     hgrie Sep 2020: in readinputTwobody():
c        -- Removed NP12p = 3rd variable in integration-grid
c           line of input file. It is actually never used for anything!
c           That also reduced number of arguments of readinputTwobody()!
c        -- Changed read of angle parameters such that input file does not need to
c           contain 5 numbers when only 2 are needed for LL
c     hgrie Aug 2020: added readinputCommonComments()
c     hgrie Aug 2020: corected couting of independent MEs:
c      now ...+2*mod(twoSnucl+1,2) for correct number -- was 2*(mod(twoSnucl,2)-1)
c      
c     hgrie June 2018: renamed "parity" to "symmetry
c     -- see notes in usesymmetry+writeoutput.densities.f
c      
c     hgrie May 2018: used to be part of 3HeCompt/common
c     now part of common-densities, backward compatibility deliberately broken
c     no changes yet
c      
c     hgrie May 2018: more extensive description of changes in main.*.f
c                     rewritten such that magnetic quantum numbers are "2xMz" etc
c      
c     -------------------------------------------------------------------
c     LEGACY: Relevant notes on original common/redinput.f:
c     this combined the 3 subroutines (with different ordering in input.dat) into one
c      
c     Written by D.P.-11/97
c     
c     Bruno Strandberg
c     ********************************************************************
c     Rev: modified 4 June 2014 from Deepshikha's 3He code
c     1. Got rid of densitytype, deuttype
c     2. Got rid of palpha,pbeta,nalpha,nbeta
c     3. Got rid of dg1p,dg2p,dg3p,dg4p,dg1n,dg2n,dg3n,dg4n
c     4. Got rid of firsttime - remove from other parts of code!
c     5. Got rid of ampfile,amp1Bfile,amp2Bfile
c     6. Got rid of ampUnitno,amp1BUnitno,amp2BUnitno
c     7. Created thetaLow, thetaHigh
c     8. Changed interval-->thetaInterval
c     9. Added Oepsilon3 option to calcstring parsing
c     10. Added whichbody option
c     11. Added variables to control quadrature settings
c     ********************************************************************
c     hgrie 20 June 2014: modified for use of LebedevLaikov or Gaussian
c     integration for theta & phi separately,
c     for solid angle integral in (12) system
c     
c     hgrie Oct 2014: added "nosymmetry" as option _not_ to use symmetry.
c     DEFAULT: use symmetry (does not have to be called)
c     
c     added "cartesian" as option for Cartesian basis of
c     photon polarisations (original choice by D.)
c     DEFAULT: spherical (i.e. cartesian=.false.)
c     
c     This routine reads in the values of the external photon momentum &
c     scattering angle in the lab frame. It takes them from the I/O
c     unit denoted by inUnitno, where they are stored in MeV and degrees.
c     
c     *******************************************************************
c     hgrie Feb 2016: added OQ4 and Odelta4 options,
c     but activated only for twobody!     
c     ********************************************************************
c
c     hgrie May 2017: implemented that j12max (max total ang mom in (12) subsystem)
c                     can be set in input file.
c                     Defaults are j12max=2 for onebody and j12max=1 for twobody.
c                     That is enough for convergence on the <1% level in amplitudes.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine ReadinputCommon(Elow,Ehigh,Einterval,frame,
     &     thetaLow,thetaHigh,thetaInterval,
     &     outfile,descriptors,densityFileName,inUnitno,
     &     nucleus,Anucl,twoSnucl,Mnucl,
     &     twoMzplimit,cartesian,verbosity)
c      
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
c     
c     
c     VARIABLES PASSED OUT:
c     
      real*8,intent(out)        :: Elow,Ehigh,Einterval
      integer,intent(out)       :: frame
      real*8,intent(out)        :: thetaLow,thetaHigh,thetaInterval
      character*500,intent(out) :: outfile
      character*200,intent(out) :: descriptors ! additional descriptors of calculation for outputfilename
      character*500,intent(out) :: densityFileName

      character*3, intent(out)  :: nucleus ! name of nucleus to be considered, for output file name
      integer,intent(out) :: Anucl               ! mass number of target nucleus
      integer,intent(out) :: twoSnucl            ! spin of target nucleus
      real*8,intent(out)  :: Mnucl               ! mass of target nucleus
      
      integer,intent(out) :: twoMzplimit         ! variable for using symmetry to cut down amplitude calculations
      logical,intent(out) :: cartesian           ! hgrie Oct 2014: for output in Cartesian basis of photon polarisations
      
      integer,intent(out) :: verbosity
c     
c     ---------------------------------------------------------------
c     BS: Ultimately variable def's in Readme?
c     ---------------------------------------------------------------
c     
c     Elow, Ehigh - beam-in low, beam-in high
c     Einterval   - step to move from Elow to Ehigh
c     
c     frame       - which frame the results are to be presented for
c     1 = lab. frame
c     2 = c.m. frame
c     
c     thetaLow, thetaHigh - low and high limits of theta angle
c     thetaInterval       - step to move from thetaLow-->thetaHigh       
c     
c     outfile - name of file to which output is to be written
c     *******************************************************************
c     
c     VARIABLES PASSED IN:
c     
      integer,intent(in) :: inUnitno ! I/O unit containing all the input information
c     
c     *******************************************************************
c     
c     INTERNAL VARIABLES:
c     
      character*500 calcstring ! holds info on type of calculation being done hgrie Aug 2020: increased from 80
      character*500 stringtolower ! function to convert string to all-lowercase, for comparisons -- defined at file end
      character*5 halfinteger     ! function to divide integer by 2 to string representing half-integer-- defined at end
c
      character*500 dummy
      
c     *******************************************************************
c     
c     set logical variables to FALSE, other parameters to "0": defines default values!
      descriptors = "-"
      
c     Read in the input file line-by-line, write output info to terminal
c     
      read (inUnitno,*) Elow,Ehigh,Einterval
      write (*,*) '   Lowest beam energy           = ', Elow, 'MeV'
      write (*,*) '   Highest beam energy          = ', Ehigh, 'MeV'
      write (*,*) '   Energy step from low to high = ', Einterval, 'MeV'
      read (inUnitno,*) thetaLow,thetaHigh,thetaInterval
      write (*,*) '   Lowest theta                 = ', thetaLow, 'deg'
      write (*,*) '   Highest theta                = ', thetaHigh, 'deg'
      write (*,*) '   Theta step from low to high  = ', thetaInterval, 'deg'

      read (inUnitno,'(A500)') outfile
      outfile = TRIM(outfile)
      write (*,*) 'First take on output filename (placeholders to be replaced below): '
      write (*,*) '      ',TRIM(outfile)
      
c     density filename and determination of target nucleus
      read (inUnitno, *) densityFileName
      write (*,*) 'First take on density filename (placeholders to be replaced below): '
      write (*,*) '      ',TRIM(densityFileName)
      
c     determine target nucleus from density filename, covering 3HE, 3he, 3He,...
c     Mnucl = 0. ! initialise Mnucl so that we can later ask if it was adapted by target determination
      if (index(stringtolower(densityFileName),'3he').ne.0) then
         nucleus = "3He"
         Anucl = 3
         twoSnucl = 1           ! 2*target spin
         Mnucl = M3He
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: 3He, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
      else if (index(stringtolower(densityFileName),'3h').ne.0) then
         nucleus = "3H"
         Anucl = 3
         twoSnucl = 1           ! 2*target spin
c         Mnucl = M3H
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: 3H, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
      else if ((index(stringtolower(densityFileName),'deuteron').ne.0).or.(index(stringtolower(densityFileName),'2h').ne.0)) then
         nucleus = "2H"
         Anucl = 2
         twoSnucl = 2           ! 2*target spin
         Mnucl = Md
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: Deuteron, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
      else if (index(stringtolower(densityFileName),'4he').ne.0) then
         nucleus = "4He"
         Anucl = 4
         twoSnucl = 0           ! 2*target spin
         Mnucl = M4He
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: 4He, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
      else if (index(stringtolower(densityFileName),'6li').ne.0) then
         nucleus = "6Li"
         Anucl = 6
         twoSnucl = 2           ! 2*target spin
         Mnucl = M6Li
         write (*,'(A,I2,A,A,A,G13.7,A)') '   Target Nucleus: 6Li, Anucl = ',Anucl,
     &        ', spin Snucl = ', halfinteger(twoSnucl),', mass Mnucl = ',Mnucl,'MeV'
c     I could also implement a routine "search string for A=7 and Z=3"...
      else
         stop "*** ERROR: Could not determine target nucleus from density filename: Abort."
      end if

      if (Mnucl.eq.0.) stop "*** ERROR: Target Nucleus Mass not defined: Abort."
      
c     finally, read string of common parameters
      read (inUnitno,'(A80)') calcstring

c      set up different levels of verbosity for debugging
      if (index(calcstring,'verbose').eq.0) then
         verbosity = 0
      else if (index(calcstring,'verbose4').ne.0) then 
         verbosity = 4
         write (*,*) '********** Verbose Mode 4 **********'
      else if (index(calcstring,'verbose3').ne.0) then
         verbosity = 3
         write (*,*) '********** Verbose Mode 3 **********'
      else if (index(calcstring,'verbose2').ne.0) then
         verbosity = 2
         write (*,*) '********** Verbose Mode 2 **********'
      else 
         verbosity = 1
c     print1Namps=.true.
         write (*,*) '********** Verbose Mode 1 **********'
      end if
c     
c     Determine frame
c     
      if (index(calcstring,'lab').ne. 0) then
         frame=lab
         write (*,*) 'Input and output in laboratory frame.'
      else if (index(calcstring,'cm').ne. 0) then
         frame=cm
         write (*,*) 'Input and output in centre-of-mass frame.'
      else 
         write (*,*) '*** ERROR: unknown frame.'
         stop
      end if
c     
c     Determine if using symmetry
c     
      if (index(calcstring,'nosymmetry').eq.0) then
         twoMzplimit=0
         dummy = ""
         if(MOD(twoSnucl,2).eq.0) then
            dummy = "(and Mz(in)>0 for Mz(out)=0, and Resultxx, Resultyy for Mz(in)=Mz(out) with Resultxy=Resultyx=0)"
         end if
         write (*,'(A,I3,A,A,A,I3,A)') " Calculate",
     &        2*(twoSnucl+1)**2,
     &        " amplitudes with Mz(out)>0", TRIM(dummy),", find other",
     &        2*(twoSnucl+1)**2," by Symmetry."
      else
         twoMzplimit=-twoSnucl
         write (*,'(A,I3,A)') "Calculate all",4*(twoSnucl+1)**2," amplitudes independently, not using Symmetry."
      end if
c     
c     Determine if using Cartesian basis of photon polarisation (default:spherical)
c     
      if (index(calcstring,'cartesian').ne.0) then
         cartesian=.true.
         write (*,*) "Output in Cartesian basis of photon polarisations."
      else
         cartesian=.false.
         write (*,*) "Output in spherical basis of photon polarisations."
      end if

      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018
c     Now the routine which reads input specific for onebody

      subroutine ReadinputOnebody(inUnitno,calctype,variedA,descriptors,
c---- Variable to control Feynman quadrature settings------------------------
     &     Nx,
     &     verbosity)
c     
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
c     
c     VARIABLES PASSED INOUT:
      character*200,intent(inout) :: descriptors ! additional descriptors of calculation for outputfilename
      
c     VARIABLES PASSED OUT:
c     
      integer,intent(out) :: calctype
      integer,intent(out) :: variedA
      integer,intent(out) :: Nx                ! number of quadratures for loop integration in one-body amp.
      
c     ---------------------------------------------------------------
c     BS: Ultimately variable def's in Readme?
c     ---------------------------------------------------------------
c     
c     calctype    - which calculation to do
c     hgrie 19 Oct 2014: modified to VaryAp or VaryAn
c     to indicate if p or n amp varied
c     hgrie 19 Oct 2014: changed IA assignation since unused. VaryA was 5;
c     now VaryAp is 0, VaryAn is 1. Leaves numbers >=5 open for future.
c     -1: uninitialised value.
c     0=VaryAp, vary a proton amplitude defined by variedA, all other amps 0
c     1=VaryAp, vary a neutron amplitude defined by variedA, all other amps 0
c     2=OQ2=Odelta0, the O(q²) calculation Chi PT calculation;
c     3=OQ3=Odelta2, the full O(q³) Chi PT calculation;
c     4=Oepsilon3=Odelta3, calculation with delta diagrams.
c     
c     BS:----
c     variedA to indicate which A is varied by calctype=VaryA
c     *******************************************************************
c     
c     VARIABLES PASSED IN:
c     
      integer,intent(in) :: inUnitno ! I/O unit containing all the input information
      integer,intent(in) :: verbosity
c     
c     *******************************************************************
c     
c     INTERNAL VARIABLES:
c     
      character*500 calcstring ! used to hold info. on type of calculation being done
c     
c     *******************************************************************
c     
      if (verbosity.eq.1000) continue ! keep for future use
      if (descriptors.ne.'a') continue ! keep for future use
      
c     Read in the input file line-by-line, write output info to terminal
c     
      read (inUnitno,'(A80)') calcstring
c     Read in Feynman quadrature settings------------------------------------
      read (inUnitno,*) Nx
c-----------------------------------------------------------------------------------------

c     Determine the calculation to be done

      if ((index(calcstring,'OQ2').ne. 0).or.(index(calcstring,'Odelta0').ne. 0)) then 
         calctype=OQ2
         write (*,*) 'O(e²delta⁰)=O(Q²) calculation.'
      else if ((index(calcstring,'OQ3').ne. 0).or.(index(calcstring,'Odelta2').ne. 0) ) then
         calctype=OQ3
         write (*,*) 'O(e²delta²)=O(Q³) calculation.'
      else if ((index(calcstring,'Oepsilon3').ne. 0).or.(index(calcstring,'Odelta3').ne. 0) ) then
         calctype=Oepsilon3
         write (*,*) 'O(e²delta³)=O(epsilon³) calculation.'
         write (*,*) '  Delta and DeltaPi parts of static scalar polarisabilities are subtracted.'
         write (*,*) "  Delta parameters:"
         write (*,*) "    Delta mass                  Mdelta    = ",Mdelta," MeV"
         write (*,*) "      => Delta-N mass splitting Delta     = ",delta," MeV"
         write (*,*) "    PionNDelta coupling         gPiNDelta = ",gpind
         write (*,*) "    DeltaGamma coupl. (nonrel.) b1        = ",b1
      else if ((index(calcstring,'OQ4').ne. 0).or.(index(calcstring,'Odelta4').ne. 0) ) then
         
c     hgrie note Feb 2017: when implemented, need to use different LECs c_i, for onebody δ⁴ vs Q⁴

         if (index(calcstring,'OQ4').ne. 0) then
            calctype=OQ4
            write (*,*) 'O(Q⁴) calculation.'
            write (*,*) '   -- NOT YET IMLPEMENTED FOR ONEBODY. -- Exiting.'
            stop
         else
            calctype=Odelta4
            write (*,*) 'O(e²delta⁴) calculation.'
            write (*,*) '   -- NOT YET IMLPEMENTED FOR ONEBODY. -- Exiting.'
            stop
         end if
c     variation of proton amplitudes         
      else if (index(calcstring,'VaryA1p').ne. 0) then
         calctype=VaryAp
         variedA=1
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA2p').ne. 0) then
         calctype=VaryAp
         variedA=2
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA3p').ne. 0) then
         calctype=VaryAp
         variedA=3
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA4p').ne. 0) then
         calctype=VaryAp
         variedA=4
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA5p').ne. 0) then
         calctype=VaryAp
         variedA=5
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA6p').ne. 0) then
         calctype=VaryAp
         variedA=6
         write (*,'(A,I1,A)') 'Varying proton Amplitude A', variedA ,'p, all other amplitudes set to 0.'
c     variation of neutron amplitudes         
      else if (index(calcstring,'VaryA1n').ne. 0) then
         calctype=VaryAn
         variedA=1
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA2n').ne. 0) then
         calctype=VaryAn
         variedA=2
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA3n').ne. 0) then
         calctype=VaryAn
         variedA=3
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA4n').ne. 0) then
         calctype=VaryAn
         variedA=4
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA5n').ne. 0) then
         calctype=VaryAn
         variedA=5
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else if (index(calcstring,'VaryA6n').ne. 0) then
         calctype=VaryAn
         variedA=6
         write (*,'(A,I1,A)') 'Varying neutron Amplitude A', variedA ,'n, all other amplitudes set to 0.'
      else
c     hgrie 19 Oct 2014: if calctype not yet specified, things went wrong. 
         write (*,*) '*** ERROR: Calculation type unknown. -- Exiting.'
         stop
      end if
c     
c     write quadrature variables to terminal
c     
      write (*,*) "No integrals in onebody."
c     onebody only: number of points of Feynman parameter integration
      write (*,180) Nx

c     end hgrie mod         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     formats
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 180  format (1X,"Number of quadratures in Feynman parameter integral of single-N amplitude = ",I2)
      
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018
c     Now the routine which reads input specific for twobody

      subroutine ReadinputTwobody(inUnitno,calctype,descriptors,
c---- Variables to control radial quadrature settings------------------------
     &     NP12A,NP12B,P12A,P12B,P12C,
c---- Variables to control angular quadrature settings------------------------
     &     AngularType12,Nanggrid12,
     &     Nordth12,Nordphi12,
     &     NthBins12,NphiBins12,
     &     j12max,verbosity)
c     
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
c     
c     VARIABLES PASSED INOUT:
      character*200,intent(inout) :: descriptors ! additional descriptors of calculation for outputfilename
      
c     VARIABLES PASSED OUT:
c     
      integer,intent(out) :: calctype
      integer,intent(out) :: AngularType12,Nordth12,Nordphi12,Nanggrid12
      integer,intent(out) :: NthBins12,NphiBins12
      integer,intent(out) :: j12max            !hgrie May 2017: maximum total ang mom in (12) system 
      
c     BS: Added quadrature settings variables------------------------
      integer,intent(out) :: NP12A,NP12B
      real*8,intent(out)  ::  P12A,P12B,P12C
c     ---------------------------------------------------------------
c     BS: Ultimately variable def's in Readme?
c     ---------------------------------------------------------------
c     
c     calctype    - which calculation to do
c     2=OQ2=Odelta0, the O(q²) calculation Chi PT calculation;
c     3=OQ3=Odelta2, the full O(q³) Chi PT calculation;
c     4=Oepsilon3=Odelta3, calculation with delta diagrams.
c     
c     NNincl - set to true if the NN cut contribution is to be included -- NEITHER USED NOR IMPLEMENTED YET
c     
c     Nordth12,NthBins12 - define quadrature distribution for theta integral
c     NordPhi12,NphiBins12 - define quadrature distribution for phi integral
c     
c     *******************************************************************
c     
c     VARIABLES PASSED IN:
c     
      integer,intent(in) :: inUnitno ! I/O unit containing all the input information
      integer,intent(in) :: verbosity
c     
c     *******************************************************************
c     
c     INTERNAL VARIABLES:
c     
      character*500 calcstring ! holds info on type of calculation being done
      character*1  dummy ! a variable only for converting he integer j12max to a character string
      
c     *******************************************************************
c     
      if (verbosity.eq.1000) continue ! keep for future use
      
c     Read in the input file line-by-line, write output info to terminal
c     
      read (inUnitno,'(A80)') calcstring
      read (inUnitno,*) NP12A,NP12B
      read (inUnitno,*) P12A,P12B,P12C
c     Read in angular quadrature settings------------------------------------
c     hgrie Sep 2020: first determine angular type, then read only those parameters needed for that type
      read (inUnitno,*) AngularType12
      backspace inUnitno                ! reset to start of line
      if (AngularType12.eq.1) then      ! Gauss
         read (inUnitno,*) AngularType12,Nordth12,Nthbins12,Nordphi12,Nphibins12
      else if (AngularType12.eq.2) then ! LebedevLaikov
         read (inUnitno,*) AngularType12,Nordth12
      else                              ! none of the two: continue -- error message will be produced below
         continue
      end if
c----------------------------------------------------------------------------------------
c     Determine the calculation to be done
      
      if ((index(calcstring,'OQ2').ne. 0).or.(index(calcstring,'Odelta0').ne. 0)) then 
         write (*,*) 'Twobody has no O(e²delta⁰)=O(Q²) MEC contribution.'
         write (*,*) '   -- Use at least O(e²delta²)=O(Q³). -- Exiting.'
         stop
      else if ((index(calcstring,'OQ3').ne. 0).or.(index(calcstring,'Odelta2').ne. 0) ) then
         calctype=OQ3
         write (*,*) 'O(e²delta²)=O(Q³) calculation. (12) subsystem is (pn) since only charged MECs.'
      else if ((index(calcstring,'Oepsilon3').ne. 0).or.(index(calcstring,'Odelta3').ne. 0) ) then
         write (*,*) 'Twobody has no O(e²delta³)=O(epsilon³) MEC contribution.'
         write (*,*) '   -- Use O(e²delta²)=O(Q³) or O(e²delta⁴)=O(Q⁴) instead. -- Exiting.'
         stop
      else if ((index(calcstring,'OQ4').ne. 0).or.(index(calcstring,'Odelta4').ne. 0) ) then
c     hgrie note Feb 2017: once implemented, need to use different LECs c_i, for onebody δ⁴ vs Q⁴
         write (*,*) 'O(e²delta⁴)=O(Q⁴) MECs identical [no Δ(1232)]. (12) subsystem is (pn) since only charged MECs.'
         calctype=Odelta4           ! again: MECs for Δ-less and Δ-ful are identical, so need no switch.
      else
c     hgrie 19 Oct 2014: if calctype not yet specified, things went wrong. 
         write (*,*) '*** ERROR: Calculation type unknown. -- Exiting.'
         stop
      end if
c
c     Now start on parameters of remaining integration
      write (*,*) "Integration parameters of (12) subsystem:"
      
c     Determine maximum total angular momentum in (12) subsystem from input file
c     
      if (index(calcstring,'j12max=').eq.0) then
         j12max = 1             ! twobody converges faster
         write(*,'(A,I4)') " Total angular momentum in (12) subsystem not specified -- using default j12max =",j12max
      else if (index(calcstring,'j12max=').ne.0) then
         if (index(calcstring,'j12max=1').ne.0) then
            j12max = 1
         else  if (index(calcstring,'j12max=0').ne.0) then
            j12max = 0
         else  if (index(calcstring,'j12max=2').ne.0) then
            j12max = 2
         else  if (index(calcstring,'j12max=3').ne.0) then
            j12max = 3
         else  if (index(calcstring,'j12max=4').ne.0) then
            j12max = 4
         else  if (index(calcstring,'j12max=5').ne.0) then
            j12max = 5
         else
            write(*,*) "*** ERROR: Input attempted to set j12max to value which is not 0,1,2,3,4 or 5. -- Exiting."
            stop
         end if    
         write(*,'(A,I4)') " Total angular momentum in (12) subsystem set to j12max = ",j12max
         write(dummy,'(I1)') j12max
         descriptors = trim(descriptors) // "j12max=" // dummy // "-"
      end if
c     
c     write quadrature variables to terminal
c        radial integration routines & parameters used    

      write (*,*) "Radial Integration in (12) subsystem:"
      write (*,*) "   total of NP12 = NP12A+NP12B = ",NP12A+NP12B," points"
      write (*,126) NP12A
      write (*,127) NP12A/2,P12A
      write (*,128) NP12A/2,P12A,P12B
      write (*,129) NP12B,P12B,P12C
      write (*,*) "   NP12p == ( NP12A+NP12B = ",NP12A+NP12B,") cannot be dialled."
c        angular integration routines & parameters used 
      if (AngularType12.eq.1) then
         write (*,*) "Angular Integration in (12) subsystem:"
         write (*,*) "   by Gauss-Legendre for theta & phi separately"
         write (*,132) Nordth12, Nthbins12, Nordth12*Nthbins12
         write (*,140) Nordphi12, Nphibins12, Nordphi12* Nphibins12
         write (*,142) (Nordth12*NthBins12)*(Nordphi12*NphiBins12)
         write (*,*) "           Nanggrid12 not used."
      else if (AngularType12.eq.2) then
         write (*,*) "Angular Integration in (12) subsystem"
         write (*,*) "   by Lebedev-Laikov for theta & phi combined:"
         Nanggrid12 = Nordth12
         if (Nanggrid12.gt.Nangmax) then
            write (*,*) '*** ERROR: Nanggrid12 (set to Nordth12) >',Nangmax,' too large -- Exiting.'
            stop
         end if
         write (*,144) Nanggrid12
         write (*,*) "   Nthbins12, Nordphi12, Nphibins12 have no meaning."
      else
         write (*,*) "*** ERROR: Illegal AngularType12 -- Exiting."
         stop
      end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     formats
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 126  format (' ','   1 bin with hyperbolic map: (NP12A = ', I2,') points:')
 127  format (' ','      (NP12A/2 = ', I2,') points in p12 interval [0; (P12A = ',F10.4,')] fm^-1')
 128  format (' ','      (NP12A/2 = ', I2,') points in p12 interval [(P12A = ',F10.4,');(P12B = ',F10.4,')] fm^-1')
 129  format (' ','   1 bin with linear map: (NP12B=', I2,') points in p12 interval [(P12B = ',F10.4,');(P12C = ',F10.4,')] fm^-1')

 132  format (' ','   theta: (Nordth12 = ',I4,') points per (NthBins12 = ', I2,') bins = ',I16,' points.')
 140  format (' ','   phi:   (Nordphi12  = ',I4,') points per (NphiBins12  = ', I2,') bins = ',I16,' points.')
 142  format (' ','   ==> size of solid angle grid = (Nordth12*NthBins12)*(Nordphi12*NphiBins12) = ',I16,' points.')
 144  format (' ','   ==> preliminary size of solid angle grid Nanggrid12(set to Nordth12) = 'I4,' points.')

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie Aug 2020
c     Now the routine which reads last line as comment and appends to descriptors

      subroutine ReadinputCommonComments(descriptors,inUnit,verbosity)     
      
c     VARIABLES PASSED INOUT:
      character*200,intent(inout) :: descriptors ! additional descriptors of calculation for outputfilename
      
c     VARIABLES PASSED OUT:
      
c     VARIABLES PASSED IN:
      integer,intent(in)     :: inUnit
      integer,intent(in)     :: verbosity
c     
c     INTERNAL VARIABLES:
      character*200 comments
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      if (verbosity.eq.1000) continue ! keep for future use
      
      read (inUnit,*) comments
      write (*,*) "ADDITIONAL ",trim(comments)
c     add comments to end of descriptors
      descriptors = trim(descriptors) // "." // comments(11:) ! strips out "COMMENTS:_"

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     function to convert string to all-lowercase, hgrie May 2018
c     from https://groups.google.com/forum/#!msg/comp.lang.fortran/CKx1L2Ahkxg/HH_kMoHAffcJ
      
      function stringtolower( string ) result (new)
      character(len=500)           :: string

      character(len=len(string)) :: new
      
      integer                    :: i
      integer                    :: k

      length = len(string)
      new    = string
      do i = 1,len(string)
         k = iachar(string(i:i))
         if ( k >= iachar('A') .and. k <= iachar('Z') ) then
            k = k + iachar('a') - iachar('A')
            new(i:i) = achar(k)
         endif
      enddo
      end function stringtolower 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     function to make half-integer -- self-made hgrie May 2018

      function halfinteger(number) result (string)
      integer :: number
      character*5 string

      if (mod(number,2).eq.0) then
         write (string,'(I2)') number/2
      else
         write (string,'(I2,A)') number,"/2"
      end if
      end function halfinteger
