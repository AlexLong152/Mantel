c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2017, revised May 2018 (see below)
c     based on Andreas Nogga's "template" files common-densities/2Ndensity-module/testcompdens.F90 in May 2017/2018.
c                                           and common-densities/2Ndensity-module/CompDens.F90 in May 2017/2018.
c
c     This file adapted from main.twobody.f, plus changes.
c          added read of density matrix,
c          eliminated spectator(3)-integrations and calls, and calls to wave function
c          densityFileName now used to set name of input density file
c          split setquad setting up angular integrations into separate routines
c                 for (12) and spectator (3) system
c          here, (12) integration only -- spectator (3) integration provided by 2Ndensity
c     j12max can now be specified in input file: default is j12max=2 for onebody and j12max=1 for twobody,
c     as in previous runs (where they were hardwired in code).
c     These values give MEs which are converged to better than 0.7% -- see documentation/.
c
c     twoSmax/twoMz dependence: only via array size of Result()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NOTE ON UNITS hgrie Nov 2023
c      
c     The "mantle" code has base unit fm (but NOT the output file/Result(), see below!)
c     p12, p12p: momenta in fm^-1
c     rho      : 2N density in fm^3 (quantum numbers per volume momentum space)
c
c     k        : photon omentum/energy still given in MeV
c
c     However, in finalstatesums.twobodyvia2Ndensity.f, the call
c              call Calculate2BIntegralI2(...,p12*HC,P12MAG(ip12p)*HC,...)
c     converts the momenta from fm^-1 to MeV.
c     That subroutine is defined in calculate2BI2.f.
c     Therefore, that routine and the subsequent "kernel" parts of the code use base unit MeV:
c     calculate2BI2.f
c     2Bspinisospintrans.f : part of "kernel"
c     spintricks.f         : part of "kernel"
c     spintricksasy.f      : part of "kernel"
c
c     Therefore, the twobody "kernel diagrams" are all using base unit of MeV.
c
c     Output units: Code constructed such that if kernel is given in units of MeV^-n, then Result() output is in MeV^(3-n).
c     ===> SEE detailed description and examples in 2Bspinisospintrans.f .
c
c     Routines which compute vectors, like calcmomenta.f, simply use the same base unit in and out, i.e. are "unit neutral".
c      
c     The integral done is symbolically:
c
c     (A choose 2) * Σ ∫ dp12 p12² ∫ dp12p p12p²/(2π)³ (HC)³ * rho(p12,p12p) * Σ ∫ dΩ12 ∫ dΩ12p * <projected kernel onto orbital angmoms via Clebsches>
c     # NN pairs       ========================       ======   =============   ================   =====================================================
c                                fm^-6                MeVfm³        fm³           no units                  base unit MeV ==> MeV^-n
c     ......................................................................   ................
c            Σ_(mt12,j12,s12,l12,m12)                  in main.twobody.f       ms: spin proj in (12)
c            Σ_(mt12p=mt12,j12p,s12p,l12p,m12p,Mzp,Mz) in finalstatesums.f     Σ_(msp,ms) and
c            integrations done in                                              ang int done in
c            finalstatesums.twobodyvia2Ndensity.f                              calculate2BI2.f
c     -----------------------------------------------------------------------------------------   -----------------------------------------------------   
c                                   "mantle" code                                                                   "kernel" code
c
c     ===> Overall units of output are MeV^(3-n) for kernel with base units MeV^-n !!!!
c 
c     The (HC)³/(2π)³ above (programmed in finalstatesums.twobodyvia2Ndensity.f 's "f=...")
c     converts between the fm units of the "mantle" and the MeV units of the "kernel".
c     It ALSO includes ONE of the Fourier volumes 1/(2π)³. There is no second Fourier volume (killed by phase space). 
c     This guarantees that onebody and twobody have the same size and can simply be summed to get the total amplitude:
c
c     amplitude = onebody + twobody, without any relative factors (provided both provide output in same base units).
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c     Does output of multiple angles & energies to same output file work? May need adjusting output file name.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c
c     hgrie Oct 2022: *HUGE CHANGE* inside twobodyfinalstatesumsvia2Ndensity():
c           Defined twobody ME to INCLUDE the factor 1/(2π)³ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c           so that final amplitudes for onebody and twobody have SAME sizes.
c      
c     hgrie Sep 2020: in readinputTwobody():
c           removed NP12p = 3rd variable in integration-grid
c           line of input file. It is actually never used for anything!
c           That also reduced number of arguments of readinputTwobody()!
c     hgrie Aug/Sep 2020: rewrote makedensityfilename() to deal with extracting densities from a .gz or downloading from server
c     hgrie June 2018: renamed "parity" to "symmetry -- see notes in usesymmetry+*.f
c       
c     hgrie May 2018: decluttered files: remove obsolete variables, unify look, documentation,...
c     hgrie May 2018: new subroutines to read input, weeded out unused variables, new input.dat format. 
c     hgrie May 2018: rewritten to accommodate hdf5 format for input files of 2N density rho
c     hgrie May 2018:
c           All quantum numbers which start with "two" run over integers
c                  Examples: 
c                     twoMz,twoMzp: magnetic quantum numbers of in/out target nucleus, times 2.
c                     twoSnucl: 2 x spin of target nucleus
c           For all such variables, run over 2xQM values.
c                  Examples:
c                     in do-loops: twoMz runs from +twoSnucl to -twoSnucl, in steps of -2
c                     in array: Result(extQnum,twoMzp,twoMz) runs over   (1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
c
c     This numbering agrees with Andreas' assignments of magnetic quantum numbers.
c     In ARRAYs, this means a slight waste of space since we leave many array entries unused.
c     Example S=1/2: 2x2=4 entries needed, but array is (-1:1,-1:1): 3x3=9  : 225% of needed size
c     Example S=1:   3x3=9 entries needed, but array is (-2:2,-2:2): 5x5=25 : 280% of needed size
c
c     Still, our arrays are for small nuclei -- we do not really waste a lot. Not a time/storage issue.
c
c     hgrie May 2018: outsourced symmetry+output into sub routine outputroutine(), identical for onebody and twobody
c      
c     Implemented symmetry for arbitrary nucleon spin:
c     Use Mzp>=0, and for Mzp=0, run only over Mz>=0
c     -- that's still 2 more than necessary since ME(+0->+0) = ME(-0->-0) and ME(+0->-0) = ME(-0->+0)
c     but it's good enough, saving lots of CPU time.
c     see manuscript "Compton Densities Approach" pp11-12
c      
c     In May 2018, Andreas replaced fkltt() by a more sophisticated and parallel routine initclebsch() in
c     common-densities/2Ndensity-module/clebsch.F .
c      
c     hgrie June 2017: implemented run over several energies and angles into same output file
c                      implemented: when input file contains a 2Ndensity filename which contains "XXX"
c                                   and "YYY" strings, then these are automatically replaced by
c                                   XXX => energy of run (in numeric format used by Andreas)
c                                   YYY => angle of run (in numeric format used by Andreas)
c                      That reduces error-proneness.
c      
c     modified by hgrie June 2014:
c     calculate nucleon amplitudes outside fewbody loops
c     use symmetry to calculate only amplitudes with 3He out-spin +1/2
c     modified hgrie 20 June 2014: add option to use LebedevLaikov or Gaussian
c     integration for theta & phi separately,
c     for solid angle integral in (12) system 
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

      PROGRAM twobodydensitymain
      
      USE CompDens              ! needs module CompDens.mod
      
      IMPLICIT NONE
c**********************************************************************
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
      include '../common-densities/constants.def'
c**********************************************************************
c     
c     program argument management
      integer narg              ! number of arguments
      character*500 inputfile   ! argument string is written to this
c     
c*********************************************************************
      integer NP12A,NP12B,NP12
      real*8  P12A,P12B,P12C
      real*8 P12MAG(Npmax),AP12MAG(Npmax)
c     
c**********************************************************************
c     
c     Input information:
c     
c     inUnitno-unit where input information is stored
c     outUnitno-unit where output is written
c     kgamma-photon momentum in cm frame in units of MeV
c     Nangles-number of angles at which calculation is to be done
c     thetacm-photon scattering angle in cm frame, in radians
c     calctype-which calculation to do
c     
c     outfile-name of output file
c     
      integer inUnitno,outUnitno
      
      real*8 Egamma,kgamma,thetaL,thetacm,Elow,Ehigh,Einterval
      
      real*8 thetaLow,thetaHigh,thetaInterval
      integer calctype,Nangles,Nenergy,ienergy,j ! number of energies/angles; index for energies/angles
      character*200 descriptors  ! additional descriptors of calculation

      integer,parameter :: variedA = 0 ! no variedA since no 1N amplitude, but define here for makeoutputfilename()
      character*3 nucleus ! name of nucleus to be considered, for output file name
      integer Anucl             ! target nucleus mass number
      integer twoSnucl          ! 2 x target nucleus spin
      real*8 Mnucl               ! mass of target nucleus
      
      character*500 outfile
      character*500 densityFileName,originaldensityFileName ! second for multiple energies or angles
      
c*********************************************************************
c     Quadrature variables:
c     
c     Nordth,Nthbins-number of quadratures/bin and number of bins for 
c     theta integration
c     Nordphi,Nphibins-number of quadratures/bin and number of bins for phi integration
c     Nth-total number of theta quadratures
c     Nphi12-total number of phi quadratures
c     pq,wp-radial quadratures and weights, set up on [0,infty]
c     thq,wth-theta quadratures and weights, set up on [0,PI]
c     phi12,wphi-phi quadratures and weights, set up on [0,2 PI]
c     
      integer Nordth12,Nordphi12,Nthbins12,Nphibins12

      integer Nth12,Nphi12
      
      real*8 th12(Nangmax),phi12(Nangmax)
      
      integer AngularType12,Nanggrid12
      real*8 angweight12(Nangmax,Nangmax)
c     
c----------------------------------------------------------------------
c     
c     Momentum variables:
c     
c     pp-magnitude of final-state relative three-momentum vector,
c     in IA=p + 1/2(k - k') as a vector sum. In IA the
c     kinematics are 
c     
c                 \    /
c                  \  /   
c     k' + k/2 + p  \/        -k/2 + p
c     ---------------------------------------
c     
c     
c     
c     ---------------------------------------
c     - k/2 - p
c     
      real*8 k,kth,kphi,kp,kpth,kpphi,Qk,Qkth,Qkphi
      real*8 t,omega
c     
c**********************************************************************
c     
      integer m12,mt12          ! projections of total ang mom & isospin of (12) subsystem: automatically integers
      
      integer i,ip12
      integer l12,s12,j12,t12   ! orb ang mom, spin, total ang mom, isospin of (12): automatically integers
      
      integer j12max            ! max total ang mom in (12) subsystem -- =1 suffices for 1% accuracy.
      
c     projections of target nucleus' in-spin, out-spin
      integer twoMz,twoMzp      !  -- these two not use right now
      
      integer extQnum, extQnumlimit      ! counter and number of combined external quantum numbers of in and out state
      
      integer symmetry          ! whether and which specific package Usesymmetry() should be used for process -- NOT YET IMPLEMENTED!!!

      real*8 frac

      complex*16, allocatable :: Result(:,:,:) ! extQnum from 1 to extQnumlimit; twoMzp from -twoSnucl to twoSnucl, stepsize 2; twoMz from -twoSnucl to twoSnucl, stepsize 2; rest blank.
      
      integer verbosity         ! verbosity index for stdout hgrie June 2014

      integer test              ! a generic integer for testing i/o
      logical testtf            ! a generic logical for testing i/o
c     if density file generated from a .gz, delete that temporary file after each energy/angle
c     if downloaded and .gz, also delete the download.
c     0: do not delete; 1: delete un-gz'd file; 2: delete downloaded and un-gz'd file 
      integer rmDensityFileLater 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end OF VARIABLE DECLARATIONS, BEGINNING OF CODING
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) 
      write(*,*) "================================================================================"
      write(*,*) "Twobody Contributions to Few-Nucleon Processes Calculated Via 2N-Density Matrix"
      write(*,*) "================================================================================"
      write(*,*) "   Mantle Code Version 1.0"
      write(*,*) "      Alexander Long/hgrie starting November 2023   "
      write(*,*) "      based on Compton density code: D. Phillips/A. Nogga/hgrie starting August 2020"
      write(*,*) "      based on 3He Compton codes: D. Phillips/A. Nogga/hgrie starting May 2018"
      write(*,*) "                                  D. Phillips/B. Strasberg/hgrie starting June 2014"
      write(*,*) "                                  with A. Margaryan 2016-17, modifying codes by D. Shukla 2007/8"
c**********************************************************************
c     Reading the input file from command line
c**********************************************************************

c     get the number of arguments
      narg=command_argument_count()

c     if you have 1 argument, write it to inputfile, otherwise stop 
      if (narg.eq.1) then
         call get_command_argument(1, inputfile)
      else
         write(*,*) "*** ERROR: Pass one input file as argument!"
         stop
      end if
c     
c     
c**********************************************************************
c     Report kernel process and version
      call KernelGreeting(verbosity)
      
c**********************************************************************
c     Reading in data from the input file
c**********************************************************************
      inUnitno=13
      outUnitno=10

      open(unit=inUnitno, file= inputfile, status= 'OLD',iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open input file!!! Aborting."
      
      call ReadinputCommon(Elow,Ehigh,Einterval,
     &     thetaLow,thetaHigh,thetaInterval,
     &     outfile,descriptors,densityFileName,inUnitno,
     &     nucleus,Anucl,twoSnucl,Mnucl,extQnumlimit,
     &     symmetry,verbosity)
c      
      call ReadinputTwobody(inUnitno,calctype,descriptors,
c---- Variables to control radial quadrature settings------------------------
     &     NP12A,NP12B,P12A,P12B,P12C,
c---- Variables to control angular quadrature settings------------------------
     &     AngularType12,Nanggrid12,
     &     Nordth12,Nordphi12,
     &     NthBins12,NphiBins12,
     &     j12max)
      
      call ReadinputCommonComments(descriptors,inUnitno,verbosity)
      
      close(unit=inUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close input file!!! Aborting."
c
      call makeoutputfilename(outfile,calctype,nucleus,descriptors,densityFileName,variedA,
     &     Elow,Ehigh,Einterval,thetaLow,thetaHigh,thetaInterval,verbosity)
      
c**********************************************************************
c     FINAL PRELIMINARIES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     hgrie June 2017: keep original filename: needed for replacements of energy & angle later 
      originaldensityFileName = densityFileName

      thetaL=0.d0
      thetacm=0.d0
c     BS: new Nangles algebra due to changed input form
      Nangles=int((thetaHigh-thetaLow)/thetaInterval)+1
      Nenergy=int((Ehigh-Elow)/Einterval)+1
      
c**********************************************************************
c     (12) integration set-up
c          spectator (3) integration is provided by 2Ndensity, so nothing to do. 
c     define total number of integration points for (12) mom magnitude
      NP12 = NP12A+NP12B
c
c     Set up radial quadratures for (12) integration
      call TRNS(NP12A,NP12B,NP12,P12A,P12B,P12C,P12MAG,AP12MAG)
      
c     Set up angular quadratures for (12) integration
      call Setquad12(th12,Nth12,phi12,Nphi12,
     &     Nordth12,Nthbins12,Nordphi12,Nphibins12,
     &     AngularType12,angweight12,Nanggrid12,1 !verbosity
     &     )

      write (*,*) "***************************** END OF INITIALISATION *****************************"
      
c**********************************************************************
      open(unit=outUnitno, file=outfile,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not open output file!!! Aborting."
c**********************************************************************
c     Loop over Energies
c**********************************************************************
      do j=1,Nenergy
         Egamma=Elow+Einterval*(j-1)
         ienergy=int(Egamma)
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
         write(*,*)
         write(*,*) "Incoming photon energy (rounded) = ", ienergy, " MeV"
         write(*,*) "Number of Angles Nangles         = ", Nangles
c**********************************************************************
c     Loop over angles
c**********************************************************************
         do i=1,Nangles
            thetacm=(thetaLow+real(i-1)*thetaInterval)*Pi/180.d0
            if (thetacm.eq.0.0d0) then
               thetacm=1.0d0*Pi/180.d0
               write(*,*) "   Replaced input angle 0 deg with 1 deg."
            end if   
            kgamma=Egamma      
            write (outUnitno,*) "cm ","omega = ",Egamma,"thetacm = ",thetacm*180.0/Pi
            write(*,*)
            write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            write(*,*) 'Calculating amps: theta =',thetacm*180.0/Pi,'deg: angle #',i
c**********************************************************************
            call calcphotonmomenta(k,kth,kphi,t,kp,kpth,kpphi,omega,
     &           Qk,Qkth,Qkphi,kgamma,thetacm,verbosity)
c**********************************************************************
c     be a good boy and initialise everything to 0, overwriting entries from previous ω/θ
c**********************************************************************
            allocate(Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl))
            Result=c0
c**********************************************************************
c     hgrie June 2017: create name of 1Ndensity file for given energy and angle, unpack it
c     define correct formats for energy and angle
c     hgrie May 2018: outsourced into subroutine common-densities/makedensityfilename.f
            densityFileName = originaldensityFileName
            call makedensityfilename(densityFileName,Egamma,thetacm,rmDensityFileLater,verbosity)
c**********************************************************************
c     hgrie May 2017: read 2Ndensity
            call read2Ndensity(densityFileName,Anucl,twoSnucl,omega,thetacm,j12max,P12MAG,AP12MAG,NP12,verbosity)
c**********************************************************************      
c     hgrie Aug/Sep 2020: delete the local .h5 file if one was generated from .gz
            if (rmDensityFileLater.gt.0) then
               call EXECUTE_COMMAND_LINE("rm "//densityFileName, WAIT=.True., EXITSTAT=test )
               if (test.ne.0) stop "*** ERROR: Could not remove .h5 file created from .gz"
               write(*,*) "   Removed .h5 file unzipped from .gz."
               if (rmDensityFileLater.ge.2) then
                  INQUIRE(FILE=TRIM(densityFileName)//".gz", EXIST=testtf)
                  if ( testtf ) then
                     call EXECUTE_COMMAND_LINE("rm "//TRIM(densityFileName)//".gz", WAIT=.True., EXITSTAT=test )
                     if (test.ne.0) stop "*** ERROR: Could not remove .h5 file created from .gz"
                     write(*,*) "   Removed .h5.gz file downloaded."
                  end if
               end if
            end if ! rmDensityFileLater
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c**********************************************************************
            write(*,*) "*********Now convoluting 2N helicity amplitudes with 2N density matrix.*********"  
            do mt12=0,0                ! only charged pion exchange at OQ4 => (12) subsystem is (pn) 
               do j12=0,j12max                      ! total ang mom (12); usually Jmax=1 for 1% convergence
                  do s12=0,1                        ! spin (12)
                     do l12=abs(j12-s12),j12+s12    ! angular mom. (12)
                        t12=(1-(-1)**(l12+s12+1))/2 ! isospin (12)
                        do m12=-j12,j12             ! spin projection (12)
                           do ip12=1,NP12           ! integration over momentum magnitude (12)
                              call twobodyfinalstatesumsvia2Ndensity(
     &                             Result,
     &                             Anucl,twoSnucl,extQnumlimit,j12,m12,l12,s12,t12,mt12,
     &                             k,thetacm,
     &                             ip12,P12MAG(ip12),AP12MAG(ip12),     
     &                             p12MAG,AP12MAG,NP12,                                 
     &                             th12,phi12,Nth12,Nphi12,j12max,
     &                             AngularType12,angweight12,calctype,symmetry,verbosity)   
                           end do !ip12
                        end do    !m12
                     end do       !l12
                  end do          !s12
               end do             !j12
            end do                !mt12
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: output to file, to stdout(if wanted) and to mathematica friendly format (if wanted)
            call outputroutine(outUnitno,twoSnucl,extQnumlimit,
     &           Result,verbosity)
            
c     be a good boy and deallocate arrays. Compilers do that automatically for simple programs. Better safe than sorry.
            deallocate (Result, STAT=test ) ! test becomes nonzero if this fails
            if (test .ne. 0) stop "*** ERROR: Arrays Result(): Deallocation error. Abort."
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         end do                 ! Nangles
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
         write(*,*) "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"
      end do                    ! Nenergy
      close(outUnitno,iostat=test)
      if (test .ne. 0) stop "*** ERROR: Could not close output file!!!"
      
      write (*,*) '*** Wrote output to file: ',TRIM(outfile)
      
      stop
      end PROGRAM

      
