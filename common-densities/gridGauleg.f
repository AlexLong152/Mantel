cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Part of MANTLE code for One/Twobody Contributions to Few-Nucleon Processes Calculated Via 1N/2N-Density Matrix
c     NEW Nov 2023: v1.0 Alexander Long/hgrie 
c               Based on Compton density code v2.0: D. Phillips/A. Nogga/hgrie starting August 2020/hgrie Oct 2022
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CONTAINS SUBROUTINES:
c              AnglePtsWts  : Generate Nbins evenly spaced bins in interval (xlow,xhigh) and fill each bin with Nord gaussian integ pts/wts
c              RadialPtsWts : Generate Radial integration pts and wts in MeV (same units as c):
c                             divide p = (0,infty) into Nbins according to spline or tangent mapping with Np pts and wts each.
c              gridGauleg   : Fill each bin in a grid with Nord gaussian integration points(&wts)
c     File contains subroutines for dividing an interval into bins and filling each bin with gaussian pts/wts.
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
      subroutine AnglePtsWts(Nord,Nbins,NXMX,xlow,xhigh, xX,dxX,Nx,verbosity)
c--------------------------------------------------------------------
c     Generate Nbins evenly spaced bins in interval (xlow,xhigh) and
c     fill each bin with Nord gaussian integ pts/wts
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT VARIABLES:

      real*8,intent(out) :: xX(NXMX),dxX(NXMX)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INPUT VARIABLES:

      integer,intent(in) :: Nord,Nbins,NXMX,Nx 
      real*8,intent(in)  :: xlow,xhigh
      
      integer,intent(in) :: verbosity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOCAL VARIABLES:
      
      integer NBINSMX
      parameter(NBINSMX=64)     ! common to AnglePtsWts & RadialPtsWts
      integer i
      real*8 binWalls(NBINSMX)
      common/binWork/binWalls
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      do i=1,Nbins+1
         binWalls(i) = xlow + (xhigh-xlow)*(i-1.)/Nbins
      end do
      call gridGauleg(binWalls,Nbins+1,Nord,NXMX,.false.,0.d0,0.d0,
     &     xX,dxX,Nx,verbosity)
c     do i=1,Nx
c     write(6,*) xX(i),dxX(i)
c     end do
c     
      if (verbosity.eq.1000) continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine RadialPtsWts(Nord,Nbins,NRADIALMX,PI, RadialType,
     &     t0,pbar,tdel, c_tail,  pX,dpX,Np,verbosity)
c-------------------------------------------------------------------
c   Generate Radial integration pts and wts in MeV (same units as c).
c     Divide p = (0,infty) into Nbins according to spline or tangent
c     mapping with Np pts and wts each.
c INPUT:
c	Nbins,Np................# bins, # integ pts per bin
c	NRADIALMX...............Maximum size of Nbins*Np
c	PI......................physical constant, 3.14...
c	binWalls................array work space (no input)
c	RadialType..............1=splMapbins, 2...10=tanmapbins, 
c				11=splmappts, >11=tanMapPts
c	t0,pbar,tdel............spline/knot parameters
c _tail..................tan mapping parameter (~1000MeV for deut)
c
c OUTPUT:
c	pX,wtpX.................integ pts & wts
c	Np......................Nbins*Np = total # of integ pts & wts
c--------------------------------------------------------------------
      implicit none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer NBINSMX
      parameter(NBINSMX=64)     ! common to AnglePtsWts & RadialPtsWts
      integer Nord,Nbins,NRADIALMX,RadialType, Np, i
      integer,intent(in) :: verbosity
      double precision PI, binWalls(NBINSMX), t0,pbar,tdel, c_tail
      double precision pX(NRADIALMX),dpX(NRADIALMX), xi,cosx,sinx
      logical tailTF
      common/binWork/binWalls
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      if (Nbins.gt.NBINSMX) then
         write(6,*) 'ERROR: too many bins. Nbins = ',Nbins
         stop
      endif
c     
c-----map bins into (0,infty), then fill bins with gauss pts/wts
      if (RadialType.lt.10) then
c     
         if (RadialType.eq.1) then ! old spl scheme
c---------map (-1,1) bins to (0,infty) by spline mapping. Fill w 
c     gauss pts/wts
            binWalls(1) = t0
            do i=2,Nbins+1
               xi = -cos( (2.*I-3.)/(2.*Nbins)*PI ) ! xi \in (-1,1)
               binWalls(i) = tdel + pbar * dsqrt((1.+xi)/(1.-xi))
            end do
            tailTF = .false.	! binWall(Nbins+1) = finite
         else                   !if (RadialType.eq.3) then ! old wally/nkd scheme
c---------map (0,1) bins to (0,infty) by tan mapping and fill w 
c     gauss pts/wts
            do i=1,Nbins
               xi = (i-1.)/Nbins ! xi \in (0,1)
               binWalls(i) = c_tail * tan(PI/2.*xi)
c     write(6,*) ' binwalls', binWalls(i)
            end do
            tailTF = .true.     ! binWall(Nbins+1) = infty
         endif
         call gridGauleg(binWalls,Nbins+1,Nord,NRADIALMX,tailTF,
     &        c_tail,PI,pX,dpX,Np,verbosity)
c     
c-----fill bins(0,1) with gauss pts/wts, then map into (0,infty)
      else
         do i=1,Nbins+1
            binWalls(i) = (i-1.)/Nbins
         end do
         tailTF = .false.
         call gridGauleg(binWalls,Nbins+1,Nord,NRADIALMX,tailTF,
     &        c_tail,PI, pX,dpX,Np,verbosity)
c     
         if (RadialType.eq.11) then ! new spl scheme
c---------map (0,1) pts&wts to (0,infty) by spline mapping
c     bins here do NOT correspond to the spl mapping above!
            do i=1,Np
               cosx = cos(pX(i))
               sinx = sin(pX(i))
               pX(i)  = pbar * dsqrt((1.-cosx)/(1.+cosx))
               dpX(i) = pbar/2. /(1.-cosx)**.5/(1.+cosx)**1.5 *sinx*PI*
     &              dpX(i)
            end do
         else                   ! new wally/nkd scheme
c---------map (0,1) pts&wts to (0,infty) by tan mapping
            do i=1,Np
               xi = pX(i)
               pX(i)  = c_tail * tan(PI/2.*xi)
               dpX(i) = c_tail * PI/2.* dpX(i) / cos(PI/2.*xi)**2
            end do
         endif
c     
      endif
c     
c     do i=1,Np
c     write(6,*) pX(i),dpX(i)
c     end do
c     
      if (verbosity.eq.1000) continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gridGauleg(binWalls,NbinWalls,Nord,NMX,tailTF,
     &     c_tail,PI,xX,dxX,NxPts,verbosity)
c--------------------------------------------------------------------
c     Fill each bin in a grid with Nord gaussian integration points(&wts).
c--------------------------------------------------------------------
c INPUT:
c	binWalls........array defining walls of bins
c	NbinWalls.......# bin walls = #bins+1
c	Nord............# of gaussian pts/wts to put in each bin
c	NMX.............maximum value of NxPts=Nord*Nbins
c	tailTF..........true -> use tan mapping for last bin; 
c                             lastwall=infty
c	c_tail..........parameter controling tan mapping
c	PI..............3.1414... (used in tan mapping)
c OUTPUT:
c	xX......array of Nord gaussian pts in each of NbinWalls-1 
c		bins defined by (binwalls(1),...,binwalls(NbinWalls))
c	dxX.....gaussian wts corresponding to xX
c	NxPts...Nord*(NbinWalls-1) = # pts in xX
c--------------------------------------------------------------------
      implicit none
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer NxPts, NbinWalls,Nord,NMX,  iNord,iGrid,i,NORDMX
      double precision xX(NMX),dxX(NMX), binWalls(NbinWalls),
     &     c_tail,PI,xxx
      integer,intent(in) :: verbosity
      logical tailTF
      parameter(NORDMX=100)
      double precision PtsGauss(NORDMX),WtsGauss(NORDMX)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      if (Nord.gt.NORDMX) then
         write(6,*) 'ERROR: Nord > NORDMX'
         stop
      elseif (Nord*(NbinWalls-1).gt.NMX) then
         write(6,*) 'ERROR: Nord*Nbins > NMX; nord,Nbins= ', nord,
     &        nbinwalls-1
         stop
      endif
c     
c-----Generate standard gaussian pts and wts between 0 and 1
      call gauleg(0.D0,1.D0,PtsGauss,WtsGauss,Nord,verbosity)
c     
c-----Fill each grid interval with Nord gaussian integration pts &wts
c-----pts & wts for each bin (except possibly last)
      do iGrid=1,NbinWalls-1
         if ((.not. TailTF) .or. iGrid.ne.NbinWalls-1) then
            do iNord=1,Nord
               i= (iGrid-1)*Nord + iNord
               xX(i) = binWalls(iGrid)
     &              +(binWalls(iGrid+1)-binWalls(iGrid)) * 
     &              PtsGauss(iNord)
               dxX(i) = (binWalls(iGrid+1)-binWalls(iGrid)) * 
     &              WtsGauss(iNord)
            end do
         else                   ! tail, last bin ends at infty
c     
c     Use tangent mapping (0-infty) for last bin
c     May later want to input c_tail as variable
c     c_tail = 1000.
c     
            do inord=1,Nord
               i= (iGrid-1)*Nord + iNord
               xxx = pi/2.d0 * PtsGauss(iNord) ! 0<xx<pi/2
               xX(i) = binWalls(iGrid) + c_tail*tan(xxx)
               dxX(i)= PI/2.d0 *c_tail*WtsGauss(iNord)/(cos(xxx)**2)
            end do
         endif
      end do
      NxPts= Nord * (NbinWalls-1)
c     
      if (verbosity.eq.1000) continue
      return
      end
c--------------------------------------------------------------------
