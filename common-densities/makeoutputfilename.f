c     hgrie Oct 2022: v2.0 fewbody-Compton
c     hgrie Nov 2023: added DATE: if specified in output file template, replaced by -date=yyyymmdd-hhmmss.nnn (milliseconds)
c     hgrie Aug 2020: added DENSITY, ORDER, NUCLEUS
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018: create name of output file for given energy and angle: XXX (ω), YYY (θ)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine makeoutputfilename(outfile,calctype,nucleus,descriptors,densityFileName,variedA,
     &     Elow,Ehigh,Einterval,thetaLow,thetaHigh,thetaInterval,verbosity)
c**********************************************************************
      IMPLICIT NONE
c**********************************************************************
      include '../common-densities/calctype.def'
c**********************************************************************
c     input variables
      character*500,intent(inout) :: outfile         ! name of output file
      character*500, intent(in)   :: densityFileName ! name of density file
      character*3,intent(in) ::  nucleus             ! name of nucleus to be considered, for output file name
      integer,intent(in) :: calctype,variedA
      character*200,intent(in) :: descriptors        ! additional descriptors of calculation
      real*8,intent(in)  :: Elow,Ehigh,Einterval
      real*8,intent(in)  :: thetaLow,thetaHigh,thetaInterval
      integer,intent(in) :: verbosity                !verbosity index for stdout hgrie June 2014
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     intrinsic variables

      integer               :: dummyint

      character*7,parameter :: densitytext = "DENSITY" ! placeholder for density name in input file
      character*7,parameter :: nucleustext = "NUCLEUS" ! placeholder for target nucleus name in input file
      character*5,parameter :: ordertext = "ORDER"     ! placeholder for order in output file
      character*3,parameter :: energytext = "XXX"      ! placeholder for energy/range in output file
      character*3,parameter :: angletext = "YYY"       ! placeholder for angle/range in output file
      character*4,parameter :: datetext = "DATE"       ! placeholder for productiondate and time in output file
      character*10          :: orderreplacement
      character*13          :: energyreplacement,anglereplacement
      character*8           :: dated
      character*10          :: timed
      character*24          :: datereplacement
      character*500         :: densityreplacement
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
c     replace ORDER in _original_ filename by order -- if none such, then just proceed
      if (calctype.eq.Odelta0) then
         orderreplacement = "Odelta0"
      else if (calctype.eq.Odelta2) then
         orderreplacement = "Odelta2"
      else if (calctype.eq.Odelta3) then
         orderreplacement = "Odelta3"
      else if (calctype.eq.Odelta4) then
         orderreplacement = "Odelta4"
      else if (calctype.eq.OQ4) then
         orderreplacement = "OQ4"
      else if (calctype.eq.varyAp) then
         write(orderreplacement,'(A,I1,A)') "VaryA",variedA,"p"
      else if (calctype.eq.varyAn) then
         write(orderreplacement,'(A,I1,A)') "VaryA",variedA,"n"
      end if
c     add descriptors to output filename
      do
         dummyint = INDEX(outfile,ordertext(:LEN_TRIM(ordertext))) ; if (dummyint == 0) EXIT
         outfile = outfile(:dummyint-1) // trim(orderreplacement) // trim(descriptors)
     &        // outfile(dummyint+LEN_TRIM(ordertext):)
      end do

c     replace DENSITY in  _original_ filename by name of potential in density filename
      do
         dummyint = INDEX(outfile,densitytext(:LEN_TRIM(densitytext))) ; if (dummyint == 0) EXIT
         if ( index(densityFileName,"-om=").ne.0 ) then
            densityreplacement =
     &           densityFileName(index(densityFileName,"-dens")+14+LEN_TRIM(nucleus):index(densityFileName,"-om=")) //
     &           densityFileName(index(densityFileName,"YYY-")+4:index(densityFileName,"-rho")-1)
         else if ( index(densityFileName,"-omega=").ne.0 ) then
            densityreplacement =
     &           densityFileName(index(densityFileName,"-dens")+14+LEN_TRIM(nucleus):index(densityFileName,"-omega=")) //
     &           densityFileName(index(densityFileName,"YYY-")+4:index(densityFileName,"-rho")-1)
         else
            write(*,*) "WARNING: Could not replace DENSITY string in output filename -- no >om=< or >omega=< in density filename."
            write(densityreplacement,"(I10)") time()
            write(*,*) "      ...so replacement will be: ", densityreplacement
         end if   
c         write(*,*) densityreplacement
         outfile =outfile(:dummyint-1) // densityreplacement(:LEN_TRIM(densityreplacement))
     &        // outfile(dummyint+LEN_TRIM(densitytext):)
      end do
      
c     replace XXX in _original_ filename by energy/range -- if none such, then just proceed
      if (Elow.eq.Ehigh) then
         write(energyreplacement,'(I0.3)') NINT(Elow)
      else
         write(energyreplacement,'(I0.3,A,I0.3,A,I0.3)')
     &        NINT(Elow),"to",NINT(Ehigh),"in",NINT(Einterval)
c         write(*,*) energyreplacement
      end if
      do
         dummyint = INDEX(outfile,energytext(:LEN_TRIM(energytext))) ; if (dummyint == 0) EXIT
         outfile = outfile(:dummyint-1) // energyreplacement(:LEN_TRIM(energyreplacement))
     &        // outfile(dummyint+LEN_TRIM(energytext):)
      end do
      
c     replace YYY in _original_ filename by angle/range -- if none such, then just proceed
      if (thetaLow.eq.thetaHigh) then
         write(anglereplacement,'(I0.3)') NINT(thetaLow)
      else
         write(anglereplacement,'(I0.3,A,I0.3,A,I0.3)')
     &        NINT(thetaLow),"to",NINT(thetaHigh),"in",NINT(thetaInterval)
c         write(*,*) anglereplacement
      end if
      do
         dummyint = INDEX(outfile,angletext(:LEN_TRIM(angletext))) ; if (dummyint == 0) EXIT
         outfile =outfile(:dummyint-1) // anglereplacement(:LEN_TRIM(anglereplacement))
     &        // outfile(dummyint+LEN_TRIM(angletext):)
      end do
      
c     replace DDD in _original_ filename by date and time -- if none such, then just proceed
      call DATE_AND_TIME(dated,timed)
      write(datereplacement,'(4A)') "date=",dated,"-",timed
c     write(*,*) datereplacement
      do
         dummyint = INDEX(outfile,datetext(:LEN_TRIM(datetext))) ; if (dummyint == 0) EXIT
         outfile =outfile(:dummyint-1) // datereplacement(:LEN_TRIM(datereplacement))
     &        // outfile(dummyint+LEN_TRIM(datetext):)
      end do
      
c     replace NUCLEUS in filename by nucleus name derived from input filename
      do
         dummyint = INDEX(outfile,nucleustext(:LEN_TRIM(nucleustext))) ; if (dummyint == 0) EXIT
         outfile = outfile(:dummyint-1) // nucleus(:LEN_TRIM(nucleus))
     &        // outfile(dummyint+LEN_TRIM(nucleustext):)
      end do

      outfile = TRIM(outfile)
cccccccccccccccccccccccccccccccccccccccccccc      
      write (*,*) 'Write output to file: ',TRIM(outfile)
      
      if (verbosity.eq.1000) continue
      
      return
      
      end

      
