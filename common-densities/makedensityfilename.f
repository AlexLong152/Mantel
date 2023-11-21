c     hgrie Oct 2022: v2.0 fewbody-Compton
c     hgrie Aug/Sep 2020: rewrote to deal with extracting densities from a .gz,
c     if input file name ends with ".gz" or if direct filename not found,
c     and of downloading files form a server
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018: outsourced from main.*via*density.f into its own subroutine
c     hgrie June 2017: create name of 1Ndensity file for given energy and angle, unpack it
c     define correct formats for energy and angle
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c
c     hgrie Aug 2020: rewrote to deal with extracting densities from a .gz, if input file name ends with ".gz" or if direct filename not found
c     hgrie May 2018: used to be part of 3HeCompt/common
c     now part of common-densities, backward compatibility deliberately broken
c     no changes yet
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine makedensityfilename(densityFileName,Egamma,thetacm,rmDensityFileLater,verbosity)
c**********************************************************************
      IMPLICIT NONE

c**********************************************************************
      include '../common-densities/constants.def'
c**********************************************************************
c     input variables
      character*500,intent(inout) :: densityFileName ! name of density file
      real*8,intent(in)           :: Egamma, thetacm ! energy [MeV] & angle [rad] of input file: input to subroutine
      integer,intent(in)          :: verbosity       ! verbosity index for stdout hgrie June 2014
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output variables
c     if density file generated from a .gz, delete that temporary file after each energy/angle
c     if downloaded and .gz, also delete the download.
c     0: do not delete; 1: delete un-gz'd file; 2: delete downloaded and un-gz'd file 
      integer,intent(out)         :: rmDensityFileLater 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     intrinsic variables

      integer               :: dummy
      character*3,parameter :: text1 = "XXX" ! placeholder for energy of the density file in input file
      character*3,parameter :: text2 = "YYY" ! placeholder for angle of the density file in input file
      
      character*500         :: originaldensityFileName ! original name of density file
      character*8           :: replacement1,replacement2
      logical               :: testtf                  ! generic logical for testing i/o
      integer               :: test                    ! generic integer for testing i/o
      character*20          :: testline                ! generic string for testing i/o

      integer*4             :: startdownload, stopdownload ! start and stop time if density downloaded from http
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      testtf=.True.             ! pre-set test variable for density file existence
      rmDensityFileLater=0      ! pre-set test variable if density file extracted from .gz (to delete it after use)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*) "   Now looking for density file..."
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc           
c     hgrie Aug 2020: keep original filename: needed if extraction from .tgz below
      originaldensityFileName = densityFileName
c      
      write(replacement1,'(ES8.2)') Egamma
      write(replacement2,'(ES8.2)') thetacm*180.0/Pi
c     replace XXX in _original_ filename by energy -- if none such, then just proceed
      do
         dummy = INDEX(densityFileName,text1(:LEN_TRIM(text1))) ; if (dummy == 0) EXIT
         densityFileName = densityFileName(:dummy-1) // replacement1(:LEN_TRIM(replacement1))
     &        // densityFileName(dummy+LEN_TRIM(text1):)
      end do
c     replace YYY in filename by angle -- if none such, then just proceed
      do
         dummy = INDEX(densityFileName,text2(:LEN_TRIM(text2))) ; if (dummy == 0) EXIT
         densityFileName = densityFileName(:dummy-1) // replacement2(:LEN_TRIM(replacement2))
     &        // densityFileName(dummy+LEN_TRIM(text2):) 
      end do
      densityFileName = TRIM(densityFileName)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*)    "   Density file name generated from energy & angle: ", TRIM(densityFileName)

c     hgrie Sep 2020: if filename starts with "https://", then interpret as weblink and download
      if (densityFileName(:8).eq."https://") then
         write(*,*) "   URL given for density file -- attempting download from ",
     &        TRIM(densityFileName(:index(densityFileName,"/",BACK=.True.)))
         startdownload = time()
         call EXECUTE_COMMAND_LINE("curl -b license_value=license%2Bis%2Baccepted "//densityFileName//" -O",
     &        WAIT=.True.,EXITSTAT=test)
         stopdownload = time()
         write(*,*) "      Download took ",stopdownload-startdownload," s. -- Should be subtracted to get actual runtime."
         if (test.ne.0) stop "*** ERROR: Could not download density file: Abort."
         densityFileName = densityFileName(index(densityFileName,"/",BACK=.True.)+1:)
c     Check if the downloaded file is actually a .gz or .dat or .h5 file, or just an error message.
c     I do that by asking of the downloaded file contains the string "html", which is indicative of a html message.
c     Fortran does not allow to pipe a result directly from shell to programme, so I write the result of the query to a
c     temporary file, read it into fortran, and then depete that file again. Mit der Brust durchs Auge.
         call EXECUTE_COMMAND_LINE("file -b "//TRIM(densityFileName)//" > "//TRIM(densityFileName)//".tmp",
     &        WAIT=.True.,EXITSTAT=test)
         open(18,FILE=TRIM(densityFileName)//".tmp",STATUS='OLD')
         read (18,*) testline
         if (index(testline,'HTML').ne.0) then
            stop "*** ERROR: downloaded file is a html file -- likely an error message from the server."
         else if ((index(testline,'gzip').ne.0).or.(index(testline,'Hierarchical').ne.0).or.(index(testline,'ASCII').ne.0)) then
            write(*,*) "      Download successful, file type apparently supported (.gz, .dat for onebody or .h5 for twobody)."
            continue
         else
            write(*,*) "  *** Download successful, but file type not .gz/dat/h5: not supported? -- Proceed, keep fingers crossed."
         end if
         close(18)
         call EXECUTE_COMMAND_LINE("rm -rf "//TRIM(densityFileName)//".tmp",WAIT=.True.)
         rmDensityFileLater = 2 ! set flag: both downloaded and unzipped version is removed after each energy/angle
      end if
      
c     hgrie Aug 2020: if .dat or .h5 ending, the look for it. 
      if ((densityFileName(LEN_TRIM(densityFileName)-3:LEN_TRIM(densityFileName)).eq.".dat").or.
     &     (densityFileName(LEN_TRIM(densityFileName)-2:LEN_TRIM(densityFileName)).eq.".h5")) then
         INQUIRE(FILE=densityFileName, EXIST=testtf)
         if (.not.testtf) then 
            write(*,*) "   Density file does not exist directly -- trying .gz version."
            densityFileName = TRIM(densityFileName)//".gz"
         else
            write(*,*) "      Density file exists, good."
         end if
      else if (densityFileName(LEN_TRIM(densityFileName)-1:LEN_TRIM(densityFileName)).ne."gz") then
         write(*,*) "*** ERROR: Density file name has neither ending .dat nor .h5 nor *gz -- do not know what to do: Abort."
         stop "*** Try maybe download from website? Density file format in input file: http://<website>/<filename>"
      end if
c     hgrie June 2022: if ending .gz or .dat/.h5 does not exist, look for .tar.gz file and tar xzf
      if (densityFileName(LEN_TRIM(densityFileName)-6:LEN_TRIM(densityFileName)).eq.".tar.gz") then
         write(*,*) "   Density file should be .tar.gz file."
         INQUIRE(FILE=densityFileName, EXIST=testtf)
         if (.not. testtf) stop "*** ERROR: .tar.gz'ipped density file does not exist: Abort."
         write(*,*) "      .tar.gz'ipped density file exists, good. Now tar xzf for density file."
c     needs to be tarred into the original directory. tar usually extracts into current working directory,
c     so change that by extracting name of density directory from filename if filename contains "/". Not elegant, not nice.
         if (index(densityFileName,"/").eq.0) then
            call EXECUTE_COMMAND_LINE("tar xzf "//densityFileName,WAIT=.True.,EXITSTAT=test)
         else
            call EXECUTE_COMMAND_LINE("tar xzf "//densityFileName//" --directory "
     &           //densityFileName(1:index(densityFileName,"/",.true.)),WAIT=.True.,EXITSTAT=test)
         end if   
         if (test.ne.0) stop "*** ERROR: Could not tar xzf density file: Abort."
         densityFileName = densityFileName(:LEN_TRIM(densityFileName)-7)
         INQUIRE(FILE=densityFileName, EXIST=testtf)
         if (.not. testtf) stop "*** ERROR: Density file does not exist: Abort."
         write(*,*) "      tar xzf'ed density file exists, good."
         if (rmDensityFileLater.ne.2) then !only reset if it's not already set by download 
            rmDensityFileLater = 1 ! set flag: only unzipped version is removed after each energy/angle
         end if
      end if
c     hgrie June 2022: if ending .tar.ga or .gz or .dat/.h5 does not exist, look for .tgz file and tar xzf
      if (densityFileName(LEN_TRIM(densityFileName)-3:LEN_TRIM(densityFileName)).eq.".tgz") then
         write(*,*) "   Density file should be .tgz file."
         INQUIRE(FILE=densityFileName, EXIST=testtf)
         if (.not. testtf) stop "*** ERROR: .tgz'ipped density file does not exist: Abort."
         write(*,*) "      .tgz'ipped density file exists, good. Now tar xzf for density file."
c     needs to be tarred into the original directory. tar usually extracts into current working directory,
c     so change that by extracting name of density directory from filename if filename contains "/". Not elegant, not nice.
         if (index(densityFileName,"/").eq.0) then
            call EXECUTE_COMMAND_LINE("tar xzf "//densityFileName,WAIT=.True.,EXITSTAT=test)
         else
            call EXECUTE_COMMAND_LINE("tar xzf "//densityFileName//" --directory "
     &           //densityFileName(1:index(densityFileName,"/",.true.)),WAIT=.True.,EXITSTAT=test)
         end if   
         if (test.ne.0) stop "*** ERROR: Could not tar xzf density file: Abort."
         densityFileName = densityFileName(:LEN_TRIM(densityFileName)-4)
         INQUIRE(FILE=densityFileName, EXIST=testtf)
         if (.not. testtf) stop "*** ERROR: Density file does not exist: Abort."
         write(*,*) "      tar xzf'ed density file exists, good."
         if (rmDensityFileLater.ne.2) then !only reset if it's not already set by download 
            rmDensityFileLater = 1 ! set flag: only unzipped version is removed after each energy/angle
         end if
      end if
c     hgrie Aug 2020: if ending .gz or .dat/.h5 does not exist, look for .gz file and unzip
      if (densityFileName(LEN_TRIM(densityFileName)-2:LEN_TRIM(densityFileName)).eq.".gz") then
         write(*,*) "   Density file should be .gz file."
         INQUIRE(FILE=densityFileName, EXIST=testtf)
         if (.not. testtf) stop "*** ERROR: .gz'ipped density file does not exist: Abort."
         write(*,*) "      .gz'ipped density file exists, good. Now unzipping for density file."
         call EXECUTE_COMMAND_LINE("gunzip --keep "//densityFileName,WAIT=.True.,EXITSTAT=test)
         if (test.ne.0) stop "*** ERROR: Could not gunzip density file: Abort."
         densityFileName = densityFileName(:LEN_TRIM(densityFileName)-3)
         INQUIRE(FILE=densityFileName, EXIST=testtf)
         if (.not. testtf) stop "*** ERROR: Density file does not exist: Abort."
         write(*,*) "      Unzipped density file exists, good."
         if (rmDensityFileLater.ne.2) then !only reset if it's not already set by download 
            rmDensityFileLater = 1 ! set flag: only unzipped version is removed after each energy/angle
         end if
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      if (verbosity.eq.1000) continue
      
      return
      
      end

      
