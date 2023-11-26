c     hgrie Oct 2022: v2.0 fewbody-Compton
c     hgrie Aug 2020: produce mathematica-friendly output to stdout of the first 2(2S+1)² MEs.
c     These are not related by time-reversal symmetry.
c     Also includes useful sub-subroutines for that.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     hgrie Aug 2020: new
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outputtomath(Result,twoSnucl,extQnumlimit,verbosity)
c**********************************************************************
      IMPLICIT NONE
c**********************************************************************
c     input variables

      integer,intent(in)    :: extQnumlimit
      integer,intent(in)    :: twoSnucl
      complex*16,intent(in) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      
      integer,intent(in)    :: verbosity
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     intrinsic variables

      integer :: extQnum,twoMzp,twoMz
c     for mathematica-friendly output, define numbers as strings. not elegant, but works
      character(len=64) string
      character(len=4*(twoSnucl+1)**2*68) longstring
      longstring = "" ! initialise auxiliary
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      if (verbosity.eq.1000) continue
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c output ALL amplitudes when symmetry NOT invoked        
      write(*,'(A,I2,A,I2,A,I4,A)') "Mathematica-friendly output: all ",extQnumlimit,
     &     "*",(twoSnucl+1)**2," = ", extQnumlimit*(twoSnucl+1)**2," amplitudes (some related by symmetries):"
      write(*,*)    "   [Sequence counting down from max values: extqnum∈[1;extQnumlimit], Mzp∈[S;-S], Mz∈[S;-S]"
      string = ""               ! initialise
      do  twoMzp=twoSnucl,-twoSnucl,-2
         do twoMz=twoSnucl,-twoSnucl,-2
            do extQnum=1,extQnumlimit
               write(string,'(SP,"(",E24.18,",",E24.18,")")') Result(extQnum,twoMzp,twoMz)
               call ConvertComplexToMath(string)
               longstring = trim(adjustl(longstring)) // string // ","
               call StripSpaces(longstring)
            end do              ! extQnum
         end do                 ! twoMz
      end do                    ! twoMzp

      longstring = '{' // trim(adjustl(longstring)) // '}'    
      longstring = longstring(:index(longstring,",}")-1) // "}"
      write(*,*) '        ',trim(adjustl(longstring))    
      
      return
      
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convert a string which used to be a Fortran Real into a mathematica-readable string,
c     e.g. 1.2345E003 => 1.2345*10*(003)
      subroutine ConvertRealToMath(string)
      character(len=*),intent(inout) :: string
      if ( index(string,"E").ne.0 ) then
         string = trim(string(1:index(string,"E")-1) // '*10^(' // string(index(string,"E")+1:)) // ')'
      else
         string = trim(string)
      end if
      string = trim(adjustl(string))
      call StripSpaces(string)
      end ! subroutine ConvertRealToMath

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convert a string which used to be a Frortran Compex into a mathematica-readable string,
c     e.g. (1.2E3,3.4E3) => 1.2*10^3+3.4*10^3*I
      subroutine ConvertComplexToMath(string)
      character(len=*),intent(inout) :: string
      character(len=30) :: restring,imstring
      write(restring,*) trim(adjustl(string(index(string,"(")+1:index(string,",")-1)))
      write(imstring,*) trim(adjustl(string(index(string,",")+1:index(string,")")-1)))
      call ConvertRealToMath(restring)
      call ConvertRealToMath(imstring)
      string = trim(adjustl(restring // imstring // "*I"))
      call StripSpaces(string)
      end ! subroutine ConvertComplexToMath
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine StripSpaces(string)
      character(len=*),intent(inout) :: string
      integer :: stringLen 
      integer :: last, actual
      
      stringLen = len (string)
      last = 1
      actual = 1
      
      do while (actual < stringLen)
         if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
         else
            last = last + 1
            if (actual.lt.last) actual = last
         endif
      end do
      
      end ! subroutine StripSpaces

      
