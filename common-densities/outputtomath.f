c     hgrie Oct 2022: v2.0 fewbody-Compton
c     hgrie Aug 2020: produce mathematica-friendly output to stdout of the first 2(2S+1)² MEs.
c     These are not related by time-reversal symmetry.
c     Also includes useful sub-subroutines for that.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c     hgrie Aug 2020: new
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine outputtomath(Result11,Result1m1,Resultm11,Resultm1m1,twoSnucl,twoMzplimit,verbosity)
c**********************************************************************
      IMPLICIT NONE
c**********************************************************************
c     input variables

      integer,intent(in)    :: twoSnucl,twoMzplimit
      complex*16,intent(in) :: Result11(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)     
      complex*16,intent(in) :: Result1m1(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16,intent(in) :: Resultm11(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16,intent(in) :: Resultm1m1(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      
      integer,intent(in)    :: verbosity
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     intrinsic variables

      integer :: twoMzp,twoMz,twoMzlimit
c     for mathematica-friendly output, define numbers as strings. not elegant, but works
      character(len=64) string
      character(len=4*(twoSnucl+1)**2*68) longstring
      longstring = "" ! initialise auxiliary
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      if (verbosity.eq.1000) continue
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (twoMzplimit.eq.0) then ! only symmetry-related amplitudes if symmetry invoked
         write(*,'(A,I2,A)') "Mathematica-friendly output: first ",2*(twoSnucl+1)**2," amplitudes (not related by symmetries):"
         write(*,*)    "   [Sequence counting down from max values: Mzp∈[S;0], Mz∈[S;-S] unless Mzp=0: Mz∈[S;0], λp,λ=±1]"
         string = ""            ! initialise
         
         do  twoMzp=twoSnucl,0,-2
c     for Mzp=0, run only over Mz>=0 -- that's still 2 more than necessary, but good enough -- see cure below
            if (twoMzp.eq.0) then
               twoMzlimit = 0
            else
               twoMzlimit = -twoSnucl
            end if   
            do twoMz=twoSnucl,twoMzlimit,-2
c     λp=1,λ=1
               write(string,'(SP,"(",E24.18,",",E24.18,")")') Result11(twoMzp,twoMz)
               call ConvertComplexToMath(string)
               longstring = trim(adjustl(longstring)) // string // ","
               call StripSpaces(longstring)
c     λp=-1,λ=1
               write(string,'(SP,"(",E24.18,",",E24.18,")")') Result1m1(twoMzp,twoMz)
               call ConvertComplexToMath(string)
               longstring = trim(adjustl(longstring)) // string // ","
               call StripSpaces(longstring)
c     noly write the opposite helicities if neither Mzp nor Mz are zero
               if ((twoMzp.eq.0).and.(twoMz.eq.0)) then
                  continue
               else
c     λp=1,λ=-1
                  write(string,'(SP,"(",E24.18,",",E24.18,")")') Resultm11(twoMzp,twoMz)
                  call ConvertComplexToMath(string)
                  longstring = trim(adjustl(longstring)) // string // ","
                  call StripSpaces(longstring)
c     λp=-1,λ=-1
                  write(string,'(SP,"(",E24.18,",",E24.18,")")') Resultm1m1(twoMzp,twoMz)
                  call ConvertComplexToMath(string)
                  longstring = trim(adjustl(longstring)) // string // ","
                  call StripSpaces(longstring)
               end if           ! twoMzp=twoMz=0
            end do              ! twoMz
         end do                 ! twoMzp
      else                      ! output ALL amplitudes
c output ALL amplitudes when symmetry NOT invoked        
         write(*,'(A,I2,A)') "Mathematica-friendly output: all ",4*(twoSnucl+1)**2," amplitudes (some related by symmetries):"
         write(*,*)    "   [Sequence counting down from max values: Mzp∈[S;-S], Mz∈[S;-S]"
         string = ""            ! initialise
         do  twoMzp=twoSnucl,-twoSnucl,-2
            do twoMz=twoSnucl,-twoSnucl,-2
c     λp=1,λ=1
               write(string,'(SP,"(",E24.18,",",E24.18,")")') Result11(twoMzp,twoMz)
               call ConvertComplexToMath(string)
               longstring = trim(adjustl(longstring)) // string // ","
               call StripSpaces(longstring)
c     λp=-1,λ=1
               write(string,'(SP,"(",E24.18,",",E24.18,")")') Result1m1(twoMzp,twoMz)
               call ConvertComplexToMath(string)
               longstring = trim(adjustl(longstring)) // string // ","
               call StripSpaces(longstring)
c     λp=1,λ=-1
               write(string,'(SP,"(",E24.18,",",E24.18,")")') Resultm11(twoMzp,twoMz)
               call ConvertComplexToMath(string)
               longstring = trim(adjustl(longstring)) // string // ","
               call StripSpaces(longstring)
c     λp=-1,λ=-1
               write(string,'(SP,"(",E24.18,",",E24.18,")")') Resultm1m1(twoMzp,twoMz)
               call ConvertComplexToMath(string)
               longstring = trim(adjustl(longstring)) // string // ","
               call StripSpaces(longstring)
            end do              ! twoMz
         end do                 ! twoMzp
      end if                    ! twoMzplimit symmetry question

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

      
