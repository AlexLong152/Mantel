c
c
c     contains 
c              outputroutine()
c
c     use to be part of usesymmetry+writeoutput-densities.f.
c     usesymmetry() eliminated
c     function arguments changed
c      
c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c
c     hgrie May 2018:
c     merges old common/useparity.f and parity/output routines at the end of main.*.f
c     contains 
c              outputroutine()
c
c     twoSmax/twoMz dependence: via array sizes & do-loops

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c
c     hgrie Aug 2020: added call to outputtomath() which produces
c     mathematica-friendly output of all MEs.
      
c     hgrie May 2018:
c     merges old common/useparity.f and parity/output routines at the end of main.*.f
c     contains useparityroutine()
c              outputroutine()
c
c     twoSmax/twoMz dependence: via array sizes & do-loops

c     hgrie May 2018: used to be part of 3HeCompt/common
c     now part of common-densities, backward compatibility deliberately broken
c     no changes yet
c      
c     hgrie May 2018: more extensive description of changes in main.*.f
c     rewritten such that magnetic quantum numbers are "2*Mz" etc
c     changed phases of parity operation in useparityroutine() to cover 3He and deuteron,
c      
c     Implemented symmetry for arbitrary nucleon spin:
c     Use Mzp>=0, and for Mzp=0, run only over Mz>=0
c     -- that's still 2 more than necessary since ME(+0->+0) = ME(-0->-0) and ME(+0->-0) = ME(-0->+0)
c     but it's good enough, saving lots of CPU time.
c     see manuscript "Compton Densities Approach" pp11-12
c*************************************************************************************
c*************************************************************************************
c*************************************************************************************

      subroutine outputroutine(outUnitno,twoSnucl,extQnumlimit,
     &     Result,verbosity)
c     hgrie May 2018: new routines, outsourced from main.*.f
c
c**********************************************************************
c     
      implicit none
c     
c**********************************************************************
c     
c     INPUT/OUTPUT VARIABLE:
c     
c     amplitudes with photon helicities 
c     
      complex*16,intent(in) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      
      integer,intent(in) :: extQnumlimit
      integer,intent(in) :: twoSnucl
      integer,intent(in) :: outUnitno
      
      integer,intent(in) :: verbosity
      
      integer :: extQnum,twoMzp,twoMz
c     
c**********************************************************************
c     
      if (verbosity.eq.1000) continue
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output as human readable number to stdout            
      if (verbosity.ge.0) then
         do twoMzp=twoSnucl,-twoSnucl,-2
            do twoMz=twoSnucl,-twoSnucl,-2
               do extQnum=1,extQnumlimit
                  write (*,'(A,I4,A,I4,A,I4,A,F24.19,SP,F24.19," i")') !E30.19 for exponential form 0.123.....E-56
     &                 "Result(exQnum=",extQnum,",twoMzp=",twoMzp,", twoMz=",twoMz,"): ",Result(extQnum,twoMzp,twoMz)
               end do           ! extQnum
            end do              ! twoMz
         end do                 ! twoMzp
      end if                    ! verbosity
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     output to file
      do twoMzp=twoSnucl,-twoSnucl,-2
         do twoMz=twoSnucl,-twoSnucl,-2
            do extQnum=1,extQnumlimit
               write (outUnitno,*) Result(extQnum,twoMzp,twoMz)
            end do              ! extQnum
         end do                 ! twoMz
      end do                    ! twoMzp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c     hgrie Aug 2020: if so wanted, output first independent MEs also to screen in a form that can directly be pasted into mathematica
      if (verbosity.ge.0) call outputtomath(Result,twoSnucl,extQnumlimit,verbosity)

      
      return
      end                       ! outputroutine
