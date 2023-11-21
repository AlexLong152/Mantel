c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2018
c     Version 7 hgrie Feb 2013 -- compatibility to previous versions broken
c     consult readME for explanation of variables, kinematics, version history, original authors etc.       
c     
c     twoSmax/twoMz dependence: via array sizes & do-loops 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     TO do:
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c
c     hgrie Oct 2022: Clarification: When usesymmetry() is set (twoMzplimit=0) and Mzp=Mz=0,
c     then the below STILL constructs 2*(2Snucl+1)² + 2 spherical ResultXXYY, namely all 4 amplitudes for Mzp=Mz=0,
c     instead of the 2*(2Snucl+1)² that would be minimum:
c     Two are superfluous, e.g. Result11 and Result1m1,
c     from which follow Resultm1m1=Result11 and Resultm11=Result1m1.
c     We can live with that overcounting because it does not affect computation time.
c     If one wanted t ocure this, one would have to make this subroutine dependent on twoMzplimit --
c     why bother?
c     
c     hgrie May 2018: used to be part of 3HeCompt/common
c     now part of common-densities, backward compatibility deliberately broken
c     no changes yet
c     hgrie May 2018: more extensive description of changes in main.*.f
c            rewritten such that magnetic quantum numbers are "2*Mz" etc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine transcarttospher(Result11,Result1m1,Resultm11,Resultm1m1,
     &     twoSnucl,Resultxx,Resultxy,Resultyx,Resultyy,verbosity)
c     
c**********************************************************************
c     
c     Transform the Cartesian amplitudes to spherical ones
c     
c     Written by Deepshikha 9/07 for deuteron; mdified for 3He hgrie Oct 2014
c     
c***********************************************************************
      implicit none
c**********************************************************************
      include 'params.def'
      include 'calctype.def'
      include 'constants.def'
c----------------------------------------------------------------------
c     
c     Compton amplitude variables: Target spin indices are (twoMzp,twoMz)
c     
c     Resultxx,Resultxy,Resultyx,Resultyy:     Deuteron Compton amplitudes for specific photon 
c                                              polarizations in cartesian basis for photon polarization
c     Result11,Result1m1,Resultm11,Resultm1m1: same in helicities/spherical basis for photon polarization
c     
c     In each basis, name of variable is Resultλiλf, i.e. "IN" first, but the array is labelled (twoMzp,twoMz), i.e. "OUT" first!!! 
c---------------------------------------------------------------------------------------
c     These are inputs
      complex*16,intent(in)  :: Resultxx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl),Resultxy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16,intent(in)  :: Resultyx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl),Resultyy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      integer,intent(in)     :: verbosity
      integer,intent(in)     :: twoSnucl
      
      
c     These are outputs.
      complex*16,intent(out) :: Result11(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl),Result1m1(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16,intent(out) :: Resultm11(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl),Resultm1m1(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      
c     Internal variables
      integer                :: twoMz,twoMzp
c----------------------------------------------------------------------
      if (verbosity.eq.1000) continue
      
c     initialise against paranoia
      Result11=c0
      Result1m1=c0
      Resultm11=c0
      Resultm1m1=c0
      
      do twoMzp=twoSnucl,-twoSnucl,-2
         do twoMz=twoSnucl,-twoSnucl,-2
            Result11(twoMzp,twoMz)=0.5*((Resultxx(twoMzp,twoMz)+Resultyy(twoMzp,twoMz))+ci*
     &           (-1.d0*(Resultxy(twoMzp,twoMz))+Resultyx(twoMzp,twoMz)))
            Result1m1(twoMzp,twoMz)=0.5*((-1.d0*(Resultxx(twoMzp,twoMz))+Resultyy(twoMzp,twoMz))-ci*
     &           (Resultxy(twoMzp,twoMz)+Resultyx(twoMzp,twoMz)))
            Resultm11(twoMzp,twoMz)=0.5*((-1.d0*(Resultxx(twoMzp,twoMz))+Resultyy(twoMzp,twoMz))+ci*
     &           (Resultxy(twoMzp,twoMz)+Resultyx(twoMzp,twoMz)))
            Resultm1m1(twoMzp,twoMz)=0.5*((Resultxx(twoMzp,twoMz)+Resultyy(twoMzp,twoMz))+ci*
     &           (Resultxy(twoMzp,twoMz)-1.d0*(Resultyx(twoMzp,twoMz))))
         end do
      end do
      return
      end
      





 
 

