c     To contain symmetry considerations and corresponding output messages
c     Can use variable "symmetry" to impose several constraints
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function invokesymmetrytwoMzp(symmetry,twoSnucl,twoMzp,verbosity)

      implicit none

      integer,intent(in) :: symmetry          ! whether and which specific package Usesymmetry() should be used for process -- NOT YET IMPLEMENTED!!!
      integer,intent(in) :: twoSnucl,twoMzp
      
      integer,intent(in) :: verbosity

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set default: no symmetry
      invokesymmetrytwoMzp=.False.

      if ((symmetry.ge.1).and.(twoMzp.le.0.)) invokesymmetrytwoMzp = .True.

      if (verbosity.eq.1000) continue
      end function invokesymmetrytwoMzp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function invokesymmetrytwoMz(symmetry,twoSnucl,twoMzp,twoMz,verbosity)

      implicit none

      integer,intent(in) :: symmetry          ! whether and which specific package Usesymmetry() should be used for process -- NOT YET IMPLEMENTED!!!
      integer,intent(in) :: twoSnucl,twoMzp,twoMz
      
      integer,intent(in) :: verbosity

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set default: no symmetry
      invokesymmetrytwoMz=.False.

      if (verbosity.eq.1000) continue
      end function invokesymmetrytwoMz
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical function invokesymmetryextQnum(symmetry,extQnumlimit,extQnum,twoSnucl,twoMzp,twoMz,verbosity)

      implicit none

      integer,intent(in) :: symmetry          ! whether and which specific package Usesymmetry() should be used for process -- NOT YET IMPLEMENTED!!!
      integer,intent(in) :: extQnumlimit,extQnum,twoSnucl,twoMzp,twoMz
      
      integer,intent(in) :: verbosity

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set default: no symmetry
      invokesymmetryextQnum=.False.

      if ((symmetry.eq.2).and.(extQnum.eq.2)) invokesymmetryextQnum = .True.

      if (verbosity.eq.1000) continue
      end function invokesymmetryextQnum
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
