c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2017, revised May 2018 (see below):
c     2N density integration (replaces I3 integration of Deepshikha's code)
c     based on Deepshikha's twobody/finalstatesums.twobody.f
c
c     twoMz dependence: arrays and do-loops
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c     TO do:
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CHANGES:
c      
c     hgrie Oct 2022: *HUGE CHANGE* inside twobodyfinalstatesumsvia2Ndensity():
c           Defined twobody ME to INCLUDE the factor 1/(2π)³ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c           so that final amplitudes for onebody and twobody have SAME sizes.
c      
c     hgrie June 2018: renamed "parity" to "symmetry -- see notes in usesymmetry+*.f
c       
c     hgrie May 2018: more extensive description of changes in main.*.f
c                     rewritten such that magnetic quantum numbers are "2*Mz" etc
c
      subroutine twobodyfinalstatesumsvia2Ndensity(
     &     Result,
     &     Anucl,twoSnucl,extQnumlimit,j12,m12,l12,s12,t12,mt12,
     &     k,thetacm,
     &     ip12,p12,wp12,
     &     P12MAG,AP12MAG,NP12,
     &     th12,phi12,Nth12,Nphi12,j12max,
     &     AngularType12,angweight12,calctype,symmetry,verbosity)
c     
c     Note structure of loops is almost the same as in one-body case, only real difference is in computation of I2, and
c     fact that p_{12}' integral now runs over full range [0,infty). Partly for this reason, the order of the p_{12}'
c     and p_3' loops has been interchanged, c.f. the one-body version of this routine.
c     
c     hgrie 20 June 2014: modified for use of LebedevLaikov or Gaussian
c     integration for theta & phi separately, for solid angle integral in (12) system 

      USE CompDens ! needs module CompDens.mod
c     
      implicit NONE
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
      integer,intent(in) :: j12,m12,l12,s12,t12,mt12 ! are automatically integers, so do NOT multiply by 2, unlike for Mz=>twoMz
      integer,intent(in) :: j12max
      real*8,intent(in)  :: k,thetacm
      real*8,intent(in)  :: p12,wp12
      real*8,intent(in)  :: P12MAG(Npmax),AP12MAG(Npmax)
      integer,intent(in) :: NP12
      real*8,intent(in)  :: th12(Nangmax),phi12(Nangmax)
      integer,intent(in) :: Nth12,Nphi12
      integer,intent(in) :: Anucl,twoSnucl
      integer,intent(in) :: extQnumlimit
      
      integer extQnum      ! counter of combined external quantum numbers of in and out state
c
      complex*16,intent(out) :: Result(1:extQnumlimit,-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
c      
      integer alpha2N,alpha2Np,rindx
      
      integer mt12p,j12p,s12p,l12p,t12p,m12p ! quantum #s of (12) system -- integer already, so no factor 2
      integer ip12,ip12p
      complex*16 Int2B(1:extQnumlimit)
      complex*16 fact
      
c      complex*16 Int2Bx,Int2By,Int2Bpx,Int2Bpy ! for STUMP, see below
c      complex*16 Int3x, Int3y, Int3px, Int3py  ! for STUMP, see below
c      complex*16 factx, facty, factpx, factpy  ! for STUMP, see below
      
      integer twoMz,twoMzp

      integer,intent(in) :: AngularType12
      real*8,intent(in)  :: angweight12(Nangmax,Nangmax)
      integer,intent(in) :: calctype,symmetry,verbosity

c      logical invokesymmetrytwoMzp  ! function to invoke symmetry/-ies of amplitudes, dependent on process
c      logical invokesymmetrytwoMz   ! function to invoke symmetry/-ies of amplitudes, dependent on process
c      logical invokesymmetryextQnum ! function to invoke symmetry/-ies of amplitudes, dependent on process

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mt12p=mt12                              ! isospin projection in (12), fixes charge of spectator(s)    
      do j12p=0,j12max                        ! total ang mom (12); usually j12max=1 for 1% convergence
         do s12p=0,1                          ! spin (12) subsystem
            do l12p=abs(j12p-s12p),j12p+s12p  ! orbital angular mom. (12)
               t12p=(1-(-1)**(l12p+s12p+1))/2 ! isospin (12) subsystem
               do m12p=-j12p,j12p             ! total ang mom projection of out-(12)
c                  
c     Angular-momentum sums are implemented exactly as in one-body version of this routine
c     
                  do ip12p=1,NP12 ! mag of momentum (12) subsystem
c
                     call Calculate2BIntegralI2(Int2B,
c     &                    Int2Bx,Int2By,Int2Bpx,Int2Bpy, ! for STUMP, see below
     &                    extQnumlimit,
     &                    j12p,m12p,l12p,s12p,t12p,mt12p,
     &                    j12,m12,l12,s12,t12,mt12,
     &                    p12*HC,P12MAG(ip12p)*HC,th12,
     &                    phi12,Nth12,Nphi12,thetacm,k,
     &                    AngularType12,angweight12,calctype,verbosity)

                     do twoMzp=twoSnucl,-twoSnucl,-2
c     check if ME to be skipped and reconstructed via a symmetry of the amplitudes
c                        if (invokesymmetrytwoMzp(symmetry,twoSnucl,twoMzp,verbosity)) exit
                        
                        do twoMz=twoSnucl,-twoSnucl,-2
c     check if ME to be skipped and reconstructed via a symmetry of the amplitudes
c                           if (invokesymmetrytwoMz(symmetry,twoSnucl,twoMzp,twoMz,verbosity)) exit
c                           
c     now call density
c                           write(*,*) l12,s12,j12,mt12,m12,twoMz
c                           write(*,*) l12p,s12p,j12p,mt12p,m12p,twoMzp
                           alpha2N = get2Nchannum(l12,s12,j12,mt12,m12,twoMz)
                           alpha2Np = get2Nchannum(l12p,s12p,j12p,mt12p,m12p,twoMzp)
                           rindx=rhoindx(alpha2N,alpha2Np)
c                           write(*,*) "rindx = ",rindx," ; alpha2N  = ",alpha2N," ; alpha2Np = ",alpha2Np
c                           write(*,*) "             ρ12 = ",rho(ip12,ip12,rindx)
c                           write(*,'(3(A,I5),A,E15.8,SP,E16.8," I")')
c     &                          "    ρ(rindx=",rindx,",ip12=",ip12,",ip12p=",ip12p,") = ",rho(ip12,ip12p,rindx)
                           
c     Next line: integration volume of p12 and p12p magnitudes and 2Ndensity ρ
c     hgrie May 2018/July 2020:
c     original factor 3.d0 in 3He-code is replaced by number of nucleon pairs inside nucleus -- see 3He-densities paper
c     hgrie Oct 2022: include here factor 1/(2π)³ so that final amplitudes for onebody and twobody have SAME sizes.
c     multiplication by HC**3.d0 transforms ME units from fm^-3 to MeV^3
                           fact=Anucl*(Anucl-1)/2*p12**2*wp12*P12MAG(ip12p)**2*AP12MAG(ip12p)*
     &                          rho(ip12,ip12p,rindx)*HC**3.d0/(2*Pi)**3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                           
c     hgrie Nov 2023: Following is a STUMP from the Compton code, used there only for OQ4 -- NOT YET IMPLEMENTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     I leave this here because maybe some of this can be recycled later for boost corrections or so?
c                           Int3x=0.d0  !quick fix to get NOthing at OQ3 -- needs changing at OQ4!
c                           Int3y=0.d0  !quick fix to get NOthing at OQ3 -- needs changing at OQ4!
c                           Int3px=0.d0 !quick fix to get NOthing at OQ3 -- needs changing at OQ4!
c                           Int3py=0.d0 !quick fix to get NOthing at OQ3 -- needs changing at OQ4!
                           
c                           factx=fact*Int3x
c                           facty=fact*Int3y                                    
c                           factpx=fact*Int3px
c                           factpy=fact*Int3py
c                           
c                           write(*,*) fact, factx, facty, factpx, factpy
c                           write(*,*) Int2Bxx, Int2Bxy, Int2Byx, Int2Byy 
c                           write(*,*) Int2Bpx, Int2Bpy, Int2Bx, Int2By
c                           Resultxx(twoMzp,twoMz) = Resultxx(twoMzp,twoMz) + fact*Int2B(1) + factx*Int2Bpx + factpx*Int2Bx ! associated with xx in Compton
c                           Resultxy(twoMzp,twoMz) = Resultxy(twoMzp,twoMz) + fact*Int2B(2) + factx*Int2Bpy + factpy*Int2Bx ! associated with xy in Compton
c                           Resultyx(twoMzp,twoMz) = Resultyx(twoMzp,twoMz) + fact*Int2B(3) + facty*Int2Bpx + factpx*Int2By ! associated with yx in Compton
c                           Resultyy(twoMzp,twoMz) = Resultyy(twoMzp,twoMz) + fact*Int2B(4) + facty*Int2Bpy + factpy*Int2By ! associated with yy in Compton
c     END OF STUMP    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
                           do extQnum=1,extQnumlimit
c     check if ME to be skipped and reconstructed via a symmetry of the amplitudes
c                              if (invokesymmetryextQnum(symmetry,extQnumlimit,extQnum,twoSnucl,twoMzp,twoMz,verbosity)) exit
                              Result(extQnum,twoMzp,twoMz) = Result(extQnum,twoMzp,twoMz) + fact*Int2B(extQnum) 
                           end do !extQnum
                        end do    !twoMz
                     end do       !twoMzp
                  end do          !ip12p
               end do             !m12p
            end do                !l12p
         end do                   !s12p
      end do                      !j12p
      
      if (verbosity.eq.1000) continue
      return
      end
