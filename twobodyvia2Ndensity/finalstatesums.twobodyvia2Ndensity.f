c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
c     hgrie May 2017, revised May 2018 (see below):
c     2N density integration (replaces I3 integration of Deepshikha's code)
c     based on Deepshikha's twobody/finalstatesums.twobody.f
c
c     twoMz, twoMzplimit dependence: arrays and do-loops
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
     &     Resultxx,Resultxy,Resultyx,Resultyy,
     &     Anucl,twoSnucl,twoMzplimit,j12,m12,l12,s12,t12,mt12,
     &     k,thetacm,
     &     ip12,p12,wp12,
     &     P12MAG,AP12MAG,NP12,
     &     th12,phi12,Nth12,Nphi12,j12max,
     &     AngularType12,angweight12,calctype,verbosity)
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
      integer,intent(in) :: Anucl,twoSnucl,twoMzplimit
c
      complex*16,intent(out) :: Resultxx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl),Resultxy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
      complex*16,intent(out) :: Resultyx(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl),Resultyy(-twoSnucl:twoSnucl,-twoSnucl:twoSnucl)
c      
      integer alpha2N,alpha2Np,rindx
      
      integer mt12p,j12p,s12p,l12p,t12p,m12p ! quantum #s of (12) system -- integer already, so no factor 2
      integer ip12,ip12p
      complex*16 Int2Bxx,Int2Bxy,Int2Byx,Int2Byy,Int3
      complex*16 Int2Bx,Int2By,Int2Bpx,Int2Bpy
      complex*16 Int3x, Int3y, Int3px, Int3py
      complex*16 fact, factx, facty, factpx,factpy,f 
      integer twoMz,twoMzp

      integer,intent(in) :: AngularType12
      real*8,intent(in)  :: angweight12(Nangmax,Nangmax)
      integer,intent(in) :: calctype,verbosity

      integer twoMzlimit  ! for symmetry calculation: Mzp>=0 *and* for Mzp=0, only Mz>=0, else Mz between +Snucl and -Snucl

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
                     call Calculate2BIntegralI2(Int2Bxx,Int2Bxy,Int2Byx,Int2Byy,
     &                    Int2Bx,Int2By,Int2Bpx,Int2Bpy,
     &                    j12p,m12p,l12p,s12p,t12p,mt12p,
     &                    j12,m12,l12,s12,t12,mt12,
     &                    p12*HC,P12MAG(ip12p)*HC,th12,
     &                    phi12,Nth12,Nphi12,thetacm,k,
     &                    AngularType12,angweight12,calctype,verbosity)

c     let twoMzprime only run over half of MEs when symmetry is used
                     do twoMzp=twoSnucl,twoMzplimit,-2
c         for Mzp=0, run only over Mz>=0 -- that's still 2 more than necessary, but good enough -- will be cured below
                        if ((twoMzp.eq.0).and.(twoMzplimit.eq.0)) then
                           twoMzlimit = 0
                        else
                           twoMzlimit = -twoSnucl
                        end if   
                        do twoMz=twoSnucl,twoMzlimit,-2
c     now call density
                           
c                           write(*,*) l12,s12,j12,mt12,m12,twoMz
c                           write(*,*) l12p,s12p,j12p,mt12p,m12p,twoMzp
                           alpha2N = get2Nchannum(l12,s12,j12,mt12,m12,twoMz)
                           alpha2Np = get2Nchannum(l12p,s12p,j12p,mt12p,m12p,twoMzp)
                           rindx=rhoindx(alpha2N,alpha2Np)
c                           write(*,*) "rindx = ",rindx," ; alpha2N  = ",alpha2N," ; alpha2Np = ",alpha2Np
c                           write(*,*) "             ρ12 = ",rho(ip12,ip12,rindx)
                           Int3=rho(ip12,ip12p,rindx)
c                           write(*,'(3(A,I5),A,E15.8,SP,E16.8," I")')
c     &                          "    ρ(rindx=",rindx,",ip12=",ip12,",ip12p=",ip12p,") = ",Int3
                           
                           Int3x=0.d0 !quick fix to get NOthing at OQ3 -- needs changing at OQ4!
                           Int3y=0.d0 !quick fix to get NOthing at OQ3 -- needs changing at OQ4!
                           Int3px=0.d0 !quick fix to get NOthing at OQ3 -- needs changing at OQ4!
                           Int3py=0.d0 !quick fix to get NOthing at OQ3 -- needs changing at OQ4!
c     multiplication by HC**3.d0 transforms ME units from fm^-3 to MeV^3
c     hgrie May 2018/July 2020:
c     original factor 3.d0 in 3He-code is replaced by number of nucleon pairs inside nucleus -- see 3He-densities paper
c     hgrie Oct 2022: include here factor 1/(2π)³ so that final amplitudes for onebody and twobody have SAME sizes.
                           f=Anucl*(Anucl-1)/2*p12**2*wp12*P12MAG(ip12p)**2*AP12MAG(ip12p)*HC**3.d0/(2*Pi)**3
c     next line replaces Int3 by 2Ndensity ρ                                    
                           fact=f*Int3
c     following only for OQ4 -- NOT YET IMPLEMENTED !!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           factx=f*Int3x
                           facty=f*Int3y                                    
                           factpx=f*Int3px
                           factpy=f*Int3py
                           
c                           write(*,*) fact, factx, facty, factpx, factpy
c                           write(*,*) Int2Bxx, Int2Bxy, Int2Byx, Int2Byy 
c                           write(*,*) Int2Bpx, Int2Bpy, Int2Bx, Int2By
                           
                           Resultxx(twoMzp,twoMz)=Resultxx(twoMzp,twoMz)+fact*Int2Bxx+factx*Int2Bpx+factpx*Int2Bx
                           Resultyy(twoMzp,twoMz)=Resultyy(twoMzp,twoMz)+fact*Int2Byy+facty*Int2Bpy+factpy*Int2By
c hgrie Aug 2020: now cure: for Mzp=Mz=0, only calculate xx and yy, since xy and yx must be zero               
                           if ((twoMzplimit.eq.0).and.(twoMzp.eq.0).and.(twoMz.eq.0)) then
                              continue
                           else
                              Resultxy(twoMzp,twoMz)=Resultxy(twoMzp,twoMz)+fact*Int2Bxy+factx*Int2Bpy+factpy*Int2Bx
                              Resultyx(twoMzp,twoMz)=Resultyx(twoMzp,twoMz)+fact*Int2Byx+facty*Int2Bpx+factpx*Int2By
                           end if
c      end cure
                        end do  !twoMz
                     end do     !twoMzp
                  end do        !ip12p
               end do           !m12p
            end do              !l12p
         end do                 !s12p
      end do                    !j12p
      if (verbosity.eq.1000) continue
      return
      end
