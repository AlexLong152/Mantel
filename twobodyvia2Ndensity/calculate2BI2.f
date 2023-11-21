c     hgrie Oct 2022: v2.0 fewbody-Compton
c     new Aug 2020, based on 3He density codes with the following datings/changes:
      subroutine Calculate2BIntegralI2(Int2Bxx,Int2Bxy,Int2Byx,Int2Byy,
     &     Int2Bx,Int2By,Int2Bpx,Int2Bpy,
     &     j12p,m12p,l12p,s12p,t12p,mt12p,j12,m12,
     &     l12,s12,t12,mt12,p12,p12p,th12,phi12,Nth12,Nphi12,
     &     thetacm,k,
     &     AngularType12,angweight12,calctype,verbosity)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hgrie May 2018: used to be part of 3HeCompt/twobody/
c     now part of twobodyvia2Ndensity/, backward compatibility deliberately broken
c     change: now "USE clebsch" refers to routine in common-densities/2Ndensity-modules/
c     that's strictly speaking only a change in Makefile and not in this file,
c     but it does break bitwise compatibility of output files (relative change is <10^-5)
c      
c     twoSmax/twoMz dependence: none, only on quantum numbers of (12) subsystem
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     hgrie 20 June 2014: modified for use of LebedevLaikov or Gaussian
c     integration for theta & phi separately,
c     for solid angle integral in (12) system 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      USE clebsch
      implicit none
c     
      include '../common-densities/constants.def'
      include '../common-densities/params.def'
      include '../common-densities/calctype.def'
c     
      integer,intent(in) :: m12p,m12,j12p,s12p,l12p,j12,s12,l12,Nth12,Nphi12
      integer,intent(in) :: t12p,t12,mt12p,mt12
      
c     hgrie 20 June 2014: added for LebedevLaikov
      integer,intent(in) :: AngularType12
      real*8,intent(in) :: angweight12(Nangmax,Nangmax)
c     index limits of iphi, depending on AngularType12:
      integer imin,imax,jmin,jmax
      
      integer,intent(in) :: calctype,verbosity
c     end add hgrie
      
      real*8,intent(in) :: thetacm,k,th12(Nangmax),phi12(Nangmax)
      
      complex*16,intent(out) :: Int2Bxx,Int2Bxy,Int2Byx,Int2Byy
      complex*16,intent(out) :: Int2Bx,Int2By,Int2Bpx,Int2Bpy
      
      complex*16 Compton2Bxx(0:1,-1:1,0:1,-1:1)
      complex*16 Compton2Bxy(0:1,-1:1,0:1,-1:1)
      complex*16 Compton2Byx(0:1,-1:1,0:1,-1:1)
      complex*16 Compton2Byy(0:1,-1:1,0:1,-1:1)
      complex*16 Compton2Bx(0:1,-1:1,0:1,-1:1)
      complex*16 Compton2By(0:1,-1:1,0:1,-1:1)
      complex*16 Compton2Bpx(0:1,-1:1,0:1,-1:1)
      complex*16 Compton2Bpy(0:1,-1:1,0:1,-1:1)
c     
      integer ith,iphi,jth,jphi,msp,ms,ml12p,ml12
      complex*16 Yl12(-5:5),Yl12p(-5:5)
      complex*16 Intxx(-5:5,-5:5),Intxy(-5:5,-5:5)
      complex*16     Intyx(-5:5,-5:5),Intyy(-5:5,-5:5)
      complex*16 Intx(-5:5,-5:5),Inty(-5:5,-5:5)
      complex*16     Intpx(-5:5,-5:5),Intpy(-5:5,-5:5)
      complex*16 Yl12pstar
      real*8 cgcp,cgc,p12x,p12y,p12z,p12px,p12py,p12pz,p12,p12p
c     
      if (verbosity.eq.1000) continue ! unused variable, kept for future use
c     
      if ((l12p .gt. 5) .or. (l12 .gt. 5)) then
         goto 100
      endif 
      Int2Bxx=c0
      Int2Bxy=c0
      Int2Byx=c0
      Int2Byy=c0
      Compton2Bxx=c0
      Compton2Bxy=c0
      Compton2Byx=c0
      Compton2Byy=c0
      Int2Bx=c0
      Int2By=c0
      Int2Bpx=c0
      Int2Bpy=c0
      Compton2Bx=c0
      Compton2By=c0
      Compton2Bpx=c0
      Compton2Bpy=c0
c     
      call initclebsch                !Initializing the factorial array
c     Loop  to sum over ms12 and ms12p (called ms and msp here). 
c     The value of ms and msp together with m12 & m12p determine ml12 and ml12p. 
      do msp=-s12p,s12p,1
         ml12p=m12p-msp
         do ms=-s12,s12,1
            ml12=m12-ms
c     Initializing to zero
            cgc=0.d0
            cgcp=0.d0
            Yl12=c0
            Yl12p=c0
            Yl12pstar=c0
            Intxx=c0
            Intxy=c0
            Intyx=c0
            Intyy=c0
            Intx=c0
            Inty=c0
            Intpx=c0
            Intpy=c0
            if ((abs(ml12p) .le. l12p) .and. (abs(ml12) .le. l12)) then
               do ith=1,Nth12
c     c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c     hgrie 20 June 2014: pick theta& phi summation parameters following AngularType12
c     for LebedevLaikov, only sum over diagonal elements of angweight12 (all others are zero)
                  if (AngularType12.eq.1) then !Gaussian in theta and phi separately
                     imin=1
                     imax=Nphi12
                  else if (AngularType12.eq.2) then !LebedevLaikov
                     imin=ith
                     imax=ith
                  else
                     write(*,*) "*** ERROR: Something went wrong with imin/imax in Calculate2BI2. -- Exiting."
                     stop
                  end if
c     c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                  do iphi=imin,imax
                     call calculatepvector(p12x,p12y,p12z,p12,
     &                    th12(ith),phi12(iphi),verbosity)
                     call getsphericalharmonics(Yl12,l12,th12(ith),phi12(iphi))
                     do jth=1,Nth12
c     c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c     hgrie 20 June 2014: pick theta& phi summation parameters following AngularType12
c     for LebedevLaikov, only sum over diagonal elements of angweight12 (all others are zero)
                        if (AngularType12.eq.1) then !Gaussian in theta and phi separately
                           jmin=1
                           jmax=Nphi12
                        else if (AngularType12.eq.2) then !LebedevLaikov
                           jmin=jth
                           jmax=jth
                        else
                           write(*,*) "*** ERROR: Something went wrong with imin/imax in Calculate2BI2. -- Exiting."
                           stop
                        end if
                        do jphi=jmin,jmax
                           call calculatepvector(p12px,p12py,p12pz,p12p,
     &                          th12(jth),phi12(jphi),verbosity)
                           call getsphericalharmonics(Yl12p,l12p,th12(jth),phi12(jphi))
                           Yl12pstar=Real(Yl12p(ml12p))-ci*Imag(Yl12p(ml12p))
                           call Calc2Bspinisospintrans(Compton2Bxx,Compton2Bxy,
     &                          Compton2Byx,Compton2Byy,
     &                          Compton2Bx,Compton2By,
     &                          Compton2Bpx,Compton2Bpy,
     &                          t12,mt12,t12p,mt12p,l12,
     &                          s12,l12p,s12p,thetacm,k,p12x,p12y,p12z,
     &                          p12px,p12py,p12pz,calctype,verbosity)
                           Intxx(ml12p,ml12)=Intxx(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
     &                          Compton2Bxx(s12p,msp,s12,ms)
                           Intxy(ml12p,ml12)=Intxy(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
     &                          Compton2Bxy(s12p,msp,s12,ms)
                           Intyx(ml12p,ml12)=Intyx(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
     &                          Compton2Byx(s12p,msp,s12,ms)
                           Intyy(ml12p,ml12)=Intyy(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
     &                          Compton2Byy(s12p,msp,s12,ms)
                           Intx(ml12p,ml12)=Intx(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
     &                          Compton2Bx(s12p,msp,s12,ms)
                           Inty(ml12p,ml12)=Inty(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
     &                          Compton2By(s12p,msp,s12,ms)
                           Intpx(ml12p,ml12)=Intpx(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
     &                          Compton2Bpx(s12p,msp,s12,ms)
                           Intpy(ml12p,ml12)=Intpy(ml12p,ml12)+Yl12pstar*Yl12(ml12)*
     &                          angweight12(ith,iphi)*angweight12(jth,jphi)*
     &                          Compton2Bpy(s12p,msp,s12,ms)
                        end do  ! jphi
                     end do     ! jth
                  end do        ! iphi
               end do           ! ith
            end if
            cgcp=CG(2*l12p,2*s12p,2*j12p,2*ml12p,2*msp,2*m12p)
            cgc=CG(2*l12,2*s12,2*j12,2*ml12,2*ms,2*m12)
            Int2Bxx=Int2Bxx+Intxx(ml12p,ml12)*cgc*cgcp
            Int2Bxy=Int2Bxy+Intxy(ml12p,ml12)*cgc*cgcp
            Int2Byx=Int2Byx+Intyx(ml12p,ml12)*cgc*cgcp
            Int2Byy=Int2Byy+Intyy(ml12p,ml12)*cgc*cgcp
            Int2Bx=Int2Bx+Intx(ml12p,ml12)*cgc*cgcp
            Int2By=Int2By+Inty(ml12p,ml12)*cgc*cgcp
            Int2Bpx=Int2Bpx+Intpx(ml12p,ml12)*cgc*cgcp
            Int2Bpy=Int2Bpy+Intpy(ml12p,ml12)*cgc*cgcp
         end do                 !ms12
      end do                    !ms12p
      if (verbosity.eq.1000) continue
 100  return
      end
