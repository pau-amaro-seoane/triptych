c       Copyright (c) 2001, 2002
c       by James Lombardi, Vassar College, Poughkeepsie, NY
c       and Jessica Sawyer Warren, Rutgers University, NJ.
c    This file is part of Make Me A Star.
c
c    Make Me A Star is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 2 of the License, or
c    (at your option) any later version.
c
c    Make Me A Star is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with Make Me A Star; if not, write to the Free Software
c    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
*************************************************************************
* Make Me A Star:
* This file contains the main subroutine makemeastar
* version 1.4
* August 9, 2002
* James Lombardi (lombardi@vassar.edu) and Jessica Sawyer Warren
*************************************************************************
      SUBROUTINE makemeastar(peri,
     $     mProfile1,rProfile1,PProfile1,rhoProfile1,chemicalProfiles1,
     $     mProfile2,rProfile2,PProfile2,rhoProfile2,chemicalProfiles2,
     $     numShells1,numShells2,
     &     verboseness,numberOfChemicalElements,
     $     NumOutLines,Mremnant,Premnant,
     &     Dremnant,Rremnant,ELremnant,jremnant)
**********************************************************************
* Input:
*     periastron = r_p / (R_1+R_2), where r_p is the periastron separation of
*        the (parabolic) orbits, and R_1 and R_2 are the radii of the parent
*        stars when isolated.  Therefore, r_p=0 for a headon collision and
*        increases as the impact parameter increases.
*     verboseness = an integer that sets how talkative the routines should be
*        verboseness=0 : no output
*        verboseness=1 : limited output
*        verboseness=2 : detailed output
*     el = an integer stating how many chemical abundances will be treated
*        (information on the parent abundance profiles must of course be
*         provided in parentfile(1) and parentfile(2))
*     NumOutLines = desired number of rows in the outputted remnant profiles
* Output:
*     Mremnant = an array of length NumOutLines which gives enclosed mass
*     Premnant = an array of length NumOutLines for pressure profile
*     Dremnant = an array of length NumOutLines for density profile
*     Rremnant = an array of length NumOutLines for radius profile
*     ELremnant = an array of size (NumOutLines x numberOfChemicalElements)
*        for chemical composition profiles
*     jremnant = an array of length NumOutLines for the specific angular
*        momentum profile
**********************************************************************
* This subroutine reads in data containing enclosed mass M,          *
* pressure P, and density D (in                                      *
* cgs units), and chemical abundances of a particular star.  It will *
* calculate the "entropy" A from P and D, and call                   *
* subroutines SORT, PROFILE, ENERGY, LOSEMASS, and SHOCK.            *
*     ~J.E.S. and J.C.L                                              *
**********************************************************************

      IMPLICIT NONE

* params.h contains declarations and values for N, NROWS, pi, G,
* Msun, Rsun, Mto,& Rto.  USE CGS UNITS!!!

      INCLUDE 'params.h'
      integer verbosity,NP,NumOutLines,verboseness
      double precision peri
      INTEGER Jpmax(N),I,K,C,Pmax(N),Jr,Jrmax
      DOUBLE PRECISION mProfile1(NROWS),rProfile1(NROWS),
     $     PProfile1(NROWS),rhoProfile1(NROWS)
      DOUBLE PRECISION mProfile2(NROWS),rProfile2(NROWS),
     $     PProfile2(NROWS),rhoProfile2(NROWS)
      DOUBLE PRECISION chemicalProfiles1(NROWS,elmax),
     $     chemicalProfiles2(NROWS,elmax)
      integer numShells1,numShells2
      DOUBLE PRECISION Mp(NROWS,N),Ap(NROWS,N),Rp(NROWS,N),Pp(NROWS,N),
     &     Dp(NROWS,N),Aptry(N),TotE(N),W,As(NROWS,N)
      DOUBLE PRECISION ELp(NROWS,elmax,N),
     &     ELm(NROWS,elmax,N)
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),ELr(NROWS,elmax),P(NROWS),
     &     D(NROWS),R(NROWS)
      DOUBLE PRECISION Mremnant(NROWS),
     &     ELremnant(elmax,NROWS),
     &     Premnant(NROWS), Dremnant(NROWS),
     &     Rremnant(NROWS),jremnant(NROWS)
      DOUBLE PRECISION INTlow,INThigh,NEWINT,DELNRG,zeroin3,
     &    periastron,TotErem,Jtot,TotJrem,MLsfrac
      COMMON/NRG/TotE,TotErem
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      COMMON/VARIABLES/P,D,R
      COMMON/PERI/periastron
      COMMON/MASSLOSS/MLsfrac
      COMMON/verbiage/verbosity
      EXTERNAL DELNRG
      logical alreadyfudged
      integer numberOfChemicalElements
      integer throwaway,Jrtransition
      common/numel/el
      double precision scaledJtot,Jtotin,unscaledJtotout(NROWS)
      double precision cs,jcs
      double precision jtryin(NROWS),jtryout(NROWS),k2,k3,diff
      character*16 OUTFN
      
* duplicate (or nearly duplicate) the values of some of the input variables so
* that they can be passed through common blocks:
      NP=NumOutLines-1
      verbosity=verboseness
      el=numberOfChemicalElements


      if(mProfile1(numShells1).ge.mProfile2(numShells2)) then
         Jpmax(1)=numShells1
         DO I=1,Jpmax(1)
            Mp(I,1)=mProfile1(I)
            Rp(I,1)=rProfile1(I)
            Pp(I,1)=PProfile1(I)
            Dp(I,1)=rhoProfile1(I)
            do C=1,el
               Elp(I,C,1)=chemicalProfiles1(I,C)
            enddo
         enddo
         Jpmax(2)=numShells2
         DO I=1,Jpmax(2)
            Mp(I,2)=mProfile2(I)
            Rp(I,2)=rProfile2(I)
            Pp(I,2)=PProfile2(I)
            Dp(I,2)=rhoProfile2(I)
            do C=1,el
               Elp(I,C,2)=chemicalProfiles2(I,C)
            enddo
         enddo
      else
C Switch stars if more massive star was second...
         Jpmax(1)=numShells2
         DO I=1,Jpmax(1)
            Mp(I,1)=mProfile2(I)
            Rp(I,1)=rProfile2(I)
            Pp(I,1)=PProfile2(I)
            Dp(I,1)=rhoProfile2(I)
            do C=1,el
               Elp(I,C,1)=chemicalProfiles2(I,C)
            enddo
         enddo
         Jpmax(2)=numShells1
         DO I=1,Jpmax(2)
            Mp(I,2)=mProfile1(I)
            Rp(I,2)=rProfile1(I)
            Pp(I,2)=PProfile1(I)
            Dp(I,2)=rhoProfile1(I)
            do C=1,el
               Elp(I,C,2)=chemicalProfiles1(I,C)
            enddo
         enddo
      endif

      if(verbosity.ge.2) then
         write(6,*) 'Here are the values of the constants being used:'
         write(6,*) 'c_1=',const1
         write(6,*) 'c_2=',const2
         write(6,*) 'c_3=',const3
         write(6,*) 'c_4=',const4
         write(6,*) 'c_5=',const5
         write(6,*) 'c_6=',const6
         write(6,*) 'c_7=',const7
         write(6,*) 'c_8=',const8
         write(6,*) 'c_9=',const9
         write(6,*) 'c_10=',const10
         write(6,*)
      endif

      if(el.GT.elmax) then
         write(6,*)'ERROR: Increase elmax in params.h to at least',el
         stop
      endif
c      if(NumOutLines+1.GT.NROWS) then
c         write(6,*)'ERROR: Increase NROWS in params.h to at least',
c     $        NumOutLines+1
c         write(6,*)'You may also have to change NROWS, or its'
c         write(6,*)'equivalent, in your driver program.'
c         write(6,*)'Current value of NROWS=',NROWS
c         stop
c      endif

      if(verbosity.ge.1) write(6,'(1X,a,g12.4)')
     &     'Entering MAKEMEASTAR w/ periastron parameter r_p/(R_1+R_2)='
     &     ,peri

* Now we have to calculate A.

      DO K=1,N
         if(verbosity.ge.2) WRITE(6,*)

         alreadyfudged=.false.
         DO I=1,Jpmax(K)
            Aptry(K)=Pp(I,K)/Dp(I,K)**(5.d0/3.d0)

* We want to fake data for which the "entropy" A does not increase.  Therefore
* the else part of the IF statement fakes bad data.  In particular,
* the "entropy" is forced to increase ever so slightly.
            IF(Aptry(K).GT.Ap(I-1,K) .OR. I.EQ.1) THEN
               Ap(I,K)=Aptry(K)
            else
               Ap(I,K)=1.00001d0*Ap(I-1,K)
C               Ap(I,K)=1.001d0*Ap(I-1,K)
               if(.not. alreadyfudged) then
                  if(verbosity.ge.2) then
                     write(6,*)
     &                 'Fudging A profile outside m/M_solar=',
     &                 Mp(I,K)/Msun,' in star',K
                     write(6,*)
                  endif
                  alreadyfudged=.true.
               endif
            endif
         ENDDO

         if(verbosity.ge.1) then
            write(6,'(1X,a,i2,a,g10.2,a,g11.3,a)')
     $           'Parent star',K,': Mass=',Mp(Jpmax(K),K)/Msun,
     $           'Radius=',Rp(Jpmax(K),K)/Rsun,'(solar units)' 
            if(verbosity.ge.2) WRITE(6,*) 'Jpmax(',K,')=',Jpmax(K)
         endif

         CALL ENERGY(Mp(1,K),Pp(1,K),Dp(1,K),Rp(1,K),
     &        Jpmax(K),K)
         
      ENDDO

      if(Mp(Jpmax(2),2).gt.Mp(Jpmax(1),1))
     $     pause 'The more massive star is supposed to be first...'

c     First find mass loss for the headon collision case between these two
c     parent stars:
      periastron=0.d0
      CALL LOSEMASS

c Determine the intercept b_1 for the shock heating for the headon case
      INTlow=-0.9d0
      INThigh=0.9d0
c      call zbrak(DELNRG,INTlow,INThigh,3,xb1,xb2,nb)
C      write(6,*) 'Back from zbrak: number of bracketing pairs was',nb
c      if(nb.ne.1 .and. verbosity.ge.2) then
c         write(6,*)
c         write(6,*) 'The number of bracketing pairs was',nb
c         if(nb.eq.0) stop
c         write(6,*) 'We are going to find the model with the most'
c         write(6,*) 'shock heating.  (The models with less shock'
c         write(6,*) 'heating most likely were not having their'
c         write(6,*) 'energies calculated accurately, which you'
c         write(6,*) 'can confirm by checking above to see if the'
c         write(6,*) 'virial of the remnant was ever abnormally'
c         write(6,*) 'large.)'
c      endif
c
c      NEWINT=ZBRENT3(DELNRG,xb1(nb),xb2(nb),1.d-4)

      NEWINT=zeroin3(INTlow,INThigh,DELNRG,1.d-4)

      if(verbosity.ge.2) write(6,*)
     &     'The headon collision gave an intercept b_1=',NEWINT

      if(peri.ne.0.d0) then
         periastron=peri
         if(verbosity.ge.2) then
            write(6,*)
            write(6,*)'WE CAN NOW TREAT THE ACTUAL PERIASTRON CASE...'
         endif
         CALL LOSEMASS
         NEWINT=NEWINT-const4*periastron*
     &        LOG10(Mp(Jpmax(1),1)/Mp(Jpmax(2),2))
         if(verbosity.ge.2) then
            write(6,*)'For this periastron separation, b_1=',NEWINT
         endif
      endif

      if(verbosity.ge.2) WRITE(6,*)
      if(verbosity.ge.2) WRITE(6,*) 'Starting final SHOCK'
      CALL SHOCK(NEWINT)
      if(verbosity.ge.2) then
         DO K=1,N
            write(6,*)'central shocked A for star',K,' =',As(1,K)
            write(6,*)'surface shocked A for star',K,' =',As(Jpmax(K),K)
            write(6,*)'shocked A for star',K,' at last bound mass=',
     &           As(Pmax(K),K)
            write(6,*)'mass lost from star',K,'=',Mp(Jpmax(K),K)
     $           -Mp(Pmax(K),K)
         ENDDO
         write(6,*)'fraction of ejecta from star 2=',
     $        (Mp(Jpmax(2),2)-Mp(Pmax(2),2))/
     $        (Mp(Jpmax(1),1)-Mp(Pmax(1),1)
     $        +Mp(Jpmax(2),2)-Mp(Pmax(2),2))
      endif

      if(el.eq.0 .and. .not. studymassloss) then
         if(verbosity.ge.2) WRITE(6,*) 'SKIPPING THE MIXING'
      else
         DO K=1,N
            if(verbosity.ge.2) then
               WRITE(6,*)
               WRITE(6,*) 'Calculating mixing for star',K
            endif
            CALL MIX(K)
         ENDDO
      endif


      if(makeshockfiles) then
         do K=1,2
            WRITE(OUTFN,102) K
 102        FORMAT('shock',I1,'.sbe')
            if(verbosity.ge.2) write(6,*)'WRITING SHOCK FILE...',OUTFN
            open(19,file=OUTFN)
            Do I=1,Jpmax(K)
               write(19,*) Mp(I,K),dlog10(As(I,K)-Ap(I,K)),
     $              dlog10(As(I,K)/Ap(I,K))
            enddo
            close(19)
            if(verbosity.ge.2) write(6,*)'...DONE'
         enddo
      endif
      
* Pmax(K) comes from the LOSEMASS subroutine.

      if(verbosity.ge.2) then
         WRITE(6,*)
         WRITE(6,*) 'Starting final SORT'
      endif
      CALL SORT
      if(verbosity.ge.2) then
         WRITE(6,*) 'surface A for remnant=',Ar(Jrmax)
         WRITE(6,*) 'central A for remnant=',Ar(1)
         WRITE(6,*) 'row max for remnant (Jrmax)=',Jrmax
      endif
      
      if(verbosity.ge.2) then
         WRITE(6,*)
         WRITE(6,*) 'Starting final PROFILE'
      endif
      CALL PROFILE
      
      W=0.d0
      DO K=1,N
         W=TotE(K)+W
      ENDDO

      if(verbosity.ge.2) then      
         WRITE(6,*)
         WRITE(6,*) 'Total energy of parent stars is',W
      endif
      if(verbosity.ge.1) then
         WRITE(6,'(a,f7.3)')' Total mass of remnant (in solar units) is'
     &        ,Mr(Jrmax)/Msun
c         DO K=1,N
c            write(6,*) 'Radius of parent star ',K,' (in solar units)='
c     &           ,Rp(Jpmax(K),K)/Rsun
c            write(6,*) 'Mass of parent star ',K,' (in solar units)='
c     &           ,Mp(Jpmax(K),K)/Msun
c         ENDDO
      endif
      Jtot=Mp(Jpmax(1),1)*Mp(Jpmax(2),2)*
     &     (2*G*periastron*(Rp(Jpmax(1),1)+Rp(Jpmax(2),2))/
     &     (Mp(Jpmax(1),1)+Mp(Jpmax(2),2))
     &     )**0.5d0
      if(verbosity.ge.2) write(6,*)
     &     'Total angular momentum of parent stars is (cgs)',Jtot
      TotJrem=Jtot*(1.d0-const9*MLsfrac)
      if(verbosity.ge.2) then
         write(6,*)'f_L * J_tot=',MLsfrac*Jtot
         write(6,*)'Total angular momentum of remnant is (cgs)',TotJrem
         write(6,*) 'J_remnant/M_remnant = ',TotJrem/Mr(Jrmax)
         write(6,*)'Total system energy E_tot is (cgs)',TotE(1)+TotE(2)
         write(6,*)'Total energy of remnant is (cgs)',TotErem
         write(6,*)'f_L * E_tot=',MLsfrac*(TotE(1)+TotE(2))
      endif

      throwaway=0
      DO Jr=1,Jrmax
         if(P(Jr).lt.0.d0 .and. Jr.lt.Jrmax) then
            throwaway=throwaway+1
         else
            Mremnant(Jr-throwaway)=Mr(Jr)
C            Aremnant(Jr-throwaway)=Ar(Jr)
            Premnant(Jr-throwaway)=dabs(P(Jr))
            Dremnant(Jr-throwaway)=D(Jr)
            Rremnant(Jr-throwaway)=R(Jr)
            do C=1,el
*     Note that the order of the arguments is different in ELremnant than it
*     is in Elr.
               Elremnant(C,Jr-throwaway)=ELr(Jr,C)
            enddo
         endif
      ENDDO

      if(verbosity.ge.2) then
         write(6,*)'Jrmax=',Jrmax
         write(6,*) 'throwing away',throwaway,' lines.'
      endif
      NumOutLines=Jrmax-throwaway

      Do Jr=1,NumOutLines
c jtryout will keep track of j in the outer layers, without yet including
c the contribution due to the k_3 term:
         jtryout(Jr)=(G*Mremnant(Jr)*Rremnant(Jr))**0.5d0*
     $        (Mremnant(Jr)/Mremnant(NumOutLines))**(1.d0/3.d0)

         cs=(5.d0/3.d0*Premnant(Jr)/Dremnant(Jr))**0.5d0
         jcs=cs*Rremnant(Jr)
c jtryin will keep track of j in the inner layers:
         jtryin(Jr)=const10*jcs*
     $        (Mremnant(Jr)/Mremnant(NumOutLines))**(1.d0/3.d0)
      enddo

c unscaledJtotout keeps track of the total angular momentum for the shells
c outside (and including) the shell in question, without the k3 contribution:
      unscaledJtotout(NumOutLines)=jtryout(NumOutLines)*0.5d0*
     $     (Mremnant(NumOutLines)-Mremnant(NumOutLines-1))
      Do Jr=NumOutLines-1,2,-1
         unscaledJtotout(Jr)=unscaledJtotout(Jr+1)+
     $        jtryout(Jr)*0.5d0*(Mremnant(Jr+1)-Mremnant(Jr-1))
      enddo
      unscaledJtotout(1)=unscaledJtotout(2)+
     $     jtryout(1)*0.5d0*(Mremnant(1)+Mremnant(2))
         
c Jtotin keeps track of the total angular momentum for the shells
c inside (and including) the shell in question:
      Jtotin=jtryin(1)*0.5d0*(Mremnant(1)+Mremnant(2))
      do Jr=2,NumOutLines
C get k_2 by matching slope of j...
         k2=(jtryin(Jr)-jtryin(Jr-1))/(jtryout(Jr)-jtryout(Jr-1))

C get k_3 by matching j itself
         k3=0.5d0*(
     $        jtryin(Jr-1)+jtryin(Jr)
     $        -k2*(jtryout(Jr-1)+jtryout(Jr))
     $        )

         scaledJtot=Jtotin+k2*unscaledJtotout(Jr)+
     $        k3*(
     $        Mremnant(NumOutLines)-0.5d0*(Mremnant(Jr-1)+Mremnant(Jr))
     $        )

c         write(6,*) Jr, scaledJtot,TotJrem,
c     $        0.5d0*(Mremnant(Jr)+Mremnant(Jr-1))/Mremnant(NumOutLines)

         if ((scaledJtot-TotJrem)*diff .le. 0.d0 .and. Jr.gt.2) then
            if(dabs(scaledJtot-TotJrem).lt.dabs(diff)) then
               Jrtransition=Jr
            else
               Jrtransition=Jr-1
            endif
            if(verbosity.ge.2) then
               write(6,*) 'Transition occurs near m/M_r=',
     $              0.5d0*
     $              (Mremnant(Jrtransition)+Mremnant(Jrtransition-1))
     $              /Mremnant(NumOutLines)
            endif

C tweak k2 so that get exactly the desired total angular momentum
            k2=(TotJrem-Jtotin-k3*(
     $        Mremnant(NumOutLines)-0.5d0*(Mremnant(Jr-1)+Mremnant(Jr))
     $        ))/unscaledJtotout(Jr)
            goto 234
         endif

C Figure out Jtotin for the next iteration....
         if(Jr.lt.NumOutLines) then
            Jtotin=Jtotin+0.5d0*jtryin(Jr)*
     $           (Mremnant(Jr+1)-Mremnant(Jr-1))
         else
            Jtotin=Jtotin+0.5d0*jtryin(Jr)*
     $           (Mremnant(Jr)-Mremnant(Jr-1))            
         endif

         diff=scaledJtot-TotJrem

      enddo

      if (verbosity.ge.2) then
         write(6,*) 'Using only the form for the outer layer'
      endif
      Jrtransition=1
      k3=0.d0
      k2=TotJrem/unscaledJtotout(1)

 234  continue
      if(verbosity.ge.2) then
         write(6,*)'k2=',k2,' k3=',k3
      endif

      Do Jr=1,Jrtransition-1
         jremnant(Jr)=jtryin(Jr)
      enddo
      Do Jr=Jrtransition,NumOutLines
         jremnant(Jr)=k2*jtryout(Jr)+k3
      enddo
      
      scaledJtot=jremnant(1)*0.5d0*(Mremnant(1)+Mremnant(2))
      do Jr=2, NumOutLines-1
         scaledJtot=scaledJtot+jremnant(Jr)*0.5d0*
     $        (Mremnant(Jr+1)-Mremnant(Jr-1))
      enddo
      scaledJtot=scaledJtot+jremnant(NumOutLines)*0.5d0*
     $     (Mremnant(NumOutLines)-Mremnant(NumOutLines-1))
      if(verbosity.ge.2) then
         write(6,*)'Total J=',scaledJtot,'=',TotJrem
         write(6,*) 'DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      endif
      END
      
*************************************************************
*************************************************************
*************************************************************
*************************************************************
*************************************************************
*************************************************************
*************************************************************


      FUNCTION DELNRG(INTERCEPT)

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER K,verbosity
      DOUBLE PRECISION W,TotE(N),TotErem,INTERCEPT,DELNRG,MLsfrac
      COMMON/NRG/TotE,TotErem
      COMMON/MASSLOSS/MLsfrac
      COMMON/verbiage/verbosity

      if(verbosity.ge.2) then
         WRITE(6,*)
         WRITE(6,*) '***minimizing to find intercept***'
      endif

      CALL SHOCK(INTERCEPT)
      CALL SORT
      CALL PROFILE

      W=0.d0
      DO K=1,N
         W=TotE(K)+W
      ENDDO

      DELNRG=TotErem/W-1.d0-const7*MLsfrac
      if(verbosity.ge.2)
     & write(6,*)'delta energy=',DELNRG,'(we want this to become small)'

      RETURN
      END
      SUBROUTINE ENERGY(M,P,D,R,Jmax,STAR)

*************************************************************
* This routine will calculate the total energy of a star.   *
*************************************************************

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER J,Jmax,STAR,T,verbosity
      DOUBLE PRECISION P(NROWS),D(NROWS),M(NROWS),R(NROWS),
     &     Ei,Eg,delM,TotE(N),TotErem
      COMMON/NRG/TotE,TotErem
      COMMON/verbiage/verbosity

* A star's total energy consists of its gravitational energy, Eg,
* and its internal energy, Ei.  These are calculated for the star's
* innermost point and then numerically integrated outward.

      Ei=3.d0/2.d0*P(1)/D(1)*M(1)
      Eg=-G*(4.d0*D(1)*pi)**2.d0/15.d0*R(1)**5.d0

      T=2
      IF(R(1).EQ.0.d0) THEN
         Ei=Ei+3.d0/2.d0*P(2)/D(2)*M(2)
         Eg=Eg-G*(4.d0*D(2)*pi)**2.d0/15.d0*R(2)**5.d0
         T=3
      ENDIF

      DO J=T,Jmax
         delM=M(J)-M(J-1)
         Eg=Eg-G*delM*(M(J)/R(J)+M(J-1)/R(J-1))*0.5d0
         Ei=Ei+3.d0/4.d0*(P(J)/D(J)+P(J-1)/D(J-1))*delM
      ENDDO

      IF(STAR.GT.N) THEN
         if(verbosity.ge.2)
     &        WRITE(6,*) 'Total energy for remnant =',Eg+Ei
         TotErem=Eg+Ei
      ELSE
         if(verbosity.ge.2)
     &        WRITE(6,*) 'Total energy for star',STAR,' =',Eg+Ei
      ENDIF

* If all goes well, the Virial Theorem should hold and 2Ei+Eg=0.

      if(verbosity.ge.2) WRITE(6,*) '2Ei+Eg=',2.d0*Ei+Eg
      
      TotE(STAR)=Eg+Ei

      if(dabs((2.d0*Ei+Eg)/TotE(STAR)).gt.1d-2) then
         write(6,*)
         write(6,*) 'WARNING ... WARNING ... WARNING'
         write(6,*) 'virial/(total energy)=',(2.d0*Ei+Eg)/TotE(STAR)
         write(6,*) 'which seems to be quite far away from 0.'
         write(6,*) 'This ratio would be zero for a star in'
         write(6,*) 'hydrostatic equilibrium.'
         write(6,*) 'Internal energy=',Ei
         write(6,*) 'Gravitational potential energy=',Eg
         if(STAR.GT.N) then
            write(6,*) 'The remnant seems to be far from equilibrium,'
            write(6,*)'which is OK as long as this is not the last call'
            write(6,*) 'to PROFILE....  If it is the last call,'
            write(6,*) 'you should probably request larger arrays'
            write(6,*) 'for your remnant profiles.  If this does not'
            write(6,*) 'help, consider emailing lombardi@vassar.edu'
         else
            write(6,*) 'Parent star',STAR,' seems to be far from'
            write(6,*) 'equilibrium.  Are you sure you the input'
            write(6,*) 'data file for this star has a reasonably'
            write(6,*) 'large number of rows in it?'
c            pause
         endif
      endif
      
      END
      SUBROUTINE LOSEMASS

******************************************************************
* This routine will estimate the mass loss for each parent star. *
* Any data beyond what is "lost" will be ignored.                *
******************************************************************

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jpmax(N),Pmax(N),I50(N),verbosity,I86(N),I95(N),K
      DOUBLE PRECISION Mp(NROWS,N),Rp(NROWS,N),
     &     ELp(NROWS,elmax,N),
     &     MLsfrac,Dp(NROWS,N),Ap(NROWS,N),MLoseTotal,
     &     periastron,As(NROWS,N),ELm(NROWS,elmax,N),
     &     Rp50(N),Rp86(N),Rp95(N)
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/MASSLOSS/MLsfrac
      COMMON/PERI/periastron
      COMMON/verbiage/verbosity
      COMMON/index/ I50
      SAVE Rp86, Rp95

      if(verbosity.ge.2) then
         WRITE(6,*)
         WRITE(6,*) 'Entering LOSEMASS FOR PERIASTRON=',periastron
         if(periastron.eq.0.d0)
     &    write(6,*)'EVEN IF YOU DID NOT REQUEST A HEADON COLLISION',
     &        'WE STILL HAVE TO WORK ON THE HEADON CASE FOR A WHILE.'
      endif

      IF(N.NE.2) THEN
         write(6,*) 'WARNING: MASS LOSS DOES NOT WORK UNLESS THERE ARE',
     &            ' TWO STARS'
         STOP
      ENDIF

* There should be no need to find the 50% radius if periastron is non-zero,
* because that should have already been done when treating the periastron=0
* case:
      if(periastron.eq.0.d0) then
         call geti50(I50,I86,I95)
         DO K=1,N
            Rp50(K)=((Mp(I50(K),K)/Mp(Jpmax(K),K)-.5d0)*Rp(I50(K)-1,K)+
     $           (.5d0-Mp(I50(K)-1,K)/Mp(Jpmax(K),K))*Rp(I50(K),K))/
     $           (Mp(I50(K),K)-Mp(I50(K)-1,K))*Mp(Jpmax(K),K)
            Rp86(K)=((Mp(I86(K),K)/Mp(Jpmax(K),K)-.86d0)*Rp(I86(K)-1,K)+
     $           (.86d0-Mp(I86(K)-1,K)/Mp(Jpmax(K),K))*Rp(I86(K),K))/
     $           (Mp(I86(K),K)-Mp(I86(K)-1,K))*Mp(Jpmax(K),K)
            Rp95(K)=((Mp(I95(K),K)/Mp(Jpmax(K),K)-.95d0)*Rp(I95(K)-1,K)+
     $           (.95d0-Mp(I95(K)-1,K)/Mp(Jpmax(K),K))*Rp(I95(K),K))/
     $           (Mp(I95(K),K)-Mp(I95(K)-1,K))*Mp(Jpmax(K),K)
        ENDDO
         if(verbosity.ge.2) then
            write(6,*)'R_{0.5} are',Rp50(1)/Rsun,' and',
     $           Rp50(2)/Rsun,' solar radii *** using these numbers'
            write(6,*)'R_{0.86} are',Rp86(1)/Rsun,' and',
     $           Rp86(2)/Rsun,' solar radii *** using these numbers'
            write(6,*)'R_{0.95} are',Rp95(1)/Rsun,' and',
     $           Rp95(2)/Rsun,' solar radii'
         endif
      endif

* Now we calculate the mass lost.  It depends on the 
* reduced mass, the radius of the star in question, and the radii at
* a mass fractions of 50% and 95% in both stars.  The coefficient
* const1 is determined 
* by calculating this information for a number of cases and averaging.

      MLoseTotal=
     &        const1*Mp(Jpmax(1),1)*Mp(Jpmax(2),2)*
     &        (Rp86(1)+Rp86(2))/
     &        ((Mp(Jpmax(1),1)+Mp(Jpmax(2),2))*
     &        (Rp(I50(1),1)+Rp(I50(2),2) +
     &         const2*periastron*(Rp(Jpmax(1),1)+Rp(Jpmax(2),2)) ))
      MLsfrac=MLoseTotal/(Mp(Jpmax(1),1)+Mp(Jpmax(2),2))

      if(verbosity.ge.1) write(6,'(a,f5.3,a,f5.3)')
     &' Mass loss fraction=',MLsfrac,' for periastron=',periastron
      if(verbosity.ge.2) then
         write(6,*)
     &', corresponding to a remnant mass of',(1.d0-MLsfrac)*
     &        (Mp(Jpmax(1),1)+Mp(Jpmax(2),2))/Msun,' Msun'
c         write(6,*)
c     &', corresponding to a table mass of',(1.d0-MLsfrac)*
c     &        ( NINT(100.d0*Mp(Jpmax(1),1)/Msun)/100.d0
c     $         +NINT(100.d0*Mp(Jpmax(2),2)/Msun)/100.d0 ),' Msun'
      endif
      RETURN
      END

******************************************************************
******************************************************************
******************************************************************
      subroutine geti50(I50,I86,I95)

******************************************************************
* This routine finds the locations in the parent inside of which *
* 50%, 86% and 95% of the total mass is enclosed                 *
******************************************************************
      
      implicit none
      INCLUDE 'params.h'
      INTEGER I,K,Jpmax(N),Pmax(N),I50(N),I86(N),I95(N)
      DOUBLE PRECISION Mp(NROWS,N),Rp(NROWS,N),
     &     ELp(NROWS,elmax,N),
     &     Dp(NROWS,N),Ap(NROWS,N),
     &     As(NROWS,N),ELm(NROWS,elmax,N)
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/NEWENT/Jpmax,Ap,Dp,As

      DO K=1,N
         DO I=1,NROWS
            IF(Mp(I,K).GE.50.d-2*Mp(Jpmax(K),K)) GOTO 47
         ENDDO
 47      I50(K)=I
         DO I=I50(K)+1,NROWS
            IF(Mp(I,K).GE.86.d-2*Mp(Jpmax(K),K)) GOTO 48
         ENDDO
 48      I86(K)=I
         DO I=I86(K)+1,NROWS
            IF(Mp(I,K).GE.95.d-2*Mp(Jpmax(K),K)) GOTO 49
         ENDDO
 49      I95(K)=I
      ENDDO

      END

      SUBROUTINE MIX(K)

*****************************************************************
* This routine will model the hydrodynamic mixing of the        *
* chemical compositions                                         *
* in the parent stars, after shock heating has occurred and     *
* before mass loss takes place.                                 *
*****************************************************************

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER C,K,L,I,J,Jpmax(N),Max,Pmax(N),verbosity
c     integer halfway
      DOUBLE PRECISION ELp(NROWS,elmax,N),As(NROWS,N),Mp(NROWS,N),
     &     DENOM(NROWS),CONTRIB,dM(NROWS),alpha,
     &     F(NROWS,NROWS),ELm(NROWS,elmax,N),ELmtot,
     &     ELptot,Dp(NROWS,N),Ap(NROWS,N),Rp(NROWS,N),dlogAdMavg
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/verbiage/verbosity
      common/numel/el
      double precision FBOUND(NROWS)
      character*16 OUTFN

      Max=Jpmax(K)

      DO I=2,Max-1
         dM(I)=0.5d0*(Mp(I+1,K)-Mp(I-1,K))
      ENDDO

      dM(1)=0.5d0*(Mp(2,K)+Mp(1,K))
      dM(Max)=0.5d0*(Mp(Max,K)-Mp(Max-1,K))

      dlogAdMavg=dlog(As(Pmax(K),K))-dlog(As(1,K))

      do I=1,Max
         FBOUND(I)=0.d0
      enddo

      alpha=-const8*dlogAdMavg**2.d0

      if(verbosity.ge.2) then
         write(6,*)'Mp(Max,K)=',Mp(Max,K)
         write(6,*)'dlogAdMavg=',dlogAdMavg,', alpha from paper=',-alpha
c Note that the alpha in the paper has the opposite sign from the alpha in
c this code
      endif

      DO I=1,Max
         DENOM(I)=0.d0
         DO J=1,Max
            F(J,I)=exp(alpha*((Mp(J,K)-Mp(I,K))/Mp(Max,K))**2.d0)
     &            +exp(alpha*((Mp(J,K)+Mp(I,K))/Mp(Max,K))**2.d0)
     &            +exp(alpha*((2.d0*Mp(Max,K)-Mp(J,K)
     &                                 -Mp(I,K))/Mp(Max,K))**2.d0)
            DENOM(I)=DENOM(I)+dM(J)*F(J,I)
         ENDDO
      ENDDO

      DO I=1,Max
         DENOM(I)=dM(I)/DENOM(I)
      ENDDO


      if(studymassloss) then
         WRITE(OUTFN,101) K
 101     FORMAT('fbound',I1,'.sbe')
         open(19,file=OUTFN)
         DO I=1,Max
            DO J=1,Max
C     F(J,I)*DENOM(J)*dM(I) is the mass that went from shell J to shell I
               IF(I.LE.Pmax(K))
     $              FBOUND(J)=FBOUND(J)+F(J,I)*DENOM(J)*dM(I)
            ENDDO
         ENDDO
         
         DO I=1,Max
            FBOUND(I)=FBOUND(I)/dM(I)
            write(19,*) Mp(I,K)/Msun,FBOUND(I)
         enddo
         close(19)
      endif

      DO C=1,el
         ELptot=0.d0
         ELmtot=0.d0
         DO L=1,Max
            CONTRIB=0.d0
            DO I=1,Max
               CONTRIB=CONTRIB+ELp(I,C,K)*F(I,L)*DENOM(I)
C Elp(I,C,K)*F(I,L)*DENOM(I)*dM(L) is the amount of mass of species C that
C went from shell I to shell L in star K
            ENDDO
            ELm(L,C,K)=CONTRIB
            ELptot=ELptot+ELp(L,C,K)*dM(L)
            ELmtot=ELmtot+ELm(L,C,K)*dM(L)
         ENDDO
         if(verbosity.ge.2) WRITE(6,*) 'Element:',C

         IF((ELptot-ELmtot)/ELptot.GT.1.d-6) THEN
            WRITE(6,*) 'WARNING!!  Mass may NOT be conserved!'
            WRITE(6,*) 'ELptot=',ELptot,' ELmtot=',ELmtot
         ENDIF

      ENDDO

      RETURN
      END

      SUBROUTINE PROFILE

********************************************************************
* Given as input the remnant file made by SORT with M, A, chemical *
* data, this subroutine will find values for the remnant's         *
* pressure,  density, and radius.                                  *
********************************************************************

* All units are cgs.

* Define variables: zeroin2 is a root-finding function
* Pclow & Pchigh are lower and upper limits for the guess of the central
* pressure Pcntr

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jr,Jrmax,STAR,NP,verbosity
      DOUBLE PRECISION Pcntr,Pclow,Pchigh,Psurf
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),P(NROWS),
     &     D(NROWS),R(NROWS),ELr(NROWS,elmax)
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      COMMON/VARIABLES/P,D,R
      COMMON/verbiage/verbosity
      DOUBLE PRECISION PRESSURE,zeroin2
      EXTERNAL PRESSURE

* Define the limits for zeroin2:
* These limits would have to be adjusted if you are creating remnants
* with either exceedingly low or exceedingly high central pressures...
      Pclow=0.0000000000001d0*G*Msun**2.d0/(Rsun**4.d0)
      Pchigh=3000.d0*G*Msun**2.d0/(Rsun**4.d0)

      Pcntr=zeroin2(Pclow,Pchigh,PRESSURE,1.d-4)
      if(verbosity.ge.2) WRITE(6,*) 'Central pressure=',Pcntr
   
*** DO NOT REMOVE THE FOLLOWING ASSIGNMENT STATEMENT - IT IS NEEDED!!***
*** ONE FINAL CALL TO PRESSURE IS NEEDED TO GET CORRECT REMNANT PROFILES!!***
      Psurf=PRESSURE(Pcntr)

      if(verbosity.ge.2) WRITE(6,*) 'Surface pressure=',Psurf
      IF(D(Jrmax).LE.0.d0) THEN
         P(Jrmax)=0.d0
         D(Jrmax)=1.d-30
         if(verbosity.ge.2) WRITE(6,*) 'Surface pressure set to 0'
      ENDIF

      Jr=1
      STAR=N+1
      CALL ENERGY(Mr(Jr),P(Jr),D(Jr),R(Jr),Jrmax,STAR)

      END

**************************************************************
      FUNCTION PRESSURE(Pctry)

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jr,Jrmax,T,NP,Jr2
      DOUBLE PRECISION dM,Pctry,Rup
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),P(NROWS),D(NROWS),
     &     R(NROWS),ELr(NROWS,elmax)
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      COMMON/VARIABLES/P,D,R
      DOUBLE PRECISION PRESSURE,RADIUS,zeroin
	double precision Rlow
      EXTERNAL RADIUS

      T=2

      P(1)=Pctry
      D(1)=(dabs(Pctry)/Ar(1))**(3.d0/5.d0)
      R(1)=(3.d0*Mr(1)/(4.d0*pi*D(1)))**(1.d0/3.d0)
      
      IF(R(1).EQ.0.d0) THEN
         dM=Mr(2)-Mr(1)
         R(2)=(3.d0*dM/(4.d0*pi*D(1)))**(1.d0/3.d0)
         P(2)=Pctry-2.d0*pi/3.d0*G*R(2)**2.d0*D(1)**2.d0
         D(2)=(dabs(P(2))/Ar(2))**(3.d0/5.d0)
         T=3
      ENDIF

      DO Jr=T,Jrmax
         dM=Mr(Jr)-Mr(Jr-1)
         Rup=R(Jr-1)+10.d0*dM/(4.d0*pi*(R(Jr-1)**2.d0)*
     &        dabs(D(Jr-1)))
         
c     Here's the thing in RADIUS that we want to avoid becoming negative:
c     PJr=P(Jr-1)-G/(8.d0*pi)*(Mr(Jr-1)/R(Jr-1)**4.d0+
c     &     Mr(Jr)/RJr**4.d0)*(Mr(Jr)-Mr(Jr-1))

         if(P(Jr-1)-G/(8.d0*pi)*(Mr(Jr-1)/R(Jr-1)**4.d0
     &        )*dM.ge.0.d0) then
C     Any problem with Rlow or Rup could be avoided...
C     
C     Let's check if there is a problem with using Rlow=R(Jr-1)...
            if(P(Jr-1)-G/(8.d0*pi)*(Mr(Jr-1)/R(Jr-1)**4.d0+
     &           Mr(Jr)/R(Jr-1)**4.d0)*dM.le.0.d0) then
               Rlow=1.000001d0*( (P(Jr-1)*8.d0*pi/(G*dM)
     &              -Mr(Jr-1)/R(Jr-1)**4.d0)
     &              /Mr(Jr)
     &              )**(-0.25d0)
               Rup=2.d0*Rlow
            else
C     There's no problem with the default Rlow, and so Rup must also be OK
               Rlow=R(Jr-1)
            endif
            
            R(Jr)=zeroin(Rlow,Rup,RADIUS,1.d-5)
         else
C     It's hopeless for fixing Rlow or Rup problems, fudge and move on.
C     This model is not going to be used anyway...
            R(Jr)=Rup
            do Jr2=Jr,Jrmax
               P(Jr2)=0.d0
               D(Jr2)=0.d0
            enddo
C     Set surface pressure to some vaguely appropriate negative value that is
C     large in magnitude, so that the algorithm knows it has guessed a bad
C     central pressure on this attempt.
            PRESSURE=P(Jr-1)-G/(8.d0*pi)*(Mr(Jr-1)/R(Jr-1)**4.d0+
     &           Mr(Jr)/R(Jr)**4.d0)*dM
     $           +4.d14*(Mr(Jr)/Mr(Jrmax)-1.d0) -1.d4
            RETURN
         endif

         P(Jr)=P(Jr-1)-G/(8.d0*pi)*(Mr(Jr-1)/R(Jr-1)**4.d0+
     &        Mr(Jr)/R(Jr)**4.d0)*dM
         D(Jr)=(dabs(P(Jr))/Ar(Jr))**(3.d0/5.d0)
      ENDDO

      PRESSURE=P(Jrmax)

      RETURN
      END


*****************************************************************
      FUNCTION RADIUS(RJr)

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jr,Jrmax,NP
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),P(NROWS),R(NROWS)
      DOUBLE PRECISION PJr,RJr,RADIUS,ELr(NROWS,elmax),D(NROWS)
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      COMMON/VARIABLES/P,D,R

      PJr=P(Jr-1)-G/(8.d0*pi)*(Mr(Jr-1)/R(Jr-1)**4.d0+
     &     Mr(Jr)/RJr**4.d0)*(Mr(Jr)-Mr(Jr-1))

      RADIUS=-RJr+R(Jr-1)+(Mr(Jr)-Mr(Jr-1))/(8.d0*pi)*
     &     ((Ar(Jr-1)/P(Jr-1))**(3.d0/5.d0)/
     &     (R(Jr-1)**2.d0)+(Ar(Jr)/PJr)**(3.d0/5.d0)/(RJr**2.d0))

      RETURN
      END






      SUBROUTINE SHOCK(INTERCEPT)

*********************************************************************
* This routine will take the parent star entropy profiles and shock *
* them, creating a new entropy profile for SORT to accept and use   *
* to calculate the remnant's entropy.                               *
*********************************************************************

      IMPLICIT NONE
      INCLUDE 'params.h'

* Mto= turn-off mass in solar units
* Rto= turn-off radius in solar units
* (The use of Mto and Rto is simply because of the chosen units for the
* intercept.)
      DOUBLE PRECISION Mto,Rto
      PARAMETER(Mto=0.8d0*Msun, Rto=0.955d0*Rsun)

      INTEGER Jpmax(N),I,Imax,K,Pmax(N),I1,I2,verbosity
      DOUBLE PRECISION Ap(NROWS,N),As(NROWS,N),Dp(NROWS,N),INTERCEPT,
     &     Mp(NROWS,N),Rp(NROWS,N),ELp(NROWS,elmax,N),INTERCEPT2,
     &     MLsfrac,MLoseTotal,massoutsideofI1,massoutsideofI2,
     &     periastron,ELm(NROWS,elmax,N)
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/MASSLOSS/MLsfrac
      COMMON/PERI/periastron
      COMMON/verbiage/verbosity
      double precision RHS

* The equation for the shocked entropy, As, comes from a plot of
* Log(A-Ainit) vs. Log(Arho^5/3).  The value INTERCEPT is the y-intercept
* of the linear function fitted to the plot of the data.  The variable
* "const3" is the slope of that same function.

      DO K=1,N
         Imax=Jpmax(K)
         IF(K.EQ.1) THEN
            DO I=1,Imax
               RHS=INTERCEPT + dlog10(G*Mto**(1.d0/3.d0)*Rto)
     &              +const3*
     $              dlog10(Rto**4/Mto**2/G*Ap(I,K)*Dp(I,K)**(5.d0/3.d0))
               As(I,K)=Ap(I,K)+10.d0**RHS
c               As(I,K)=Ap(I,K)+10.d0**INTERCEPT*
c     &          G*Mto**(1.d0/3.d0)*Rto*
c     &          (Rto**4/Mto**2/G*Ap(I,K)*Dp(I,K)**(5.d0/3.d0))**const3
            ENDDO
         ELSE
            INTERCEPT2=INTERCEPT+((const4+const5)*periastron-const6)*
     &           LOG10(Mp(Jpmax(1),1)/Mp(Jpmax(K),K))
            DO I=1,Imax
               RHS=INTERCEPT2 + dlog10(G*Mto**(1.d0/3.d0)*Rto)
     &              +const3*
     $              dlog10(Rto**4/Mto**2/G*Ap(I,K)*Dp(I,K)**(5.d0/3.d0))
c               As(I,K)=Ap(I,K)+10.d0**INTERCEPT2*
c     &          G*Mto**(1.d0/3.d0)*Rto*
c     &          (Rto**4/Mto**2/G*Ap(I,K)*Dp(I,K)**(5.d0/3.d0))**const3
               As(I,K)=Ap(I,K)+10.d0**RHS
            ENDDO
            if(verbosity.ge.2) then
               write(6,*) 'Intercept for star 2=',INTERCEPT2
               write(6,*) 'Intercept for star 1=',INTERCEPT
            endif
         ENDIF
      ENDDO

      MLoseTotal=MLsfrac*(Mp(Jpmax(1),1)+Mp(Jpmax(2),2))
      I1=Jpmax(1)
      massoutsideofI1=0.d0
      massoutsideofI2=0.d0
      Do while (massoutsideofI1+massoutsideofI2 .lt. MLoseTotal)
         I1=I1-1
         massoutsideofI1=massoutsideofI1 + Mp(I1+1,1) - Mp(I1,1)

         I2=Jpmax(2)
         massoutsideofI2=0.d0
         do while(As(I2,2).gt.As(I1,1))
           I2=I2-1
           massoutsideofI2=massoutsideofI2 + Mp(I2+1,2) - Mp(I2,2)
         enddo
      enddo

      Pmax(1)=I1
      Pmax(2)=I2

      if(verbosity.ge.2) then
         write(6,*)'Pmax(1)=',Pmax(1),' Pmax(2)=',Pmax(2)
         write(6,*)'actual mass loss fraction=',
     &        (massoutsideofI1+massoutsideofI2)/
     &        (Mp(Jpmax(1),1)+Mp(Jpmax(2),2)),
     &        ' for r_p=',periastron
      endif

      RETURN
      END

      
      SUBROUTINE SORT

******************************************************************
* This subroutine will take the parent star data                 *
* of M, A and chemical abundances.  It will create similar arrays*
* for the remnant, sorted by increasing A.  In order             *
* to keep the remnant profile smooth, we'll account for the      *
* profile derivatives.                                           *
******************************************************************

* Define the variables:

* M is the mass enclosed in a star, A is the entropic variable at
* M, and EL is the chemical abundance at M.  Ji is the Jith row in 
* the data file for a star.

      IMPLICIT NONE
      INCLUDE 'params.h'
      INTEGER Jr,Jp(N),Jrmax,K,C,Pmax(N),Jpmax(N),NP
      DOUBLE PRECISION Mp(NROWS,N),ELp(NROWS,elmax,N),As(NROWS,N)
      DOUBLE PRECISION Mr(NROWS),Ar(NROWS),ELr(NROWS,elmax),Rp(NROWS,N)
      DOUBLE PRECISION dAr,dAsdMp(N),dArdMr,dArp(N),Mtmp,dAdMtmp,
     &     Armax,MofAr(N),Chemtmp(elmax),dArtmp,Chem(elmax,N),
     &     ELm(NROWS,elmax,N),
     &     Dp(NROWS,N),Ap(NROWS,N)
      COMMON/NEWINFO/ELp,Mp,Rp,Pmax,ELm
      COMMON/NEWENT/Jpmax,Ap,Dp,As
      COMMON/REMNANT/Mr,Ar,ELr,Jrmax,Jr,NP
      LOGICAL NOTDONE
      common/numel/el

* The maximum entropy in the remnant should be the highest in any of
* the parents:

      Armax=0.d0
      DO K=1,N
         Armax=max(Armax,As(Pmax(K),K))
      ENDDO

c      write(6,*)'Armax=',Armax,' Pmax(1)=',Pmax(1),' Pmax(2)=',Pmax(2)

* Initialize the rows:

      DO K=1,N
         Jp(K)=1
      ENDDO
      Jr=1

* Initialize the entropy - the first value of Ar (at the lowest Mr)
* will be the lowest entropy in any of the parent stars.

      Ar(1)=1.d30
      DO K=1,N
         Ar(1)=min(Ar(1),As(1,K))
      ENDDO

C      dlogAr=(DLOG10(Armax)-DLOG10(Ar(1)))/NP

* Define the logical variable that will let the program leave the
* DO loop at the right time, and initialize dAr.

      NOTDONE=.TRUE.

* Begin the DO loop.

      DO WHILE (NOTDONE)

* If Ar gets to be the maximum, then we're done. (hence the IF statement
* below).

         IF(Ar(Jr).GE.0.999999d0*Armax) NOTDONE=.FALSE.

* We must do a loop over each parent star.

         DO K=1,N
            DO WHILE (As(Jp(K),K).LE.Ar(Jr) .AND. Jp(K).LE.Pmax(K))
               Jp(K)=Jp(K)+1
            ENDDO
            IF(As(Pmax(K),K).LE.Ar(Jr) .AND. 
     &           As(Pmax(K),K).GE.Ar(Jr)-1.d-06) THEN
               Jp(K)=Pmax(K)
            ENDIF

            IF(Jp(K).GE.2 .AND. Jp(K).LE.Pmax(K)) THEN
               dAsdMp(K)=(As(Jp(K),K)-As(Jp(K)-1,K))/
     &              (Mp(Jp(K),K)-Mp(Jp(K)-1,K))
               MofAr(K)=(Ar(Jr)-As(Jp(K)-1,K))*Mp(Jp(K),K)/(As(Jp(K),K)-
     &              As(Jp(K)-1,K))+(As(Jp(K),K)-Ar(Jr))*Mp(Jp(K)-1,K)/
     &              (As(Jp(K),K)-As(Jp(K)-1,K))
               dArp(K)=As(Jp(K),K)-As(Jp(K)-1,K)
               DO C=1,el
                  Chem(C,K)=(Ar(Jr)-As(Jp(K)-1,K))*ELm(Jp(K),C,K)/
     &                 (As(Jp(K),K)-As(Jp(K)-1,K))+(As(Jp(K),K)-Ar(Jr))
     &                 *ELm(Jp(K)-1,C,K)/(As(Jp(K),K)-As(Jp(K)-1,K))
               ENDDO
            ELSE IF(Jp(K).EQ.1) THEN
               dAsdMp(K)=1.d30
               DO C=1,el
                  Chem(C,K)=0.d0
               ENDDO
               dArp(K)=1.d30
               MofAr(K)=0.d0
            ELSE IF(Jp(K).GT.Pmax(K)) THEN
               dAsdMp(K)=1.d30
               DO C=1,el
                  Chem(C,K)=0.d0
               ENDDO
               dArp(K)=1.d30
               MofAr(K)=Mp(Pmax(K),K)
            ENDIF
         ENDDO

         dAdMtmp=0.d0
         DO C=1,el
            Chemtmp(C)=0.d0
         ENDDO
         Mtmp=0.d0
         dArtmp=1.d30
         DO K=1,N
            dAdMtmp=1.d0/dAsdMp(K)+dAdMtmp
            DO C=1,el
               Chemtmp(C)=Chem(C,K)/dAsdMp(K)+Chemtmp(C)
            ENDDO
            Mtmp=MofAr(K)+Mtmp
            dArtmp=1.d0*min(dArp(K),dArtmp)
         ENDDO

         dArdMr=dAdMtmp**(-1.d0)
         DO C=1,el
            ELr(Jr,C)=dArdMr*Chemtmp(C)
         ENDDO
         Mr(Jr)=Mtmp

         IF(NP.lt.0) then
C N.B.!  The NP here is 1 less than the NP in the driver program...
C Use the smallest reasonable dAr if NP<0 : 
            dAr=dArtmp/abs(NP)
         else

            dAr=dArdMr*(Mp(Pmax(1),1)+Mp(Pmax(2),2))/(NP-3.5d0)

C These lines of data will have the *logarithm*
C of the entropic variable A spaced regularly...
c         dAr=10.d0**(DLOG10(Ar(Jr))+dlogAr)-Ar(Jr)

         endif

C The following line of code forces more mesh points in the envelope:
         IF(Ar(Jr)+2.d0*dAr.GE.0.999999d0*Armax) dAr=dAr/100.d0

         IF(Ar(Jr)+dAr.GE.0.999999d0*Armax) THEN
            dAr=Armax-Ar(Jr)
         ENDIF       

         Jr=Jr+1
         Ar(Jr)=Ar(Jr-1)+dAr

         IF(Jr.GT.NROWS) THEN
            write(6,*) 'WARNING!  You may need to increase NROWS in'
            write(6,*) 'the header file.  NROWS=',NROWS
            STOP
         ENDIF

      ENDDO

      Jrmax=Jr-1

      END

*************************************************************
* zeroin minimizes the radius                               *
*************************************************************
      double precision function zeroin(ax,bx,f,tol)
      double precision ax,bx,f,tol
C
C This routine is a very slightly modified version of a routine
C from the book Computer Methods for Mathematical
C Computations, by George Forsythe, Mike Malcolm, and Cleve Moler.
C
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0d0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs,dsign
c
c  compute eps, the relative machine precision
c


c      eps = 1.0d0
c   10 eps = eps/2.0d0
c      tol1 = 1.0d0 + eps
c      if (tol1 .gt. 1.0d0) go to 10

      eps=3.d-8

c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
      if(fa*fb.gt.0.d0)then
         write(6,*)
         write(6,*)'You must be pushing the limits of this code.'
         write(6,*)'Try letting Rup be even larger [increase the'
         write(6,*)'size of the term being added to R(Jr-1)]'
         write(6,*)'and see if that helps you get a model...'
         write(6,*)'Also consider emailing lombardi@vassar.edu'
         pause 'root must be bracketed for zeroin'
      endif
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
      xm = .5d0*(c - b)
      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = dabs(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 zeroin = b
      return
      end

*************************************************************
* zeroin2 minimizes the pressure                            *
*************************************************************
      double precision function zeroin2(ax,bx,f,tol)
      double precision ax,bx,f,tol
C
C This routine is a very slightly modified version of a routine
C from the book Computer Methods for Mathematical
C Computations, by George Forsythe, Mike Malcolm, and Cleve Moler.
C
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0d0)
c
c
c  output..
c
c  zeroin2 abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin2  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs,dsign
c
c  compute eps, the relative machine precision
c


c      eps = 1.0d0
c   10 eps = eps/2.0d0
c      tol1 = 1.0d0 + eps
c      if (tol1 .gt. 1.0d0) go to 10

      eps=3.d-8

c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
      if(fa*fb.gt.0.d0)then
         write(6,*)
         write(6,*)'You must be pushing the limits of this code.'
         write(6,*)'Try expanding the range of Pclow and Pchigh'
         write(6,*)'and see if that helps you get a model...'
         write(6,*)'Also consider emailing lombardi@vassar.edu'
         pause 'root must be bracketed for zeroin2'
      endif
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
      xm = .5d0*(c - b)
      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = dabs(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 zeroin2 = b
      return
      end

*************************************************************
* zeroin3 used to minimize delta Energy                     *
*************************************************************
      double precision function zeroin3(ax,bx,f,tol)
      double precision ax,bx,f,tol
C
C This routine is a very slightly modified version of a routine
C from the book Computer Methods for Mathematical
C Computations, by George Forsythe, Mike Malcolm, and Cleve Moler.
C
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0d0)
c
c
c  output..
c
c  zeroin3 abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin3  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs,dsign
c
c  compute eps, the relative machine precision
c


c      eps = 1.0d0
c   10 eps = eps/2.0d0
c      tol1 = 1.0d0 + eps
c      if (tol1 .gt. 1.0d0) go to 10

      eps=3.d-8

c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
      if(fa*fb.gt.0.d0)then
         write(6,*)
         write(6,*)
         write(6,*)'You must be pushing the limits of this code.'
         write(6,*)'Try expanding the range of INTlow and INThigh'
         write(6,*)'and see if that helps you get a model...'
         write(6,*)'Also consider emailing lombardi@vassar.edu'
         pause 'root must be bracketed for zeroin3'
      endif
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*dabs(b) + 0.5d0*tol
      xm = .5d0*(c - b)
      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = dabs(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 zeroin3 = b
      return
      end
