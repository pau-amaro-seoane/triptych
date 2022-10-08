************************************************************************
*      MODEST_EVOLUTION.F 
*
*      THIS IS STELLAR EVOLUTION MODULE #2 (stored in sev2)
*      [comprised of modest_evolution.f, modest_star.h, libstr.a]
*      Version 0.3 January 9, 2003
*            Written by Jarrod Hurle some minor modifications by James Lombardi
*
*      Requires libstr.a - SSE library of stellar evolution 
*      functions described in the paper 
*      "Comprehensive Analytic Formulae for Stellar Evolution as a 
*       Function of Mass and Metallicity", 
*       Hurley J.R., Pols O.R. and Tout C.A. 2000, MNRAS, 315, 543. 
*
*      The function EvolveStar calls the SSE library functions and 
*      evolves the star for a user specified interval controlled 
*      by dtmax (change in time/Myr), dMmax (change in mass/Msun), 
*      dRmax (change in radius/Rsun), dYmax (change in helium content), 
*      and dZmax (change in metallicity). 
*      EvolveStar recommends values of dtmax, dMmax, and dRmax, for 
*      the next step and these may or may not be adopted by the user. 
*      The following functions are supplied for the purpose of 
*      extracting stellar data for star i: 
*
*          getTime(i)   - current age of star (Myr)
*          getMass(i)   - current mass (Msun)
*          getRadius(i) - current radius (Rsun)
*          getY(i)      - helium fraction (0->1)
*          getZ(i)      - metallicity (0->1)
*          getMStime(i) - MS lifetime of the star (Myr)
*          getLum(i)    - luminosity (Lsun)
*          getMc(i)     - core mass (Msun)
*          getRc(i)     - core radius (Rsun)
*          getSpin(i)   - stellar rotation (1/yr)
*          getType(i)   - stellar type (0->15, e.g. MS star = 1.0)
*
*      Stellar variables required by the stellar evolution module are 
*      declared in modest_star.h which is provided along with the module. 
*
*      The subroutine one_star provides the link between the main 
*      program and the stellar evolution module. It controls the 
*      calls to EvolveStar and associated routines (CreateStar, 
*      print_star, etc.). one_star is not technically part of the 
*      stellar evolution module. 
* 
************************************************************************
*
      SUBROUTINE one_star(m_init,r_init,Y_init,Z_init,outfile)
*
      implicit none
      include "modest_common.h"
      integer i,idstar,iunit,iteration
      real*8 t,tmax,kwold,kw
      character*70 outfile
      real*8 time,mass,radius,y,z,loglum,teff
      real*8 oldtime,oldmass,oldradius,oldy,oldz,oldloglum,oldteff
*
* Allocate a reference id for the star and initialize stellar 
* quantities required by EvolveStar. 
*
      idstar = CreateStar(m_init,r_init,Y_init,Z_init)
      if(idstar.lt.0)then
         WRITE(6,*)' ERROR: star index not found'
         RETURN
      endif
*
* Set the start and end-point of the evolution and initialize 
* the evolution timestep. 
*
      dtmax = .5d0
      t = 0.d0
      tmax = 15658.3d0
*
* Open the output file. 
*
      iunit=50
      if(iunit.ne.6)then
         OPEN(iunit,file=outfile)
         WRITE(6,*)'... WRITING EVOLUTION DATA TO FILE'
      endif
      WRITE(iunit,*)
     $'# 1:Age [Myr] 2:M [M_sun] 3:R [R_sun] 4:Y 5:Z 6:log_10 L [L_sun]'
     $,' 7:T_eff [K] 8:type'
*
* Evolution loop. 
* Note that EvolveStar returns recommend values of dtmax, dMmax, and 
* dRmax for the next step but of these only dtmax is adopted at present. 
* Note also that it is not necessary for the entire evolution of the 
* star to be performed in a single loop - as is the case in this example. 
* 
      iteration=0
      kw=0
 10   iteration=iteration+1
      dMmax = 0.1d0*getMass(idstar)
      dRmax = .5d0*getRadius(idstar)
      dYmax = 0.1d0
      dZmax = 0.1d0
      dtmax = MIN(dtmax,tmax-age(idstar))
*
      oldtime=time
      oldmass=mass
      oldradius=radius
      oldy=y
      oldz=z
      oldloglum=loglum
      oldteff=teff
      kwold=kw

      t = EvolveStar(idstar,dtmax,dmmax,drmax,dymax,dzmax)

      time=getTime(idstar)
      mass=getMass(idstar)
      radius=getRadius(idstar)
      y=getY(idstar)
      z=getZ(idstar)
      loglum=dlog10(getLum(idstar))
      teff=getTeff(idstar)
      kw=getType(idstar)

      if(iteration.gt.1 .and. kwold.ne.kw) then
         WRITE(iunit,600) oldtime,oldmass,oldradius,oldy,oldz,
     $     oldloglum,oldteff,kwold
      endif
 600  format(3g13.8,1x,2g13.7,1x,2g13.7,1g7.2)


      if(kwold.ne.kw .and. kw.ge.11) tmax=2.d0*time

      if(mod(iteration,10).eq.0 .or. kwold.ne.kw .or. kw.ge.7) then
         WRITE(iunit,600) time,mass,radius,y,z,
     $     loglum,teff,kw
      endif
*
      if(t.lt.tmax) goto 10
*
* Close the output file and free-up the particle id number. 
*
      if(iunit.ne.6) close(iunit)
      i = finished_star(idstar)
*
      RETURN
      END
***
c      SUBROUTINE print_star(idstar,iunit,printlasttoo)
c      include "modest_common.h"
c      integer idstar,iunit
c      SAVE time,mass,radius,y,z,loglum,teff
c      logical printlasttoo
c
c
c      WRITE(iunit,600) time,mass,radius,y,z,loglum,teff,type
c      RETURN
c      END
***
      block data block
      include "modest_common.h"
      data is_used/nmax*0/
      end
***
      FUNCTION first_unused_star()
      include "modest_common.h"
      integer i
      do i=1,nmax
         if(is_used(i).eq.0)then
            first_unused_star = i
            RETURN
         endif
      enddo
      first_unused_star = -1
      RETURN
      END
***
      FUNCTION finished_star(idstar)
      include "modest_common.h"
      integer idstar
      if(idstar.lt.0.or.idstar.gt.nmax)then
         finished_star = -1
      else
         if(is_used(idstar).eq.0)then
            finished_star = -1
         else
            finished_star = idstar
            is_used(idstar) = 0
         endif
      endif
      RETURN
      END
***
*
* "Black Box" starts here. 
*
***
      FUNCTION CreateStar(m_init,r_init,Y_init,Z_init) 
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
*
      idstar = first_unused_star()
      if(idstar.lt.0)then
         CreateStar= -1
         RETURN
      endif
      CreateStar = idstar
*
      age(idstar) = 0.d0
      zmass0(idstar) = m_init
      zmass(idstar) = zmass0(idstar)
      zinit(idstar) = Z_init
      yinit(idstar) = Y_init
      is_used(idstar) = 1
      radius(idstar) = r_init
*
      RETURN
      END
***
      FUNCTION getTime(idstar)
      include "modest_common.h"
      integer idstar
      getTime = age(idstar)
      RETURN
      END
***
      FUNCTION getMass(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      getMass = zmass(idstar)
      RETURN
      END
***
      FUNCTION getRadius(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      getRadius = radius(idstar)
      RETURN
      END
***
      FUNCTION getY(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      getY = ycurr(idstar)
      RETURN
      END
***
      FUNCTION getZ(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      getZ = zcurr(idstar)
      RETURN
      END
***
      FUNCTION getMStime(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      getMStime = ms_lifetime(idstar)
      RETURN
      END
***
      FUNCTION getLum(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      getLum = zlum(idstar)
      RETURN
      END
***
      FUNCTION getMc(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      getMc = zmassc(idstar)
      RETURN
      END
***
      FUNCTION getRc(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      getRc = radc(idstar)
      RETURN
      END
***
      FUNCTION getSpin(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      getSpin = spin(idstar)
      RETURN
      END
***
      FUNCTION getType(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      
      getType = DBLE(kstar(idstar))
      RETURN
      END
***
      FUNCTION getTeff(idstar)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      double precision Rsun, Lsun, sigma
      PARAMETER(Rsun=6.9598d10, Lsun=3.826d33, sigma=5.67051d-5)
      getTeff = (zlum(idstar)*Lsun/
     $     ( 16.d0*atan(1.d0)*(radius(idstar)*Rsun)**2*sigma )
     $          )**0.25d0
      RETURN
      END
***
      FUNCTION EvolveStar(i,dtmax,dmmax,drmax,dymax,dzmax)
c
c-------------------------------------------------------------c
c
c     Evolves a single star with mass loss.
c     The timestep is not constant but determined by certain criteria.
c
c     Note: no chemical evolution in SSE at this time 
c           (Uses simple evolution of Y and Z provided by 
c            Module #1; dymax & dzmax are ignored).
c
c     Written by Jarrod Hurley 05/10/02 at AMNH, NY. 
c
c     Requires the following variables for star i via common 
c     from the main routine: 
c        age(i)    = time to which star was last evolved (Myr)
c        epoch(i)  = additional time variable used for rejuvenated 
c                    or remnant stars (Myr) 
c        ms_lifetime(i) = lifetime on MS (Myr)
c        standard_timestep(i) = recommended evolution timestep (Myr)
c        zmass0(i) = ZAMS mass (Msun)
c        zmass(i)  = current mass (Msun)
c        zmassc(i) = current core mass (Msun)
c        zinit(i)  = ZAMS metallicity (0.02 is solar)
c        zcurr(i)  = current metallicity 
c        yinit(i)  = ZAMS helium content (X = 0.76 - 3*Z)
c        ycurr(i)  = current helium content 
c        spin(i)   = spin frequency (1/yr). 
c        kstar(i)  = index of current stellar type (0-15, see below). 
c     Also saves the following useful quantities in the same common: 
c        radius(i) = current stellar radius (Rsun)
c        radc(i)   = current core radius (Rsun)
c        zlum(i)   = current stellar luminosity (Lsun)
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
      implicit none
      include "modest_common.h"
      include "modest_star.h"
*
      integer i,j,it,kw,kwold,nv
      parameter(nv=50000)
*
      real*8 mass,mt,z,aj,tm,tn
      real*8 tphys,tphysf,tphys2,tmold,tbgold
      real*8 tscls(20),lums(10),GB(10),zpars(20)
      real*8 r,lum,mc,rc,menv,renv
      real*8 ospin,jspin,djt,djmb,k2,k3
      parameter(k3=0.21d0)
      real*8 m0,mt2,mc1,ajhold,rm0
      real*8 dt,dtm,dtr,dtm0,dr,dtdr,dms,dml,rl
      real*8 tiny,alpha2
      parameter(tiny=1.0d-14,alpha2=0.09d0)
      real*8 pts1,pts2,pts3
      parameter(pts1=0.0005,pts2=0.01,pts3=0.02)
C      parameter(pts1=0.05,pts2=0.01,pts3=0.02)
      real*8 mlwind,rzamsf,vrotf
      external mlwind,rzamsf,vrotf
      real*8 neta,bwind
      common /value1/ neta,bwind
*
      tphys = age(i)
      tphysf = tphys + dtmax
      dtm = 0.d0
      mass = zmass0(i)
      mt = zmass(i)
      z = zinit(i)
      CALL zcnsts(z,zpars)
*
      if(tphys.lt.tiny)then
         kstar(i) = 1
         epoch(i) = 0.d0
         zmassc(i) = 0.d0
         radius(i) = rzamsf(mt)
         radc(i) = 0.d0
         spin(i) = 45.35d0*vrotf(mt)/radius(i)
         zcurr(i) = zinit(i)
         ycurr(i) = yinit(i)
      endif
      mc = zmassc(i)
      r = radius(i)
      ospin = spin(i)
      aj = tphys - epoch(i)
      kw = kstar(i)
      k2 = 0.15d0
*
      CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
      CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &            r,lum,kw,mc,rc,menv,renv,k2)
      jspin = ospin*(k2*r*r*(mt-mc)+k3*rc*rc*mc)
*
      mc1 = 0.d0
      rl = 0.d0
      neta = 0.5d0
      bwind = 0.d0
*
      do 10 , j = 1,nv
*
* Base new time scale for changes in radius & mass on stellar type.
*
         if(kw.le.1)then
            dtm = pts1*tm
            dtr = tm - aj
         elseif(kw.eq.2)then
            dtm = pts3*(tscls(1) - tm)
            dtr = tscls(1) - aj
         elseif(kw.eq.3)then
            if(aj.lt.tscls(6))then
               dtm = pts2*(tscls(4) - aj)
            else
               dtm = pts2*(tscls(5) - aj)
            endif
            dtr = MIN(tscls(2),tn) - aj
         elseif(kw.eq.4)then
            dtm = pts2*tscls(3)
            dtr = MIN(tn,tscls(2) + tscls(3)) - aj
         elseif(kw.eq.5)then
            if(aj.lt.tscls(9))then
               dtm = pts2*(tscls(7) - aj)
            else
               dtm = pts2*(tscls(8) - aj)
            endif
            dtr = MIN(tn,tscls(13)) - aj
         elseif(kw.eq.6)then
            if(aj.lt.tscls(12))then
               dtm = pts2*(tscls(10) - aj)
            else
               dtm = pts2*(tscls(11) - aj)
            endif
            dtm = MIN(dtm,0.005d0)
            dtr = tn - aj
         elseif(kw.eq.7)then
            dtm = pts3*tm
            dtr = tm - aj
         elseif(kw.eq.8.or.kw.eq.9)then
            if(aj.lt.tscls(6))then
               dtm = pts2*(tscls(4) - aj)
            else
               dtm = pts2*(tscls(5) - aj)
            endif
            dtr = tn - aj
         else
            dtm = MAX(0.1d0,aj*10.d0)
            dtr = dtm
         endif
         dtm0 = dtm
         dtm = MIN(dtm,dtr)
*
* Choose minimum of time-scale and remaining interval (> 100 yrs).
*
         dtm = MAX(dtm,1.0d-07*aj)
         dtm = MIN(dtm,tphysf-tphys)
*
         rm0 = r
         ajhold = aj
         tmold = tm
*
* Calculate mass loss.
*
         if(kw.lt.10)then
            dt = 1.0d+06*dtm
            dms = mlwind(kw,lum,r,mt,mc,rl,z)*dt
            dml = mt - MAX(mc,zmass(i)-dmmax)
            if(dml.lt.dms)then
               dtm = (dml/dms)*dtm
               dms = dml
            endif
         else
            dms = 0.d0
         endif
*
* Limit to 1% mass loss.
*
         if(dms.gt.0.01d0*mt)then
            dtm = 0.01d0*mt*dtm/dms
            dms = 0.01d0*mt
         endif
*
* Calculate the rate of angular momentum loss due to magnetic braking
* and/or mass loss.
*
         if(dtm.gt.tiny)then
            djt = (2.d0/3.d0)*(dms/(1.0d+06*dtm))*r*r*ospin
            if(mt.gt.0.35d0.and.kw.lt.10)then
               djmb = 5.83d-16*menv*(r*ospin)**3/mt
               djt = djt + djmb
            endif
         endif
*
* Update mass and time and reset epoch for a MS (and possibly a HG) star.
*
         if(dms.gt.tiny)then
            mt = mt - dms
            if(kw.le.2.or.kw.eq.7)then
               m0 = mass
               mc1 = mc
               mass = mt
               tbgold = tscls(1)
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               if(kw.eq.2)then
                  if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                     mass = m0
                  else
                     epoch(i) = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                            (tbgold - tmold)
                     epoch(i) = tphys - epoch(i)
                  endif
               else
                  epoch(i) = tphys - ajhold*tm/tmold
               endif
            endif
         endif
         tphys2 = tphys
         tphys = tphys + dtm
*
* Find the landmark luminosities and timescales as well as setting
* the GB parameters.
*
         aj = tphys - epoch(i)
         CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
*
* Find the current radius, luminosity, core mass and stellar type
* given the initial mass, current mass, metallicity and age
*
         kwold = kw
         CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &               r,lum,kw,mc,rc,menv,renv,k2)
*
* Check that the radius has not exceeded its allowed change 
* (unless a type change has occurred). 
*
         if(kw.eq.kwold)then
            mt2 = mt + dms
            dml = dms/dtm
            it = 0
 20         dr = r - radius(i)
            if(ABS(dr).gt.drmax)then
               it = it + 1
               if(it.eq.20.and.kw.eq.4) goto 30
               if(it.ge.30) goto 30
               dr = r - rm0
               if(dr.lt.tiny) goto 30
               dtdr = dtm/ABS(dr)
               dtm = alpha2*MAX(rm0,r)*dtdr
               if(it.ge.20) dtm = 0.5d0*dtm
               if(dtm.lt.1.0d-07*aj) goto 30
               dms = dtm*dml
               mt = mt2 - dms
               if(kw.le.2.or.kw.eq.7)then
                  mass = mt
                  CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
                  if(kw.eq.2)then
                     if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                        mass = m0
                     else
                        epoch(i) = tm + (tscls(1) - tm)*(ajhold-tmold)/
     &                               (tbgold - tmold)
                        epoch(i) = tphys2 - epoch(i)
                     endif
                  else
                     epoch(i) = tphys2 - ajhold*tm/tmold
                  endif
               endif
               tphys = tphys2 + dtm
               aj = tphys - epoch(i)
               mc = mc1
               CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
               CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                     r,lum,kw,mc,rc,menv,renv,k2)
               goto 20
            endif
 30         continue
         endif
*
* Adjust the spin of the star and reset epoch.
*
         jspin = MAX(1.0d-10,jspin - djt*1.0d+06*dtm)
         ospin = jspin/(k2*r*r*(mt-mc)+k3*rc*rc*mc)
         epoch(i) = tphys - aj
*
* Force new NS or BH to have a one second period. 
* 
         if(kw.ne.kwold.and.(kw.eq.13.or.kw.eq.14))then
            ospin = 2.0d+08
            jspin = k3*rc*rc*mc*ospin
         endif
*
* Check exit conditions. 
*
         if(tphys.ge.(tphysf-tiny)) goto 90
         if(kw.eq.15) goto 90
         if(mt.le.(zmass(i)-dmmax+tiny)) goto 90
         dr = ABS(r-radius(i)) + tiny
         if(dr.ge.drmax) goto 90
*
 10   continue
*
 90   continue
*
* Store stellar quantities in common arrays. 
*
      age(i) = tphys
      ms_lifetime(i) = tm
      standard_timestep(i) = dtm0
      zmass0(i) = mass
      zmass(i) = mt
      zmassc(i) = mc
      radius(i) = r
      radc(i) = rc
      zlum(i) = lum
      spin(i) = ospin
      kstar(i) = kw
*
* Calculate current Y and Z. 
* [SSE does not follow chemical evolution at present so this is simply 
*  a fudge for now ... the same fudge as in triptych0.1]
*
      CALL calc_y(i,tscls(2)+tscls(3),tscls(14))
      CALL calc_z(i,tscls(2)+tscls(3),tscls(14))
*
* Recommend values for next iteration. 
*
      dtmax = dtm0
      dmmax = 0.1d0*mt
      drmax = r
      EvolveStar = tphys
*
      RETURN
      END
***
      SUBROUTINE calc_y(i,tendrg,tend)
      include "modest_common.h"
      include "modest_star.h"
      integer i
      real*8 tendrg,tend
      real*8 aj,q,q1,q2,dt1,dt2
*
      if(kstar(i).ge.7)then
         ycurr(i) = MIN(yinit(i)+0.8d0,0.6d0-zinit(i))
      else
         aj = age(i) - epoch(i)
         if(kstar(i).le.1)then
            dt1 = aj
            dt2 = ms_lifetime(i)
            q1 = yinit(i)
            q2 = MIN(yinit(i)+0.1d0,1.d0-zinit(i))
            q = q1 + (dt1/dt2)*(q2-q1)
         elseif(aj.le.tendrg)then
            dt1 = aj - ms_lifetime(i)
            dt2 = tendrg - ms_lifetime(i)
            q1 = MIN(yinit(i)+0.1d0,1.d0-zinit(i))
            q2 = MIN(yinit(i)+0.4d0,0.8d0-zinit(i))
            q = q1 + (dt1/dt2)*(q2-q1)
         else
            dt1 = aj - tendrg
            dt2 = tend - tendrg
            q1 = MIN(yinit(i)+0.4d0,0.8d0-zinit(i))
            q2 = MIN(yinit(i)+0.8d0,0.6d0-zinit(i))
            q = q1 + (dt1/dt2)*(q2-q1)
         endif
         ycurr(i) = q
      endif
*
      RETURN
      END
***
      SUBROUTINE calc_z(i,tendrg,tend)
      include "modest_common.h"
      include "modest_star.h"
      integer i
      real*8 tendrg,tend
      real*8 aj,q,q1,q2,dt1,dt2
*
      if(kstar(i).ge.7)then
         zcurr(i) = MIN(1.d0,zinit(i)+0.4d0)
      elseif(kstar(i).le.1)then
         zcurr(i) = zinit(i)
      else
         aj = age(i) - epoch(i)
         if(aj.le.tendrg)then
            dt1 = aj - ms_lifetime(i)
            dt2 = tendrg - ms_lifetime(i)
            q1 = zinit(i)
            q2 = MIN(1.d0,zinit(i)+0.2d0)
            q = q1 + (dt1/dt2)*(q2-q1)
         else
            dt1 = aj - tendrg
            dt2 = tend - tendrg
            q1 = MIN(1.d0,zinit(i)+0.2d0)
            q2 = MIN(1.d0,zinit(i)+0.4d0)
            q = q1 + (dt1/dt2)*(q2-q1)
         endif
         zcurr(i) = MAX(zcurr(i),q)
      endif
*
      RETURN
      END
*
************************************************************************
