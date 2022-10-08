C------------------------------------------------------------
C        MODEST_EVOLUTION.F
C
C       Piet Hut and Jun Makino
C       Version 0.0 July 1, 2001
C
C       THIS IS STELLAR EVOLUTION MODULE #1 (stored in sev1)
C       [comprised of modest_evolution.f, modest_star.h]
C       Minor alterations to original by Jarrod Hurley
C       Version 0.1 October 23, 2002
C
C------------------------------------------------------------
C
C  Each physical quantity Q that describes the state of a star
C  can be encapsulated as an object of class `stellar_quantity',
C  where Q can stand for mass, radius, etc.  Such an object
C  contains the full information about the evolution of Q over
C  the life time of a star.  The value of a quantity Q at time t
C  can be obtained by invoking its member function Q.at(t) .
C 
C  In our current toy-model implementation, Q(t) is described by
C  a piece-wise linear function with one discontinuity, specified
C  by the following six numbers:
C 
C    t_endms   time at which the star leaves the main sequence
C    t_endrg   time at which the star leaves the red giant branch.
C    f_init    value of Q at birth of star (at time t=0)
C    f_endms   value of Q at time t=t_endms
C    f_endrg   value of Q just before time t=t_endrg
C    f_remnant value of Q just after time t=t_endrg
C 
C                                       ....... f_endrg
C                                      /|
C                                     / |
C                                    /  |
C                                   /   |
C   ^                              /    |
C   |                             /     |
C   Q                            /      |
C                          _____/.......|...... f_endms
C                _____-----             |
C      _____-----.......................|...... f_init
C                                       |------ f_remnant
C      +-----------------------+--------+------ 0
C      0                       t_endms  t_endrg
C                 t -->
C 
C  Note that we stick-figure version of stellar evolution mimics
C  only low-mass stars, while leaving out completely the horizontal
C  branch part of the evolution.  In addition, we will make the
C  following simplified assumption: t_endrg = 1.1 t_endms.
C  With no need to specify t_endrg separately, we thus are left
C  with only five independent values.
C
        
      function interpolate(t0,t1,t,f0,f1)
      real*8 interpolate,t0,t1,t,f0,f1
      interpolate= f0+(f1-f0)*(t-t0)/(t1-t0)
      end

      subroutine setup_stellar_quantity(idstar,idfunc,
     $     mstime,f0,f1,f2,f3)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar,idfunc
      real*8 mstime,f0,f1,f2,f3
      ttable(1,idfunc,idstar)= mstime
      ttable(2,idfunc,idstar)= mstime*1.1d0
      ftable(1,idfunc,idstar)= f0
      ftable(2,idfunc,idstar)= f1
      ftable(3,idfunc,idstar)= f2
      ftable(4,idfunc,idstar)= f3
      end

      function stellar_quantity(idstar,idfunc,t)
      real*8 stellar_quantity
      include "modest_common.h"
      include "modest_star.h"
      integer idstar,idfunc
      real*8 t,value,interpolate
      if (t .lt. ttable(1,idfunc,idstar)) then
         value = interpolate(0d0, ttable(1,idfunc,idstar),
     $        t, ftable(1,idfunc,idstar), ftable(2,idfunc,idstar))
      elseif(t .lt. ttable(2,idfunc,idstar)) then
         value = interpolate(ttable(1,idfunc,idstar),
     $        ttable(2,idfunc,idstar),
     $        t, ftable(2,idfunc,idstar), ftable(3,idfunc,idstar))
      else
         value = ftable(4,idfunc,idstar)
      endif
      stellar_quantity = value
      end

      
      block data block
      include "modest_common.h"
      data is_used/nmax*0/
      end

      function first_unused_star()
      include "modest_common.h"
      integer i
      do i=1,nmax
         if (is_used(i) .eq. 0) then
            first_unused_star = i
            return
         endif
      enddo
      first_unused_star = -1
      end
      

      function CreateStar(m_init, Y_init, Z_init)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      real*8 m,r,z0,z1,z2,z3,y0,y1,y2,y3,mstime
      idstar = first_unused_star()
      if (idstar .lt. 0) then
         CreateStar= -1
         return
      endif
      CreateStar=idstar
      m = m_init
      mstime =1d4 /(m*m*m)
      ms_lifetime(idstar) = mstime
      standard_timestep(idstar) = mstime*1d-5
      call setup_stellar_quantity(idstar,1,
     $     mstime,m,m, m*0.5D0,m*0.25d0)
      r = m
      call setup_stellar_quantity(idstar,2,mstime,
     $     r,r,r*100d0,0.01d0/m)
      z0 = Z_init
      z1 = Z_init
      z2 = min(1.0d0,z0+0.2d0)
      z3 = min(1.0d0,z0+0.4d0)
      call setup_stellar_quantity(idstar,4,mstime,z0,z1,z2,z3)
      y0 = Y_init
      y1 = min(y0+0.1d0,1d0-z1)
      y2 = min(y0+0.4d0,1d0-z2)
      y3 = min(y0+0.8d0,1d0-z3)
      call setup_stellar_quantity(idstar,3,mstime,y0,y1,y2,y3)
      age(idstar) = 0
      discontinuity_flag(idstar) = 0
      end
      
      function EvolveStar(idstar,dtmax,dMmax, dRmax, dYmax, dZmax)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      real*8 old_age,last_age,max_age
      real*8 m0,r0,y0,z0
      integer dmax_exceeded
      real*8 t0,t1,newt
      real*8 stellar_quantity
      if (discontinuity_flag(idstar) .eq. -1)then
         discontinuity_flag(idstar) = 1
         age(idstar) = age(idstar)
     $        +2*discon_dt*ms_lifetime(idstar)
         EvolveStar = age(idstar)
         return
      endif
      old_age = age(idstar)
      if (discontinuity_flag(idstar) .eq. 1) then
         discontinuity_flag(idstar) =0
      endif
      m0 = stellar_quantity(idstar,1,age)
      r0 = stellar_quantity(idstar,2,age)
      y0 = stellar_quantity(idstar,3,age)
      z0 = stellar_quantity(idstar,4,age)
      last_age = age(idstar)
      max_age = old_age + dtmax
 100  if (dmax_exceeded(idstar,age(idstar),old_age,
     $     dMmax,dRmax,dYmax,dZmax)
     $     .eq. 0 .and. age(idstar) .lt. max_age) then
         last_age = age(idstar)
         age(idstar) = age(idstar)+standard_timestep(idstar)
         if (age(idstar) .gt. max_age ) age(idstar) = max_age
         goto 100
      endif
      if (dmax_exceeded(idstar,age(idstar),old_age,
     $     dMmax,dRmax,dYmax,dZmax)
     $     .ne. 0) then
         t0 = last_age
         t1 = age(idstar)
 200     if ((t1-t0) .gt. discon_dt*ms_lifetime(idstar))then
            newt = (t1+t0)/2.0d0
            if(dmax_exceeded(idstar,newt,old_age,
     $           dMmax,dRmax,dYmax,dZmax) .ne. 0) then
               t1 = newt
            else
               t0 = newt
            endif
            goto 200
         endif
         if(dmax_exceeded(idstar,t1,t0,
     $        dMmax,dRmax,dYmax,dZmax) .ne. 0) then
            discontinuity_flag(idstar) = -1
         endif
         age(idstar) = t0
      endif
      EvolveStar = age(idstar)
      end
      
      function dmax_exceeded(idstar, new_age,old_age,
     $     dMmax,dRmax,dYmax,dZmax)
      integer dmax_exceeded,idstar
      real*8 stellar_quantity
      real*8 new_age,old_age,dMmax,dRmax,dYmax,dZmax
      if (abs(stellar_quantity(idstar,1,new_age)
     $     -stellar_quantity(idstar,1,old_age)) .ge. dMmax .or.
     $     abs(stellar_quantity(idstar,2,new_age)
     $     -stellar_quantity(idstar,2,old_age)) .ge. dRmax .or.
     $     abs(stellar_quantity(idstar,3,new_age)
     $     -stellar_quantity(idstar,3,old_age)) .ge. dYmax .or.
     $     abs(stellar_quantity(idstar,4,new_age)
     $     -stellar_quantity(idstar,4,old_age)) .ge. dZmax)then
          dmax_exceeded = 1
       else
          dmax_exceeded = 0
       endif
       end
      function get_discontinuity_flag(idstar)
      integer get_discontinuity_flag
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      get_discontinuity_flag = discontinuity_flag(idstar)
      end

      function finished_star(idstar)
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
      end

      function getMass(idstar)
      include "modest_common.h"
      integer idstar
      real*8 stellar_quantity
      getMass = stellar_quantity(idstar,1,age(idstar))
      end
      
      function getRadius(idstar)
      include "modest_common.h"
      integer idstar
      real*8 stellar_quantity
      getRadius = stellar_quantity(idstar,2,age(idstar))
      end
      
      function getTime(idstar)
      include "modest_common.h"
      integer idstar
      getTime = age(idstar)
      end
      
      function getY(idstar)
      include "modest_common.h"
      integer idstar
      real*8 stellar_quantity
      getY = stellar_quantity(idstar,3,age(idstar))
      end
      
      
      function getZ(idstar)
      include "modest_common.h"
      integer idstar
      real*8 stellar_quantity
      getZ = stellar_quantity(idstar,4,age(idstar))
      end

      subroutine print_star(idstar,iunit)
      integer idstar,iunit
      include "modest_common.h"
      include "modest_star.h"
      write(iunit,600)getTime(idstar),getMass(idstar),getRadius(idstar),
     $     getY(idstar),getZ(idstar)
 600  format(3g15.8,2g15.7)
      end

      subroutine one_star(m_init,Y_init,Z_init,outfile)
      include "modest_common.h"
      include "modest_star.h"
      integer idstar
      integer discon_reached,countdown,get_discontinuity_flag
      real*8 t,new_time
      character*512 outfile
      idstar = CreateStar(m_init,Y_init,Z_init)
      dtmax = 1d0
      discon_reached = 0
      t = 0
      countdown = 3
      iunit=50
      if(iunit.ne.6) then
         open(iunit,file=outfile)
         write(6,*) '... WRITING EVOLUTION DATA TO FILE'
      endif
      write(iunit,*) 'Age [Myr], M [M_sun], R [R_sun], Y, Z'
 100  if (countdown .gt. 0) then
         dMmax = 0.1d0*getMass(idstar)
         dRmax = 1d0*getRadius(idstar)
         dYmax = 0.1d0
         dZmax = 0.1d0
         new_time = EvolveStar(idstar,dtmax,dMmax,dRmax,dYmax,dZmax)
         call print_star(idstar,iunit)
         if ((new_time -t).gt. dtmax - 1d-10) dtmax = dtmax*2
         t = new_time
         if (get_discontinuity_flag(idstar) .ne. 0) discon_reached = 1
         if (discon_reached .ne. 0 ) then
            countdown  = countdown - 1
         endif
         goto 100
      endif
      if(iunit.ne.6) close(iunit)
      i = finished_star(idstar)
      end
