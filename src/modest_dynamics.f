      integer function dynamicsdriver(mass1,rad1,mass2,rad2,vinf,
     $     peri,sep0,outfile)
C     ************************************************************
C     Driver for the stellar dynamics
C     ************************************************************

      INCLUDE 'dynamics.h'

      double precision totalEnergy
      double precision mass1,rad1,mass2,rad2,vinf,peri,sep0
      integer numproducts,dynamics
      double precision vinfinity
      common /vel/ vinfinity
      character*70 outfile,outputfile
      common /outputfilename/ outputfile

C     ************************************************************

      outputfile=outfile

      dump = .TRUE.

      m1=mass1
      r1=rad1
      m2=mass2
      r2=rad2
      vinfinity=vinf
      rmin=peri*(r1+r2)
      r0=sep0*(r1+r2)
      
      if(verbosity.ge.2)
     $     write(6,*)'About to call dynamics with vinfinity=',
     $     vinfinity,'cm/s'
      dynamicsdriver=dynamics(rmin)

      end

C     **************************************************************
      integer FUNCTION dynamics(rminl)
C     Input: periastron separation rminl in cgs units
C     Output: final total orbital energy
      INCLUDE 'dynamics.h'
      character*70 outputfile
      common /outputfilename/ outputfile

      double precision rdot, rminl
      double precision thetadot, costheta, sintheta, 
     $     k, rdot2
      double precision x1, y1, x2, y2, vx1, vy1, vx2, vy2
      double precision args(8)
      integer i, nok, nbad
      external rkqs, derivs     
      double precision theta
      double precision h1
      double precision tunit,vunit,lunit,eunit
      double precision vinfinity
      common /vel/ vinfinity
      integer collided,odeint
C     ***************************************************************

      tunit=sqrt(runit**3/g/munit)
      vunit=runit/tunit
      lunit=munit*vunit*runit
      eunit=g*munit**2/runit

C     Calculate bounds for r0
      k = g*m1*m2
      mu = m1*m2/(m1+m2)

      e0=vinfinity**2*rminl/(g*(m1+m2))+1.d0

      alpha = rminl*(1.d0+e0)
C     Equation (8.40) of Marion and Thornton
      ltot = sqrt(alpha*mu*k)
C     Calculate the other quantities from what is given
C     (This fails if e=0 because it's encountered in the denominator)

C     Equation (8.10) of Marion and Thornton
      thetadot = ltot/mu/(r0*r0)      

C     Equation (8.40) of Marion and Thornton
C      etot = (e0*e0-1.d0)*mu*k*k/(ltot*ltot)*0.5d0 
      etot=0.5d0*mu*vinfinity**2
      
C     Equation (8.41) of Marion and Thornton
      costheta = (alpha/r0-1.d0)/e0
      theta=-acos(costheta)
      sintheta = sin(theta)

C     The minus signs on the position and velocity component equations
C     have been chosen so that the separation vector r equals r2-r1,
C     *not* r1-r2 as in Marion and Thornton.  This was done so that the
C     code would give the same initial conditions as twostars.f when the
C     eccentricity is 1
      x1 = -m2/(m1+m2)*r0*costheta
      y1 = -m2/(m1+m2)*r0*sintheta
      x2 = m1/(m1+m2)*r0*costheta
      y2 = m1/(m1+m2)*r0*sintheta

C     Differentiating eqn. 8.41
      rdot2 = alpha/(1+e0*costheta)**2*e0*sintheta*thetadot
C     Another expr. for rdot (using eqn. 8.14) and compare
      rdot = -sqrt((etot - 0.5d0*ltot**2/mu/r0**2 + k/r0)/0.5d0/mu)

      if(dump) then
         outhandle=20
         open(outhandle,file=outputfile)
      endif

      vx1 = m2/(m1+m2)*(r0*sintheta*thetadot - rdot*costheta)
      vy1 = m2/(m1+m2)*(-r0*costheta*thetadot - rdot*sintheta)
      vx2 = m1/(m1+m2)*(-r0*sintheta*thetadot + rdot*costheta)
      vy2 = m1/(m1+m2)*(r0*costheta*thetadot + rdot*sintheta)

      if(verbosity.ge.2) then
         write(6,*) 'rdot(from alpha): ',rdot2,'  rdot(from E): '
     $        ,rdot
         write(6,*)
         write(6,*) 'NEW COLLISION: e0=',e0
         write(6,*) 'etot=',etot,'(cgs) =',etot/eunit
         write(6,*) 'ltot=',ltot,'(cgs) =',ltot/lunit
         write(6,*) 'alpha', alpha,'(cgs) =',alpha/runit
         write(6,*) 'k=',k
         write(6,*) 'r_min/runit: ', rminl/runit,
     $        '  r_max: ', alpha/(1.d0-e0)
         write(6,*) 'Setting r0: ', r0
         write(6,*) 'r1: ',r1/runit,' r2: ',r2/runit,
     $        ' r1+r2: ',(r1+r2)/runit
         write(6,*) 'm1: ',m1/munit,'  m2: ',m2/munit,'   mu: ',mu/munit 
         write(6,*) 'theta=',theta,'cos: ',costheta,' sin: ',sintheta
         write(6,*) 'Position of star 1: ',x1/runit,y1/runit
         write(6,*) 'Position of star 2: ',x2/runit,y2/runit
         write(6,*) 'Velocity of star 1: ',vx1/vunit,vy1/vunit
         write(6,*) 'Velocity of star 2: ',vx2/vunit,vy2/vunit
         write(6,*) 'c.o.m. position=',(m1*x1+m2*x2)/(m1+m2),
     $        (m1*y1+m2*y2)/(m1+m2)
         write(6,*) 'c.o.m. velocity=',(m1*vx1+m2*vx2)/(m1+m2),
     $        (m1*vy1+m2*vy2)/(m1+m2)
      endif

C     Runge-Kutta integration

      args(1) = x1
      args(2) = y1
      args(3) = vx1
      args(4) = vy1
      args(5) = x2
      args(6) = y2
      args(7) = vx2
      args(8) = vy2

      r = sqrt((x1-x2)**2 + (y1-y2)**2)
c      etot = (m1*(args(3)**2+args(4)**2) +
c     $     m2*(args(7)**2+args(8)**2))/2.d0 - g*m1*m2/r
c      write(6,*) 'etot(from args array): ', etot
      etot = 0.5d0*(m1*(vx1**2+vy1**2) +
     $     m2*(vx2**2+vy2**2)) - g*m1*m2/r
      ltot = m1*(x1*vy1 - y1*vx1) +
     $     m2*(x2*vy2 - y2*vx2)
      e = sqrt(1.d0 + 2.d0*etot*ltot**2/mu/k**2)
      alpha=ltot**2/(mu*k)

      if(verbosity.ge.2) then
         write(6,*) 'etot(from x,y,vx,vy): ', etot
         write(6,*) 'l,e,alpha:',ltot,e,alpha
      endif

      h1=((vx1-vx2)**2+(vy1-vy2)**2)**0.5d0/(g*mu/r**2)*acc
      dxsav = h1
      kmax = 1500

      dynamics=odeint(args,8,0.d0,acc,h1,0.d0,
     $     nok,nbad,derivs,rkqs)

C     Print results to file, also calculating the energy and ang. momentum
      if(dump) then
         write(outhandle,'(a128)')
     $ '# 1: t [hour]    2:r/(R1+R2)   3:e  4:alpha [R_sun]  5:l/lunit'
     $,' 6:E/eunit    7:x1 [R_sun]  8:y1 [R_sun]  9:x2 [R_sun]  10:y2'
     $,' [R_sun]    11:theta0'
         do i=1, kount
            call getorbitalelements(i)
            if(vxcm**2+vycm**2.lt.acc*vunit**2 .and.
     $           (xcm/runit)**2+(ycm/runit)**2.lt.acc) then
               write(outhandle,'(2e12.4,f8.4,8e12.4)')
     $              xp(i)/3600,r/(r1+r2),e,
     $              alpha/runit,ltotint/lunit,etot/eunit,
     $              yp(1,i)/runit,yp(2,i)/runit, 
     $              yp(5,i)/runit,yp(6,i)/runit, 
     $              theta0
            else
               if(verbosity.ge.2)
     $              write(6,*) 'STELLAR DYNAMICS: BAD DATA AT t=',
     $              xp(i)/tunit,vxcm**2+vycm**2,acc*vunit**2,
     $              xcm**2+ycm**2,acc*runit**2
            endif
         enddo
         if(outhandle.ne.6) then
            write(6,*)'... WRITING DYNAMICS DATA TO FILE'
            close(outhandle)
         endif
      else
         call getorbitalelements(kount)
      endif

      if(verbosity.ge.2) then
         write(6,*) 'For rminl=',rminl,
     $        ' dynamics returns etot/eunit=', etot/eunit
      endif

      return 

      end


      SUBROUTINE getorbitalelements(i)
      INCLUDE 'dynamics.h'
      integer i
      double precision x1, y1, x2, y2, vx1, vy1, vx2, vy2
      double precision k, costheta,sintheta,theta,dotproduct,
     $     cosdeltatheta,deltatheta

      k = g*m1*m2
      mu = m1*m2/(m1+m2)

      x1=yp(1,i)
      y1=yp(2,i)
      vx1=yp(3,i)
      vy1=yp(4,i)
      x2=yp(5,i)
      y2=yp(6,i)
      vx2=yp(7,i)
      vy2=yp(8,i)
      r = sqrt((x1-x2)**2 + (y1-y2)**2)
      
C     Find the center of mass position and velocity:
      xcm=(m1*x1+m2*x2)/(m1+m2)
      ycm=(m1*y1+m2*y2)/(m1+m2)
      vxcm=(m1*vx1+m2*vx2)/(m1+m2)
      vycm=(m1*vy1+m2*vy2)/(m1+m2)
      
C     Calculate total orbital energy and angular momentum (with respect to
C     the center of mass, which may be drifting slightly because of numerical
c     integration errors):
      etot = 0.5d0*(m1*((vx1-vxcm)**2+(vy1-vycm)**2) +
     $     m2*((vx2-vxcm)**2+(vy2-vycm)**2)) - g*m1*m2/r
      ltotint = m1*((x1-xcm)*(vy1-vycm) - (y1-ycm)*(vx1-vxcm)) +
     $     m2*((x2-xcm)*(vy2-vycm) - (y2-ycm)*(vx2-vxcm))
      
C     Calculate the orbital elements alpha (the semi-latus rectum)
c     and e (the eccentricity)
      alpha=ltotint**2/(mu*k)
      e = sqrt(1.d0 + 2.d0*etot*alpha/k)
      
C     Calculate the angular position theta and the orbital element theta0
C     theta0 is the angle theta at which 
C     the separation r would be a minimum were the stars to continue on their
c     current paths (but without drag forces). Some "tricks" are done to make
C     sure theta and theta0 are always between -pi and +pi.
      costheta=(x2-x1)/r
      sintheta=(y2-y1)/r
      theta=sign(acos(costheta),sintheta)
      dotproduct=(x2-x1)*(vx2-vx1)+(y2-y1)*(vy2-vy1)
      cosdeltatheta=(alpha/r-1.d0)/e
      if(dabs(cosdeltatheta).gt.1.d0) then
C         write(6,*)'cos(deltatheta)=',cosdeltatheta,'??? reset to 1'
         cosdeltatheta=sign(1.d0,cosdeltatheta)
      endif
      deltatheta=sign(acos(cosdeltatheta),dotproduct)
      theta0=theta-deltatheta
      theta0=theta0-int(theta0/pi)*2.d0*pi
      
      return 

      end


      SUBROUTINE derivs(t,y,dydt)

C     ROUTINE TO SET DIFF EQS.

      INCLUDE 'dynamics.h'

      DOUBLE PRECISION t,y(*),dydt(*)
      double precision x1,y1,x2,y2
      double precision xforce,yforce

C y(1)= x1
C y(2)= y1
C y(3)= vx1
C y(4)= vy1
C y(5)= x2
C y(6)= y2
C y(7)= vx2
C y(8)= vy2

      x1=y(1)
      y1=y(2)
      x2=y(5)
      y2=y(6)
      r = sqrt((x1-x2)**2 + (y1-y2)**2)

c     if(r.le.r1+r2) then
c            The stars are overlapping:
c     else
C            The stars aren't touching:
c     endif

      xforce=g*m1*m2*(x2-x1)/r**3
      yforce=g*m1*m2*(y2-y1)/r**3

      dydt(1)=y(3)
      dydt(2)=y(4)
      dydt(3)=xforce/m1
      dydt(4)=yforce/m1
      dydt(5)=y(7)
      dydt(6)=y(8)
      dydt(7)=-xforce/m2
      dydt(8)=-yforce/m2
      return
      END
