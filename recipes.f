      integer function odeint(ystart,nvar,x1,eps,h1,hmin,nok,nbad,
     $     derivs,rkqs)
C     THIS function SOLVES DIFF EQS. USING RKQS

C     DECLARING VARIABLES

      implicit none
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      DOUBLE PRECISION eps,h1,hmin,x1,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=1500,TINY=1.d-30)
      INTEGER i,kmax,kount,nstp
      DOUBLE PRECISION dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),
     &xp,y(NMAX),yp,yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav
      COMMON /output/ xp(KMAXX),yp(NMAX,KMAXX)

      double precision m1,m2,r1,r2,r
      common /params/ m1,m2,r1,r2

C     SETTING VARIABLES

      x=x1
      h=h1
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.d0*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif

        r = sqrt((y(1)-y(5))**2 + (y(2)-y(6))**2)

        if(r.le.r1+r2)then
           do 14 i=1,nvar
              ystart(i)=y(i)
 14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          odeint=1
          return
        endif
        if(abs(hnext).lt.hmin) pause
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      write(6,*) 'too many steps in odeint.'
      write(6,*) 'YOUR STARS NEVER COLLIDED!'
c      write(6,*) 'y(i)=',(y(i),i=1,nvar)
c      write(6,*) 'dydx(i)=',(dydx(i),i=1,nvar)
c      write(6,*) 't=',x
c      write(6,*)' hdid=',hdid,' hnext=',hnext
c      write(6,*)'eps=',eps
c      write(6,*)'yscal(i)=',(yscal(i),i=1,nvar)

      odeint=2
      return

      END
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      double precision h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      double precision ytemp(NMAX)
CU    USES derivs
      INTEGER i
      double precision ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),
     *ak6(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,
     *B21=.2d0,B31=3.d0/40.d0,B32=9.d0/40.d0,B41=.3d0,
     *B42=-.9d0,B43=1.2d0,B51=-11.d0/54.d0,B52=2.5d0,
     *B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,
     *B62=175.d0/512.d0,B63=575.d0/13824.d0,B64=44275.d0/110592.d0,
     *B65=253.d0/4096.d0,C1=37.d0/378.d0,
     *C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,
     *DC1=C1-2825.d0/27648.d0,DC3=C3-18575.d0/48384.d0,
     *DC4=C4-13525.d0/55296.d0,DC5=-277.d0/14336.d0,DC6=C6-.25d0)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs, rkck
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER i
      double precision errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),
     &SAFETY,PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,
     &ERRCON=1.89d-12)
      integer merged
      common /state/ merged
      
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.d0
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1d0*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x)then
           write(6,*)'stepsize underflow in rkqs (stars have merged?)'
           write(6,*) 'y(i)=',(y(i),i=1,n)
           write(6,*) 'dydx(i)=',(dydx(i),i=1,n)
           write(6,*) 't=',x
           write(6,*)'htry=',htry,' hdid=',hdid,' hnext=',hnext
           write(6,*)'eps=',eps
           write(6,*)'yerr(i)=',(yerr(i),i=1,n)
           write(6,*)'yscal(i)=',(yscal(i),i=1,n)
           merged=1
           return
        endif
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.d0*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
