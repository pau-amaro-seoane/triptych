      implicit none

      double precision g,m1,m2,ltot,tstart,tend,r1,r2,tendtilde
      double precision acc, dxsav, rmin, theta0, r
      double precision alpha,mu,e0,e,r0,ltotint,
     &   etot
      integer nsteps, NMAX, KMAXX, MAXSTP, kmax, kount
      double precision xp,yp
      double precision pi, INF
      integer outhandle
      common/files/ outhandle
      
      logical dump

C     acc--the accuracy for the Runge-Kutta integration, 
      
      parameter (tstart = 0.d0, tendtilde = 40.d0)
      parameter (acc = 0.0000001d0, nsteps = 10000)
C      parameter (acc = 0.000001d0, nsteps = 1000)
      PARAMETER (NMAX=50, KMAXX=1500, MAXSTP=20000)
      parameter (pi=3.1415926536d0)
      parameter (g = 6.6725985d-08, INF = 1.d+30)

      common /params/ m1,m2,r1,r2,alpha,mu,e0,r0,
     $                rmin,e,theta0,r,ltotint,etot
      common /output/ xp(KMAXX),yp(NMAX,KMAXX)
      COMMON /path/ kmax,kount,dxsav
      common /flags/ dump

      double precision munit,runit
      parameter(munit=1.989d33,runit=6.9598d10)
      integer verbosity
C set verbosity of DYNAMICS routines here:
      parameter(verbosity=1)
      double precision xcm,ycm,vxcm,vycm
      common /com/ xcm,ycm,vxcm,vycm
