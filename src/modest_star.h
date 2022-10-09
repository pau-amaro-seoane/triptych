C------------------------------------------------------------
C       modest_star.h 
C
C       Common for stellar evolution related variables. 
C
C------------------------------------------------------------
C
      integer kstar(nmax)
C
      real*8 epoch(nmax)
      real*8 ms_lifetime(nmax),standard_timestep(nmax)
      real*8 zmass0(nmax),zmass(nmax),zmassc(nmax)
      real*8 zinit(nmax),yinit(nmax)
      real*8 zcurr(nmax),ycurr(nmax)
      real*8 spin(nmax),radius(nmax),radc(nmax),zlum(nmax)
C
      common /stellar/ epoch,ms_lifetime,standard_timestep,
     &                 zmass0,zmass,zmassc,zinit,yinit,
     &                 zcurr,ycurr,spin,radius,radc,zlum,kstar
C
C------------------------------------------------------------
