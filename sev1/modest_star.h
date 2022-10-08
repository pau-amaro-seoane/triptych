C------------------------------------------------------------
C       modest_star.h  
C
C       Common for stellar evolution related variables.
C
C       Piet Hut and Jun Makino
C       Version 0.0 July 1, 2001
C------------------------------------------------------------
C      
      integer nfmax
      parameter (nfmax = 4)
      integer discontinuity_flag(nmax)
C
      real*8 discon_dt
      parameter (discon_dt = 1d-8)
      real*8 ftable(4,nfmax,nmax),ttable(2,nfmax,nmax)
      real*8 ms_lifetime(nmax)
      real*8 standard_timestep(nmax)
      real*8 discontinuity_time(nmax)

      common /stellar/ ftable,ttable,ms_lifetime,standard_timestep,
     $   discontinuity_time,discontinuity_flag
C
C------------------------------------------------------------
