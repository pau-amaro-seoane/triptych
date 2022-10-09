C Version 0.3
C 1/9/03
      integer nmax
      parameter(nmax=100000)
      integer is_used(nmax)
      integer CreateStar,first_unused_star,finished_star
      real*8 age(nmax)
      real*8 moverall,roverall,Yoverall,Zoverall
      real*8 dtmax,dMmax,dRmax,dYmax,dZmax
      real*8 getTime,getMass,getRadius,getY,getZ
      real*8 getMStime,getLum,getMc,getRc,getSpin,getType
      real*8 getTeff
      real*8 EvolveStar
C
      common /modest/ age,is_used
