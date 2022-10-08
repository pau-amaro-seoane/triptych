      program triptych
C     version 0.3
C     primary code developer: James Lombardi (Vassar)
C     1/9/03
      implicit none
      integer NROWS,Jr,C,verbosity,elmax,NP,el,N
      parameter(NROWS=4000,elmax=14,N=2)
      DOUBLE PRECISION Mp(NROWS,N),Rp(NROWS,N),Pp(NROWS,N),
     &     Dp(NROWS,N), ELp(NROWS,elmax,N)
      integer numShells(N)
      double precision Mproduct(NROWS),Pproduct(NROWS),
     &     Dproduct(NROWS),Rproduct(NROWS),ELproduct(elmax,NROWS),
     &     jproduct(NROWS)
      double precision periastron
      CHARACTER*70 parentfile(N), outfile(3)
      double precision mtest,X_init,Y_init,Z_init, m_init
      double precision vinfinity,separation0
      integer numproducts,dynamicsdriver

      verbosity=1

      if(verbosity.ge.1) then
         write(6,*)'--------------------'
         write(6,*)'WELCOME TO TRIPTYCH!'
         write(6,*)'--------------------'
         write(6,*)
      endif

      if(verbosity.ge.2) write(6,*)'verbosity=',verbosity
      
***********************************************************************
* Get information from user about what stars she would like to collide:
***********************************************************************
      if (verbosity.ge.1) then
       write(6,*)'INITIAL CONDITIONS'
       write(6,*)'------------------'
      endif
      if (verbosity.ge.2) then
       write(6,*)'First we need to choose which two stars to collide.'
       write(6,*)'Included with this software package are three'
       write(6,*)'sample data files m04.mdl, m06.mdl and m08.mdl that'
       write(6,*)'model 0.4 solar mass, 0.6 solar mass and 0.8 solar'
       write(6,*)'mass population II main sequence stars, respectively.'
       write(6,*)
      endif
      write(6,*)'Enter data file name for the first parent star:'
      if(verbosity.ge.2) then
       write(6,*)'(Some systems, including SGIs, require that you type'
       write(6,*)'single quotes around the filename.)'
      endif
      read *, parentfile(1)
      write(6,'(1X,a71)') parentfile(1)
      write(6,*) 'Enter data file name for the second parent star:'
      read *, parentfile(2)
      write(6,'(1X,a71)') parentfile(2)
C      write(6,*)'parentfiles=',parentfile

      el=2
      if(verbosity.ge.2) then
         write(6,*)
         write(6,*)'Number of chemical profiles being followed=',el
      endif

      call readparentfiles(2,parentfile,verbosity,el,
     $     Mp,Rp,Pp,Dp,ELp,numShells)

      write(6,*)
      write(6,*)'Enter the relative velocity at infinity, in km/s:'
      read *, vinfinity
      write(6,'(1X,a,g15.7)')'vinfinity=',vinfinity
C convert to cm/s:
      vinfinity=vinfinity*1d5

      write(6,*)
      write(6,*)'Enter the periastron separation, normalized to the sum'
      write(6,*)'of the parent star radii (for a headon collision, use'
      write(6,*)'zero here):'
      read *, periastron
      write(6,'(1X,a,g15.7)') 'periastron/(R_1+R_2)=',periastron

      write(6,*) 
      write(6,*)'Enter the initial separation, normalized to the sum'
      write(6,*)'of the parent star radii:'
      read *, separation0
      write(6,'(1X,a,g15.7)') 'separation0/(R_1+R_2)=',separation0

      write(6,*)
      write(6,*) 'Enter output data file names for the dynamics,'
      write(6,*) 'hydrodynamics and evolution respectively.  (One'
      write(6,*) 'file name per line.)'
      read *, outfile(1)
      read *, outfile(2)
      read *, outfile(3)
      if(verbosity.ge.2) write(6,'(1X,a,3a71)')
     $     'outfiles=',(outfile(C),C=1,3)

***********************************************************************
* Now that we know what initial conditions to consider, make the call
* to the dynamics routine:
***********************************************************************
      write(6,*)
      write(6,*) 'EPISODE I: STELLAR DYNAMICS'
      write(6,*) '---------------------------'
      numproducts=dynamicsdriver(Mp(numShells(1),1),Rp(numShells(1),1),
     $                    Mp(numShells(2),2),Rp(numShells(2),2),
     $     vinfinity,periastron,separation0,outfile(1))

      if(numproducts.ne.1) then
         write(6,*) 'numproducts=',numproducts
         stop
      endif

      write(6,*)
      write(6,*) 'EPISODE II: STELLAR HYDRODYNAMICS'
      write(6,*) '---------------------------------'
      NP=0
c      write(6,*)'NP=',NP
      call makemeastar(periastron,
     $     Mp(1,1),Rp(1,1),Pp(1,1),Dp(1,1),ELp(1,1,1),
     $     Mp(1,2),Rp(1,2),Pp(1,2),Dp(1,2),ELp(1,1,2),
     $     numShells(1),numShells(2),
     $     verbosity,
     &     el,NP,Mproduct,Pproduct,Dproduct,Rproduct,
     &     ELproduct,jproduct)
      if(verbosity.ge.2) then
         write(6,*) 'Back from the makemeastar subroutine...'
         write(6,*)
         write(6,*) 'NP has been set to',NP
         write(6,*)
      endif
       

***********************************************************************
* Create a data file for the collision product model:
***********************************************************************
      if(verbosity.ge.1) then
         write(6,'(1X,a,g12.4,a)')
     $             'Collision product radius (before thermal relaxation)
     $=', Rproduct(NP)/6.9598d10,'R_sun'
         write(6,*)'... WRITING COLLISION PRODUCT MODEL TO FILE'
      endif
      OPEN(50,FILE=outfile(2))
      write(50,*)'# 1:m [M_sun] 2:P [dyne/cm^2] 3:rho [g/cm^3]'
     $           '4:r [R_sun] 5:X 6:Z 7:j [cm^2/s]' 
      do Jr=1,NP
         if(Mproduct(Jr).gt.0.95d0*Mproduct(NP).or.mod(Jr,5).eq.1)
     $        WRITE(50,54) Mproduct(Jr)/1.989d33,Pproduct(Jr),
     $        Dproduct(Jr),Rproduct(Jr)/6.9598d10,
     &        (real(ELproduct(C,Jr)),C=1,el),jproduct(Jr)
      enddo
      CLOSE(50)
 54   FORMAT(E14.8,3(1X,E14.8),2F8.5,99(1X,E10.4))

      mtest=0.d0
      X_init=0.5d0*(Mproduct(2)+Mproduct(1))*Elproduct(1,1)
      Z_init=0.5d0*(Mproduct(2)+Mproduct(1))*Elproduct(2,1)
      mtest=mtest+0.5d0*(Mproduct(2)+Mproduct(1))
      Do Jr=2,NP-1
         X_init=X_init+
     $        Elproduct(1,Jr)*0.5d0*(Mproduct(Jr+1)-Mproduct(Jr-1))
         Z_init=Z_init+
     $        Elproduct(2,Jr)*0.5d0*(Mproduct(Jr+1)-Mproduct(Jr-1))
         mtest=mtest+0.5d0*(Mproduct(Jr+1)-Mproduct(Jr-1))
      enddo
      X_init=X_init+
     $     Elproduct(1,NP)*0.5d0*(Mproduct(NP)-Mproduct(NP-1))
      Z_init=Z_init+
     $     Elproduct(2,NP)*0.5d0*(Mproduct(NP)-Mproduct(NP-1))
      mtest=mtest+0.5d0*(Mproduct(NP)-Mproduct(NP-1))
      X_init=X_init/mtest
      Z_init=Z_init/mtest
      Y_init=1.d0-X_init-Z_init
      if(abs(mtest/Mproduct(NP)-1.d0).gt.1.d-5) then
         write(6,*)'average chemical abundances are likely wrong...'
         stop
      endif

      m_init=Mproduct(NP)/1.989d33
      if(verbosity.ge.2) then
         write(6,*) 'Collision Product mass=',m_init
         write(6,*) 'Average Y=',Y_init
         write(6,*) 'Average Z=',Z_init
      endif

      write(6,*)
      write(6,*) 'EPISODE III: STELLAR EVOLUTION'
      write(6,*) '------------------------------'
      call one_star(m_init,Rproduct(NP)/6.9598d10,
     $     Y_init,Z_init,outfile(3))

      end

**********************************************************************
**********************************************************************
      SUBROUTINE readparentfiles(N,parentfile,verbosity,el,
     $     Mp,Rp,Pp,Dp,ELp,numShells)
      implicit none
      integer K,N,verbosity,NROWS,JJ,I,el,numShells(N),C,elmax
      CHARACTER*70 parentfile(N)
      parameter(NROWS=4000,elmax=14)
      DOUBLE PRECISION Mp(NROWS,N),Rp(NROWS,N),Pp(NROWS,N),
     &     Dp(NROWS,N), ELp(NROWS,elmax,N)
**********************************************************************
* Input:
*     parentfile = ARRAY of character strings for input file names, so ...
*        parentfile(1)=Name of data file for more massive parent star (star 1)
*        parentfile(2)=Name of data file for less massive parent star (star 2)
**********************************************************************

* read in the input files for each parent star K:
      DO K=1,N
         if(verbosity.ge.2) then
            WRITE(6,*)
            WRITE(6,*) 'Input: ',parentfile(K)
         endif
         
         OPEN(26,FILE=parentfile(K))
         
* I refers to rows in the input data file.  JJ refers to rows where the
* enclosed mass M is increasing.   We will throw away any rows where this does
* not occur.  (Some input files may have truncated the value for the mass at
* a small enough number of digits that two rows appear to have the same
* enclosed mass M.)
         JJ=0
         DO I=1,NROWS
            READ(26,*,END=30) Mp(I,K),Rp(I,K),Pp(I,K),
     &           Dp(I,K),(ELp(I,C,K), C=1,el)
            IF(Mp(I,K).GT.Mp(JJ,K) .OR. JJ.EQ.0) THEN
               JJ=JJ+1
               Mp(JJ,K)=Mp(I,K)
               Pp(JJ,K)=Pp(I,K)
               Dp(JJ,K)=Dp(I,K)
               DO C=1,el
                  ELp(JJ,C,K)=ELp(I,C,K)
               ENDDO
               Rp(JJ,K)=Rp(I,K)
            ENDIF
         ENDDO
            
 30      numShells(K)=JJ
         CLOSE(26)

      ENDDO

      return
      end
