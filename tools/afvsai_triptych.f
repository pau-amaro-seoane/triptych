      program afvsai_triptych
c     This code assumes that triptych collided two identical stars and calls
c     the hydro output file product.dat.
c
c     This code reads in product.dat line by line. On each line, it
c     (1) reads in the values of m, P, and rho
c     (2) calculates A=P/rho^(5/3).
c     (3) finds the two lines in the .mdl file that have mini values that straddle m/2, the first of these
c     two lines gives "last" quantities in the code below
c     (4) interpolates between these lines to get Piniuse and ainiuse
c     (5) prints to a file a line of data that contains log10(A-Aini) and log10(Pini) for plotting later.

c     We expect to find a slope of -1.1 on a log-log plot, since that is what
c     triptych/mmas is using under the hood.

      implicit none
      integer i
      real*8 m,p,rho
      real*8 mini,pini,rhoini,dummy,aini
c     Determine conversion factors from code units to cgs
      real*8 munit
      parameter(munit=1.9891d33) ! number of g in unit of mass.  1.9891d33 is for solar mass
      real*8 minilast,pinilast,ainilast,ainiuse,piniuse

      open(13,file='product.dat')
c     The first line of product.dat is a header line that says enclosed mass m
c     is in the first column, pressure P is in the second column, and density
c     rho is in the third column.
      read(13,*)

      open(14,file='m08.mdl')
c     The triptych README file says that in the
c     *.mdl files these quantities (call them mini, Pini and rhoini) are in
c     the 1st, 3rd, and 4th columns, respectively.
      read(14,*) mini,dummy,pini,rhoini
      aini=pini/rhoini**(5d0/3d0)
      minilast=mini
      ainilast=aini
      pinilast=pini

      open(15,file='loglog.dat')
      write(15,*)'log(A-A_init)         log(P_init) [both in cgs units]'
      do i=1,9999
         read(13,*,end=99) m,p,rho
c     Convert enclosed mass to cgs units
         m=m*munit
        
c         print *, m, p, rho
         do while (mini.lt.m/2)
            minilast=mini
            ainilast=aini
            pinilast=pini
            read(14,*) mini,dummy,pini,rhoini
            aini=pini/rhoini**(5d0/3d0)
         enddo
         ainiuse=((mini-m/2)*ainilast+(m/2-minilast)*aini)/
     $        (mini-minilast)
         piniuse=((mini-m/2)*pinilast+(m/2-minilast)*pini)/
     $        (mini-minilast)
         write(15,*) log10(p/rho**(5d0/3d0)-ainiuse),log10(piniuse)
      enddo
 99   continue
      close(13)
      close(14)
      close(15)
      end
