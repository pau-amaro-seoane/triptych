c       Copyright (c) 2001
c       by James Lombardi, Vassar College, Poughkeepsie, NY
c       and Jessica Sawyer Warren, Rutgers University, NJ.
c       All rights reserved.
c
c       Redistribution and use in source and binary forms are permitted
c       provided that the above copyright notice and this paragraph are
c       duplicated in all such forms and that any documentation,
c       advertising materials, and other materials related to such
c       distribution and use acknowledge that the software was developed
c       by the authors named above.
c
c       THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
c       IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
c       WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
c
*************************************************************************
* Make Me A Star:
* This is a parameter file params.h
* version 1.2
* November 1, 2001
* James Lombardi (lombardi@vassar.edu) and Jessica Sawyer Warren
*************************************************************************
* The variables const1 through const10 correspond to the constants c_{1}
* through c_{10} that appear in Lombardi, Sawyer Warren, Rasio, Sills and
* and Warren (2002).  This paper, along with the latest version of this
* software package, is available for download from
* http://vassun.vassar.edu/~lombardi/mmas/
      double precision const1,const2,const3,const4,const5,
     &   const6,const7,const8,const9,const10,const11
      parameter(const1=0.157d0,const2=1.8d0,const3=-1.1d0,
     &	const4=0.5d0,const5=5.d0,const6=2.5d0,const7=2.5d0,
     &  const8=5.d0,const9=2.d0,const10=.6d0)

*****Other useful parameters follow....
* N= # of parent stars
* el= # of chemical elements to be treated
* NROWS= large number bigger than # of rows in any input or output data file
* NP= # of desired rows in output data file that models the merger remnant
* Msun= mass of sun in cgs units
* Rsun= radius of sun in cgs units
*****ALL VALUES SHOULD BE MANIPULATED WITH CGS UNITS!!!!!*****
      INTEGER NROWS,N,elmax,el
      PARAMETER(NROWS=4000, N=2, elmax=14)
      DOUBLE PRECISION pi,G,Msun,Rsun
      PARAMETER(Msun=1.989d33, Rsun=6.9598d10,
     &  pi=3.1415926536d0, G=6.67259d-8)

*************************************************************************
* The following two variables control whether files are generated that are
* useful for studying mass loss and shock heating distribution.  
      logical studymassloss,makeshockfiles
      PARAMETER(studymassloss=.false.,makeshockfiles=.false.)

