NAME = triptych
# Stellar Evolution Module 1 (cp modest_star.h modest_evolution.f from sev1)
#LFLAGS = 
#LFLGS = 
# Stellar Evolution Module 2 (cp modest_star.h modest_evolution.f from sev2)
#LFLAGS = sev2/libstrifc.a
#LFLGS = -Lsev2 -lstrifc

#FC= gfortran -w # If on linux
FC= egfortran -w # If on OpenBSD
FFILES = triptych.f modest_dynamics.f makemeastar.f modest_evolution.f recipes.f zcnsts.f zfuncs.f star.f hrdiag.f mlwind.f mrenv.f

OBJFILES = $(FFILES:.f=.o)

$(NAME):$(OBJFILES) $(LFLAGS)
	$(FC) -o ../bin/$(NAME) -O $(OBJFILES) $(LFLGS)

triptych.o : triptych.f
modest_dynamics.o : modest_dynamics.f dynamics.h
makemeastar.o : makemeastar.f params.h
modest_evolution.o : modest_evolution.f modest_common.h modest_star.h
recipes.o : recipes.f
