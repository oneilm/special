PROJECT=int2

OBJSUF=o
FC=gfortran -w -c
FLINK=gfortran -w -o $(PROJECT)

.PHONY: $(PROJECT) clean list

.f.$(OBJSUF):
		$(FC) $<

vpath %.f .:../utils

FSRCS   =   qlegefuns_dr.f qlegefuns.f  prini.f legeexps.f

OBJS    =  $(FSRCS:.f=.$(OBJSUF)) 

$(PROJECT): $(OBJS)
			rm -f $(PROJECT)
			$(FLINK) $(OBJS)
			./$(PROJECT)

clean:
		rm -f $(OBJS)
		rm -f $(PROJECT)
