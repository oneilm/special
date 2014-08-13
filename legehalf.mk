PROJECT=int2

OBJSUF=o
FC=f90  -w -c
FLINK=f90  -w -o $(PROJECT)

.PHONY: $(PROJECT) clean list

.f.$(OBJSUF):
		$(FC) $<

vpath %.f .:../bin

FSRCS   =   legehalf_dr.f legehalf.f  prini.f

OBJS    =  $(FSRCS:.f=.$(OBJSUF)) 

$(PROJECT): $(OBJS)
			rm -f $(PROJECT)
			$(FLINK) $(OBJS)
			./$(PROJECT)

clean:
		rm -f $(OBJS)
		rm -f $(PROJECT)
