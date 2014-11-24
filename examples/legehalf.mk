
HOST = osx-gfortran
##HOST = osx-intel
PROJECT = int2

ifeq ($(HOST),osx-gfortran)
  FC = gfortran -c -w
  FFLAGS = -O3
  FLINK = gfortran -w -o $(PROJECT) -framework accelerate
endif

ifeq ($(HOST),osx-intel)
  FC = ifort -c -w
  FFLAGS = -O3
  FLINK = ifort -w -mkl -o $(PROJECT)
endif


.PHONY: all clean list


SOURCES = legehalf_dr.f ../src/legehalf.f  ../../utils/prini.f

OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(PROJECT)
	$(FLINK) $(OBJECTS)
	./$(PROJECT)

clean:
	rm -f $(OBJECTS)
	rm -f $(PROJECT)

list: $(SOURCES)
	$(warning Requires:  $^)



