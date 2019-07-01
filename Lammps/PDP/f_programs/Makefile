CDIR=.
COMPILER = ifort
COMPILER_MPI = mpif90 #mpiifort
CFLAGS   = -O3 -vec-report
#CFLAGS   = -check
LDFLAGS =

###################################
TARGET1 = poly
OBJECT1 = poly.o
TARGET2 = dp_poly
OBJECT2 = dp_poly.o

SOURCES  = $(wildcard *.f90)
SRC_MODULES = $(filter class%, $(SOURCES))
#SRC_OBJECTS_ALL = $(filter-out class%, $(SOURCES))
#SRC_OBJECTS = $(filter-out main%, $(SRC_OBJECTS_ALL))

MODULES  = $(SRC_MODULES:.f90=.o)
#OBJECTS  = $(SRC_OBJECTS:.f90=.o)

all: $(TARGET1) $(TARGET2) 

$(TARGET1): $(MODULES) $(OBJECT1)
	$(COMPILER) -o $@ $^ $(LDFLAGS)

$(TARGET2): $(MODULES) $(OBJECT2)
	$(COMPILER) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(COMPILER_MPI) $(CFLAGS) -o $@ -c $<

clean:
	@echo Cleaning files ...
	rm -f $(OBJECTS) $(MODULES) $(TARGET) *~
