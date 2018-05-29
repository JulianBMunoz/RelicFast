IDIR = ./include
ODIR = ./objects
SRCDIR = ./source

#remove this line to cancel parallelization (it's only used in collapse.cpp)
parallel = -fopenmp

CC=gcc

CFLAGS=-I $(IDIR) -O4 -ffast-math -lstdc++ $(parallel) -lm

DEPS = auxiliar.h bias.h boltzmann_Calls.h collapse.h common.h pressure.h
OBJ = auxiliar.o bias.o boltzmann_Calls.o collapse.o pressure.o RelicFast.o

_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))


%.o: $(SRCDIR)/%.cpp $(IDIR)/%.h
	$(CC) -c -o $(ODIR)/$@ $< $(CFLAGS)

relicfast: $(OBJ)
	$(CC) -o $@ $(_OBJ) $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
