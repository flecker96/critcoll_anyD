CF = gfortran
FLAGS = -O -Wall -mcmodel=large
.f.o:
	$(CF) -c $(FLAGS) -o $*.o $*.f 

.KEEP_STATE:

FOUR = fourier/four1.o fourier/fourdiff2.o fourier/inhom.o fourier/lamint.o \
	fourier/halve.o fourier/double.o fourier/mik.o fourier/fourinterpol.o \
	fourier/four_to_three.o fourier/three_to_four.o fourier/fourdiff1.o \
	fourier/fourdiff1smoothe.o

CRITINNER =  inner/shootlin.o \
	$(FOUR) solve/lubksb.o solve/ludcmp.o solve/mprove.o \
	inner/derivslin.o inner/yfromfields.o inner/fieldsfromy.o \
	inner/leftdatalin.o inner/rightdatalin.o inner/detofgam.o\
	shoot/pack2.o shoot/unpack2.o inner/fieldsfromylin.o \
	shoot/implicit_steplin.o inner/critmain_inner.o


shootlin_inner.exe: $(CRITINNER)
	$(CF) $(CRITINNER) -o $@ $(FLAGS)

clean:
	rm fort.* *~ *.exe *.out *.junk */*.o */*~ *~ *.o


