CF = gfortran
#CF = pgf77
FLAGS = -O -Wall
.f.o:
	$(CF) -c $(FLAGS) -o $*.o $*.f 

.KEEP_STATE:

FOUR = fourier/four1.o fourier/fourdiff2.o fourier/inhom.o fourier/lamint.o \
	fourier/halve.o fourier/double.o fourier/mik.o fourier/fourinterpol.o \
	fourier/four_to_three.o fourier/three_to_four.o fourier/fourdiff1.o \
	fourier/fourdiff1smoothe.o fourier/changeresolution.o

BG_TO_FILE = inner/bg_to_file.o inner/shoot_to_file.o \
	shoot/implicit_step.o \
	$(FOUR) \
	inner/derivs.o inner/yfromfields.o inner/fieldsfromy.o \
	inner/leftdata.o inner/rightdata.o\
	shoot/pack2.o shoot/unpack2.o


bg_to_file.exe: $(BG_TO_FILE)
	$(CF) $(BG_TO_FILE) -o $@ $(FLAGS)

clean:
	rm fort.* *~ *.exe *.out *.junk */*.o */*~ *~ *.o


