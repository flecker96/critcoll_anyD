CF = gfortran
FLAGS = -O -Wall
.f.o:
	$(CF) -c $(FLAGS) -o $*.o $*.f 

.KEEP_STATE:

FOUR = fourier/four1.o fourier/fourdiff2.o fourier/inhom.o fourier/lamint.o \
	fourier/halve.o fourier/double.o fourier/mik.o fourier/fourinterpol.o \
	fourier/four_to_three.o fourier/three_to_four.o fourier/fourdiff1.o \
	fourier/fourdiff1smoothe.o

SHOOTINNER = inner/shootmain_inner.o inner/shoot1new.o \
	shoot/implicit_step.o shoot/irk2_stepper.o \
	shoot/gen_grid.o \
	$(FOUR) solve/lubksb.o solve/ludcmp.o solve/mprove.o \
	inner/derivs.o inner/yfromfields.o inner/fieldsfromy.o \
	inner/leftdata.o inner/rightdata.o inner/mismatch.o\
	shoot/pack2.o shoot/unpack2.o

N_conversion = N_conversion/128to512.o $(FOUR)

shoot_inner.exe: $(SHOOTINNER)
	$(CF) $(SHOOTINNER) -o $@ $(FLAGS)

128to512.exe: $(N_conversion)
	$(CF) $(N_conversion) -o $@ $(FLAGS)

clean:
	rm fort.* *~ *.exe *.out *.junk */*.o */*~ *~ *.o


