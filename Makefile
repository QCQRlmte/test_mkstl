SHELL=/bin/sh
PROGRAM=./mkwall.exe
OBJ=mkwall.o findwall.o vectoranalysis.o

.SUFFIXES: .f .o

.f.o:
	ifort -c -O3 -convert little_endian -r8 -i4  $<

$(PROGRAM): $(OBJ)
	ifort $(OBJ)  -o $(PROGRAM)

clean:
	\rm -rf *~ *.o *.l *.exe

