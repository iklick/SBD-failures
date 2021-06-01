
include ../make.var

all: lib

lib:
	$(CC) $(CFLAGS) -I../include $(INCS) -c SBD.c SBD_qgraph.c SBD_arpack.c SBD_utils.c SBD_lanczos.c SBD_nauty.c
	ar -rcs libSBD.a SBD.o SBD_qgraph.o SBD_arpack.o SBD_utils.o SBD_lanczos.o SBD_nauty.o
	
clean:
	rm -f *.o *.a
	
