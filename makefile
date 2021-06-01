
include ./make.var

all: lib

doc:
	doxygen

lib:
	cd lib && make all
	
ex:
	cd src && make all
	$(CC) -o ./bin/runSBD ./src/runSBD.o ./lib/libSBD.a $(LIBS)
	./bin/runSBD ./examples/example1
	
clean:	
	cd bin && make clean
	cd doc && make clean
	cd lib && make clean
	cd examples && make clean
	
.PHONY: lib doc
