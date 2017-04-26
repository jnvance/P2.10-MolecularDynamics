SOURCE=simplemd.cpp Vector.cpp Random.cpp
MPICC=mpic++

all: simplemd.x

simplemd.x: $(SOURCE)
	$(MPICC) -g $(SOURCE) -O -o simplemd.x -D_STRIDE

clean:
	rm -rf simplemd.x *.o *.dSYM
	
run: simplemd.x
	cd ../input; /usr/bin/time mpirun -np 6 ../src/simplemd.x in_parallel

purge:
	cd ../input;
