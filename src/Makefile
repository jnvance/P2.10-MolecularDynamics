SOURCE=simplemd.cpp Vector.cpp Random.cpp
MPICC=mpic++

all: simplemd.x simplemd_ncl.x

simplemd.x: $(SOURCE)
	$(MPICC) -g $(SOURCE) -O -o simplemd.x -D_CELL_LIST

simplemd_ncl.x: $(SOURCE)
	$(MPICC) -g $(SOURCE) -O -o simplemd_ncl.x -DNN_CELL_LIST

clean:
	rm -rf simplemd.x *.o *.dSYM

run: simplemd.x
	cd ../input; /usr/bin/time mpirun -np 8 ../src/simplemd.x *.params

purge:
	cd ../input;
