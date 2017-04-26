SOURCE=simplemd.cpp Vector.cpp Random.cpp

all: simplemd.x

simplemd.x: $(SOURCE)
	g++-6 -g $(SOURCE) -O -o simplemd.x

clean:
	rm -rf simplemd.x *.o *.dSYM
	
run: simplemd.x
	cd ../input; ../serial/simplemd.x in
