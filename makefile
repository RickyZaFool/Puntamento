CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`

all: CelestialMap

CelestialMap: CelestialMap.cpp
	$(CC) $(CFLAGS) -o CelestialMap CelestialMap.cpp $(CFLAGSROOT) $(LIBSROOT)
clean:
	rm *.o