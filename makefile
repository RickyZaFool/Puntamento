CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`

all: CelestialMapV3

CelestialMapV3: CelestialMapV3.cpp
	$(CC) $(CFLAGS) -o CelestialMapV3 CelestialMapV3.cpp $(CFLAGSROOT) $(LIBSROOT)

Windows: CelestialMapV3.cpp
	$(CC) $(CFLAGS) -o CelestialMapV3.exe CelestialMapV3.cpp $(CFLAGSROOT) $(LIBSROOT)

clean:
	rm *.o
