CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`

all: CelestialMapV4

CelestialMapV4: CelestialMapV4.cpp
	$(CC) $(CFLAGS) -o CelestialMapV4 CelestialMapV4.cpp $(CFLAGSROOT) $(LIBSROOT)

Windows: CelestialMapV4.cpp
	$(CC) $(CFLAGS) -o CelestialMapV4.exe CelestialMapV4.cpp $(CFLAGSROOT) $(LIBSROOT)

clean:
	rm *.o
