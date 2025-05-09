CC         = g++ 
CFLAGS     = --std=c++11 -g -Wall
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`

all: CelestialMapV2

CelestialMapV2: CelestialMapV2.cpp
	$(CC) $(CFLAGS) -o CelestialMapV2.exe CelestialMapV2.cpp $(CFLAGSROOT) $(LIBSROOT)
clean:
	rm *.o
