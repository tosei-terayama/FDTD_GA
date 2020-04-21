LIBS = -L. -lnrlmsise

all: main libnrlmsise.a

main: main.o libnrlmsise.a
	g++ -o $@ main.o $(LIBS)

%.o: %.cpp
	g++ -c $< -std=c++1y

%.o: %.c
	g++ -c $< -std=c++1y

LIBOBJS = nrlmsise-00.o nrlmsise-00_data.o
libnrlmsise.a: $(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

clean:
	rm -rf *.o libnrlmsise.a main

