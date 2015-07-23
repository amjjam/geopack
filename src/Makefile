GEOPACK_OBJ=geopack.o geopack2008dp.o
FFLAGS=-fPIC -g
FC=gfortran
CPPFLAGS=-fPIC -g
CPP=g++


build: libgeopack.so

install: build  ../lib/libgeopack.so

libgeopack.so: $(GEOPACK_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o$@ $^

../lib/libgeopack.so: libgeopack.so
	- mkdir ../lib
	cp libgeopack.so ../lib

clean:
	-rm $(GEOPACK_OBJ)

uninstall:
	- rm ../lib/libgeopack.so
	- rmdir ../lib/
