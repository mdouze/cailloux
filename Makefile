
CC=g++ -std=c++11
CFLAGS=


all: _Cgeometry.so

ifeq ($(shell uname),Darwin)
  # mac flags
  SHARED_LDFLAGS=-Wl,-F. -bundle -undefined dynamic_lookup -g
else 
  # linux flags
  SHARED_LDFLAGS=-shared -g
endif



geometry2D.o : geometry2D.cpp
	$(CC) $(CFLAGS) -c geometry2D.cpp 


Cgeometry_wrap.cxx : Cgeometry.swig geometry2D.h
	swig -c++ -python $<


Cgeometry_wrap.o: Cgeometry_wrap.cxx 
	$(CC) -c $(CFLAGS)  $< -o $@  \
            -I $(shell python -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc())" )

#             -I $(shell python -c "import numpy; print numpy.get_include()")

_Cgeometry.so: geometry2D.o Cgeometry_wrap.o
	$(CC) -o $@ -g $(SHARED_LDFLAGS) $^ 


test: _Cgeometry.so
	python -m unittest test_geometry.py
