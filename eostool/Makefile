CC=/opt/spack/var/spack/environments/myenv/.spack-env/view/bin/g++
CFLAGS=-lgsl -lgslcblas -lm -std=c++17 -shared -fPIC
OBJ=param_eos.so
SOURCE=param_eos.cpp
CPP_INC=$(CMAKE_PREFIX_PATH)/include
CPP_LIB=$(CMAKE_PREFIX_PATH)/lib
LIB_GSL=/opt/spack/opt/spack/darwin-monterey-m1/apple-clang-13.1.6/gsl-2.7.1-3osqju4o5mct3e6u7hlotuhtqz3d7atm/lib
INC_BOOST=/opt/spack/opt/spack/darwin-monterey-m1/apple-clang-13.1.6/boost-1.78.0-fzzmdhxtkw6msmd6larufyol27twvr6s/include
INC_GSL=/opt/spack/opt/spack/darwin-monterey-m1/apple-clang-13.1.6/gsl-2.7.1-3osqju4o5mct3e6u7hlotuhtqz3d7atm/include

param_eos:
	$(CC) $(SOURCE) -I$(INC_BOOST) -I$(INC_GSL) -L$(LIB_GSL) $(CFLAGS) -o $(OBJ)

all: param_eos

clean:
	rm $(OBJ)
