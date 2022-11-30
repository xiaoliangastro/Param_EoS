# execute this file 'bash file parameter'
if [ $1 == 'dynamic' ]
then
    g++ param_eos.cpp -lgsl -lgslcblas -lm -std=c++17 -shared -fPIC -o param_eos.so
fi
