# COVID-19
Code for modelling the effect of social deprivation measures on COVID-19 using a simple SIR model.

C++ code uses routines from the GNU Scientific Library (GSL), spefically the functions contained in the Ordinary Differential Equation library (gsl_odeiv2.h).

C++ code can be compiled using the g++ compiler.
to compile run: "g++ -std=c++11  main.cpp -L/usr/local/lib -lgsl -lgslcblas -lm -o covid-19"
