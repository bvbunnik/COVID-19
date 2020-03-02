# COVID-19
Code for modelling the effect of social distancing measures on COVID-19 using a simple SIR model.

C++ code uses routines from the GNU Scientific Library (GSL) (https://www.gnu.org/software/gsl/), specifically the functions contained in the Ordinary Differential Equation library (gsl_odeiv2.h).

C++ code has been tested to work with the g++ compiler (version 7.4.0).
To compile run: `g++ -std=c++11  main.cpp -lgsl -lgslcblas -lm -o covid-19`.
