/**
@file main.cpp
*/
/** \mainpage COVID-19 SIR model
 * Author: Bram van Bunnik
 * E-mail: bram.vanbunnik@ed.ac.uk
 * Simulate the effects of SDM during an outbreak of COVID-19.\n
 * Date: 01-03-2020
 *
 */


#include <iostream>
#include <stdio.h>
#include <fstream>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <math.h>
#include <string>


using namespace std;

typedef struct model_params{
    double beta;
    double gamma;
} model_params;


int odes(double t, const double y[], double f[], void *params);
int scenario1(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min, double R0);
int scenario3(int simtime, int onset_tprime, int duration_tprime, ofstream& output, double beta, double gamma, double beta_min, double scen);

int vary_intervention_start();
int vary_R0();
int vary_beta_red();

double gentime(double T2, double R0){
  double G = T2 * ((R0-1)/log(2));
  return(G);
}


int odes(double t, const double y[], double f[], void *params)
{
    (void)(t);
    model_params *p = static_cast<model_params *>(params);
    f[0] = -p->beta * y[0] * y[1];
    f[1] = p->beta * y[0] *y[1] - p->gamma*y[1];
    f[2] = p->gamma*y[1];
    f[3] = p->beta*y[0]*y[1];
    return GSL_SUCCESS;
}

int main()
{
    vary_intervention_start();
    vary_R0();
    vary_beta_red();
    return 0;
}

int vary_intervention_start()
{
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string output_dir, filename;
    ofstream output;

    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    gamma = 1/gentime(T2,R0);
    simtime = 365;
   
    intervention_duration = (7*12);
    duration_cycle = 14;
    output_dir = "/media/bram/DATA/covid-19/output/";
    
    filename = output_dir + "scenario1_vary_startday.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    double beta_min=(beta*0.625);
    for (intervention_start_day = 1; intervention_start_day<=100; ++intervention_start_day){
        scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_min,R0);
    }
    output.close();
    

    filename = output_dir + "scenario3_vary_startday.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    beta_min=(beta*0.625)/2.0;
    for (intervention_start_day = 1; intervention_start_day<=100; ++intervention_start_day){
        scenario3(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_min,R0);
    }
    output.close();
    return 0;
}


int vary_R0(){
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string output_dir, filename;
    ofstream output;
    T2 = 6.0;

    gamma = 1/gentime(T2,2);
    simtime = 365;

    intervention_start_day = 49;
    intervention_duration = (7*12);
    duration_cycle = 14;
    output_dir = "/media/bram/DATA/covid-19/output/";

    filename = output_dir + "scenario1_varying_R0_startday.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    for (double R0=1.2; R0<=3.05;R0+=0.1){
        beta = R0*gamma;
        double beta_red = (beta*0.625);
        scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_red, R0);
    }
    output.close();


    intervention_start_day = 29;
    filename = output_dir + "scenario3_varying_R0_startday.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    for (R0=1.2; R0<=3.05; R0+=0.1){
        beta = R0*gamma;
        double beta_red = (beta*0.625)/2.0;
        scenario3(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_red, R0);
    }
    output.close();
    return 0;
}

int vary_beta_red(){
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day, intervention_duration;
    string output_dir, filename;
    ofstream output;

    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    gamma = 1/gentime(T2,R0);
    simtime = 365;
    intervention_start_day = 49;
    intervention_duration = (7*12);
    output_dir = "/media/bram/DATA/covid-19/output/";
    filename = output_dir + "scenario1_varying_beta.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,beta_scen\n";
    for (double beta_red = 0.05; beta_red<=0.5; beta_red+=0.005){
        scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta*(1.0-beta_red),beta_red);
    }
    output.close();


    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    gamma = 1/gentime(T2,R0);
    simtime = 365;
    intervention_start_day = 29;
    intervention_duration = (7*12);
    output_dir = "/media/bram/DATA/covid-19/output/";
    filename = output_dir + "scenario3_varying_beta.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,beta_scen\n";
    for (double beta_red = 0.05; beta_red<=0.5; beta_red+=0.005){
        scenario3(simtime, intervention_start_day, intervention_duration, output, beta, gamma, (beta*(1.0-beta_red))/2.0,beta_red);
    }
    output.close();
    return 0;
}


int scenario1(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min, double R0)
{
    model_params pars = {beta, gamma};
    gsl_odeiv2_system sys = {odes, nullptr, 4, &pars};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    double y[4] = {0.9999, 0.0001, 0.0, 0.0};
    double t = 0.0;
    for (int i = 0; i<onset_tprime; ++i){
        double ti = t + 1.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",1," << onset_tprime << "," << duration_tprime << "," << R0 <<"\n";
    }
    for (int i = onset_tprime; i<(onset_tprime+duration_tprime); ++i){
        double ti = t+1.0;
        double beta1;
        beta1=beta_min;
        pars = {beta1, gamma};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",1," << onset_tprime << "," << duration_tprime << "," << R0 << "\n";
    }
    for (int i = (onset_tprime+duration_tprime); i<=simtime; ++i){
        double ti = t+1.0;
        double beta1 = beta;
        pars = {beta1, gamma};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",1," << onset_tprime << "," << duration_tprime << "," << R0 <<"\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}


int scenario3(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min, double scen)
{
    double beta_decrease;
    beta_decrease = (beta-beta_min)/duration_tprime;
    model_params pars = {beta, mu};
    gsl_odeiv2_system sys = {odes, nullptr, 4, &pars};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    double y[4] = {0.9999, 0.0001, 0, 0};
    double t = 0.0;
    for (int i = 0; i<onset_tprime; ++i){
        double ti = t+1.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",3," << onset_tprime << "," << duration_tprime << "," << scen << "\n";
    }
    int a = 0;
    for (int i = onset_tprime; i<onset_tprime+duration_tprime; ++i){
        double ti = t+1.0;
        double beta1 = beta-beta_decrease*a;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",3," <<onset_tprime << "," << duration_tprime << "," << scen << "\n";
        ++a;
    }
    for (int i = onset_tprime+duration_tprime; i<=simtime; ++i){
        double ti = t+1.0;
        pars = {beta, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",3, " << onset_tprime << "," << duration_tprime << "," << scen << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}
