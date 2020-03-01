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

using namespace std;

typedef struct model_params{
    double beta;
    double mu;
} model_params;


int odes(double t, const double y[], double f[], void *params);
int scenario0(int simtime, ofstream &output, double beta, double mu);
int scenario1(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min=0);
int scenario2(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min=0);
int scenario3(int simtime, int onset_tprime, int duration_tprime, ofstream& output, double beta, double mu, double beta_min=0);
int scenario3_rampup(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min);
int scenario4(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min=0);
int scenario5(int simtime, int onset_tprime, int duration_cycle, ofstream &output, double beta, double mu, double beta_min=0);

int baseline_sims();
int shortened_extended_sims();
int increased_R0();
int decreased_R0();
int less_reduced_beta1();
int less_reduced_beta2();
int rampup();


double gentime(double T2, double R0){
  double G = T2 * ((R0-1)/log(2));
  return(G);
}


int odes(double t, const double y[], double f[], void *params)
{
    (void)(t);
    model_params *p = static_cast<model_params *>(params);
    f[0] = -p->beta * y[0] * y[1];
    f[1] = p->beta * y[0] *y[1] - p->mu*y[1];
    f[2] = p->mu*y[1];
    f[3] = p->beta*y[0]*y[1];
    return GSL_SUCCESS;
}

int main()
{
    baseline_sims();
    shortened_extended_sims();
    increased_R0();
    decreased_R0();
    less_reduced_beta1();
    less_reduced_beta2();
    rampup();
    return 0;
}

int baseline_sims()
{
    double R0, T2, beta, mu;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string filename;
    ofstream output;

    //Baseline, R0=2, T2=6, reduced R0=0.5, t'=41, duration = 12 weeks
    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    mu = 1/gentime(T2,R0);
    simtime = 365;
    intervention_start_day = 41;
    intervention_duration = (7*12);
    duration_cycle = 14;

    filename = "/home/bram/Documents/COVID-19/output/final/output_all_scenarios_baseline.csv";
    output.open(filename, ios::out);

    //Set file headers
    output << "t,S,I,R,C,beta,scen\n";
    //Run simulations
    scenario0(simtime, output, beta, mu);
    scenario1(simtime, intervention_start_day, intervention_duration, output, beta, mu, 0.1444);
    scenario2(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.05775);
    scenario3(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.05775);
    scenario4(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.05775);
    scenario5(simtime, intervention_start_day, duration_cycle, output, beta, mu ,0.05775);
    output.close();
    return 0;
}

int shortened_extended_sims()
{
    double R0, T2, beta, mu;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string filename;
    ofstream output;

    //Shortened / Extended intervantion, R0=2, T2=6, reduced R0=0.5, t'=41, duration = 10 weeks and duration = 14 weeks
    // Initially for S3 only
    //Shortened:
    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    mu = 1/gentime(T2,R0);
    simtime = 365;
    intervention_start_day = 41;
    //changed duration:
    intervention_duration = (7*10);
    //Todo: recalculate cycle length for S5
    duration_cycle = 14;

    filename = "/home/bram/Documents/COVID-19/output/final/output_scenario3_shortened_baseline.csv";
    output.open(filename, ios::out);
    //Set file headers
    output << "t,S,I,R,C,beta,scen\n";
    //Run simulations
    //scenario0(simtime, output, beta, mu);
    //scenario1(simtime, intervention_start_day, intervention_duration, output, beta, mu, 0.12705);
    //scenario2(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.0321);
    scenario3(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.0321);
    //scenario4(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.0321);
    //scenario5(simtime, intervention_start_day, duration_cycle, output, beta, mu ,0.0321);
    output.close();


    //Extended:
    intervention_duration = (7*14);
    //Todo: recalculate cycle length for S5
    duration_cycle = 14;

    filename = "/home/bram/Documents/COVID-19/output/final/output_scenario3_extended_baseline.csv";
    output.open(filename, ios::out);
    //Set file headers
    output << "t,S,I,R,C,beta,scen\n";
    //Run simulations
    //scenario0(simtime, output, beta, mu);
    //scenario1(simtime, intervention_start_day, intervention_duration, output, beta, mu, 0.15675);
    //scenario2(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.083);
    scenario3(simtime, intervention_start_day, intervention_duration, output, beta, mu, 0.083);
    //scenario4(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.083);
    //scenario5(simtime, intervention_start_day, duration_cycle, output, beta, mu, 0.083);
    output.close();
    return 0;
}


int increased_R0()
{
    double R0, T2, beta, mu;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string filename;
    ofstream output;

    //R0=3, T2=6, reduced R0=0.5, t'=41, duration = 12 weeks
    R0 = 3;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    mu = 1/gentime(T2,R0);
    simtime = 365;
    intervention_start_day = 41;
    intervention_duration = (7*12);
    duration_cycle = 14;

    filename = "/home/bram/Documents/COVID-19/output/final/output_all_scenarios_increased_R0.csv";
    output.open(filename, ios::out);

    //Set file headers
    output << "t,S,I,R,C,beta,scen\n";
    //Run simulations
    scenario0(simtime, output, beta, mu);
    scenario1(simtime, intervention_start_day, intervention_duration, output, beta, mu, 0.1001);
    scenario2(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.0288);
    scenario3(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.0288);
    scenario4(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.0288);
    scenario5(simtime, intervention_start_day, duration_cycle, output, beta, mu ,0.0288);
    output.close();
    return 0;
}

int decreased_R0()
{
    double R0, T2, beta, mu;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string filename;
    ofstream output;

    //R0=1.5, T2=6, reduced R0=0.5, t'=41, duration = 12 weeks
    R0 = 1.5;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    mu = 1/gentime(T2,R0);
    simtime = 365;
    intervention_start_day = 41;
    intervention_duration = (7*12);
    duration_cycle = 14;

    filename = "/home/bram/Documents/COVID-19/output/final/output_all_scenarios_decreased_R0.csv";
    output.open(filename, ios::out);

    //Set file headers
    output << "t,S,I,R,C,beta,scen\n";
    //Run simulations
    scenario0(simtime, output, beta, mu);
    scenario1(simtime, intervention_start_day, intervention_duration, output, beta, mu, 0.289);
    scenario2(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.116);
    scenario3(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.116);
    scenario4(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.116);
    scenario5(simtime, intervention_start_day, duration_cycle, output, beta, mu ,0.116);
    output.close();
    return 0;
}

int less_reduced_beta1()
{
    double R0, T2, beta, mu;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string filename;
    ofstream output;

    //R0=2, T2=6, reduced R0=1.0, t'=41, duration = 12 weeks
    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    mu = 1/gentime(T2,R0);
    simtime = 365;
    intervention_start_day = 41;
    intervention_duration = (7*12);
    duration_cycle = 14;

    filename = "/home/bram/Documents/COVID-19/output/final/output_all_scenarios_less_reduced_beta.csv";
    output.open(filename, ios::out);

    //Set file headers
    output << "t,S,I,R,C,beta,scen\n";
    //Run simulations
    scenario0(simtime, output, beta, mu);
    scenario1(simtime, intervention_start_day, intervention_duration, output, beta, mu, 0.1732);
    scenario2(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.116);
    scenario3(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.116);
    scenario4(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.116);
    scenario5(simtime, intervention_start_day, duration_cycle, output, beta, mu ,0.116);
    output.close();
    return 0;
}

int less_reduced_beta2()
{
    double R0, T2, beta, mu;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string filename;
    ofstream output;

    //R0=2, T2=6, reduced R0=1.5, t'=41, duration = 12 weeks
    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    mu = 1/gentime(T2,R0);
    simtime = 365;
    intervention_start_day = 41;
    intervention_duration = (7*12);
    duration_cycle = 14;

    filename = "/home/bram/Documents/COVID-19/output/final/output_all_scenarios_less_reduced_beta2.csv";
    output.open(filename, ios::out);

    //Set file headers
    output << "t,S,I,R,C,beta,scen\n";
    //Run simulations
    scenario0(simtime, output, beta, mu);
    scenario1(simtime, intervention_start_day, intervention_duration, output, beta, mu, 0.202);
    scenario2(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.1733);
    scenario3(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.1733);
    scenario4(simtime, intervention_start_day, intervention_duration, output, beta, mu ,0.1733);
    scenario5(simtime, intervention_start_day, duration_cycle, output, beta, mu ,0.1733);
    output.close();
    return 0;
}

int rampup()
{
    double R0, T2, beta, mu;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string filename;
    ofstream output;

    //R0=2, T2=6, reduced R0=0.5, t'=41, duration = 12 weeks
    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    mu = 1/gentime(T2,R0);
    simtime = 365;
    intervention_start_day = 41;
    intervention_duration = (7*12);
    duration_cycle = 14;

    filename = "/home/bram/Documents/COVID-19/output/final/output_scenario3_rampup.csv";
    output.open(filename, ios::out);

    //Set file headers
    output << "t,S,I,R,C,beta,scen\n";
    //Run simulations
    scenario3_rampup(simtime, intervention_start_day, intervention_duration, output, beta, mu, 0.05775);
    output.close();
    return 0;
}

int scenario0(int simtime,  ofstream &output, double beta, double mu)
{
    model_params pars = {beta, mu};
    gsl_odeiv2_system sys = {odes, nullptr, 4, &pars};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    double y[4] = {0.9999, 0.0001, 0,0};
    double t = 0.0;
    for (int i = 0; i<=simtime; ++i){
        double ti = i + 1;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",0" << "\n";
    }
    return 0;
}

int scenario1(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min)
{
    model_params pars = {beta, mu};
    gsl_odeiv2_system sys = {odes, nullptr, 4, &pars};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    double y[4] = {0.9999, 0.0001, 0, 0};
    double t = 0.0;
    for (int i = 0; i<onset_tprime; ++i){
        double ti = t + 1.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",1" << "\n";
    }
    for (int i = onset_tprime; i<(onset_tprime+duration_tprime); ++i){
        double ti = t+1.0;
        double beta1;
        beta1=beta_min;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",1" << "\n";
    }
    for (int i = (onset_tprime+duration_tprime); i<=simtime; ++i){
        double ti = t+1.0;
        double beta1 = beta;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",1" << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}


int scenario2(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min)
{
    double beta_decrease;
    beta_decrease = (beta-beta_min)/duration_tprime;
    model_params pars = {beta, mu};
    gsl_odeiv2_system sys = {odes, nullptr, 4, &pars};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    double y[4] = {0.9999, 0.0001, 0,0};
    double t = 0.0;
    for (int i = 0; i<onset_tprime; ++i){
        double ti = t+1.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",2" << "\n";
    }
    int a = 0;
    for (int i = onset_tprime; i<onset_tprime+duration_tprime; ++i){
        double ti = t+1.0;
        double beta1 = beta_min+beta_decrease*a;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",2" << "\n";
        ++a;
    }
    for (int i = (onset_tprime+duration_tprime); i<=simtime; ++i){
        double ti = t+1.0;
        double beta1 = beta;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",2" << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}

int scenario3(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min)
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",3"  << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",3" << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",3" << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}

int scenario3_rampup(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min)
{
    double beta_decrease;
    beta_decrease = (beta-beta_min)/4.0;
    int duration = (duration_tprime/4);
    model_params pars = {beta, mu};
    gsl_odeiv2_system sys = {odes, nullptr, 4, &pars};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    double y[4] = {0.9999, 0.0001, 0,0};
    double t = 0.0;
    for (int i = 0; i<(onset_tprime+round(0.5*duration)); ++i){
        double ti = t+1.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << "," << onset_tprime << "\n";
    }
    int a = 0;
    int onset = (onset_tprime+round(0.5*duration));
    for (int i = onset; i<onset+duration; ++i){
        double ti = t+1.0;
        double beta1 = beta-beta_decrease;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << "," << onset_tprime << "\n";
        ++a;
    }
    onset = onset+(duration_tprime/4);
    for (int i = onset; i<onset+duration; ++i){
        double ti = t+1.0;
        double beta1 = beta-2*beta_decrease;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << "," << onset_tprime << "\n";
        ++a;
    }
    onset = onset+(duration_tprime/4);
    for (int i = onset; i<onset+duration; ++i){
        double ti = t+1.0;
        double beta1 = beta-3*beta_decrease;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << "," << onset_tprime << "\n";
        ++a;
    }
    onset = onset+(duration_tprime/4);
    for (int i = onset; i<(onset+round(0.5*duration)); ++i){
        double ti = t+1.0;
        double beta1 = beta-4*beta_decrease;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << "," << onset_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << "," << onset_tprime << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}


int scenario4(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min)
{
    double beta_change;
    beta_change = (beta-beta_min)/(duration_tprime/2);
    double beta_temp=0;
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",4"  << "\n";
    }
    int a = 0;
    for (int i = onset_tprime; i<onset_tprime+(duration_tprime/2); ++i){
        double ti = t+1.0;
        double beta1 = beta-beta_change*a;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",4" << "\n";
        ++a;
        beta_temp = beta1;
    }
    a = 0;
    for (int i = onset_tprime+(duration_tprime/2); i<onset_tprime+duration_tprime; ++i){
        double ti = t+1.0;
        double beta1 = beta_temp+beta_change*a;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",4" << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",4" << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}

int scenario5(int simtime, int onset_tprime, int duration_cycle, ofstream &output, double beta, double mu, double beta_min)
{
    double beta_low;
    beta_low = beta_min;
    int duration = duration_cycle;
    model_params pars = {beta, mu};
    gsl_odeiv2_system sys = {odes, nullptr, 4, &pars};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    double y[4] = {0.9999, 0.0001, 0,0};
    double t = 0.0;
    for (int i = 0; i<(onset_tprime+7); ++i){
        double ti = t+1.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",5"  << "\n";
    }
    int a = 0;
    int onset = (onset_tprime+7); //t'+1
    for (int i = onset; i<onset+duration; ++i){
        double ti = t+1.0;
        double beta1 = beta_low;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5"  << "\n";
        ++a;
    }
    onset = onset+duration; //t'+3
    for (int i = onset; i<onset+duration; ++i){
        double ti = t+1.0;
        double beta1 = beta;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5" << "\n";
        ++a;
    }
    onset = onset+duration; //t'+5
    for (int i = onset; i<onset+duration; ++i){
        double ti = t+1.0;
        double beta1 = beta_low;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5" << "\n";
        ++a;
    }
    onset = onset+duration; //t'+7
    for (int i = onset; i<(onset+duration); ++i){
        double ti = t+1.0;
        double beta1 = beta;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5" << "\n";
        ++a;
    }
    onset = onset+duration; //t'+9
    for (int i = onset; i<(onset+duration); ++i){
        double ti = t+1.0;
        double beta1 = beta_low;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5" << "\n";
        ++a;
    }
    onset = onset+duration; //t'+11
    for (int i = onset; i<(onset+duration); ++i){
        double ti = t+1.0;
        double beta1 = beta;
        pars = {beta1, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5" << "\n";
        ++a;
    }
    for (int i = onset+duration; i<=simtime; ++i){
        double ti = t+1.0;
        pars = {beta, mu};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",5" << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}
