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
int scenario0(int simtime, ofstream &output, double beta, double gamma);
int scenario1(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min, double R0);
int scenario2(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min);
int scenario3(int simtime, int onset_tprime, int duration_tprime, ofstream& output, double beta, double gamma, double beta_min);
int scenario3_rampup(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min);
int scenario4(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min);
int scenario5(int simtime, int onset_tprime, int duration_cycle, ofstream &output, double beta, double gamma, double beta_min);

double calc_beta_min(double beta, double fraction, int old_duration, int new_duration, bool linear=1);

int vary_intervention_start();
int vary_intervention_duration();
int vary_R0();

double gentime(double T2, double R0){
  double G = T2 * ((R0-1)/log(2));
  return(G);
}

double calc_beta_min_R0(double beta, double fraction, int old_duration, int new_duration, bool linear)
{
    double auc = old_duration*fraction*beta/2.0;
    double new_fraction = auc/(new_duration*beta/2);
    if (!linear){
        new_fraction = auc/(new_duration*beta);
    }
    double beta_min = (1-new_fraction)*beta;
    return beta_min;
}

double calc_beta_min(double beta, double fraction, int old_duration, int new_duration, bool linear)
{
    double auc = old_duration*fraction*beta/2.0;
    double new_fraction = auc/(new_duration*beta/2);
    if (!linear){
    	new_fraction = auc/(new_duration*beta);
    }
    double beta_min = (1-new_fraction)*beta;
    return beta_min;
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
    //vary_intervention_start();
    //vary_intervention_duration();
    vary_R0();
    return 0;
}

int vary_intervention_start()
{
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string output_dir, filename;
    ofstream output;

    //Baseline, R0=2, T2=6, reduced R0=0.5, t'=41, duration = 12 weeks
    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    gamma = 1/gentime(T2,R0);
    simtime = 365;
   
    intervention_duration = (7*12);
    //calc_beta_min(beta, fraction, old duration, new_duration);
    duration_cycle = 14;
    output_dir = "/mnt/d/covid-19/output/";
    
    
    filename = output_dir + "scenario1-5_varying_startday.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day\n";
    for (intervention_start_day = 26; intervention_start_day<=56; ++intervention_start_day){
    	//Set file headers
    	    	//Run simulations
        scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, 0.1444,R0);
    	scenario2(simtime, intervention_start_day, intervention_duration, output, beta, gamma,0.05775);
    	scenario3(simtime, intervention_start_day, intervention_duration, output, beta, gamma,0.05775);
    	scenario4(simtime, intervention_start_day, intervention_duration, output, beta, gamma,0.05775);
    	scenario5(simtime, intervention_start_day, duration_cycle, output, beta, gamma, 0.05775);
    }
    output.close();
    return 0;
}

int vary_intervention_duration()
{
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string output_dir, filename;
    ofstream output;

    //Baseline, R0=2, T2=6, reduced R0=0.5, t'=41, duration = 12 weeks
    R0 = 2;
    T2 = 6;
    beta = R0*(1/gentime(T2,R0));
    gamma = 1/gentime(T2,R0);
    simtime = 365;
	intervention_start_day = 41;
    //intervention_duration = (7*12);
    //calc_beta_min(beta, fraction, old duration, new_duration);
    duration_cycle = 14;
    output_dir = "/mnt/d/covid-19/output/";
    
    filename = output_dir + "scenario1-5_varying_duration.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration\n";
    for (intervention_duration = 8*7; intervention_duration<=16*7; ++intervention_duration){
        scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, 0.1444,R0);
    	scenario2(simtime, intervention_start_day, intervention_duration, output, beta, gamma,0.05775);
    	scenario3(simtime, intervention_start_day, intervention_duration, output, beta, gamma,0.05775);
    	scenario4(simtime, intervention_start_day, intervention_duration, output, beta, gamma,0.05775);
    	//scenario5(simtime, intervention_start_day, duration_cycle, output, beta, gamma, 0.05775);
    }
    output.close();
    return 0;
}

int vary_R0(){
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day, intervention_duration, duration_cycle;
    string output_dir, filename;
    ofstream output;
    T2 = 6;
    gamma = 1/gentime(T2,R0);
    simtime = 365;
    //intervention_start_day = 41;
    intervention_duration = (7*12);
    //calc_beta_min(beta, fraction, old duration, new_duration);
    duration_cycle = 14;
    output_dir = "/mnt/d/covid-19/output/";

    filename = output_dir + "scenario1-5_varying_R0_startday.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    for (R0=1.5; R0<=3.0; R0+=0.5){
        beta = R0*(1/gentime(T2,R0));
        double beta_red = beta*0.625;
        for (intervention_start_day = 26; intervention_start_day<=56; ++intervention_start_day){
            scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_red, R0);
        }
    }
    output.close();
    intervention_start_day = 41;
    filename = output_dir + "scenario1-5_varying_R0_duration.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    for (R0=1.5; R0<=3.0; R0+=0.5){
        beta = R0*(1/gentime(T2,R0));
        double beta_red = beta*(0.625);
        for (intervention_duration = 8*7; intervention_duration<=16*7; ++intervention_duration){
            scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_red, R0);
        }
    }
    output.close();
    return 0;
}

int scenario0(int simtime, ofstream &output, double beta, double mu)
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",2," << onset_tprime << "," << duration_tprime <<"\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",2," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",2," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",3," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",3," <<onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",3, " << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << "," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << "," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << "," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << "," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << "," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << "," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",4," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",4," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",4," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",4," << onset_tprime << "," << duration_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",5," << onset_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",5," << onset_tprime << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}
