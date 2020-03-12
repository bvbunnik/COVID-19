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
int scenario0(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min, double R0);
int scenario1(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min, double R0);
int scenario1_mult_int(int simtime, int onset_tprime, int onset_tprime1, int duration_tprime, int duration_tprime1, ofstream &output, double beta, double gamma, double beta_min);
int scenario2(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min, double R0);
int scenario3(int simtime, int onset_tprime, int duration_tprime, ofstream& output, double beta, double gamma, double beta_min, double R0);
int scenario4(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double gamma, double beta_min, double R0);
int scenario5(int simtime, int onset_tprime, int duration_cycle, ofstream &output, double beta, double gamma, double beta_min, double R0);
int vary_intervention_start();
int vary_R0();
int vary_beta_red();
int vary_T2();
int multiple_intervention();


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
    //vary_intervention_start();
    //vary_R0();
    //vary_beta_red();
    //vary_T2();
    multiple_intervention();
    return 0;
}


int multiple_intervention()
{
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day1, intervention_duration, intervention_duration1, intervention_duration2;
    string output_dir, filename;
    ofstream output;

    R0 = 2.4;
    T2 = 4.6;
    //beta = R0*(1/gentime(T2,R0));
    gamma = 1/gentime(T2,R0);
    //reduced beta due to background reduction (37.5% reduction)
    beta = gamma*1.5;
    simtime = 364*2;

    //intervention_duration1 = (7*12);
    output_dir = "/media/bram/DATA/covid-19/output/";

    filename = output_dir + "scenario1_multiple_interventions_extended_range.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day1,start_day2\n";
    //beta_min is further reduced by 40% to get to a total reduction of 75% of original R0
    double beta_min=beta*0.4;
    for (intervention_start_day1 = 0; intervention_start_day1<=560; ++intervention_start_day1){
        for (int intervention_start_day2 = intervention_start_day1+intervention_duration; intervention_start_day2 <= simtime-intervention_duration; ++intervention_start_day2){
            scenario1_mult_int(simtime, intervention_start_day1, intervention_start_day2, intervention_duration, intervention_duration, output, beta, gamma, beta_min);
        }
    }
    for (intervention_start_day1 = 0; intervention_start_day1<=560; ++intervention_start_day1){
        for (intervention_duration1 = 4*7; intervention_duration1<=20*7; intervention_duration1+=7){
            for (int intervention_start_day2 = intervention_start_day1+intervention_duration1; intervention_start_day2 <= simtime-intervention_duration1; ++intervention_start_day2){
                for (intervention_duration2 = 4*7; intervention_duration2<=20*7; intervention_duration2+=7){
                    scenario1_mult_int(simtime, intervention_start_day1, intervention_start_day2, intervention_duration1, intervention_duration2, output, beta, gamma, beta_min);
                }
            }
        }
    }
    return 0;
}


int vary_intervention_start()
{
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day, intervention_duration;
    string output_dir, filename;
    ofstream output;

    R0 = 2.4;
    T2 = 4.6;
    //beta = R0*(1/gentime(T2,R0));
    gamma = 1/gentime(T2,R0);
    //reduced beta due to background reduction (37.5% reduction)
    beta = gamma*1.5;
    simtime = 364*2;
   
    intervention_duration = (7*12);
    output_dir = "/media/bram/DATA/covid-19/output/";
    
    filename = output_dir + "scenario1_new_paras.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    double beta_min=beta*0.4;

    for (intervention_start_day = 1; intervention_start_day<=simtime-intervention_duration; ++intervention_start_day){
        scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_min,R0);
    }
//    int duration_cycle = 14;
//    beta_min=(beta*(1.0-0.375))/2.0;
//    for (intervention_start_day = 1; intervention_start_day<=100; ++intervention_start_day){
//        scenario2(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_min,R0);
//        scenario3(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_min,R0);
//        scenario4(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_min,R0);
//        scenario5(simtime, intervention_start_day,duration_cycle, output, beta, gamma, beta_min,R0);
//    }
    output.close();
    return 0;
}


int vary_R0(){
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day, intervention_duration;
    string output_dir, filename;
    ofstream output;
    T2 = 6.0;

    gamma = 1/gentime(T2,2);
    simtime = 365;

    intervention_start_day = 49;
    intervention_duration = (7*12);
    output_dir = "/media/bram/DATA/covid-19/output/";

    filename = output_dir + "scenario1_varying_R0_startday.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    for (double R0=1.2; R0<=3.05;R0+=0.1){
        beta = R0*gamma;
        double beta_red = (beta*(1-0.375));
        scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta_red, R0);
    }
    output.close();


    intervention_start_day = 29;
    filename = output_dir + "scenario3_varying_R0_startday.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    for (R0=1.2; R0<=3.05; R0+=0.1){
        beta = R0*gamma;
        double beta_red = (beta*(1-0.375))/2.0;
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

int vary_T2(){
    double R0, T2, beta, gamma;
    int simtime, intervention_start_day, intervention_duration;
    string output_dir, filename;
    ofstream output;

    R0 = 2.4;
    simtime = 365;
    intervention_start_day = 31;
    intervention_duration = (7*12);
    int duration_cycle = 14;
    double beta_red = 0.375;
    output_dir = "/media/bram/DATA/covid-19/output/";
    filename = output_dir + "scenario1-5_new_RWC_vary_T2.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,T2\n";
    for (T2=3.0; T2<=8.05; T2+=0.1) {
        beta = R0*(1.0/gentime(T2,R0));
        gamma = 1.0/gentime(T2,R0);
        scenario0(simtime, 0, intervention_duration, output, beta, gamma, beta*(1.0-beta_red),T2);
        scenario1(simtime, 31, intervention_duration, output, beta, gamma, beta*(1.0-beta_red),T2);
        scenario2(simtime, 45, intervention_duration, output, beta, gamma, beta*(1.0-beta_red)/2.0,T2);
        scenario3(simtime, 14, intervention_duration, output, beta, gamma, beta*(1.0-beta_red)/2.0,T2);
        scenario4(simtime, 29, intervention_duration, output, beta, gamma, beta*(1.0-beta_red)/2.0,T2);
        scenario5(simtime, 32, duration_cycle, output, beta, gamma, beta*(1.0-beta_red)/2.0,T2);
    }
    output.close();
    return 0;
}


int vary_start_duration()
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
    double beta_red=0.375;
    intervention_duration = (7*12);
    //calc_beta_min(beta, fraction, old duration, new_duration);
    duration_cycle = 14;
    output_dir = "/media/bram/DATA/covid-19/output/";

    filename = output_dir + "scenario1_startday+duration.csv";
    output.open(filename, ios::out);
    output << "t,S,I,R,C,beta,scen,start_day,duration,R0\n";
    for (intervention_start_day = 1; intervention_start_day<=100; ++intervention_start_day){
        for (intervention_duration=1; intervention_duration<=16*7;++intervention_duration){
            scenario1(simtime, intervention_start_day, intervention_duration, output, beta, gamma, beta*(1.0-beta_red),R0);
        }
    }
    output.close();
    return 0;
}



int scenario0(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min, double R0)
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",0" << ",0,0," << R0 << "\n";
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


int scenario1_mult_int(int simtime, int onset_tprime, int onset_tprime1, int duration_tprime, int duration_tprime1, ofstream &output, double beta, double gamma, double beta_min)
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",1," << onset_tprime << "," << onset_tprime1 << "," << duration_tprime << "," << duration_tprime1 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",1," << onset_tprime << "," << onset_tprime1 << "," << duration_tprime << "," << duration_tprime1 << "\n";
    }
    for (int i = (onset_tprime+duration_tprime); i<onset_tprime1; ++i){
        double ti = t+1.0;
        double beta1 = beta;
        pars = {beta1, gamma};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",1," << onset_tprime << "," << onset_tprime1 << "," << duration_tprime << "," << duration_tprime1 <<"\n";
    }
    for (int i = onset_tprime1; i<(onset_tprime1+duration_tprime1); ++i){
        double ti = t+1.0;
        double beta1;
        beta1=beta_min;
        pars = {beta1, gamma};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",1," << onset_tprime << "," << onset_tprime1 << "," << duration_tprime << "," << duration_tprime1 << "\n";
    }
    for (int i = (onset_tprime1+duration_tprime1); i<=simtime; ++i){
        double ti = t+1.0;
        double beta1 = beta;
        pars = {beta1, gamma};
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            return 1;
        }
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",1," << onset_tprime << "," << onset_tprime1 << "," << duration_tprime << "," << duration_tprime1 << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}



int scenario2(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min, double R0)
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",2," << onset_tprime << "," << duration_tprime << "," << R0 <<"\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",2," << onset_tprime << "," << duration_tprime << "," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",2," << onset_tprime << "," << duration_tprime << "," << R0 << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}

int scenario3(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min, double R0)
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",3," << onset_tprime << "," << duration_tprime << "," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",3," <<onset_tprime << "," << duration_tprime << "," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",3, " << onset_tprime << "," << duration_tprime << "," << R0 << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}



int scenario4(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_min, double R0)
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",4," << onset_tprime << "," << duration_tprime << "," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",4," << onset_tprime << "," << duration_tprime << "," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",4," << onset_tprime << "," << duration_tprime << "," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",4," << onset_tprime << "," << duration_tprime << "," << R0 << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}

int scenario5(int simtime, int onset_tprime, int duration_cycle, ofstream &output, double beta, double mu, double beta_min, double R0)
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",5," << onset_tprime << ", ," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << ", ," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << ", ," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << ", ," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << ", ," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << ", ," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta1 << ",5," << onset_tprime << ", ," << R0 << "\n";
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
        output << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << beta << ",5," << onset_tprime << ", ," << R0 << "\n";
    }
    gsl_odeiv2_driver_free(d);
    return 0;
}
