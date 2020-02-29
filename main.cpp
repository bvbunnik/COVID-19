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
int scenario0(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu);
int scenario1(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction);
int scenario2(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction);
int scenario3(int simtime, int onset_tprime, int duration_tprime, ofstream& output, double beta, double mu, double beta_fraction);
int scenario3_rampup(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu);
int scenario4(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction);
int scenario5(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction);


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
    double R0 = 2;
    double T2 = 6;
    double beta = R0*(1/gentime(T2,R0));
    double mu = 1/gentime(T2,R0);
    int simtime = 365;
    int intervention_start_day = 41;
    int intervention_duration = (7*12);
    double beta_fraction = 1.0/2.0;
    string filename = "/home/bram/Documents/Alex/COVID-19/output/output_all_scenarios_R0=2_beta=0.5.csv";
    ofstream output(filename, ios::out);
    output << "t,S,I,R,C,beta,scen\n";
    scenario0(simtime, intervention_start_day, intervention_duration, output, beta, mu);
    scenario1(simtime, intervention_start_day, intervention_duration, output, beta, mu, beta_fraction);
    scenario2(simtime, intervention_start_day, intervention_duration, output, beta, mu, beta_fraction);
    scenario3(simtime, intervention_start_day, intervention_duration, output, beta, mu, beta_fraction);
    scenario4(simtime, intervention_start_day, intervention_duration, output, beta, mu, beta_fraction);
    scenario5(simtime, intervention_start_day, intervention_duration, output, beta, mu, beta_fraction);

    output.close();
    return 0;
}

int scenario0(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu)
{
    model_params pars = {beta, mu};
    gsl_odeiv2_system sys = {odes, nullptr, 4, &pars};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    double y[4] = {0.9999, 0.0001, 0,0};
    double t = 0.0;
    ulong t1 = 364;
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

int scenario1(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction)
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
        double beta1 = (1-(beta_fraction+(0.5*beta_fraction)))*beta;
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


int scenario2(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction)
{
    double beta_decrease = (beta-beta_fraction*beta)/duration_tprime;
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
        double beta1 = beta_fraction*beta+beta_decrease*a;
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

int scenario3(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction)
{
    double beta_decrease = (beta-beta_fraction*beta)/duration_tprime;
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

int scenario3_rampup(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction)
{
    double beta_decrease = (beta-beta_fraction*beta)/4;
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


int scenario4(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction)
{
    double beta_change = (beta-beta_fraction*beta)/(duration_tprime/2);
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

int scenario5(int simtime, int onset_tprime, int duration_tprime, ofstream &output, double beta, double mu, double beta_fraction)
{
    double beta_low = beta_fraction*beta;
    int duration = 14;
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
