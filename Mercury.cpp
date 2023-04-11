#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/*
Constants
*/

const double M = 1.989e30; //Mass of sun in KG
const double m = 3.285e23; //Mass of mercury in KG
const double c = 3e8;
const double mu = M*m/(M+m); //Reduced mass
const double G = 6.67430e-11; //Gravitational constant
const double L = 9.1e38; //Momentum 

/*
Define ODEs below, only first order
*/

double theta_ode(double r)
{
    return L/(mu*pow(r,2.0));
}

double r_newtonian(double r)
{
    double v_dev = -G*M*m/pow(r,2.0);
    return (1/mu)*(v_dev + pow(L,2.0)/(mu*pow(r,3.0)));
}

double r_relativ(double r)
{
    double v_dev = -G*M*m/(pow(r,2.0));
    double v_rel = -3.0*(G/(pow(c,2.0)))*((M+m)*pow(L,2.0)/(mu*pow(r,4.0)));
    return (1/mu)*(v_dev + v_rel + pow(L,2.0)/(mu*pow(r,3.0)));
}

/*
Define methods below
*/

double euler(double step_size, double r0, double v0, double theta0)
{
    double theta = theta0 + step_size*theta_ode(r0)
     - pow(step_size,2.0)*v0/(3.0*mu*pow(r0,3.0));
    return theta;
}

void verlet_newtonian(double step_size, double *returnArray1, 
double *returnArray2, double *returnArray3)
{
    double r0n = returnArray1[1];
    double r0n_1 = returnArray2[1];
    double r = 2.0*r0n - r0n_1 + pow(step_size,2.0)*r_newtonian(r0n);
    double v = (r - r0n_1)/2.0/step_size;
    returnArray3[0] = r;
    returnArray3[1] = v;
}

void verlet_relativ(double step_size, double *returnArray1, 
double *returnArray2, double *returnArray3)
{
    double r0n = returnArray1[1];
    double r0n_1 = returnArray2[1];
    double r = 2.0*r0n - r0n_1 + pow(step_size,2.0)*r_relativ(r0n);
    double v = (r - r0n_1)/2.0/step_size;
    returnArray3[0] = r;
    returnArray3[1] = v;
}

/*
Main function below
*/

int main() 
{   
    /*
    Relativistic procession of perihelion of Mercury
    */
    double step_size = 60.0*60.0; //1 hour
    int nsteps = 24*160; //160 days
    
    double r0 = 7.5e10 ; //Orbit radius of mercury
    double v0 = -1000; //m/s
    double t0 = 0.0;
    double theta0 = 0.0;
    
    //Number of steps
    int row = nsteps;
    //Number of variables (t,r,v,theta)
    int column = 4;
    
    double values_newt[row][column] = {0.0};
    values_newt[0][0] = t0;
    values_newt[0][1] = r0;
    values_newt[0][2] = v0;
    values_newt[0][3] = theta0; 
    
    values_newt[1][0] = t0 + step_size;
    values_newt[1][1] = r0 + v0*step_size + 
    (pow(step_size,2.0)/2)*r_newtonian(r0);
    values_newt[1][2] = v0 + step_size*r_newtonian(r0);
    values_newt[1][3] = theta0 + step_size*theta_ode(r0); 
    
    double values_rel[row][column] = {0.0};
    values_rel[0][0] = t0;
    values_rel[0][1] = r0;
    values_rel[0][2] = v0;
    values_rel[0][3] = theta0;
    
    values_rel[1][0] = t0 + step_size;
    values_rel[1][1] = r0 + v0*step_size + 
    (pow(step_size,2.0)/2)*r_relativ(r0);
    values_rel[1][2] = v0 + step_size*r_relativ(r0);
    values_rel[1][3] = theta0 + step_size*theta_ode(r0);   
    
    double energy[row][2] = {0.0};
    energy[0][0] = t0;
    energy[0][1] = (mu/2)*(pow(v0,2.0) + pow(r0,2.0)*pow(theta_ode(r0),2.0));
    
    energy[1][0] = t0 + step_size;
    energy[1][1] = (mu/2)*(pow(values_rel[1][2],2.0) + 
    pow(r0,2.0)*pow(values_rel[1][3],2.0));
    
    for (int i = 2; i < nsteps; i++)
    {
        double step_values_newt_new[2] = {0.0};
        double step_values_newt_n[4] = {values_newt[i-1][0],values_newt[i-1][1],
        values_newt[i-1][2],values_newt[i-1][3]};
        double step_values_newt_n_1[4] = {values_newt[i-2][0],values_newt[i-2][1],
        values_newt[i-2][2],values_newt[i-2][3]};
        
        double step_values_rel_new[2] = {0.0};
        double step_values_rel_n[4] = {values_rel[i-1][0],values_rel[i-1][1],
        values_rel[i-1][2],values_rel[i-1][3]};
        double step_values_rel_n_1[4] = {values_rel[i-2][0],values_rel[i-2][1],
        values_rel[i-2][2],values_rel[i-2][3]};
        
        double new_theta_newt = euler(step_size, step_values_newt_n[1], 
        step_values_newt_n[2], step_values_newt_n[3]);
        verlet_newtonian(step_size, step_values_newt_n, step_values_newt_n_1,
        step_values_newt_new);
        
        double new_theta_rel = euler(step_size, step_values_rel_n[1], 
        step_values_rel_n[2], step_values_rel_n[3]);
        verlet_relativ(step_size, step_values_rel_n, step_values_rel_n_1,
        step_values_rel_new);
        
        values_newt[i][0] = step_values_newt_n[0] + step_size;
        values_newt[i][1] = step_values_newt_new[0];
        values_newt[i][2] = step_values_newt_new[1];
        values_newt[i][3] = new_theta_newt;
        
        energy[i][0] = step_values_rel_n[0] + step_size;
        energy[i][1] = (mu/2)*(pow(step_values_rel_new[1],2.0) + 
        pow(step_values_rel_new[0],2.0)
        *pow(theta_ode(step_values_rel_new[0]),2.0));
        
        values_rel[i][0] = step_values_rel_n[0] + step_size;
        values_rel[i][1] = step_values_rel_new[0];
        values_rel[i][2] = step_values_rel_new[1];
        values_rel[i][3] = new_theta_rel;
    }
    
    ofstream myfile ("Newtonian.txt");
    if (myfile.is_open())
    {
        for (int i = 0; i < row; i++) 
        {
            for (int j = 0; j < column; j++)
            {
                myfile << to_string(values_newt[i][j]) + " ";
            }
            myfile << "\n";
        }
    myfile.close();
    }
    
    ofstream myfile2 ("Relativistic.txt");
    if (myfile2.is_open())
    {
        for (int i = 0; i < row; i++) 
        {
            for (int j = 0; j < column; j++)
            {
                myfile2 << to_string(values_rel[i][j]) + " ";
            }
            myfile2 << "\n";
        }
    myfile2.close();
    }
    
    ofstream myfile3 ("Energy.txt");
    if (myfile3.is_open())
    {
        for (int i = 0; i < row; i++) 
        {
            for (int j = 0; j < 2; j++)
            {
                myfile3 << to_string(energy[i][j]) + " ";
            }
            myfile3 << "\n";
        }
    myfile3.close();
    }
    return 0;
    
}
