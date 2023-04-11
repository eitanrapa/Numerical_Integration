#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/*
Define ODEs below, only first order
*/

double ODE1 (double t, double y)
{
    double y_prime;
    y_prime = 2 - exp(-4*t) - 2*y;
    return y_prime;
}

/*
Define methods below
*/

void euler(double step_size, double *returnArray)
{
    double t0 = returnArray[0];
    double y0 = returnArray[1];
    double y = y0 + step_size*ODE1(t0,y0);
    returnArray[1] = y;
}

void runge_kutta(double step_size, double *returnArray)
{
    double t0 = returnArray[0];
    double y0 = returnArray[1];
    double k1 = ODE1(t0,y0);
    double k2 = ODE1(t0 + (step_size/2), y0 + step_size*(k1/2));
    double k3 = ODE1(t0 + (step_size/2), y0 + step_size*(k2/2));
    double k4 = ODE1(t0 + step_size, y0 + (step_size*k3));
    double y = y0 + (1.0/6.0)*step_size*(k1 + 2*k2 + 2*k3 + k4);
    returnArray[1] = y;
}

void verlet(double step_size, double *returnArray)
{
    double t0 = returnArray[0];
    double y0 = returnArray[1];
    double v0 = returnArray[2];
    double v_step = v0 + (step_size/2)*ODE1(t0,y0);
    double y = y0 + step_size*v_step;
    double v = v_step + (step_size/2)*ODE1(t0 + step_size, y);
    returnArray[1] = y;
    returnArray[2] = v;
}

/*
Main function below
*/

int main() 
{   
    /*
    Check in problem
    */
    double step_size = 0.1;
    int nsteps = 11;
    double t0 = 0.0 ;
    double y0 = 1.0;
    double v0 = ODE1(t0,y0);
    //Number of steps
    int row = nsteps;
    //Number of variables (x,y,v)
    int column = 3;
    double values[row][column] = {0.0};
    values[0][0] = t0;
    values[0][1] = y0;
    values[0][2] = v0; 
    for (int i = 1; i < nsteps; i++)
    {
        double step_values[3] = {values[i-1][0],values[i-1][1],values[i-1][2]};
        runge_kutta(step_size, step_values);
        values[i][0] = step_values[0] + step_size;
        values[i][1] = step_values[1];
        values[i][2] = step_values[2];
    }
    
    ofstream myfile ("Output.txt");
    if (myfile.is_open())
    {
        for (int i = 0; i < row; i++) 
        {
            for (int j = 0; j < column; j++)
            {
                myfile << to_string(values[i][j]) + " ";
            }
            myfile << "\n";
        }
    myfile.close();
    }
 
    return 0;
    
}


