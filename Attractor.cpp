#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/*
Define ODEs below, only first order
*/

double strange(double *returnArray, int position)
{
    double sigma = 10.0;
    double rho = 28.0;
    double beta = 8.0/3.0;
    if (position == 0) 
    {
        return sigma*(returnArray[2]-returnArray[1]);
    }
    else if (position == 1)
    {
        return (returnArray[1]*(rho - returnArray[3]) - returnArray[2]);
    }
    else 
    {
        return (returnArray[1]*returnArray[2] - beta*returnArray[3]);
    }
}


/*
Define methods below
*/

void euler_attractor(double step_size, double *returnArray, int position)
{
    double y0 = returnArray[1 + position];
    double y = y0 + step_size*strange(returnArray, position);
    returnArray[1 + position] = y;
}

/*
Main function below
*/

int main() 
{   
    double step_size = 0.01;
    int nsteps = 1000;
    double x0 = 0.1;
    double y0 = 0.2;
    double z0 = 0.3;
    //Number of steps
    int row = nsteps;
    //Number of variables (t,x,y,z)
    int column = 4;
    double values[row][column] = {0.0};
    values[0][1] = x0;
    values[0][2] = y0; 
    values[0][3] = z0; 
    
    //Step 1
    for (int i = 1; i < nsteps; i++)
    {
        double xstep_values[4] = {values[i-1][0], values[i-1][1],
        values[i-1][2], values[i-1][3]};
        double ystep_values[4] = {values[i-1][0], values[i-1][1],
        values[i-1][2], values[i-1][3]};
        double zstep_values[4] = {values[i-1][0], values[i-1][1],
        values[i-1][2], values[i-1][3]};
        euler_attractor(step_size, xstep_values, 0);
        euler_attractor(step_size, ystep_values, 1);
        euler_attractor(step_size, zstep_values, 2);
        values[i][0] = xstep_values[0] + step_size;
        values[i][1] = xstep_values[1];
        values[i][2] = ystep_values[2];
        values[i][3] = zstep_values[3];
    }
    
    ofstream myfile ("Strange1.txt");
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
    
    //Step 2
    x0 = values[1][nsteps-1];
    y0 = values[2][nsteps-1];
    z0 = values[3][nsteps-1];
    values[row][column] = {0.0};
    values[0][1] = x0;
    values[0][2] = y0; 
    values[0][3] = z0;
    double values_pert[row][column] = {0.0};
    values_pert[0][1] = x0 + 0.00000001;
    values_pert[0][2] = y0 - 0.00000001; 
    values_pert[0][3] = z0 - 0.00000001;
    
    //Step 3
    for (int i = 1; i < nsteps; i++)
    {
        double xstep_values[4] = {values[i-1][0], values[i-1][1],
        values[i-1][2], values[i-1][3]};
        double ystep_values[4] = {values[i-1][0], values[i-1][1],
        values[i-1][2], values[i-1][3]};
        double zstep_values[4] = {values[i-1][0], values[i-1][1],
        values[i-1][2], values[i-1][3]};
        euler_attractor(step_size, xstep_values, 0);
        euler_attractor(step_size, ystep_values, 1);
        euler_attractor(step_size, zstep_values, 2);
        values[i][0] = xstep_values[0] + step_size;
        values[i][1] = xstep_values[1];
        values[i][2] = ystep_values[2];
        values[i][3] = zstep_values[3];
    }
    
    for (int i = 1; i < nsteps; i++)
    {
        double xstep_values[4] = {values_pert[i-1][0], values_pert[i-1][1],
        values_pert[i-1][2], values_pert[i-1][3]};
        double ystep_values[4] = {values_pert[i-1][0], values_pert[i-1][1],
        values_pert[i-1][2], values_pert[i-1][3]};
        double zstep_values[4] = {values_pert[i-1][0], values_pert[i-1][1],
        values_pert[i-1][2], values_pert[i-1][3]};
        euler_attractor(step_size, xstep_values, 0);
        euler_attractor(step_size, ystep_values, 1);
        euler_attractor(step_size, zstep_values, 2);
        values_pert[i][0] = xstep_values[0] + step_size;
        values_pert[i][1] = xstep_values[1];
        values_pert[i][2] = ystep_values[2];
        values_pert[i][3] = zstep_values[3];
    }
    
    double displacement[row][column] = {0.0};
    for (int i = 0; i < row; i++)
    {
        displacement[i][0] = values[i][0];
        for (int j = 1; j < column; j++)
        {
            displacement[i][j] = values_pert[i][j] - values[i][j];
        }
    }
    
    //Step 4
    
    double L_t[nsteps][2] = {0.0};
    double dt0 = sqrt(pow(displacement[0][1],2) + pow(displacement[0][2],2)
          + pow(displacement[0][3],2));
    for (int i = 0; i < nsteps; i++)
    {
        double dt = sqrt(pow(displacement[i][1],2) + pow(displacement[i][2],2)
          + pow(displacement[i][3],2));
        L_t[i][0] = displacement[i][0];
        L_t[i][1] = log(dt/dt0);
    }
     
    ofstream myfile2 ("Lyapunov.txt");
    if (myfile2.is_open())
    {
        for (int i = 0; i < row; i++) 
        {
            for (int j = 0; j < 2; j++)
            {
                myfile2 << to_string(L_t[i][j]) + " ";
            }
            myfile2 << "\n";
        }
        
    myfile2.close();
    }
    
    //Part 2 
    
    int e_steps = 9;
    double e_value = 1; 
    double n_e_matrix[e_steps][2] = {0.0};
    for (int l = 0; l < e_steps; l++)
    {
        int box_number = floor(120 / e_value);
        int box_matrix[box_number][box_number][box_number]= {0};
        
        for (int i = 0; i < row; i++)
        {
          
            int yvalue = floor(ceil(values[i][0] + 30.0)/e_value);
            int xvalue = floor(ceil(values[i][1] + 30.0)/e_value);
            int zvalue = floor(ceil(values[i][2])/e_value);
            box_matrix[xvalue][yvalue][zvalue] = 1;
        }
    
        int box_count = 0;
        for (int i = 0; i < box_number; i++)
        {
            for (int j = 0; j < box_number; j++)
            {
                for (int k = 0; k < box_number; k++)
                {
                    if (box_matrix[i][j][k] == 1)
                    {
                        box_count++;
                    }
                }
            }
        }
        double n_e = log(box_count);
        double e_e = fabs(log(e_value));
        n_e_matrix[l][1] = n_e;
        n_e_matrix[l][0] = e_e;
        e_value = e_value + 0.1;
    }
    
    ofstream myfile3 ("Strange_geometry.txt");
    if (myfile3.is_open())
    {
        for (int i = 0; i < e_steps; i++) 
        {
            for (int j = 0; j < 2; j++)
            {
                myfile3 << to_string(n_e_matrix[i][j]) + " ";
            }
            myfile3 << "\n";
        }
        
    myfile3.close();
    }
    
    return 0;
    
}
