#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
using namespace std;

class Quadrature
{
    private:
        double lowBound;
        double uppBound;
        double solution;
        double error;

        double func_to_integrate(double val)
        {
            val = exp(val);
            return val;
        }; 

        double func_to_integrate_2(double val)
        {
            val = cos(2*val) * exp(sin(2*val));
            return val;
        };

        double func_to_integrate_3(double val)
        {
            val = tanh(val);
            return val;
        };

        double func_to_integrate_4(double val)
        {
            val = (pow(val,2) + 1) / val;
            return val;
        };

        double func_to_integrate_5(double val)
        {
            val = tanh(val);
            return val;
        };

        double trapezoidal_rule(double lowBound, double uppBound)
        {
            double f0, f1, val;

            f0 = func_to_integrate(lowBound);
            f1 = func_to_integrate(uppBound);
            val = 0.5 * (uppBound - lowBound) * (f0 + f1);

            return val;
        };

        double simpson_rule_1(double lowBound, double uppBound)
        {
            double f0, f1, f2, mid, val;

            f0 = func_to_integrate(lowBound);
            mid = 0.5 * (uppBound - lowBound) + lowBound;
            f1 = func_to_integrate(mid);
            f2 = func_to_integrate(uppBound);
            val = ((uppBound - lowBound) / 6) * (f0 + 4*f1 + f2);

            return val;
        };

        double midpoint_rule(double lowBound, double uppBound)
        {
            double f0, mid, val;

            mid = 0.5 * (uppBound - lowBound) + lowBound;
            f0 = func_to_integrate(mid);
            val = (uppBound - lowBound) * f0;

            return val;
        };

        double open_newton_cotes_two_points(double lowBound, double uppBound)
        {
            double f0, f1, diff, input, val;

            diff = (uppBound - lowBound);
            input = diff / 3 + lowBound;
            f0 = func_to_integrate(input);
            input = (2 * diff) / 3 + lowBound;
            f1 = func_to_integrate(input);
            val = 0.5 * (uppBound - lowBound) * (f0 + f1);

            return val;
        };

    public:
    // Overloaded Constructors
        Quadrature(double lowerBound, double UpperBound)
        {
            lowBound = lowerBound;
            uppBound = UpperBound;
        }; 

        double integrate(int numIntervals, string method)
        {
            double val = 0.0; 
            vector<double> mesh = generate_uniform_mesh(lowBound, uppBound, numIntervals);
            for (int i = 0; i < (int)mesh.size()-1; i++)
            {
                double partial_sum;
                if (method == "open newton-cotes one point")
                {
                    partial_sum = midpoint_rule(mesh[i], mesh[i+1]);
                }else if (method == "open newton-cotes two point")
                {
                    partial_sum = open_newton_cotes_two_points(mesh[i], mesh[i+1]);
                }else if (method == "closed newton-cotes two point")
                {
                    partial_sum = trapezoidal_rule(mesh[i], mesh[i+1]);
                }else if (method == "closed newton-cotes three point")
                {
                    partial_sum = simpson_rule_1(mesh[i], mesh[i+1]);
                }
                
                val += partial_sum;
            };
            return val;
        };

        vector<double> generate_uniform_mesh(double start, double end, int numpoints)
        {
            vector<double> mesh;
            double h = ( end - start ) / numpoints;
            // Generate mesh
            for (int i = 0; i < numpoints; i++)
            {
                float val = start + (h * i);
                mesh.push_back(val);
            };
            mesh.push_back(end);

            return mesh;
        };
};

int main()
{
    double lowerbound = 0;
    double upperbound = 1;
    int numIntervals = 3;
    double solution = 0.0;

    Quadrature Quadrature_1(lowerbound, upperbound);
    solution = Quadrature_1.integrate(numIntervals, "closed newton-cotes three point");
    cout << solution;
};