#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "interpolation.h"

// Function to evaluate f(x) = 1 / (1 + x^2)
double f(double x) {
    return 1.0 / (1 + x * x);
}

int main() {
    std::ofstream output_file("Newton_output.txt");

    // Define the range and step size for the x values in the interval [-5, 5]
    std::vector<int> n_values = {2, 4, 6, 8};
    int num_points = 100;
    double x_min = -5, x_max = 5;
    
    for (int n : n_values) {
        // Generate x_i points
        std::vector<double> x_i, y_i;
        for (int i = 0; i <= n; ++i) {
            double x_val = x_min + i * (x_max - x_min) / n;
            x_i.push_back(x_val);
            y_i.push_back(f(x_val));
        }

        // Calculate divided differences
        std::vector<double> coeffs = dividedDifferences(x_i, y_i);

        // Print the results to the file for plotting
        output_file << "n = " << n << ":\n";
        for (int i = 0; i <= num_points; ++i) {
            double x = x_min + i * (x_max - x_min) / num_points;
            double y_exact = f(x);
            double y_interp = newtonInterpolation(x, x_i, coeffs);
            output_file << std::fixed << std::setprecision(5) 
                        << "x: " << x << " Exact: " << y_exact << " Interpolated: " << y_interp << "\n";
        }
        output_file << "\n";
    }

    output_file.close();   

    return 0;
}
