#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "interpolation.h"

// Function to evaluate f(x) = 1 / (1 + 25 * x^2)
double f(double x) {
    return 1.0 / (1 + 25 * x * x);
}


int main() {
    std::ofstream output_file("Chebyshev_output.txt");

    // Define the values of n as required
    std::vector<int> n_values = {5, 10, 15, 20};
    int num_points = 100;
    double x_min = -1, x_max = 1;
    
    for (int n : n_values) {
        // Generate Chebyshev nodes and corresponding function values
        std::vector<double> x_i = chebyshevNodes(n);
        std::vector<double> y_i;
        for (double x : x_i) {
            y_i.push_back(f(x));
        }

        // Calculate divided differences
        std::vector<double> coeffs = dividedDifferences(x_i, y_i);

        // Print the results for plotting
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

    return 0;
}
