#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <string>

// Function to evaluate f(x) = 1 / (1 + 25 * x^2)
double f(double x);

// Function to calculate divided differences for Newton interpolation
std::vector<double> dividedDifferences(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    std::vector<double> coeffs = y; // Start with y values
    for (int j = 1; j < n; ++j) {
        for (int i = n - 1; i >= j; --i) {
            coeffs[i] = (coeffs[i] - coeffs[i - 1]) / (x[i] - x[i - j]);
        }
    }
    return coeffs;
}

// Function to evaluate the Newton interpolation polynomial at a point x
double newtonInterpolation(double x, const std::vector<double>& x_points, const std::vector<double>& coeffs) {
    double result = coeffs[0];
    double product = 1.0;
    for (int i = 1; i < coeffs.size(); ++i) {
        product *= (x - x_points[i - 1]);
        result += coeffs[i] * product;
    }
    return result;
}

// Function to calculate Chebyshev nodes for interpolation
std::vector<double> chebyshevNodes(int n) {
    std::vector<double> nodes;
    for (int k = 0; k < n; ++k) {
        double node = cos(M_PI * (2 * k + 1) / (2 * n));
        nodes.push_back(node);
    }
    return nodes;
}


#endif // INTERPOLATION_H
