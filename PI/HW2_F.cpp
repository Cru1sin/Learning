#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

// Struct for a 2D point
struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};

// Generate heart shape points
std::vector<Point> generateHeartPoints(int m) {
    std::vector<Point> points;
    for (int i = 0; i < m; ++i) {
        double theta = 2 * M_PI * i / m;
        double x = cos(theta);
        double y = (2.0 / 3.0) * (sqrt(fabs(x)) + sin(theta));
        points.emplace_back(x, y);
    }
    return points;
}

// Function to evaluate a cubic Bézier curve at a given t
Point evaluateBezier(const std::vector<Point>& controlPoints, double t) {
    int n = controlPoints.size() - 1;
    Point result(0, 0);

    for (int i = 0; i <= n; ++i) {
        // Compute binomial coefficient
        double binom = tgamma(n + 1) / (tgamma(i + 1) * tgamma(n - i + 1));
        double coeff = binom * pow(1 - t, n - i) * pow(t, i);
        result.x += coeff * controlPoints[i].x;
        result.y += coeff * controlPoints[i].y;
    }
    return result;
}

// Approximate heart shape using piecewise cubic Bézier curves
std::vector<Point> approximateHeartWithBezier(const std::vector<Point>& heartPoints, int segments = 100) {
    std::vector<Point> bezierCurvePoints;
    int numControlPoints = heartPoints.size();

    for (int i = 0; i < numControlPoints - 3; i += 3) {
        std::vector<Point> controlPoints = {heartPoints[i], heartPoints[i + 1], heartPoints[i + 2], heartPoints[i + 3]};
        for (int j = 0; j <= segments; ++j) {
            double t = static_cast<double>(j) / segments;
            bezierCurvePoints.push_back(evaluateBezier(controlPoints, t));
        }
    }
    return bezierCurvePoints;
}

// Save points to a file
void savePointsToFile(const std::vector<Point>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& point : points) {
            file << point.x << " " << point.y << "\n";
        }
        file.close();
        std::cout << "Data saved to " << filename << "\n";
    } else {
        std::cerr << "Failed to open file: " << filename << "\n";
    }
}

int main() {
    std::vector<int> ms = {10, 40, 160}; // Values of m

    for (int m : ms) {
        std::vector<Point> heartPoints = generateHeartPoints(m);
        std::vector<Point> bezierCurvePoints = approximateHeartWithBezier(heartPoints);
        
        // Save the points to a text file
        std::string filename = "heart_" + std::to_string(m) + ".txt";
        savePointsToFile(bezierCurvePoints, filename);
    }

    return 0;
}
