#include <iostream>
#include "hermite.h"

int main() {
    std::vector<double> times = {0, 3, 5, 8, 13};
    std::vector<double> displacements = {0, 225, 383, 623, 993};
    std::vector<double> velocities = {75, 77, 80, 74, 72};

    std::vector<double> dividedDifferences;
    std::vector<double> timePoints;
    
    // 计算分差商
    computeDividedDifferences(times, displacements, velocities, dividedDifferences, timePoints);

    double t = 10;
    double predictedPosition = evaluateHermite(t, timePoints, dividedDifferences);
    double predictedSpeed = evaluateHermiteDerivative(t, timePoints, dividedDifferences);

    std::cout << "Predicted position at t = 10: " << predictedPosition << " feet\n";
    std::cout << "Predicted speed at t = 10: " << predictedSpeed << " feet/sec\n";

    // 检查是否超速
    bool speedLimitExceeded = false;
    for (double time = times.front(); time <= times.back(); time += 0.1) {
        if (evaluateHermiteDerivative(time, timePoints, dividedDifferences) > 81.0) {
            speedLimitExceeded = true;
            break;
        }
    }

    std::cout << "Did the car exceed 81 feet/sec? " << (speedLimitExceeded ? "Yes" : "No") << "\n";

    return 0;
}
