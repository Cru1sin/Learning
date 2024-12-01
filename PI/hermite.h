#ifndef HERMITE_H
#define HERMITE_H

#include <vector>

#include "hermite.h"
#include <vector>

// 计算 Hermite 插值的分差商
void computeDividedDifferences(const std::vector<double>& times, 
                               const std::vector<double>& displacements, 
                               const std::vector<double>& velocities, 
                               std::vector<double>& dividedDifferences,
                               std::vector<double>& timePoints) {
    int n = times.size();
    timePoints.resize(2 * n);
    dividedDifferences.resize(2 * n);

    for (int i = 0; i < n; ++i) {
        timePoints[2 * i] = times[i];
        timePoints[2 * i + 1] = times[i];
        dividedDifferences[2 * i] = displacements[i];
        dividedDifferences[2 * i + 1] = displacements[i];
        
        // 使用速度初始化重复节点的分差商
        if (i > 0) {
            dividedDifferences[2 * i + 1] = velocities[i];
        }
    }

    // 计算分差商，处理重复节点的情况
    for (int j = 1; j < 2 * n; ++j) {
        for (int i = 2 * n - 1; i >= j; --i) {
            if ((timePoints[i] - timePoints[i - j])==0) {
                // 处理重复节点的情况，使用速度代替位置差分
                dividedDifferences[i] = velocities[i / 2];
            } else {
                dividedDifferences[i] = 
                    (dividedDifferences[i] - dividedDifferences[i - 1]) / (timePoints[i] - timePoints[i - j]);
            }
        }
    }
}

// 在给定点 t 处评估 Hermite 插值多项式
double evaluateHermite(double t, const std::vector<double>& timePoints, 
                       const std::vector<double>& dividedDifferences) {
    double result = dividedDifferences[0];
    double term = 1.0;
    int n = timePoints.size();
    
    for (int i = 1; i < n; ++i) {
        term *= (t - timePoints[i - 1]);
        result += dividedDifferences[i] * term;
    }
    
    return result;
}

// 在给定点 t 处评估 Hermite 插值多项式的导数
double evaluateHermiteDerivative(double t, const std::vector<double>& timePoints, 
                                 const std::vector<double>& dividedDifferences) {
    double derivative = 0.0;
    int n = timePoints.size();

    for (int i = 1; i < n; ++i) {
        double termDerivative = dividedDifferences[i];
        double product = 1.0;

        for (int j = 0; j < i; ++j) {
            termDerivative *= (t - timePoints[j]);
            if (j < i - 1) {
                product *= (t - timePoints[j]);
            }
        }

        derivative += termDerivative * product;
    }

    return derivative;
}


#endif
