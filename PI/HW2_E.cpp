#include "interpolation.h"
#include <iostream>
#include <vector>

int main() {
    // 给定的数据
    std::vector<double> days = {0, 6, 10, 13, 17, 20, 28};
    std::vector<double> weightsSp1 = {6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7};
    std::vector<double> weightsSp2 = {6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89};
    
    // 计算差商表
    std::vector<double> dividedDifferencesSp1 = dividedDifferences(days, weightsSp1);
    std::vector<double> dividedDifferencesSp2 = dividedDifferences(days, weightsSp2);
    
    // 预测43天重量
    double t = 43;
    double predictedWeightSp1 = newtonInterpolation(t, days, dividedDifferencesSp1);
    double predictedWeightSp2 = newtonInterpolation(t, days, dividedDifferencesSp2);

    // 输出预测结果
    std::cout << "Predicted weight for Sp1 at day " << t << ": " << predictedWeightSp1 << " grams\n";
    std::cout << "Predicted weight for Sp2 at day " << t << ": " << predictedWeightSp2 << " grams\n";

    // 判断是否死亡
    if (predictedWeightSp1 <= 0) {
        std::cout << "Sp1 larvae might die after another 15 days.\n";
    } else {
        std::cout << "Sp1 larvae might survive after another 15 days.\n";
    }

    if (predictedWeightSp2 <= 0) {
        std::cout << "Sp2 larvae might die after another 15 days.\n";
    } else {
        std::cout << "Sp2 larvae might survive after another 15 days.\n";
    }

    return 0;
}