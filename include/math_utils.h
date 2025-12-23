#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <vector>

class MathHelper {
   public:
    static double calculateScalarProductOfVectors(
        const std::vector<double>& vec1, const std::vector<double>& vec2);

    static std::vector<double> subtractVectors(const std::vector<double>& vec1,
                                               const std::vector<double>& vec2);

    static std::vector<double> additionVectors(const std::vector<double>& vec1,
                                               const std::vector<double>& vec2);

    static std::vector<double> multiplyVectorByMultiplier(
        const std::vector<double>& vec1, double multiplier);

    static std::vector<double> calculateCrossProductOfVectors(
        const std::vector<double>& vec1, const std::vector<double>& vec2);

    static double calculateNormOfVector(const std::vector<double>& vec1);

    static std::vector<double> calculateNormalizedVector(
        const std::vector<double>& vec1);

    static double calculateOrientedTetrahedronVolume6(
        const std::vector<double>& a, const std::vector<double>& b,
        const std::vector<double>& c, const std::vector<double>& d);
};

#endif