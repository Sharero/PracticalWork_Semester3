#include "../include/math_utils.h"

#include <cmath>

double MathHelper::calculateScalarProductOfVectors(
    const std::vector<double>& vec1, const std::vector<double>& vec2) {
    double result = 0.0;

#pragma unroll 4
    for (int i = 0; i < vec1.size(); i++) {
        result += vec1[i] * vec2[i];
    }

    return result;
}

std::vector<double> MathHelper::subtractVectors(
    const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> result(vec1.size());

#pragma unroll 4
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] - vec2[i];
    }

    return result;
}

std::vector<double> MathHelper::additionVectors(
    const std::vector<double>& vec1, const std::vector<double>& vec2) {
    std::vector<double> result(vec1.size());

#pragma unroll 4
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] + vec2[i];
    }

    return result;
}

std::vector<double> MathHelper::multiplyVectorByMultiplier(
    const std::vector<double>& vec1, double multiplier) {
    std::vector<double> result(vec1.size());

#pragma unroll 4
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] * multiplier;
    }

    return result;
}

std::vector<double> MathHelper::calculateCrossProductOfVectors(
    const std::vector<double>& vec1, const std::vector<double>& vec2) {
    return {(vec1[1] * vec2[2]) - (vec1[2] * vec2[1]),
            (vec1[2] * vec2[0]) - (vec1[0] * vec2[2]),
            (vec1[0] * vec2[1]) - (vec1[1] * vec2[0])};
}

double MathHelper::calculateNormOfVector(const std::vector<double>& vec1) {
    return std::sqrt(calculateScalarProductOfVectors(vec1, vec1));
}

std::vector<double> MathHelper::calculateNormalizedVector(
    const std::vector<double>& vec1) {
    double norm_of_vector = calculateNormOfVector(vec1);
    std::vector<double> result(vec1.size());

#pragma unroll 4
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] / norm_of_vector;
    }

    return result;
}

double MathHelper::calculateOrientedTetrahedronVolume6(
    const std::vector<double>& a, const std::vector<double>& b,
    const std::vector<double>& c, const std::vector<double>& d) {
    std::vector<double> u = subtractVectors(b, a);
    std::vector<double> v = subtractVectors(c, a);
    std::vector<double> w = subtractVectors(d, a);

    return (u[0] * (v[1] * w[2] - v[2] * w[1])) -
           (u[1] * (v[0] * w[2] - v[2] * w[0])) +
           (u[2] * (v[0] * w[1] - v[1] * w[0]));
}