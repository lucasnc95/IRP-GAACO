#ifndef UTILS_HPP
#define UTILS_HPP

#include <random>

// Gerador global de números aleatórios, declarado aqui e definido em utils.cpp
extern std::mt19937 rng;

// Funções utilitárias inline para números aleatórios
inline int randint(int a, int b) {
    return std::uniform_int_distribution<int>(a, b)(rng);
}

inline double randreal() {
    return std::uniform_real_distribution<double>(0.0, 1.0)(rng);
}

#endif // UTILS_HPP