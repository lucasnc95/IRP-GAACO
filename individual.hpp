#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <vector>

struct Individual {
    std::vector<std::vector<int>> deliveries;
    double fitness;
};

#endif // INDIVIDUAL_HPP