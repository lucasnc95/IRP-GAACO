#ifndef IRP_HPP
#define IRP_HPP

#include <vector>
#include <string>

class Customer {
public:
    int id;
    double x, y;
    double invCost;
    int initialInv, minLevelInv, maxLevelInv;
    std::vector<int> demand;
};

class Depot {
public:
    int id;
    double x, y;
    double invCost;
    int initialInv, minLevelInv, maxLevelInv;
    std::vector<int> production;
};

class IRP {
public:
    int nVehicles, nDepots, nCustomers, nPeriods;
    int Vehicle_Type, Fleet_Size, Capacity;
    std::vector<Depot> depots;
    std::vector<Customer> customers;
    std::vector<std::vector<double>> costMatrix;

    void readDataFromFile(const std::string& filename);
    void readDataFromArchettiFile(const std::string& filename);
    void printData() const;
};

#endif // IRP_HPP
