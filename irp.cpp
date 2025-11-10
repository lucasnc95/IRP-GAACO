#include "irp.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>
#include <map>
#include <algorithm>
using std::vector;

void IRP::readDataFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::string aux;
    file >> aux >> nVehicles >> aux >> nDepots >> aux >> nCustomers >> aux >> nPeriods;

    getline(file, line);
    getline(file, line);
    file >> Vehicle_Type >> Fleet_Size >> Capacity;
    getline(file, line);
    getline(file, line);

    for (int i = 0; i < nDepots; ++i) {
        Depot depot;
        file >> depot.id >> depot.x >> depot.y >> depot.invCost >> depot.initialInv >> depot.minLevelInv >> depot.maxLevelInv;
        depot.production.resize(nPeriods);
        for (int j = 0; j < nPeriods; ++j)
            file >> depot.production[j];
        depots.push_back(depot);
    }

    getline(file, line);
    getline(file, line);

    for (int i = 0; i < nCustomers; ++i) {
        Customer customer;
        file >> customer.id >> customer.x >> customer.y >> customer.invCost >> customer.initialInv >> customer.minLevelInv >> customer.maxLevelInv;
        customer.demand.resize(nPeriods);
        for (int j = 0; j < nPeriods; ++j) {
            file >> customer.demand[j];
        }
        customers.push_back(customer);
    }

    getline(file, line);
    getline(file, line);
    getline(file, line);
    getline(file, line);

    costMatrix.resize(nDepots + nCustomers, std::vector<double>(nDepots + nCustomers, 0));
    int contAux = 0;
    for (int i = 0; i < nDepots + nCustomers; i++) {
        if (i > 0)
            file >> aux;

        for (int j = 0; j <= contAux; j++) {
            if (j == i)
                costMatrix[i][j] = 0;
            else {
                file >> costMatrix[i][j];
                costMatrix[j][i] = costMatrix[i][j];
            }
        }
        contAux++;
    }
}

void IRP::printData() const {
    std::cout << "=== IRP Instance Data ===\n";
    std::cout << "Vehicles: " << nVehicles
              << " Vehicle_Type: " << Vehicle_Type
              << " Fleet_Size: " << Fleet_Size
              << " Capacity: " << Capacity << "\n";
    std::cout << "Depots: " << nDepots << " Customers: " << nCustomers
              << " Periods: " << nPeriods << "\n\n";

    std::cout << "-- Depots --\n";
    for (const auto& d : depots) {
        std::cout << "Depot " << d.id
                  << " (x=" << d.x << ", y=" << d.y << ")"
                  << " invCost=" << d.invCost
                  << " initialInv=" << d.initialInv
                  << " min=" << d.minLevelInv
                  << " max=" << d.maxLevelInv
                  << "\n Production per period: ";
        for (int t = 0; t < nPeriods; ++t)
            std::cout << d.production[t] << (t + 1 < nPeriods ? ", " : "");
        std::cout << "\n";
    }
    std::cout << "\n-- Customers --\n";
    for (const auto& c : customers) {
        std::cout << "Customer " << c.id
                  << " (x=" << c.x << ", y=" << c.y << ")"
                  << " invCost=" << c.invCost
                  << " initialInv=" << c.initialInv
                  << " min=" << c.minLevelInv
                  << " max=" << c.maxLevelInv
                  << "\n Demand per period: ";
        for (int t = 0; t < nPeriods; ++t)
            std::cout << c.demand[t] << (t + 1 < nPeriods ? ", " : "");
        std::cout << "\n";
    }
    std::cout << "\n-- Cost Matrix (" << nDepots + nCustomers << " x " << nDepots + nCustomers << ") --\n";
    std::cout << std::fixed << std::setprecision(2);
    for (size_t i = 0; i < costMatrix.size(); ++i) {
        std::cout << std::setw(4) << i << ": ";
        for (size_t j = 0; j < costMatrix[i].size(); ++j) {
            std::cout << costMatrix[i][j]
                      << (j + 1 < costMatrix[i].size() ? " " : "");
        }
        std::cout << "\n";
    }
    std::cout << "=========================\n";
}

void IRP::readDataFromArchettiFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    
    depots.clear();
    customers.clear();
    costMatrix.clear();
    
    int n_locations;
    file >> n_locations >> this->nPeriods >> this->Capacity >> this->nVehicles;
    
    // <-- CORREÇÃO CRÍTICA: O número de clientes é o total de locais - 1 -->
    this->nCustomers = n_locations - 1;
    this->nDepots = 1;
    //this->nVehicles = 1; 
    this->Fleet_Size = nVehicles; 
    
    std::map<int, std::pair<double, double>> coords_map;

    // Lê os dados do depósito (supplier)
    Depot depot;
    int depot_id_from_file;
    int depot_prod_per_period;
    file >> depot_id_from_file >> depot.x >> depot.y >> depot.initialInv >> depot_prod_per_period >> depot.invCost;
    depot.id = depot_id_from_file;
    depot.minLevelInv = 0;
    depot.maxLevelInv = depot.initialInv * 100;
    depot.production.assign(this->nPeriods, depot_prod_per_period);
    depots.push_back(depot);
    coords_map[depot.id] = {depot.x, depot.y};

    // Lê os dados dos clientes (retailers)
    for (int i = 0; i < this->nCustomers; ++i) {
        Customer customer;
        int customer_demand_per_period;
        file >> customer.id >> customer.x >> customer.y 
             >> customer.initialInv >> customer.maxLevelInv >> customer.minLevelInv
             >> customer_demand_per_period >> customer.invCost;
        
        customer.demand.assign(this->nPeriods, customer_demand_per_period);
        customers.push_back(customer);
        coords_map[customer.id] = {customer.x, customer.y};
    }

    file.close();

    // --- Calcula a Matriz de Custos ---
    int matrixSize = 1 + this->nCustomers;
    costMatrix.resize(matrixSize, std::vector<double>(matrixSize));

    std::vector<int> index_to_id_map(matrixSize);
    index_to_id_map[0] = depots[0].id;
    for(int i = 0; i < nCustomers; ++i) {
        index_to_id_map[i+1] = customers[i].id;
    }

for (int i = 0; i < matrixSize; ++i) {
    for (int j = 0; j < matrixSize; ++j) {
        if (i == j) {
            costMatrix[i][j] = 0;
        } else {
            int id_i = index_to_id_map[i];
            int id_j = index_to_id_map[j];
            
            double xi = coords_map.at(id_i).first;
            double yi = coords_map.at(id_i).second;
            double xj = coords_map.at(id_j).first;
            double yj = coords_map.at(id_j).second;
            
            double dx = xi - xj;
            double dy = yi - yj;

            costMatrix[i][j] = static_cast<int>(std::sqrt(dx * dx + dy * dy) + 0.5);
        }
    }
}

}



