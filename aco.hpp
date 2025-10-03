#ifndef ACO_HPP
#define ACO_HPP

#include "irp.hpp"
#include "parameters.hpp" // <-- MUDANÇA
#include <vector>

// A struct ACO_Params foi movida para parameters.hpp

struct ACO_Result {
    double bestCost;
    std::vector<std::vector<int>> bestRoutes;
};

ACO_Result runACO_for_period(
    const IRP& irp,
    const std::vector<int>& deliveries_for_period,
    const ACO_Params& params,
    bool verbose = false
);

#endif // ACO_HPP