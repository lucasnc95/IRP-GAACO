#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP

#include "irp.hpp" // <-- MUDANÇA: Incluído para ter acesso às definições
#include <vector>

// Assinatura principal atualizada para receber mais informações
double improve_routes(
    std::vector<std::vector<int>>& routes,
    const IRP& irp,
    const std::vector<int>& deliveries_for_period,
    const std::vector<std::vector<double>>& dist,
    bool verbose = false
);

#endif // LOCAL_SEARCH_HPP