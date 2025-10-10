#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#include "irp.hpp"
#include "individual.hpp"
#include "parameters.hpp"

// Função central que AVALIA e PREENCHE os dados de um indivíduo
void evaluate_and_fill(
    Individual& ind, // Passado por referência para ser modificado
    const IRP& irp, 
    const ACO_Params& aco_params
);

// Função de verificação de factibilidade
bool check_feasibility(const Individual& ind, const IRP& irp);

#endif // EVALUATION_HPP