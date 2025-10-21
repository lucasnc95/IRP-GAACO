#include "irp.hpp"
#include "ga.hpp"
#include "parameters.hpp"
#include "utils.hpp"
#include "evaluation.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <numeric>

void print_usage(const char* prog_name) {
    std::cerr << "Uso: " << prog_name << " <arquivo_instancia> <arquivo_parametros> <semente_aleatoria> [--verbose]\n\n";
    std::cerr << "Argumentos:\n";
    std::cerr << "  <arquivo_instancia>    Caminho para o arquivo de instância IRP (.irp ou .dat).\n";
    std::cerr << "  <arquivo_parametros>   Caminho para o arquivo de configuração (.txt).\n";
    std::cerr << "  <semente_aleatoria>    Um número inteiro para a semente aleatória.\n";
    std::cerr << "  --verbose              (Opcional) Ativa logs de progresso.\n";
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }

    std::string instance_file = argv[1];
    std::string params_file = argv[2];
    unsigned int seed;
    try {
        seed = std::stoul(argv[3]);
    } catch (const std::exception& e) {
        std::cerr << "Erro: A semente aleatória '" << argv[3] << "' é inválida.\n";
        return 1;
    }
    bool verbose_mode = (argc > 4 && std::string(argv[4]) == "--verbose");
    
    rng.seed(seed);

    GA_Params ga_params;
    ACO_Params aco_params;
    load_parameters_from_file(params_file, ga_params, aco_params);
    
    IRP irp;
    std::string extension = instance_file.substr(instance_file.find_last_of(".") + 1);
    if (extension == "dat") {
        irp.readDataFromArchettiFile(instance_file);
    } else {
        irp.readDataFromFile(instance_file);
    }
    
    std::cout << "Instância lida: " << instance_file << "\n";
    std::cout << "Semente aleatória: " << seed << "\n";
    
    std::cout << "--- Parâmetros do GA ---\n";
    std::cout << " População: " << ga_params.popSize << ", Gerações: " << ga_params.nGen 
              << ", Limite Estagnação: " << ga_params.stagnation_threshold << "\n"
              << " Crossover: " << ga_params.pCrossover << ", Mutação: " << ga_params.pMutation 
              << ", Tentativas Crossover: " << ga_params.crossover_max_tries
              << ", Torneio K: " << ga_params.tournamentK << "\n";
              
    std::cout << "--- Parâmetros do ACO ---\n";
    std::cout << " Formigas: " << aco_params.nAnts << ", Iterações: " << aco_params.nIter 
              << ", P(Busca Local): " << aco_params.pLocalSearch << "\n"
              << " Alpha: " << aco_params.alpha
              << ", Beta: " << aco_params.beta << ", Rho: " << aco_params.rho 
              << ", Q: " << aco_params.Q << "\n";
    irp.printData();          

    Individual bestOverall = run_genetic_algorithm(irp, ga_params, aco_params, verbose_mode);

    std::cout << "\n=== Resumo da solução ===\n";
    std::cout << "Matriz de entregas da melhor solução:\n";
    printDeliveriesMatrix(bestOverall, irp);

    std::cout << "\n--- Análise de Custos ---\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Custo total roteamento: ................... " << bestOverall.routing_cost << "\n";
    std::cout << "  Custo total armazenagem (clientes): ...... " << bestOverall.customer_holding_cost << "\n";
    std::cout << "  Custo total armazenagem (depósito): ...... " << bestOverall.depot_holding_cost << "\n";
    std::cout << "  Penalidade (Estoque Final > 0): .......... " << bestOverall.final_inventory_penalty << "\n";
    std::cout << "  --------------------------------------------------\n";
    std::cout << "  Custo Total (Fitness Final): ............ " << bestOverall.fitness << "\n";
    std::cout << "  Solução é " << (bestOverall.is_feasible ? "FACTÍVEL" : "INFACTÍVEL") << "\n";

    std::cout << "\n--- Rotas da Melhor Solução por Período ---\n";
    for (int t = 0; t < irp.nPeriods; ++t) {
        std::cout << "Periodo " << t << ":\n";
        if (bestOverall.routes_per_period[t].empty()) {
            std::cout << "  sem entregas.\n";
        } else {
            for (size_t r = 0; r < bestOverall.routes_per_period[t].size(); ++r) {
                std::cout << "    Rota " << r << ": ";
                for (size_t i = 0; i < bestOverall.routes_per_period[t][r].size(); ++i) {
                    std::cout << bestOverall.routes_per_period[t][r][i] << (i + 1 < bestOverall.routes_per_period[t][r].size() ? " " : "");
                }
                std::cout << "\n";
            }
        }
    }
    
    exportAndPlotRoutes(irp, bestOverall, aco_params, "routes_data.txt", "plot_routes.py");
    
    std::cout << "\nExecução concluída.\n";
    return 0;
}