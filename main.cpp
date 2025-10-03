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
        std::cout << "Detectado formato de arquivo '.dat' (Archetti et al.). Lendo instância..." << std::endl;
        irp.readDataFromArchettiFile(instance_file);
    } else {
        std::cout << "Detectado formato de arquivo '.irp'. Lendo instância..." << std::endl;
        irp.readDataFromFile(instance_file);
    }
    
    std::cout << "Instância lida: " << instance_file << "\n";
    std::cout << "Semente aleatória: " << seed << "\n";
    
    std::cout << "--- Parâmetros do GA ---\n";
    std::cout << " População: " << ga_params.popSize << ", Gerações: " << ga_params.nGen 
              << ", Limite Estagnação: " << ga_params.stagnation_threshold << "\n"
              << " Crossover: " << ga_params.pCrossover << ", Mutação: " << ga_params.pMutation 
              << ", Torneio K: " << ga_params.tournamentK << "\n";
              
    std::cout << "--- Parâmetros do ACO ---\n";
    std::cout << " Formigas: " << aco_params.nAnts << ", Iterações: " << aco_params.nIter 
              << ", P(Busca Local): " << aco_params.pLocalSearch << "\n"
              << " Alpha: " << aco_params.alpha
              << ", Beta: " << aco_params.beta << ", Rho: " << aco_params.rho 
              << ", Q: " << aco_params.Q << "\n";

    Individual bestOverall = run_genetic_algorithm(irp, ga_params, aco_params, verbose_mode);

    std::cout << "\n=== Resumo da solução ===\n";
    std::cout << "Matriz de entregas da melhor solução:\n";
    printDeliveriesMatrix(bestOverall, irp);

    EvaluationResult final_result = simulate_and_evaluate(bestOverall, irp, aco_params);

    std::cout << "\n--- Análise de Custos ---\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Custo total roteamento: ................... " << final_result.routing_cost << "\n";
    std::cout << "  Custo total armazenagem (clientes): ...... " << final_result.customer_holding_cost << "\n";
    std::cout << "  Custo total armazenagem (depósito): ...... " << final_result.depot_holding_cost << "\n";
    std::cout << "  --------------------------------------------------\n";
    std::cout << "  Custo Operacional Total (Fitness Final): .. " << final_result.operational_cost() << "\n";
    std::cout << "  Solução é FACTÍVEL\n";

    std::cout << "\n--- Exemplo de Rotas por Período ---\n";
    for (int t = 0; t < irp.nPeriods; ++t) {
        const auto& del = bestOverall.deliveries[t];
        long sumD = std::accumulate(del.begin(), del.end(), 0L);
        std::cout << "Periodo " << t << ":\n";
        if (sumD == 0) {
            std::cout << "  sem entregas.\n";
        } else {
            ACO_Result res = runACO_for_period(irp, del, aco_params, false);
            for (size_t r = 0; r < res.bestRoutes.size(); ++r) {
                std::cout << "    Rota " << r << ": ";
                for (size_t i = 0; i < res.bestRoutes[r].size(); ++i) {
                    std::cout << res.bestRoutes[r][i] << (i + 1 < res.bestRoutes[r].size() ? " " : "");
                }
                std::cout << "\n";
            }
        }
    }
    
    exportAndPlotRoutes(irp, bestOverall, aco_params, "routes_data.txt", "plot_routes.py");
    
    std::cout << "\nExecução concluída.\n";
    return 0;
}