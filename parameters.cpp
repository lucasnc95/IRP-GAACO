#include "parameters.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm> // Para std::find_if
#include <cctype>    // Para std::isspace

// Função auxiliar para remover espaços em branco do início e fim de uma string
void trim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void load_parameters_from_file(const std::string& filename, GA_Params& ga_params, ACO_Params& aco_params) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "AVISO: Não foi possível abrir o arquivo de parâmetros '" << filename << "'. Usando valores padrão." << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::stringstream ss(line);
        std::string key, value;
        if (std::getline(ss, key, '=') && std::getline(ss, value)) {
            trim(key);
            trim(value);

 try {
                // Parâmetros do GA
                if (key == "popSize") ga_params.popSize = std::stoi(value);
                else if (key == "nGen") ga_params.nGen = std::stoi(value);
                else if (key == "pCrossover") ga_params.pCrossover = std::stod(value);
                else if (key == "pMutation") ga_params.pMutation = std::stod(value);
                else if (key == "tournamentK") ga_params.tournamentK = std::stoi(value);
                else if (key == "stagnation_threshold") ga_params.stagnation_threshold = std::stoi(value);
                else if (key == "crossover_max_tries") ga_params.crossover_max_tries = std::stoi(value);
                // Parâmetros do ACO
                else if (key == "nAnts") aco_params.nAnts = std::stoi(value);
                else if (key == "nIter") aco_params.nIter = std::stoi(value);
                else if (key == "alpha") aco_params.alpha = std::stod(value);
                else if (key == "beta") aco_params.beta = std::stod(value);
                else if (key == "rho") aco_params.rho = std::stod(value);
                else if (key == "Q") aco_params.Q = std::stod(value);
                else if (key == "pLocalSearch") aco_params.pLocalSearch = std::stod(value);
            } catch (const std::exception& e) {
                std::cerr << "AVISO: Erro ao ler o valor para a chave '" << key << "'. Valor '" << value << "' é inválido." << std::endl;
            }
        }
    }
    std::cout << "Parâmetros lidos de '" << filename << "'." << std::endl;
}