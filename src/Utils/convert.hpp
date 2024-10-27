#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <nlohmann/json.hpp>
#include "../Option/option.hpp"
#include "../black_scholes_model/black_scholes_model.hpp"
#include "../monte_carlo/monte_carlo.hpp"
#include "../json_helper.hpp"
#include <iostream>

Option *convert_json_to_option(nlohmann::json json);

BlackScholesModel *convert_json_to_model(nlohmann::json json);

MonteCarlo *convert_json_to_monte_carlo(std::string file_path);

#endif