#include <iostream>
#include "MonteCarlo.hpp"
#include "pcpd_helper.hpp"
#include <fstream>

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cout << "Le nombre d'arguments attendu est 3" << std::endl;
    }
    std::ifstream file(argv[2]);
    if (!file.is_open())
    {
        std::cerr << "Error opening file" << std::endl;
        exit(1);
    }
    nlohmann::json json = nlohmann::json::parse(file);
    MonteCarlo *monte_carlo = new MonteCarlo(json);
    PnlMat *data = pnl_mat_create_from_file(argv[1]);
    double prix = 0.0;
    double prix_std_dev = 0.0;
    PnlVect *delta = pnl_vect_create_from_zero(monte_carlo->option->size);
    PnlVect *delta_std_dev = pnl_vect_create_from_zero(monte_carlo->option->size);
    double erreur_couverture = 0.0;
    monte_carlo->price_delta(prix, prix_std_dev, delta, delta_std_dev);
    monte_carlo->calculPAndL(data, erreur_couverture);

    HedgingResults res(prix, prix_std_dev, erreur_couverture);
    std::cout << res << std::endl;

    delete monte_carlo;
    pnl_vect_free(&delta);
    pnl_vect_free(&delta_std_dev);
    return 0;
}