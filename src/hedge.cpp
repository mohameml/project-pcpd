#include <iostream>
#include "monte_carlo/monte_carlo.hpp"
#include "Utils/convert.hpp"
#include "monte_carlo/monte_carlo.hpp"

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cout << "Le nombre d'arguments attendu est 3" << std::endl;
    }

    // MonteCarlo *monte_carlo = convert_json_to_monte_carlo(argv[2]);
    // PnlMat *data = pnl_mat_create_from_file(argv[1]);
    // double prix;
    // double prix_std_dev;
    // double erreur_couverture;
    // monte_carlo->price(prix, prix_std_dev);
    // monte_carlo->calculPAndL(data, erreur_couverture);

    // HedgingResults res(prix, prix_std_dev, erreur_couverture);
    // std::cout << res << std::endl;

    // delete monte_carlo;

    return 0;
}