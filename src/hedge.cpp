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

    MonteCarlo *monte_carlo = convert_json_to_monte_carlo(argv[2]);
    PnlMat *data = pnl_mat_create_from_file(argv[1]);
    double prix;
    double prix_std_dev;
    PnlVect *delta = pnl_vect_create_from_zero(monte_carlo->option->option_size);
    PnlVect *delta_std_dev = pnl_vect_create_from_zero(monte_carlo->option->option_size);
    double erreur_couverture;
    monte_carlo->price(prix, prix_std_dev, delta, delta_std_dev);
    monte_carlo->calculPAndL(data, erreur_couverture);

    HedgingResults res(prix, prix_std_dev, erreur_couverture);
    std::cout << res << std::endl;

    delete monte_carlo;
    pnl_vect_free(&delta);
    pnl_vect_free(&delta_std_dev);
    return 0;
}