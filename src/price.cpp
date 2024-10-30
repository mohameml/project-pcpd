#include <iostream>
#include "monte_carlo/monte_carlo.hpp"
#include "Utils/convert.hpp"
#include "monte_carlo/monte_carlo.hpp"
#include "pnl/pnl_vector.h"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << "Le nombre d'arguments attendu est 2" << std::endl;
        exit(1);
    }

    MonteCarlo *monte_carlo = convert_json_to_monte_carlo(argv[1]);

    double prix;
    double prix_std_dev;
    PnlVect *delta = pnl_vect_create_from_zero(monte_carlo->option->option_size);
    PnlVect *delta_std_dev = pnl_vect_create_from_zero(monte_carlo->option->option_size);

    monte_carlo->price(prix, prix_std_dev, delta, delta_std_dev);

    PricingResults res(prix, prix_std_dev, delta, delta_std_dev);
    std::cout << res << std::endl;

    delete monte_carlo;
    pnl_vect_free(&delta);
    pnl_vect_free(&delta_std_dev);

    return 0;
}
