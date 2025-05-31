//
//  project_examples.cpp
//  
//
//  Created by Céline Nguyen on 14/02/2020.
//

#include <stdio.h>
#include "monte_carlo.hpp"
#include <random>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>

int main(){
    std::mt19937 G(time(NULL));
    {
    // OPTION DIGITALE
    std::cout << "Option Digitale" << std::endl;
    // Parameters:
    int T = 1;
    double S0 = 100;
    double sigma = 0.2;
    double r = 0.05;
    double min_K = 140;
    double Q = 1;
    int n = 10000;
    DigitalOption option(S0, sigma, r, min_K, T, Q);
    // Monte Carlo estimation:
    Stats MonteCarlo = option.price(G,n);
    Stats Delta = option.delta(G,n);
    Stats Gamma = option.gamma(G,n);
    Stats Vega = option.vega(G,n);
    std::cout << "Monte Carlo : " << MonteCarlo.get_MonteCarlo() << std::endl;
    std::cout << "Delta : " << Delta.get_MonteCarlo() << std::endl;
    std::cout << "Gamma : " << Gamma.get_MonteCarlo() << std::endl;
    std::cout << "Vega : " << Vega.get_MonteCarlo() << std::endl;
    // Export data:
    Export("DigitalPrice.dat", MonteCarlo);
    Export("DigitalDelta.dat", Delta);
    Export("DigitalGamma.dat", Gamma);
    Export("DigitalVega.dat", Vega);
    }
    {
    // OPTION CORRIDOR
    std::cout << "Option Corridor" << std::endl;
    // Parameters:
    int T = 1;
    double S0 = 100;
    double sigma = 0.2;
    double r = 0.1;
    double min_K = 100;
    double max_K = 110;
    double Q = 1;
    int n = 10000;
    DigitalOptionInterval option(S0, sigma, r, min_K, T, Q, max_K);
    // Monte Carlo estimation:
    Stats MonteCarlo = option.price(G,n);
    Stats Delta = option.delta(G,n);
    Stats Gamma = option.gamma(G,n);
    Stats Vega = option.vega(G,n);
    std::cout << "Monte Carlo : " << MonteCarlo.get_MonteCarlo() << std::endl;
    std::cout << "Delta : " << Delta.get_MonteCarlo() << std::endl;
    std::cout << "Gamma : " << Gamma.get_MonteCarlo() << std::endl;
    std::cout << "Vega : " << Vega.get_MonteCarlo() << std::endl;
    // Export data:
    Export("CorridorPrice.dat", MonteCarlo);
    Export("CorridorDelta.dat", Delta);
    Export("CorridorGamma.dat", Gamma);
    Export("CorridorVega.dat", Vega);
    }
    {
    // OPTION ASIATIQUE
    std::cout << "Option Asiatique" << std::endl;
    // Parameters:
    int T = 1;
    double S0 = 100;
    double sigma = 0.2;
    double r = 0.1;
    double K = 100;
    int nb = 150;
    int n = 10000;
    AsianOption option(S0, sigma, r, K, T, nb, G);
    // Monte Carlo estimation:
    Stats MonteCarlo = option.price(n);
    Stats Delta = option.delta(n);
    std::cout << "Monte Carlo : " << MonteCarlo.get_MonteCarlo() << std::endl;
    std::cout << "Delta : " << Delta.get_MonteCarlo() << std::endl;
    // Export data:
    Export("AsianPrice.dat", MonteCarlo);
    Export("AsianDelta.dat", Delta);
    }
    return 0;

}
