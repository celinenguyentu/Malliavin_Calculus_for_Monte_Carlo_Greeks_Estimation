//
//  monte_carlo.hpp
//  
//
//  Created by Céline Nguyen on 14/02/2020.
//

#ifndef monte_carlo_hpp
#define monte_carlo_hpp

#include <stdio.h>
#include <random>
#include <algorithm>
#include <cmath>
#include <utility>
#include <string>
#include <fstream>
#include <vector>

class Stats{
protected:
    std::vector<double> Values;
    std::vector<double> Mean;
    std::vector<double> Cumsum;
    std::vector<double> Var;
public:
    Stats(double x=0.): Values(), Mean(), Cumsum(), Var() {};
    std::vector<double> get_Values() const;
    std::vector<double> get_Mean() const;
    std::vector<double> get_Var() const;
    std::vector<double> get_uconfidence() const;
    std::vector<double> get_lconfidence() const;
    double get_MonteCarlo() const;
    friend Stats & operator+=(Stats &, double);
    friend std::ostream & operator << (std::ostream &, const Stats &);
    friend void Export(std::string, const Stats &);
};

template <class Statistique, class RandomVariable, class Measurement, class RNG>
void MonteCarlo(Statistique & res, RandomVariable & X, const Measurement & f, RNG & G, long unsigned int n){
    for (long unsigned i=0;i<n; i++){
        double W1 = X(G);
        double W2 = X(G);
        res += f(W1, W2);
    }
};


class Product {
protected:
    double S0;
    double sigma;
    double r;
    double K;
    int T;
public:
    Product(double S0_, double sigma_, double r_, double K_, int T_): S0(S0_), sigma(sigma_), r(r_), K(K_), T(T_) {};
    virtual double payoff(double N, double M, double theta=0.) const =0;
    virtual Stats price(std::mt19937 & G, long unsigned int n, double t=0.) const;
    virtual Stats delta(std::mt19937 & G, long unsigned int n, double t=0.) const;
    Stats gamma(std::mt19937 & G, long unsigned int n, double t=0.) const;
    Stats vega(std::mt19937 & G, long unsigned int n, double t=0.) const;
    virtual double DeltaWeight(double x, double y) const = 0;
    virtual double GammaWeight(double x, double y) const = 0;
    virtual double VegaWeight(double x, double y) const = 0;
    
};

class DigitalOption : public Product {
protected:
    double Q;
public:
    DigitalOption(double S0_, double sigma_, double r_, double K_,int T_, double Q_): Product(S0_, sigma_, r_, K_, T_), Q(Q_){};
    double payoff(double N, double M, double theta=0.) const ;
    double DeltaWeight(double x, double y) const;
    double GammaWeight(double x, double y) const;
    double VegaWeight(double x, double y) const;
};

class DigitalOptionInterval: public Product {
protected:
    double Q;
    double max_K;
public:
    DigitalOptionInterval(double S0_, double sigma_, double r_, double min_K_, int T_, double Q_, double max_K_): Product(S0_, sigma_, r_, min_K_, T_), Q(Q_), max_K(max_K_){};
    double payoff(double N, double M, double theta=0.) const ;
    double DeltaWeight(double x, double y) const;
    double GammaWeight(double x, double y) const;
    double VegaWeight(double x, double y) const;
};




/*


template <class Statistique, class RandomVariable, class Measurement, class RNG>
void MonteCarloExotic(Statistique & res, RandomVariable & X, const Measurement & f, RNG & G, long unsigned int n){
    for (long unsigned i=0;i<n; i++){
        double N = X(G);
        res += f(N, G);
    }
};



class AsianOption : public Product {
protected:
    int nb;
public:
    AsianOption(double S0_, double sigma_, double r_, double K_, int T_, int nb_): Product(S0_, sigma_, r_, K_, T_), nb(nb_) {};
    std::vector<double> Riemann(std::mt19937 & G) const;
    Stats price(std::mt19937 & G, long unsigned int n, double t=0.) const override;
    Stats delta(std::mt19937 & G, long unsigned int n, double t=0.) const override;
    double payoff(double N, double M, double theta=0.) const override;
    double DeltaWeight(double N, double approx0, double approx1, double approx2) const;
    double DeltaWeight(double x, double y) const override {return 0;};
    double GammaWeight(double x, double y) const override {return 0;};
    double VegaWeight(double x, double y) const override {return 0;};
};

*/


template <class Statistique, class Measurement>
void MonteCarloAsian(Statistique & res, const Measurement & f, long unsigned int n){
    for (long unsigned i=0;i<n; i++){
        res += f();
    }
};

class AsianOption {
protected:
    double S0;
    double sigma;
    double r;
    double K;
    int T;
    int nb;
    std::mt19937& G;
public:
    AsianOption(double S0_, double sigma_, double r_, double K_, int T_, int nb_, std::mt19937 & G_): S0(S0_), sigma(sigma_), r(r_), K(K_), T(T_), nb(nb_), G(G_){};
    double approx(int k) const;
    double payoff() const;
    Stats price(long unsigned int n, double t=0.) const;
    Stats delta(long unsigned int n, double t=0.) const;
    double DeltaWeight() const;
};

#endif /* monte_carlo_hpp */
