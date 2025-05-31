//
//  monte_carlo.cpp
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
#include <math.h>

// STATS
std::vector<double> Stats::get_Values() const {return Values;};

std::vector<double> Stats::get_Mean() const {return Mean;};

std::vector<double> Stats::get_Var() const {return Var;};

std::vector<double> Stats::get_uconfidence() const {
    std::vector<double> Uconfidence(Values.size());
    for (int i=0; i<Values.size(); i++){
        Uconfidence[i] = Mean[i]+(1.96*sqrt(Var[i]/(i+1)));
    }
    return Uconfidence;
};

std::vector<double> Stats::get_lconfidence() const {
    std::vector<double> Lconfidence(Values.size());
    for (int i=0; i<Values.size(); i++){
        Lconfidence[i] = Mean[i]-(1.96*sqrt(Var[i]/(i+1)));
    }
    return Lconfidence;
};

double Stats::get_MonteCarlo() const {return Mean.back();};

Stats & operator += (Stats & stat, double x){
    (stat.Values).push_back(x);
    int n = (stat.Values).size();
    if (n == 1){
        (stat.Mean).push_back(x);
        (stat.Cumsum).push_back(x*x);
        (stat.Var).push_back(0);
    }
    else{
        double newMean = ((stat.Mean).back()*(stat.Mean).size()+x)/n;
        (stat.Mean).push_back(newMean);
        double newCumsum = (stat.Cumsum).back()+x*x;
        (stat.Cumsum).push_back(newCumsum);
        double newVar = (stat.Cumsum).back()/n - (stat.Mean).back()*(stat.Mean).back();
        (stat.Var).push_back(newVar);
    }
    return stat;
}

std::ostream & operator << (std::ostream & flux, const Stats & stat){
    std::vector<double> uconf = stat.get_uconfidence();
    std::vector<double> lconf = stat.get_lconfidence();
    for (unsigned k=0; k<(stat.Mean).size(); k++){
        flux << k+1 << " " << (stat.Values)[k] << " " << (stat.Mean)[k] << " " << uconf[k] << " " << lconf[k] <<  std::endl;
    }
    return flux;
}

void Export(std::string s, const Stats & stat){
    std::ofstream fichier(s);
    std::vector<double> values = stat.get_Values();
    std::vector<double> mean = stat.get_Mean();
    std::vector<double> uconf = stat.get_uconfidence();
    std::vector<double> lconf = stat.get_lconfidence();
    for (int i=0; i<values.size(); i++){
        fichier << i+1 << " " << values[i] << " " << mean[i] << " " << uconf[i] << " " << lconf[i] << std::endl;
    }
    fichier.close();
};


// PRODUCT

Stats Product::price(std::mt19937 & G, long unsigned int n, double t) const {
    std::normal_distribution<double> N(0,1);
    auto Function_to_evaluate = [=](double x, double y){return exp(-r*T)*payoff(x, y);};
    Stats stat;
    MonteCarlo(stat, N, Function_to_evaluate, G, n);
    return stat;
};


Stats Product::delta(std::mt19937 & G, long unsigned int n, double t) const {
    std::normal_distribution<double> N(0,1);
    auto Function_to_evaluate = [=](double x, double y){return exp(-r*T)*payoff(x, y)*DeltaWeight(x,y);};
    Stats stat;
    MonteCarlo(stat, N, Function_to_evaluate, G, n);
    return stat;
};


Stats Product::gamma(std::mt19937 & G, long unsigned int n, double t) const {
    std::normal_distribution<double> N(0,1);
    auto Function_to_evaluate = [=](double x, double y){return exp(-r*T)*payoff(x, y)*GammaWeight(x, y);};
    Stats stat;
    MonteCarlo(stat, N, Function_to_evaluate, G, n);
    return stat;
};

Stats Product::vega(std::mt19937 & G, long unsigned int n, double t) const {
    std::normal_distribution<double> N(0,1);
    auto Function_to_evaluate = [=](double x, double y){return exp(-r*T)*payoff(x, y)*VegaWeight(x, y);};
    Stats stat;
    MonteCarlo(stat, N, Function_to_evaluate, G, n);
    return stat;
};




// DIGITAL OPTION

double DigitalOption::payoff(double N, double M, double theta) const {
    double ST = S0*exp((r-sigma*sigma/2)*T+sigma*(N+theta)*sqrt(T));
    return Q*(ST > K);
}

double DigitalOption::DeltaWeight(double x, double y) const {
    return (x*sqrt(T))/(S0*sigma*T);
}
    
double DigitalOption::GammaWeight(double x, double y) const {
    double WT = x*sqrt(T);
    double weight = ((WT*WT/sigma*T)-WT-(1/sigma))/(S0*S0*sigma*T);
    return weight;
}

double DigitalOption::VegaWeight(double x, double y) const {
    double WT = x*sqrt(T);
    double weight = (WT*WT/sigma*T)-WT-(1/sigma);
    return weight;
}

// CORRIDOR OPTION

double DigitalOptionInterval::payoff(double N, double M, double theta) const {
    double ST = S0*exp((r-sigma*sigma/2)*T+sigma*(N+theta)*sqrt(T));
    return Q*((ST >= K) & (ST<= max_K));
}

double DigitalOptionInterval::DeltaWeight(double x, double y) const {
    return (x*sqrt(T))/(S0*sigma*T);
}
    
double DigitalOptionInterval::GammaWeight(double x, double y) const {
    double WT = x*sqrt(T);
    double weight = ((WT*WT/sigma*T)-WT-(1/sigma))/(S0*S0*sigma*T);
    return weight;
}

double DigitalOptionInterval::VegaWeight(double x, double y) const {
    double WT = x*sqrt(T);
    double weight = (WT*WT/sigma*T)-WT-(1/sigma);
    return weight;
}




// ASIAN OPTION


/* ESSAI 2
std::vector<double> AsianOption::Riemann(std::mt19937 & G) const {
    std::vector<double> approx(4);
    std::normal_distribution<double> N(0,1);
    approx[0] = S0; approx[1] = S0; approx[2] = S0;
    double h = (double)T/(double)nb;
    double St = S0;
    for (int i=1; i<nb; i++){
        St = St*exp((r-sigma*sigma/2)*h+sigma*sqrt(h)*N(G));
        approx[0] += St;
        approx[1] += pow((double)i/(double)nb,1)*St;
        approx[2] += pow((double)i/(double)nb,2)*St;
    }
    approx[0]/=h;
    approx[1]/=h;
    approx[2]/=h;
    approx[3] = St*exp((r-sigma*sigma/2)*h+sigma*sqrt(h)*N(G));
    return approx;
};

Stats AsianOption::price(std::mt19937 & G, long unsigned int n, double t) const {
    std::normal_distribution<double> N(0,1);
    auto Function_to_evaluate = [=](double N, std::mt19937 & G){
        std::vector<double> approx = Riemann(G);
        return exp(-r*T)*payoff(approx[0],0);
        
    };
    Stats stat;
    MonteCarloExotic(stat, N, Function_to_evaluate, G, n);
    return stat;
}

Stats AsianOption::delta(std::mt19937 & G, long unsigned int n, double t) const {
    std::normal_distribution<double> N(0,1);
    auto Function_to_evaluate = [=](double N, std::mt19937 & G){
        std::vector<double> approx = Riemann(G);
        return exp(-r*T)*payoff(approx[0],0)*DeltaWeight(approx[3],approx[0], approx[1], approx[2])/S0;
    };
    Stats stat;
    MonteCarloExotic(stat, N, Function_to_evaluate, G,  n);
    return stat;
}

double AsianOption::payoff(double N, double M, double theta) const {
    if (N>K) return N-K;
    else return 0;
}

double AsianOption::DeltaWeight(double ST, double approx0, double approx1, double approx2) const {
    double WT = (log(ST/S0) - (r - sigma*sigma/2)*T)/sigma;
    return (approx0/approx1)*((WT/sigma + approx2/approx1)-1);
}
        
*/



double AsianOption::approx(int k) const {
    double h = (double)T/(double)nb;
    double sum = 0;
    double St = S0;
    std::normal_distribution<double> N(0,1);
    for (int i=0; i<nb; i++){
        St = St*exp((r-sigma*sigma/2)*h + sigma*sqrt(h)*N(G));
        sum += pow(((double)i/(double)nb), k) * St;
    }
    return sum*h;
};

double AsianOption::payoff() const {
    double Diff = approx(0) -K;
    if (Diff>0) return Diff;
    else return 0;
}

Stats AsianOption::price(long unsigned int n, double t) const {
    auto Function_to_evaluate = [=](){return exp(-r*T)*payoff();};
    Stats stat;
    MonteCarloAsian(stat, Function_to_evaluate, n);
    return stat;
};

Stats AsianOption::delta(long unsigned int n, double t) const {
    auto Function_to_evaluate = [=](){return exp(-r*T)*payoff()*DeltaWeight()/S0;};
    Stats stat;
    MonteCarloAsian(stat, Function_to_evaluate, n);
    return stat;
};

double AsianOption::DeltaWeight() const{
    std::normal_distribution<double> N(0,1);
    double approx0 = approx(0); double approx1 = approx(1); double approx2 = approx(2);
    return (approx0/approx1)*((sqrt(T)*N(G)/sigma + approx2/approx1)-1);
    
}


