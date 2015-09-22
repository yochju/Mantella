#pragma once

// Armadillo
#include <armadillo>

namespace mant {

  class RankingSupportVectorMachine {
  public:
    explicit RankingSupportVectorMachine(const arma::Mat<double>& parameters, arma::uword niter = 1000, double epsilon = 1.0, arma::Col<double> upperBound = -arma::ones(1), double kernelParameter = 1.0, double sigmaPow = 1.0) noexcept;
    void learn();
    
    void evaluatePoints(const arma::Mat<double>& points, const arma::Mat<double>& parameters, const arma::Col<double>& rankingDirection, double twoSigmaPow2);

    arma::Col<double>& getRankingDirection();
    
    arma::Col<double>& getFitnessForEvaluatedPoints();
    
    double getTwoSigmaPow2();

    void setMaximalNumberOfIterations(const unsigned long long maximalNumberOfIterations);

    void setUpperBound(const arma::Col<double>& upperBound);

    void setKernelParameter(const double kernelParameter);
    
    double getKernelParameter();
    
    void setIterations(const arma::uword iterations);

    arma::uword getIterations();

  protected:
    //RankSVMLearn.cpp Input Variables
    arma::Mat<double> parameters;
    arma::uword dimension; //N or nx
    arma::uword ntraining; //ntrain
    arma::uword niter = 1000;
    double epsilon = 1.0;
    arma::Col<double> upperBound; //Ci
    double kernelParameter = 1.0; //sigmaA
    double sigmaPow = 1.0;
    //invsqrtC not listed since we moved encoding outside of RankSVM

    //RankSVMLearn.cpp Output Variables
    arma::Col<double> rankingDirection; //optAlphas
    double twoSigmaPow2 = 0;

    //RankSVMLearn.cpp Variables Between Functions
    arma::uword nAlpha; //may need to be int, not sure if ntraining can be 0
    arma::Mat<double> K;
    
    //RankSVMFunc.cpp
    arma::Col<double> fitness; //Fit
    arma::Mat<double> points;
    double ntest;
  };
}