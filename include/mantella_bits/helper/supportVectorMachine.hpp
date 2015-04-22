namespace mant {

  class RankingSVM {
    public:
      explicit RankingSVM(arma::Mat<double> X_training, unsigned int nx, unsigned int ntraining, unsigned int niter, double epsilon, arma::Mat<double> invertedSquaredC, arma::Col<double> p_Ci,double sigma_A, double sigma_Pow, arma::Col<double> xmean) noexcept;
      
      
    protected:
      DistPow2_Euclidean(arma::Col<double> x1, arma::Col<double> x2);
      
      double sigmaA = 1;
      double sigmaPow = 1;
      double epsilon = 1;
      double
  };

  //
  // Implementation
  //

}