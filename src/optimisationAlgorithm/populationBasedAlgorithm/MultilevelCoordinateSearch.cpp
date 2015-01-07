#include <hop_bits/optimisationAlgorithm/populationBasedAlgorithm/MultilevelCoordinateSearch.hpp>


////COMMENTS
//step1_ is used at a lot of places where matrix/column sizes arent calculable in advance (since size is dynamic in matlab).
//some variables have been renamed, there original name is noted in the header file (if they are renamed).
//there are a LOT of unexplained numbers in matlab where it's not clear if they have to be decremented to fit arrays starting at 0...
//numberOfIterations counting should be looked at, currently just called on every matlabs "ncall" increase
namespace hop {

  MultilevelCoordinateSearch::MultilevelCoordinateSearch(const std::shared_ptr<OptimisationProblem<double>> optimisationProblem,
      const unsigned int& populationSize, arma::Mat<double> boundaries, unsigned int boxDivisions, arma::Mat<double> hess, arma::Col<arma::uword> initialPointIndex, unsigned int maxLocalSearchSteps, double localStopThreshold)
  : PopulationBasedAlgorithm<double>(optimisationProblem, populationSize), boxDivisions_(boxDivisions), boundaries_(boundaries),
  maxLocalSearchSteps_(maxLocalSearchSteps), localStopThreshold_(localStopThreshold), hess_(hess), initialPointIndex_(initialPointIndex) {

    //assigning standard values for variables. Can't do in header-file since dependend on input variable "boundaries"
    hess_ = arma::ones(boundaries.col(0).n_elem);
    boxDivisions_ = 50 * boundaries.col(0).n_elem + 10;

    //for convenience
    unsigned int numberOfDimensions = optimisationProblem_->getNumberOfDimensions();

    //length of u or v and thus length of hess should equal numberOfDimensions
    if (numberOfDimensions != boundaries.col(0).n_elem) {
      std::cout << "lower boundaries dimensions don't match!" << std::endl;
    }
    if (numberOfDimensions != boundaries.col(1).n_elem) {
      std::cout << "upper boundaries dimensions don't match!" << std::endl;
    }
    if (numberOfDimensions != hess.row(0).n_elem && optimisationProblem_->getNumberOfDimensions() != hess.col(0).n_elem) {
      std::cout << "hess dimensions don't match!" << std::endl;
    }
    //Check if bounds are malformed
    for (std::size_t i = 0; i < numberOfDimensions; i++) {
      if (boundaries(i, 0) >= boundaries(i, 1)) {
        std::cout << "boundary malformed! u: " << boundaries(i, 0) << ", v:" << boundaries(i, 1) << std::endl;
      }
    }

    //init of large arrays
    //TODO: type completely unclear right now, pdf page 6f
    isplit_ = arma::Col<arma::uword>(step1_, arma::fill::zeros);
    level_ = arma::Col<arma::uword>(step1_, arma::fill::zeros);
    ipar_ = arma::Col<arma::uword>(step1_, arma::fill::zeros);
    ichild_ = arma::Col<arma::uword>(step1_, arma::fill::zeros);
    boxBaseVertexFunctionValues_ = arma::Mat<double>(2, step1_, arma::fill::zeros);
    z_ = arma::Mat<double>(2, step1_, arma::fill::zeros);
    nogain_ = arma::Col<arma::uword>(step1_, arma::fill::zeros);

    //initialization list
    //this is the equivalent of iinit = 0 in matlab
    //TODO: this is basically initialpopulation
    x0_.col(0) = boundaries_.col(0);
    x0_.col(1) = (boundaries_.col(0) + boundaries_.col(1)) / 2.0;
    x0_.col(2) = boundaries_.col(1);
    initListEvaluations_ = arma::Mat<double>(populationSize_, numberOfDimensions);

    //TODO: for custom initialisation lists there is a check here to see if they violate the boundaries


    //l_ L and x0_ are the custom initialisation list variables
    //l_ is supposed to point to the initial point x^0 in x0_ 
    //l_ also never gets changed in matlab as far as i could see
    //L gives the amount of predefined values per dimension (basically populationSize_ with more finetuning possible)
    //TODO: mcs.m does infinity check on x0_ here, not sure if needed

    arma::Col<double> initialPoint(numberOfDimensions, arma::fill::zeros);
    for (std::size_t i = 0; i < numberOfDimensions; i++) {
      initialPoint(i) = x0_(i, initialPointIndex_(i));
    }

    bestObjectiveValue_ = optimisationProblem_->getObjectiveValue(initialPoint);
    bestParameter_ = initialPoint;
    initListEvaluations_(0, initialPointIndex_(0)) = bestObjectiveValue_;

    for (std::size_t r = 0; r < numberOfDimensions; r++) {
      bestPointIndex_(r) = initialPointIndex_(r);
      for (std::size_t c = 0; c < populationSize_; c++) {
        if (c == initialPointIndex_(r)) {
          if (r != 1) {
            initListEvaluations_(r, c) = initListEvaluations_(r - 1, bestPointIndex_(r - 1));
          }
        } else {
          initialPoint(r) = x0_(r, c);
          initListEvaluations_(r, c) = optimisationProblem_->getObjectiveValue(initialPoint);
          numberOfIterations_++;
          if (initListEvaluations_(r, c) < bestObjectiveValue_) {
            bestObjectiveValue_ = initListEvaluations_(r, c);
            bestParameter_ = initialPoint;
            bestPointIndex_(r) = c;
          }
        }
      }
      initialPoint(r) = x0_(r, bestPointIndex_(r));
    }
    //in init.m all operations are done transposed (reason unknown yet), we simply do it at the end
    //TODO: About halfway through the code this still seems to make no sense
    initListEvaluations_ = initListEvaluations_.t();

    //base vertex and opposite vertex init
    for (std::size_t i = 0; i < numberOfDimensions; i++) {
      baseVertex_(i) = x0_(i, initialPointIndex_(i));

      //if true, use u, else use v
      if (std::abs(baseVertex_(i) - boundaries.col(0)(i)) > std::abs(baseVertex_(i) - boundaries.col(1)(i))) {
        oppositeVertex_(i) = boundaries.col(0)(i);
      } else {
        oppositeVertex_(i) = boundaries.col(1)(i);
      }
    }

    //init of record list, nboxes, nbasket,nbasket0,nsweep, m, nloc, xloc
    //"flag" is not needed since it is equal to optProblem->isFinished()
    //values not listed here are defined in header
    record_ = arma::Col<arma::uword>(boxDivisions_ - 1, arma::fill::zeros);
    m_ = numberOfDimensions;
    record_(0) = 1;
    xloc_ = arma::Mat<double>(maxLocalSearchSteps_, step1_);

    //generate boxes
    initBoxes();
  }

  void MultilevelCoordinateSearch::optimiseImplementation() {
    //for convenience
    unsigned int numberOfDimensions = optimisationProblem_->getNumberOfDimensions();

    double f0min = bestObjectiveValue_;
    //TODO: find better value than step1_...

    arma::Mat<double> pointsInBasket = arma::Mat<double>(numberOfDimensions, step1_); //xmin
    arma::Col<double> pointsInBasketValue = arma::Col<double>(step1_); //fmi

    // s- the vector record is updated, and the minimal level s containing non-split boxes is computed
    unsigned int minimalLevel = startSweep();

    //TODO: mcs.m checks for "minimalLevel < boxDivisions" as a while condition.
    //Which makes sense, since we cannot calculate anything if we cannot divide further.
    while (!isFinished() && !isTerminated() && minimalLevel < boxDivisions_) {
      unsigned int par = record_(minimalLevel); //the best box at level s is the current box
      //vertex.m START
      //TODO: check if populationSize_ comparisons need to be '-1'd
      arma::Col<double> x = arma::Col<double>(numberOfDimensions, arma::datum::inf);
      arma::Col<double> y = arma::Col<double>(numberOfDimensions, arma::datum::inf);
      arma::Col<double> x1 = arma::Col<double>(numberOfDimensions, arma::datum::inf);
      arma::Col<double> x2 = arma::Col<double>(numberOfDimensions, arma::datum::inf);
      arma::Col<double> f1 = arma::Col<double>(numberOfDimensions, arma::fill::zeros);
      arma::Col<double> f2 = arma::Col<double>(numberOfDimensions, arma::fill::zeros);
      arma::Col<arma::uword> n0 = arma::Col<arma::uword>(numberOfDimensions, arma::fill::zeros);
      double fold = boxBaseVertexFunctionValues_(0, par);
      unsigned int mVertex = par;
      while (mVertex > 1) {
        double i = std::abs(isplit_(ipar_(mVertex)));
        n0(i) = n0(i) + 1;
        //matlab checks for 1, since it's index use 0
        if (ichild_(mVertex) == 0) {
          if (x(i) == arma::datum::inf || x(i) == z_(0, ipar_(mVertex))) {
            //matlab passes 2, but it's used as an index so we need to use 1
            vert1(i, 1, mVertex, x, x1, x2, f1, f2);
          } else {
            updtf(numberOfDimensions, i, fold, x1, x2, f1, f2, boxBaseVertexFunctionValues_(0, ipar_(mVertex)));
            fold = boxBaseVertexFunctionValues_(0, ipar_(mVertex));
            //matlab passes 1, but it's used as an index so we need to use 0
            vert2(i, 0, mVertex, x, x1, x2, f1, f2);
          }
          //matlab checks for 2, since it's index use 1
        } else if (ichild_(mVertex) >= 1) {
          updtf(numberOfDimensions, i, fold, x1, x2, f1, f2, boxBaseVertexFunctionValues_(0, ipar_(mVertex)));
          fold = boxBaseVertexFunctionValues_(0, ipar_(mVertex));
          if (x(i) == arma::datum::inf || x(i) == z_(1, ipar_(mVertex))) {
            //matlab passes 1, but it's used as an index so we need to use 0
            vert1(i, 0, mVertex, x, x1, x2, f1, f2);
          } else {
            //matlab passes 2, but it's used as an index so we need to use 1
            vert2(i, 1, mVertex, x, x1, x2, f1, f2);
          }
        }
        //matlab checks for 1/2, since it's index use 0/1
        //original matlab code: 1 <= ichild(m) & ichild(m) <= 2 & y(i) == Inf
        if ((ichild_(mVertex) == 0 || ichild_(mVertex) == 1) && y(i) == arma::datum::inf) {
          y(i) = splitByGoldenSectionRule(z_(0, ipar_(mVertex)), z_(1, ipar_(mVertex)), boxBaseVertexFunctionValues_(0, ipar_(mVertex)), boxBaseVertexFunctionValues_(1, ipar_(mVertex)));
        }
        //box m was generated by splitting according to the init. list
        if (ichild_(mVertex) < 0) {
          int j1 = 0;
          int j2 = 0;
          int j3 = 0;
          int k = 0;
          if (boundaries_.col(0)(i) < x0_(i, 0)) {
            //TODO: since these are indexes too they also might need to be decremented
            //the if also needs to be adjusted. same in else below.
            j1 = std::ceil(std::abs(ichild_(mVertex)) / 2.0);
            j2 = std::floor(std::abs(ichild_(mVertex)) / 2.0);
            if ((std::abs(ichild_(mVertex)) / 2.0 < j1 && j1 > 1) || j1 == populationSize_) {
              j3 = -1;
            } else {
              //TODO: this probably needs to be 0. same below. 
              //but -1 should be correct (since that is a pseudoindex to recognize initlist boxes)
              j3 = 0;
            }
          } else {
            j1 = std::floor(std::abs(ichild_(mVertex)) / 2.0) + 1;
            j2 = std::ceil(std::abs(ichild_(mVertex)) / 2.0);
            if ((std::abs(ichild_(mVertex)) / 2.0 + 1) > j1 && j1 < populationSize_) {
              j3 = 0;
            } else {
              j3 = -1;
            }
          }
          //box m was generated in the init. procedure
          if (isplit_(ipar_(mVertex)) < 0) {
            k = i;
            //box m was generated by a later split according to the init.list
            //k points to the corresponding function values  
          } else {
            //matlab passes 2, but it's used as an index so we need to use 1
            k = z_(1, ipar_(mVertex));
          }
          if (j1 != initialPointIndex_(i) || (x(i) != arma::datum::inf && x(i) != x0_(i, (initialPointIndex_(i))))) {
            updtf(numberOfDimensions, i, fold, x1, x2, f1, f2, initListEvaluations_(initialPointIndex_(i), k));
            fold = initListEvaluations_(initialPointIndex_(i), k);
          }
          if (x(i) == arma::datum::inf || x(i) == x0_(i, j1)) {
            x(i) = x0_(i, j1);
            if (x1(i) == arma::datum::inf) {
              vert3(i, j1, mVertex, k, x1, x2, f1, f2);
            } else if (x2(i) == arma::datum::inf && x1(i) != x0_(i, j1 + j3)) {
              x2(i) = x0_(i, j1 + j3);
              f2(i) = f2(i) + initListEvaluations_(j1 + j3, k);
            } else if (x2(i) == arma::datum::inf) {
              //matlab checks for 1, since it's index use 0
              if (j1 != 0 && j1 != populationSize_) {
                x2(i) = x0_(i, j1 - j3);
                f2(i) = f2(i) + initListEvaluations_(j1 - j3, k);
              } else {
                x2(i) = x0_(i, j1 + 2 * j3);
                f2(i) = f2(i) + initListEvaluations_(j1 + 2 * j3, k);
              }
            }
          } else {
            if (x1(i) == arma::datum::inf) {
              x1(i) = x0_(i, j1);
              f1(i) = f1(i) + initListEvaluations_(j1, k);
              if (x(i) != x0_(i, j1 + j3)) {
                x2(i) = x0_(i, j1 + j3);
                f2(i) = f2(i) + initListEvaluations_(j1 + j3, k);
              }
            } else if (x2(i) == arma::datum::inf) {
              if (x1(i) != x0_(i, j1)) {
                x2(i) = x0_(i, j1);
                f2(i) = f2(i) + initListEvaluations_(j1, k);
              } else if (x(i) != x0_(i, j1 + j3)) {
                x2(i) = x0_(i, j1 + j3);
                f2(i) = f2(i) + initListEvaluations_(j1 + j3, k);
              } else {
                //matlab checks for 1, since it's index use 0
                if (j1 != 0 && j1 != populationSize_) {
                  x2(i) = x0_(i, j1 - j3);
                  f2(i) = f2(i) + initListEvaluations_(j1 - j3, k);
                } else {
                  x2(i) = x0_(i, j1 + 2 * j3);
                  f2(i) = f2(i) + initListEvaluations_(j1 + 2 * j3, k);
                }
              }
            }
          }
          if (y(i) == arma::datum::inf) {
            //TODO: why is there an index check for 0 in matlab?!!?
            if (j2 == 0) {
              y(i) = boundaries_.col(0)(i);
            } else if (j2 == populationSize_) {
              y(i) = boundaries_.col(1)(i);
            } else {
              y(i) = splitByGoldenSectionRule(x0_(i, j2), x0_(i, j2 + 1), initListEvaluations_(j2, k), initListEvaluations_(j2 + 1, k));
            }
          }
        }
        mVertex = ipar_(mVertex);
      }
      for (int i = 0; i < numberOfDimensions; i++) {
        if (x(i) == arma::datum::inf) {
          x(i) = x0_(i, initialPointIndex_(i));
          vert3(i, initialPointIndex_(i), mVertex, i, x1, x2, f1, f2);
        }
        if (y(i) == arma::datum::inf) {
          y(i) = oppositeVertex_(i);
        }
      }
      //vertex.m END

      bool doSplit = false; //splt
      if (minimalLevel > 2 * numberOfDimensions * (arma::min(n0) + 1)) {
        splitByRank(par, numberOfDimensions, n0);
        doSplit = true;
      } else {
        //TODO: this if should be unnecessary in c++. else be !if
        if (nogain_(par)) {
          doSplit = false;
        } else {
          arma::Col<double> expectedGain = expectedGainOfSplit(par, numberOfDimensions, n0, x1, x, f1, f2); //e
          //index again so use 0, matlab=1
          double fexp = boxBaseVertexFunctionValues_(0, par) + arma::min(expectedGain);
          if (fexp < bestObjectiveValue_) {
            doSplit = true;
          } else {
            doSplit = false;
            nogain_(par) = 1;
          }
        }
      }
      if (doSplit) {
        int i = isplit_(par);
        level_(par) = 0;
        //index again so use 1, matlab=2
        if (z_(1, par) == arma::datum::inf) {
          m_++;
          //index again so use 1, matlab=2
          z_(1, par) = m_;
          //little different than matlab, if this returns true = isFinished() | isTerminated()
          if (splitByInitList(i, minimalLevel, par, n0, x1, x2, pointsInBasket, pointsInBasketValue)) {
            break; //should break out of major while loop
          }
        } else {
          //index again so use 0, matlab=1
          z_(0, par) = baseVertex_(i);
          //little different than matlab, if this returns true = isFinished() | isTerminated()
          if (split(i, minimalLevel, par, n0, x1, x2, pointsInBasket, pointsInBasketValue)) {
            break;
          }
        }
        //if the pre-assigned size of the `large' arrays has already been exceeded, these arrays are made larger
        if(nboxes_ > dim) {
          //TODO: are the additional elements automatically set to zero? if not, need to do that
          isplit_.resize(nboxes_+step);
          level_.resize(nboxes_+step);
          ipar_.resize(nboxes_+step);
          ichild_.resize(nboxes_+step);
          z_.resize(z_.n_rows, nboxes_+step);
          nogain_.resize(nboxes_+step);
          boxBaseVertexFunctionValues_.resize(z_.n_rows, nboxes_+step);
          dim = nboxes_+step;
        }
        if(isFinished() || isTerminated()) {
          break;
        }
      } else {//no splitting, increase the level by 1
        if(minimalLevel+1 < boxDivisions_) {
          level_(par) = minimalLevel+1;
          //index again so use 0, matlab=1
          updateRecord(par,minimalLevel+1,boxBaseVertexFunctionValues_.row(0));
        } else {
          level_(par) = 0;
          nbasket_++;
          pointsInBasket.col(nbasket_) = baseVertex_;
          //index again so use 0, matlab=1
          pointsInBasketValue(nbasket_) = boxBaseVertexFunctionValues_(0, par);
        }
      }
      //of prepare for splitting
      minimalLevel++;
      while(minimalLevel < boxDivisions_) {
        if(record_(minimalLevel) == 0) {
          minimalLevel++;
        } else {
          break;
        }
      }
      //if smax is reached, a new sweep is started 
      if(minimalLevel==boxDivisions_) {
        if(maxLocalSearchSteps_ > 0) {
          //armadillo can't return both sort-index and sorted result, so we have to do it in two lines:
          arma::Col<arma::uword> j = arma::sort_index(pointsInBasketValue.rows(nbasket0_+1,nbasket_));
          pointsInBasketValue.rows(nbasket0_+1,nbasket_) = pointsInBasketValue.elem(j);
          
          
          for(arma::uword ind : j) {
            
          }
        }
      }
    }
  }

  void MultilevelCoordinateSearch::initBoxes() {
    //for convenience
    unsigned int numberOfDimensions = optimisationProblem_->getNumberOfDimensions();

    //parameter values of box 1
    ipar_(0) = 0;
    level_(0) = 1;
    ichild_(0) = 1;
    boxBaseVertexFunctionValues_(0, 0) = initListEvaluations_(initialPointIndex_(0), 0);

    int par = 0;

    arma::Col<double> var = arma::Col<double>(step1_);

    for (std::size_t i = 0; i < numberOfDimensions; i++) {
      isplit_(par) = -i; //boxes split in the init. procedure get a negative splitting index
      int nchild = 0;
      if (x0_(i, 1) > boundaries_.col(0)(i)) {
        nboxes_++;
        nchild++;
        genBox(nboxes_, par, level_(par) + 1, -nchild, initListEvaluations_(0, i));
      }
      //TODO: in matlab this is L(i) instead of populationSize,
      //seems kind of pointless since we don't have dynamic popSize per dimension so this always or never true
      //TODO: also, this might need to be 2 instead of 3
      double v1 = 0;
      if (populationSize_ == 3) {
        v1 = boundaries_.col(1)(i);
      } else {
        v1 = x0_(i, 2);
      }
      arma::Col<double> d = quadraticPolynomialInterpolation(x0_.submat(i, 0, i, 2), initListEvaluations_.submat(0, i, 2, i));
      double xl = minimumQuadraticPolynomial(boundaries_.col(0)(i), v1, d, x0_.submat(i, 0, i, 2));
      double fl = quadraticPolynomial(xl, d, x0_.submat(i, 0, i, 2));
      double xu = minimumQuadraticPolynomial(boundaries_.col(0)(i), v1, -d, x0_.submat(i, 0, i, 2));
      double fu = quadraticPolynomial(xu, d, x0_.submat(i, 0, i, 2));

      int par1 = 0;
      int j1 = 0; //auxiliary index
      if (bestPointIndex_(i) == 0) {
        if (xl < x0_(i, 0)) {
          par1 = nboxes_; //label of the current box for the next coordinate
        } else {
          par1 = nboxes_ + 1;
          j1 = 1;
        }
      }
      for (std::size_t j = 0; j < populationSize_ - 1; j++) {
        nboxes_++;
        nchild++;
        int s = 0;
        if (initListEvaluations_(j, i) <= initListEvaluations_(j + 1, i)) {
          s = 1;
        } else {
          s = 2;
        }
        genBox(nboxes_, par, level_(par) + s, -nchild, initListEvaluations_(j, i));
        if (j >= 1) {
          if (bestPointIndex_(i) == j) {
            if (xl <= x0_(i, j)) {
              par1 = nboxes_ - 1;
              j1 = j - 1;
            } else {
              par1 = nboxes_;
              j1 = j + 1;
            }
          }
          if (j <= populationSize_ - 2) {
            d = quadraticPolynomialInterpolation(x0_.submat(i, j, i, j + 2), initListEvaluations_.submat(j, i, j + 2, i));
            double u1 = 0;
            if (j < populationSize_ - 2) {
              u1 = x0_(i, j + 2);
            } else {
              u1 = boundaries_.col(1)(i);
            }
            xl = minimumQuadraticPolynomial(x0_(i, j), u1, d, x0_.submat(i, j, i, j + 2));
            fl = std::min(quadraticPolynomial(xl, d, x0_.submat(i, j, i, j + 2)), fl);
            xu = minimumQuadraticPolynomial(x0_(i, j), u1, -d, x0_.submat(i, j, i, j + 2));
            fu = std::max(quadraticPolynomial(xu, d, x0_.submat(i, j, i, j + 2)), fu);
          }
        }
        nboxes_++;
        nchild++;
        genBox(nboxes_, par, level_(par) + 3 - s, -nchild, initListEvaluations_(j + 1, i));
      }
      if (x0_(i, populationSize_) < boundaries_.col(0)(i)) {
        nboxes_++;
        nchild++;
        genBox(nboxes_, par, level_(par) + 1, -nchild, initListEvaluations_(populationSize_, i));
      }
      if (bestPointIndex_(i) == populationSize_) {
        if (x0_(i, populationSize_) < boundaries_.col(0)(i)) {
          if (xl <= x0_(i, populationSize_)) {
            par1 = nboxes_ - 1;
            j1 = populationSize_ - 1;
          } else {
            par1 = nboxes_;
            j1 = populationSize_ + 1;
          }
        } else {
          par1 = nboxes_;
          j1 = populationSize_ - 1;
        }
      }
      var(i) = fu - fl; // the quadratic model is taken as a crude measure of the variability in the ith component
      level_(par) = 0; //box is marked as split
      par = par1;
      //TODO: no idea what this splval in this next section is for. it is never used...
      double splval = 0;
      if (j1 == 0) {
        splval = boundaries_.col(0)(i);
      } else if (j1 == populationSize_ + 1) {
        splval = boundaries_.col(1)(i);
      } else {
        if (j1 < bestPointIndex_(i)) {
          splval = splitByGoldenSectionRule(x0_(i, j1), x0_(i, bestPointIndex_(i)), initListEvaluations_(j1, i), initListEvaluations_(bestPointIndex_(i), i));
        } else {
          splval = splitByGoldenSectionRule(x0_(i, bestPointIndex_(i)), x0_(i, j1), initListEvaluations_(bestPointIndex_(i), i), initListEvaluations_(j1, i));
        }
      }
    }
    //from matlab: best function value after the init. procedure
    //doesnt make much sense to me since we never changed initListEvaluations
    bestObjectiveValue_ = initListEvaluations_(bestPointIndex_(numberOfDimensions), numberOfDimensions);
    double var0 = 0;
    for (std::size_t i = 0; i < numberOfDimensions; i++) {
      //TODO: next two lines should equal [var0,p(i)] = max(var); not sure if correct
      var0 = arma::max(var);
      variabilityRanking_(i) = var0;
      var(var0) = -1;
      bestParameter_(i) = x0_(i, bestPointIndex_(i)); //from matlab: best point after the init. procedure
    }
  }

  std::string MultilevelCoordinateSearch::to_string() const noexcept {
    return "MultilevelCoordinateSearch";
  }

  void MultilevelCoordinateSearch::genBox(int nbox, int par, int level, int nchild, double baseVertexFunctionValue) {
    ipar_(nbox) = par;
    level_(nbox) = level;
    //TODO: nchild probably needs to be decremented by 1, or all calls to genBox's nchild
    ichild_(nbox) = nchild;
    boxBaseVertexFunctionValues_(0, nbox) = baseVertexFunctionValue;
  }

  arma::Col<double> MultilevelCoordinateSearch::quadraticPolynomialInterpolation(arma::Col<double> supportPoints, arma::Col<double> functionValues) {
    arma::Col<double> d(3);
    d(0) = functionValues(0);
    d(1) = (functionValues(1) - functionValues(0)) / (supportPoints(1) - supportPoints(0));
    double f23 = (functionValues(2) - functionValues(1)) / (supportPoints(2) - supportPoints(1));
    d(2) = (f23 - d(1)) / (supportPoints(2) - supportPoints(0));
    return d;
  }

  double MultilevelCoordinateSearch::minimumQuadraticPolynomial(double a, double b, arma::Col<double> d, arma::Mat<double> x0_fragment) {
    double x = 0;
    if (d(2) == 0) {
      if (d(1) == 0) {
        x = a;
      } else {
        x = b;
      }
      return x;
    } else if (d(2) > 0) {
      double x1 = 0.5 * (x0_fragment(0)) + x0_fragment(1) - 0.5 * (d(1) / d(2));
      if (a <= x1 && x1 <= b) {
        x = x1;
        return x;
      }
    }
    if (quadraticPolynomial(a, d, x0_fragment) < quadraticPolynomial(b, d, x0_fragment)) {
      x = a;
    } else {
      x = b;
    }
    return x;
  }

  double MultilevelCoordinateSearch::quadraticPolynomial(double x, arma::Col<double> d, arma::Mat<double> x0_fragment) {
    return d(0) + d(1)*(x - x0_fragment(0)) + d(2)*(x - x0_fragment(0))*(x - x0_fragment(1));
  }

  double MultilevelCoordinateSearch::splitByGoldenSectionRule(double x1, double x2, double f1, double f2) {
    if (f1 <= f2) {
      return x1 + 0.5 * (-1 + std::sqrt(5))*(x2 - x1);
    } else {
      return x1 + 0.5 * (3 - std::sqrt(5))*(x2 - x1);
    }
  }

  void MultilevelCoordinateSearch::splitByRank(unsigned int par, unsigned int numberOfDimensions, arma::Col<arma::uword> n0) {
    //index again so use 0, matlab=1 for all 3 variables
    int isplit = 0;
    int n1 = n0(0);
    int p1 = variabilityRanking_(0);
    //matlab starts at 2 obviously
    for (int i = 1; i < numberOfDimensions; i++) {
      if (n0(i) < n1 || (n0(i) == n1 && variabilityRanking_(i) < p1)) {
        isplit = i;
        n1 = n0(i);
        p1 = variabilityRanking_(i);
      }
    }
    if (n1 > 0) {
      //index again so use 1, matlab=2
      z_(1, par) = splitBySubint(baseVertex_(isplit), oppositeVertex_(isplit));
    } else {
      z_(1, par) = arma::datum::inf;
    }
  }

  double MultilevelCoordinateSearch::splitBySubint(double x, double y) {
    double x2 = y;
    if (x == 0 && std::abs(y) > 1000) {
      x2 = std::copysign(1.0, y);
    } else if (x != 0 && std::abs(y) > 100 * std::abs(x)) {
      //TODO: c++ standardlibraries have no signum. wat. using copysign instead...
      //original matlab: x2 = 10.*sign(y)*abs(x);
      x2 = 10 * std::copysign(x, y);
    }
    return x + 2 * (x2 - x) / 3.0;
  }

  bool MultilevelCoordinateSearch::splitByInitList(unsigned int splittingIndex, unsigned int minimalLevel, unsigned int par, arma::Col<arma::uword> n0, arma::Col<double> x1, arma::Col<double> x2, arma::Mat<double> pointsInBasket, arma::Col<double> pointsInBasketValue) {
    initListEvaluations_.col(m_).zeros();
    for (int j = 0; j < populationSize_; j++) {
      if (j != initialPointIndex_(splittingIndex)) {
        //TODO: why are we writing into the baseVertex? oO
        baseVertex_(splittingIndex) = x0_(splittingIndex, j);
        initListEvaluations_.col(m_)(j) = optimisationProblem_->getObjectiveValue(baseVertex_);
        numberOfIterations_++;
        if (initListEvaluations_.col(m_)(j) < bestObjectiveValue_) {
          bestObjectiveValue_ = initListEvaluations_.col(m_)(j);
          bestParameter_ = baseVertex_;
        }
        //In matlab this if is a little different and inside the if directly before this.
        //Our stopping conditions are a little different so it's here.
        if (isFinished() || isTerminated()) {
          return true;
        }
      } else {
        //index again so use 0, matlab=1
        initListEvaluations_.col(m_)(j) = boxBaseVertexFunctionValues_(0, par);
      }
    }
    double fm = arma::min(initListEvaluations_.col(m_));
    double i1 = fm;

    //splval again, still seems to serve no function except for printing...
    double splval1 = 0;
    double splval2 = 0;
    //this shouldn't be an index, so it's left at 1
    if (i1 > 1) {
      //TODO: might need to be decremented one more
      splval1 = splitByGoldenSectionRule(x0_(splittingIndex, i1 - 1), x0_(splittingIndex, i1), initListEvaluations_.col(m_)(i1 - 1), initListEvaluations_.col(m_)(i1));
    } else {
      splval1 = boundaries_.col(0)(splittingIndex);
    }
    if (i1 < populationSize_) {
      splval2 = splitByGoldenSectionRule(x0_(splittingIndex, i1), x0_(splittingIndex, i1 + 1), initListEvaluations_.col(m_)(i1), initListEvaluations_.col(m_)(i1 + 1));
    } else {
      splval2 = boundaries_.col(1)(splittingIndex);
    }
    if (minimalLevel + 1 < boxDivisions_) {
      int nchild = 0;
      //index again so use 0, matlab=1
      if (boundaries_.col(0)(splittingIndex) < x0_(splittingIndex, 0)) {//in that case the box at the boundary gets level s + 1
        nchild++;
        nboxes_++;
        //index again so use 0, matlab=1
        genBox(nboxes_, par, minimalLevel + 1, -nchild, initListEvaluations_.col(m_)(0));
        //index again so use 0, matlab=1
        updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
      }
      for (int j = 0; j < populationSize_ - 1; j++) {
        nchild++;
        //splval again, still seems to serve no function except for printing...
        double splval = splitByGoldenSectionRule(x0_(splittingIndex, i1), x0_(splittingIndex, i1 + 1), initListEvaluations_.col(m_)(i1), initListEvaluations_.col(m_)(i1 + 1));
        int level0 = 0;
        if (initListEvaluations_.col(m_)(j) <= initListEvaluations_.col(m_)(j + 1) || minimalLevel + 2 < boxDivisions_) {
          nboxes_++;
          if (initListEvaluations_.col(m_)(j) <= initListEvaluations_.col(m_)(j)) {
            level0 = minimalLevel + 1;
          } else {
            level0 = minimalLevel + 2;
          }
          genBox(nboxes_, par, level0, -nchild, initListEvaluations_.col(m_)(j));
          //index again so use 0, matlab=1
          updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
        } else {
          baseVertex_(splittingIndex) = x0_(splittingIndex, j);
          nbasket_++;
          pointsInBasket.col(nbasket_) = baseVertex_;
          pointsInBasketValue(nbasket_) = initListEvaluations_.col(m_)(j);
        }
        nchild++;
        if (initListEvaluations_.col(m_)(j + 1) < initListEvaluations_.col(m_)(j) || minimalLevel + 2 < boxDivisions_) {
          nboxes_++;
          if (initListEvaluations_.col(m_)(j + 1) < initListEvaluations_.col(m_)(j)) {
            level0 = minimalLevel + 1;
          } else {
            level0 = minimalLevel + 2;
          }
          genBox(nboxes_, par, level0, -nchild, initListEvaluations_.col(m_)(j + 1));
          //index again so use 0, matlab=1
          updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
        } else {
          baseVertex_(splittingIndex) = x0_(splittingIndex, j + 1);
          nbasket_++;
          pointsInBasket.col(nbasket_) = baseVertex_;
          pointsInBasketValue(nbasket_) = initListEvaluations_.col(m_)(j + 1);
        }
      }
      if (x0_(splittingIndex, populationSize_) < boundaries_.col(1)(splittingIndex)) {//in that case the box at the boundary gets level s + 1
        nchild++;
        nboxes_++;
        genBox(nboxes_, par, minimalLevel + 1, -nchild, initListEvaluations_.col(m_)(populationSize_));
        //index again so use 0, matlab=1
        updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
      }
    } else {
      for (int j = 0; j < populationSize_; j++) {
        baseVertex_(splittingIndex) = x0_(splittingIndex, j);
        nbasket_++;
        pointsInBasket.col(nbasket_) = baseVertex_;
        pointsInBasketValue(nbasket_) = initListEvaluations_.col(m_)(j);
      }
    }
  }

  bool MultilevelCoordinateSearch::split(unsigned int splittingIndex, unsigned int minimalLevel, unsigned int par, arma::Col<arma::uword> n0, arma::Col<double> x1, arma::Col<double> x2, arma::Mat<double> pointsInBasket, arma::Col<double> pointsInBasketValue) {
    //index again so use 1, matlab=2
    baseVertex_(splittingIndex) = z_(1);
    //index again so use 1, matlab=2
    boxBaseVertexFunctionValues_(1, par) = optimisationProblem_->getObjectiveValue(baseVertex_);
    numberOfIterations_++;
    //index again so use 1, matlab=2
    if (boxBaseVertexFunctionValues_(1, par) < bestObjectiveValue_) {
      //index again so use 1, matlab=2
      bestObjectiveValue_ = boxBaseVertexFunctionValues_(1, par);
      bestParameter_ = baseVertex_;
      //Our stopping conditions are a little different
      if (isFinished() || isTerminated()) {
        return true;
      }
    }
    //splval again...
    //index again, all decremented by 1
    double splval = splitByGoldenSectionRule(z_(par, 0), z_(par, 1), boxBaseVertexFunctionValues_(0, par), boxBaseVertexFunctionValues_(1, par));
    if (minimalLevel + 1 < boxDivisions_) {
      //index again, all decremented by 1
      if (boxBaseVertexFunctionValues_(0, par) < boxBaseVertexFunctionValues_(1, par)) {
        nboxes_++;
        //index again so use 0, matlab=1
        genBox(nboxes_, par, minimalLevel + 1, 1, boxBaseVertexFunctionValues_(0, par));
        //index again so use 0, matlab=1
        updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
        if (minimalLevel + 2 < boxDivisions_) {
          nboxes_++;
          //index again so use 1, matlab=2
          genBox(nboxes_, par, minimalLevel + 2, 2, boxBaseVertexFunctionValues_(1, par));
          //index again so use 0, matlab=1
          updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
        } else {
          //index again so use 1, matlab=2
          baseVertex_(splittingIndex) = z_(par, 1);
          nbasket_ = nbasket_ + 1;
          pointsInBasket.col(nbasket_) = baseVertex_;
          //index again so use 1, matlab=2
          pointsInBasketValue(nbasket_) = boxBaseVertexFunctionValues_(1, par);
        }

      } else {
        if (minimalLevel + 2 < boxDivisions_) {
          nboxes_++;
          //index again so use 0, matlab=1
          genBox(nboxes_, par, minimalLevel + 2, 1, boxBaseVertexFunctionValues_(0, par));
          //index again so use 0, matlab=1
          updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
        } else {
          //index again so use 0, matlab=1
          baseVertex_(splittingIndex) = z_(par, 0);
          nbasket_ = nbasket_ + 1;
          pointsInBasket.col(nbasket_) = baseVertex_;
          //index again so use 0, matlab=1
          pointsInBasketValue(nbasket_) = boxBaseVertexFunctionValues_(0, par);
        }
        nboxes_++;
        //index again so use 1, matlab=2
        genBox(nboxes_, par, minimalLevel + 1, 2, boxBaseVertexFunctionValues_(1, par));
        //index again so use 0, matlab=1
        updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
      }

      // if the third box is larger than the smaller of the other two boxes,
      // it gets level s + 1; otherwise it gets level s + 2
      //index again so use 1, matlab=2
      if (z_(par, 1) != oppositeVertex_(splittingIndex)) {
        //index again so use 1, matlab=2
        if (std::abs(z_(par, 1) - oppositeVertex_(splittingIndex)) > std::abs(z_(par, 1))*(3 - std::sqrt(5)*0.5)) {
          nboxes_++;
          //index again so use 1, matlab=2
          genBox(nboxes_, par, minimalLevel + 1, 3, boxBaseVertexFunctionValues_(1, par));
          //index again so use 0, matlab=1
          updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
        } else {
          if (minimalLevel + 2 < boxDivisions_) {
            nboxes_++;
            //index again so use 1, matlab=2
            genBox(nboxes_, par, minimalLevel + 2, 3, boxBaseVertexFunctionValues_(1, par));
            //index again so use 0, matlab=1
            updateRecord(nboxes_, level_(nboxes_), boxBaseVertexFunctionValues_.row(0));
          } else {
            //index again so use 1, matlab=2
            baseVertex_(splittingIndex) = z_(par, 1);
            nbasket_ = nbasket_ + 1;
            pointsInBasket.col(nbasket_) = baseVertex_;
            //index again so use 1, matlab=2
            pointsInBasketValue(nbasket_) = boxBaseVertexFunctionValues_(1, par);
          }
        }
      }
    } else {
      //index again so use 0, matlab=1
      baseVertex_(splittingIndex) = z_(par, 0);
      nbasket_ = nbasket_ + 1;
      pointsInBasket.col(nbasket_) = baseVertex_;
      //index again so use 0, matlab=1
      pointsInBasketValue(nbasket_) = boxBaseVertexFunctionValues_(0, par);

      //index again so use 1, matlab=2
      baseVertex_(splittingIndex) = z_(par, 1);
      nbasket_ = nbasket_ + 1;
      pointsInBasket.col(nbasket_) = baseVertex_;
      //index again so use 1, matlab=2
      pointsInBasketValue(nbasket_) = boxBaseVertexFunctionValues_(1, par);
    }
  }

  arma::Col<double> MultilevelCoordinateSearch::expectedGainOfSplit(unsigned int par, unsigned int numberOfDimensions, arma::Col<arma::uword> n0, arma::Col<double> x1, arma::Col<double> x2, arma::Col<double> f1, arma::Col<double> f2) {
    double emin = arma::datum::inf;
    arma::Col<double> expectedGain = arma::Col<double>(numberOfDimensions);
    for (int i = 0; i < numberOfDimensions; i++) {
      if (n0(i) == 0) {
        //expected gain for splitting according to the initialization list
        expectedGain(i) = arma::min(initListEvaluations_.col(i)) - initListEvaluations_(initialPointIndex_(i), i);
        if (expectedGain(i) < emin) {
          emin = expectedGain(i);
          isplit_(par) = i;
          //index again so use 1, matlab=2
          z_(1, par) = arma::datum::inf;
        }
      } else {
        arma::Col<double> z1;
        z1 << baseVertex_(i) << x1(i) << x2(i) << arma::endr;
        //index again at "boxBaseVertexFunctionValues_(0,par)" so use 0, matlab=1
        arma::Col<double> z2;
        z2 << 0 << f1(i) - boxBaseVertexFunctionValues_(0, par) << f2(i) - boxBaseVertexFunctionValues_(0, par) << arma::endr;
        arma::Col<double> d = quadraticPolynomialInterpolation(z1, z2);
        arma::Col<double> eta = subint(baseVertex_(i), oppositeVertex_(i));
        //safeguard against splitting too close to x(i)
        double xi1 = arma::min(eta);
        double xi2 = arma::max(eta);
        double z = minimumQuadraticPolynomial(xi1, xi2, d, z1);
        expectedGain(i) = quadraticPolynomial(z, d, z1);
        if (expectedGain(i) < emin) {
          emin = expectedGain(i);
          isplit_(par) = i;
          //index again so use 1, matlab=2
          z_(1, par) = z;
        }
      }
    }
    return expectedGain;
  }

  unsigned int MultilevelCoordinateSearch::startSweep() {
    record_ = arma::Col<arma::uword>(boxDivisions_ - 1, arma::fill::zeros);
    unsigned int s = boxDivisions_;
    for (unsigned int i = 0; i < nboxes_; i++) {
      if (level_(i) > 0) {
        if (level_(i) < s) {
          s = level_(i);
        }
        if (!record_(level_(i))) {
          record_(level_(i)) = i;
        } else if (boxBaseVertexFunctionValues_(0, i) < boxBaseVertexFunctionValues_(0, record_(level_(i)))) {
          record_(level_(i)) = i;
        }
      }
    }
    return s;
  }

  void MultilevelCoordinateSearch::vert1(int updateIndex, unsigned int j, unsigned int m, arma::Col<double> x, arma::Col<double> x1, arma::Col<double> x2, arma::Col<double> f1, arma::Col<double> f2) {
    int j1 = 0;
    //matlab checks for 1, since it's index use 0
    if (j == 0) {
      //also index again
      j1 = 1;
    } else {
      j1 = 0;
    }
    x(updateIndex) = z_(j1);
    if (x1(updateIndex) == arma::datum::inf) {
      x1(updateIndex) = z_(j);
      f1(updateIndex) = f1(updateIndex) + boxBaseVertexFunctionValues_(j);
    } else if (x2(updateIndex) == arma::datum::inf && x1(updateIndex) != z_(j)) {
      x2(updateIndex) = z_(j);
      f2(updateIndex) = f2(updateIndex) + boxBaseVertexFunctionValues_(j);
    }
  }

  void MultilevelCoordinateSearch::vert2(int updateIndex, unsigned int j, unsigned int m, arma::Col<double> x, arma::Col<double> x1, arma::Col<double> x2, arma::Col<double> f1, arma::Col<double> f2) {
    int j1 = 0;
    //matlab checks for 1, since it's index use 0
    if (j == 0) {
      //also index again
      j1 = 1;
    } else {
      j1 = 0;
    }
    if (x1(updateIndex) == arma::datum::inf) {
      x1(updateIndex) = z_(j);
      f1(updateIndex) = f1(updateIndex) + boxBaseVertexFunctionValues_(j);
      if (x(updateIndex) != z_(j1)) {
        x2(updateIndex) = z_(j1);
        f2(updateIndex) = f2(updateIndex) + boxBaseVertexFunctionValues_(j1);
      }
    } else if (x2(updateIndex) == arma::datum::inf && x1(updateIndex) != z_(j)) {
      x2(updateIndex) = z_(j);
      f2(updateIndex) = f2(updateIndex) + boxBaseVertexFunctionValues_(j);
    } else if (x2(updateIndex) == arma::datum::inf) {
      x2(updateIndex) = z_(j1);
      f2(updateIndex) = f2(updateIndex) + boxBaseVertexFunctionValues_(j1);
    }
  }

  void MultilevelCoordinateSearch::vert3(int updateIndex, unsigned int j, unsigned int m, unsigned int f0Index, arma::Col<double> x1, arma::Col<double> x2, arma::Col<double> f1, arma::Col<double> f2) {
    int k1 = 0;
    int k2 = 0;
    //matlab checks for 1, since it's index use 0
    if (j == 0) {
      //also index again
      k1 = 1;
      k2 = 2;
    } else if (j == populationSize_) {
      //TODO: might need to decrement one more
      k1 = populationSize_ - 2;
      k2 = populationSize_ - 1;
    } else {
      //TODO: might need to in-/decrement one more
      k1 = j - 1;
      k2 = j + 1;
    }
    x1(updateIndex) = x0_(k1);
    x2(updateIndex) = x0_(k2);
    f1(updateIndex) = f1(updateIndex) + boxBaseVertexFunctionValues_(k1);
    f2(updateIndex) = f2(updateIndex) + boxBaseVertexFunctionValues_(k2);
  }

  void MultilevelCoordinateSearch::updtf(unsigned int numberOfDimensions, unsigned int splittingIndex, double fold, arma::Col<double> x1, arma::Col<double> x2, arma::Col<double> f1, arma::Col<double> f2, double baseVertexValueCurrentBox) {
    for (int i = 0; i < numberOfDimensions; i++) {
      if (i != splittingIndex) {
        if (x1(i) == arma::datum::inf) {
          f1(i) = f1(i) + fold - baseVertexValueCurrentBox;
        }
        if (x2(i) == arma::datum::inf) {
          f2(i) = f2(i) + fold - baseVertexValueCurrentBox;
        }
      }
    }
    //updtf.m sets fold = f here. since the inputvalue fold never gets changed, this doesn't actually belong here.
  }

  arma::Col<double> MultilevelCoordinateSearch::subint(double x, double y) {
    int x2 = y;
    int f = 1000;
    if (f * std::abs(x) < 1) {
      if (std::abs(y) > f) {
        x2 = std::copysign(1.0, y);
      }
    } else {
      if (std::abs(y) > f * std::abs(x)) {
        //TODO: c++ standardlibraries have no signum. wat. using copysign instead...
        //original matlab: x2 = 10.*sign(y)*abs(x);
        x2 = 10 * std::copysign(x, y);
      }
    }
    return arma::Col<double>(x + (x2 - x) / 10.0, x2);
  }

  void MultilevelCoordinateSearch::updateRecord(unsigned int label, int level, arma::Col<double> f) {
    if (record_.n_elem < level) {
      record_(level) = label;
    } else if (record_(level) == 0) {
      record_(level) = label;
    } else if (f(label) < f(record_(level))) {
      record_(level) = label;
    }
  }
}