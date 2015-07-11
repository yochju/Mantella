namespace mant {
  template <typename T>
  class ModifiedCuckooSearch : public PopulationBasedOptimisationAlgorithm<T> {
    public:
      explicit ModifiedCuckooSearch(
          const std::shared_ptr<OptimisationProblem<T>> optimisationProblem,
          const unsigned int populationSize) noexcept;
		  
      void setTopNestPortion(
					const double topNestPortion) noexcept;
			void setWorstNestPortion(
					const double worstNestPortion) noexcept;
      void setMaxLevyStepSize(
					const double maxLevyStepSize) noexcept;
	  
			std::string toString() const noexcept override;

    protected:
      double topNestPortion_;
			double worstNestPortion_;
      double maxLevyStepSize_;

      void optimiseImplementation() override;
			
			virtual arma::Col<T> levyFlight(
					const arma::Col<T> startingPoint,
					const double levyStepSize) noexcept;

			virtual unsigned int getRandomTopNest(
					const arma::Col<unsigned int> ranking) noexcept;
			
			virtual unsigned int getRandomNest() noexcept;
  };

  //
  //Implementation
  //

	template <typename T>
  ModifiedCuckooSearch<T>::ModifiedCuckooSearch(
      const std::shared_ptr<OptimisationProblem<T>> optimisationProblem,
      const unsigned int populationSize) noexcept
    : PopulationBasedOptimisationAlgorithm<T>(optimisationProblem, populationSize) {
    setTopNestPortion(0.25);
		setWorstNestPortion(0.75);
    setMaxLevyStepSize(1.0);
  }

  template <typename T>
  void ModifiedCuckooSearch<T>::optimiseImplementation() {
		const double goldenRatio = (1 + std::sqrt(5))/2;
    arma::Mat<T> nests = this->getRandomPopulation();

    arma::Col<double> objectiveValues(this->populationSize_);
    for (std::size_t n = 0; n < this->populationSize_; ++n) {
			++this->numberOfIterations_;
			objectiveValues(n) = this->getObjectiveValue(nests.col(n));
      this->updateBestParameter(nests.col(n), this->getSoftConstraintsValue(nests.col(n)), objectiveValues(n));
    }
    
		unsigned int generation = 0;
    while(!this->isFinished() && !this->isTerminated()) {
			++this->numberOfIterations_;
			++generation;
			arma::Col<unsigned int> ranking = arma::sort_index(objectiveValues);
			
			for(std::size_t i = 0; i < this->populationSize_; ++i) {
				if(ranking(i) > std::ceil((1 - worstNestPortion_) * this->populationSize_)) {
					++this->numberOfIterations_;
					double levyStepSize = maxLevyStepSize_ / std::sqrt(generation);	
				
					nests.col(i) = levyFlight(nests.col(i), levyStepSize);
					objectiveValues(i) = this->getObjectiveValue(nests.col(i));
					this->updateBestParameter(nests.col(i), 0.0, objectiveValues(i));
					
					if(this->isFinished() || this->isTerminated()) {
						break;
					}
				}
			}
			
			for(std::size_t i = 0; i < this->populationSize_; ++i) {
				if(ranking(i) < std::floor(topNestPortion_ * this->populationSize_)) {
					++this->numberOfIterations_;
					unsigned int randTopNest = getRandomTopNest(ranking);
					
					if(randTopNest = i) {
						double levyStepSize = maxLevyStepSize_ / std::pow(generation, 2);
						
						const arma::Col<T>& cuckooCandidate = levyFlight(nests.col(i), levyStepSize);
						const double objectiveValue = this->getObjectiveValue(cuckooCandidate);
						
						unsigned int randNest = getRandomNest();
						if(objectiveValue < objectiveValues(randNest)) {
							nests.col(randNest) = cuckooCandidate;
							objectiveValues(randNest) = objectiveValue;
							this->updateBestParameter(cuckooCandidate, 0.0, objectiveValue);
						}
					}
					else {
						//generating new Nest on the line between the two chosen Nests
						unsigned int worseNest = ranking(i) > ranking(randTopNest) ? i : randTopNest;
						const arma::Col<T>& cuckooCandidate = this->boundaryHandling(nests.col(worseNest) + arma::abs(nests.col(i) - nests.col(randTopNest)) / goldenRatio);
						const double objectiveValue = this->getObjectiveValue(cuckooCandidate);
						
						unsigned int randNest = getRandomNest();
						if(objectiveValue < objectiveValues(randNest)) {
							nests.col(randNest) = cuckooCandidate;
							objectiveValues(randNest) = objectiveValue;
							bool test = this->updateBestParameter(cuckooCandidate, 0.0, objectiveValue);
						}
					}					
				}
				if(this->isFinished() || this->isTerminated()) {
					break;
				}
			}
    }
  }
	
	template <typename T>
	arma::Col<T> ModifiedCuckooSearch<T>::levyFlight(
			const arma::Col<T> startingPoint,
			const double levyStepSize) noexcept {
		//calculates a new Nest from the starting point, using the levyStepSize and an equation for levy flights
		return this->boundaryHandling(startingPoint + levyStepSize * (arma::randn<arma::Col<T>>(this->numberOfDimensions_) * (std::pow(std::tgamma(2.5) * std::sin(arma::datum::pi * 0.75) / (std::tgamma(1.25) * 1.5 * std::pow(2, 0.25)), 2/3)) / (arma::pow(arma::abs(arma::randn<arma::Col<T>>(this->numberOfDimensions_)), 2/3))));
	}
	
	template <typename T>
	unsigned int ModifiedCuckooSearch<T>::getRandomTopNest(
			const arma::Col<unsigned int> ranking) noexcept {
		unsigned int pos = getRandomNest();
		while(ranking(pos) >= std::floor(topNestPortion_ * this->populationSize_)) {
			pos = getRandomNest();
		}
		return pos;
	}
	
	template <typename T>
	unsigned int ModifiedCuckooSearch<T>::getRandomNest() noexcept{
		return std::rand() % this->numberOfDimensions_;
	}
  template <typename T>
  void ModifiedCuckooSearch<T>::setTopNestPortion(
			const double topNestPortion) noexcept {
    topNestPortion_ = topNestPortion;
  }
	
	template <typename T>
  void ModifiedCuckooSearch<T>::setWorstNestPortion(
			const double worstNestPortion) noexcept {
    worstNestPortion_ = worstNestPortion;
  }

  template <typename T>
  void ModifiedCuckooSearch<T>::setMaxLevyStepSize(
      const double maxLevyStepSize) noexcept {
    maxLevyStepSize_ = maxLevyStepSize;
  }

  template <typename T>
  std::string ModifiedCuckooSearch<T>::toString() const noexcept {
    return "modified_cuckoo_search";
  }
}
