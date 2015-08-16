#include <mantella_bits/samplesSelection/bestNeighbourhoodSamplesSelection.hpp>

namespace mant {
  void BestNeighbourhoodSamplesSelection::selectImplementation() {
    arma::Col<double> bestParameter;
    double bestObjectiveValue = std::numeric_limits<double>::infinity();
    for (const auto& sample : samples_) {
      if(sample.second < bestObjectiveValue) {
        bestParameter = sample.first;
        bestObjectiveValue = sample.second;
      }
    }
    
    arma::Col<double> distances(samples_.size());
    arma::uword n = 0;
    for (const auto& sample : samples_) {
      distances(++n) = arma::norm(bestParameter - sample.first);
    }

    for (const auto& i : static_cast<arma::Col<arma::uword>>(static_cast<arma::Col<arma::uword>>(arma::sort_index(distances)).head(numberOfSelectedSamples_))) {
      const auto& selectedSample = std::next(std::begin(samples_), i);
      selectedSamples_.insert({selectedSample->first, selectedSample->second});
    }
  }
  
  std::string BestNeighbourhoodSamplesSelection::toString() const {
    return "nearest_to_best_samples_selection";
  }
}