/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef AthenaMonitoringKernel_HistogramFiller_HistogramFiller2D_h
#define AthenaMonitoringKernel_HistogramFiller_HistogramFiller2D_h

#include "TH2.h"

#include <boost/range/combine.hpp>
#include "AthenaMonitoringKernel/HistogramFiller.h"


namespace Monitored {
  /**
   * @brief Filler for plain 2D histogram
   */
  class HistogramFiller2D : public HistogramFiller {
  public:
    HistogramFiller2D(const HistogramDef& definition, std::shared_ptr<IHistogramProvider> provider)
      : HistogramFiller(definition, provider) {}

    virtual unsigned fill() override {
      if (m_monVariables.size() != 2) {
        return 0;
      }

      // handling of the weight
      std::function<double(size_t)> weightAccessor = [] (size_t ){ return 1.0; };  // default is always 1.0
      std::vector<double> weightVector;
      if ( m_monWeight != nullptr ) {
	weightVector = m_monWeight->getVectorRepresentation();
	weightAccessor = [&](size_t i){
	  if ( weightVector.size() == 1 )
	    return weightVector[0];
	  return weightVector[i];
	};
      }

      struct Extractor {
	std::function<double(size_t)> doublesAccessor = nullptr;
	std::function<const char*(size_t)> stringsAccessor = nullptr;
	std::vector<double> doublesVector;
	std::vector<std::string> stringsVector;
	size_t size() const {
	  return std::max( doublesVector.size(), stringsVector.size() );
	}
      };
      Extractor value1;
      Extractor value2;
      auto prepareExtractor = [&]( Extractor& ex, int index ) {
	if ( not m_monVariables[index].get().hasStringRepresentation() ) {
	  ex.doublesVector = m_monVariables[index].get().getVectorRepresentation();
	  ex.doublesAccessor = [&](size_t i){
	    if ( ex.doublesVector.size() == 1 )
	      return ex.doublesVector[0];
	    return ex.doublesVector[i];
	  };
	} else {
	  ex.stringsVector = m_monVariables[index].get().getStringVectorRepresentation();
	  ex.stringsAccessor = [&]( size_t i) {
	    if ( ex.stringsVector.size() == 1 )
	      return ex.stringsVector[0].c_str();
	    return ex.stringsVector[i].c_str();
	  };
	}
      };
      prepareExtractor( value1, 0 );
      prepareExtractor( value2, 1 );




      // rather unpleaseant code but I did not want complicate it further
      // we need to handle now 4 cases,
      // double-double, string-double, double-string and string-string in calling the fill, we always pass the weigh
      const size_t maxsize = std::max( value1.size(), value2.size() );
      if ( value1.doublesAccessor and value2.doublesAccessor )
	fill( maxsize, value1.doublesAccessor, value2.doublesAccessor, weightAccessor );
      else if ( value1.stringsAccessor and value2.doublesAccessor )
	fill( maxsize, value1.stringsAccessor, value2.doublesAccessor, weightAccessor );
      else if ( value1.doublesAccessor and value2.stringsAccessor )
	fill( maxsize, value1.doublesAccessor, value2.stringsAccessor, weightAccessor );
      else
	fill( maxsize, value1.stringsAccessor, value2.stringsAccessor, weightAccessor );
      return maxsize;
    }

  protected:
    template<typename F1, typename F2, typename F3>
    void fill( size_t n, F1 f1, F2 f2, F3 f3 ) {

      std::lock_guard<std::mutex> lock(*(this->m_mutex));
      auto histogram = this->histogram<TH2>();
      for ( size_t i = 0; i < n; ++i ) {
	histogram->Fill( f1(i), f2(i), f3(i) );
      }
    }

  };
}

#endif /* AthenaMonitoringKernel_HistogramFiller_HistogramFiller2D_h */
