/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef JIVEXML_BJETRETRIEVER_H
#define JIVEXML_BJETRETRIEVER_H

#include <string>
#include <vector>
#include <map>

#include "JiveXML/IDataRetriever.h"
#include "AthenaBaseComps/AthAlgTool.h"

class JetCollection;

namespace JiveXML{
  
  /**
   * @class BJetRetriever
   * @brief Retrieves @c Jet @c objects, filter from weight for BJet (JetEvent/Jet)
   *
   *  - @b Properties
   *    - FavouriteJetCollection
   *    - OtherJetCollections
   *    - DoWriteHLT
   *
   *  - @b Retrieved @b Data
   *    - Usual four-vector: phi, eta, et, mass, energy, px, py, pz
   *  . - Some extras for jet: flavourTagWeight, charge
   *    - No cells in AOD Jets.
   */
  class BJetRetriever : virtual public IDataRetriever,
                                   public AthAlgTool {
    
    public:
      
      /// Standard Constructor
      BJetRetriever(const std::string& type,const std::string& name,const IInterface* parent);
      
      /// Retrieve all the data
      virtual StatusCode retrieve(ToolHandle<IFormatTool> &FormatTool); 
      const DataMap getData(const JetCollection*);

      /// Return the name of the data type
      virtual std::string dataTypeName() const { return m_typeName; };

    private:
      ///The data type that is generated by this retriever
      const std::string m_typeName;

      std::string m_sgKeyFavourite;
      std::vector<std::string> m_otherKeys;
      bool m_doWriteHLT;
      float m_weightCut;
  };
}
#endif
