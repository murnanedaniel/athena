/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef JIVEXML_SEGMENTRETRIEVER_H
#define JIVEXML_SEGMENTRETRIEVER_H

#include "JiveXML/IDataRetriever.h"
#include "AthenaBaseComps/AthAlgTool.h"

namespace JiveXML {
  
  /**
   * @class SegmentRetriever
   * @brief Retrieves all @c Trk::SegmentCollection objects
   *
   *  - @b Properties
   *    - none
   *
   *  - @b Retrieved @b Data
   *    - <em>x,y,z</em>: coordinates of segment starting point
   *    - <em>eta,phi</em>: segment direction
   *    - <em>NHits</em>: number of hits associated with this segment
   *    - <em>hits</em>: ID of hits associated with this segment
   *  .
   */
  class SegmentRetriever : virtual public IDataRetriever,
                                   public AthAlgTool {
    
    public:
      
      /// Standard Constructor
      SegmentRetriever(const std::string& type,const std::string& name,const IInterface* parent);
      
      /// Retrieve all the data
      virtual StatusCode retrieve(ToolHandle<IFormatTool> &FormatTool); 

      /// Return the name of the data type
      virtual std::string dataTypeName() const { return typeName; };

    private:
      ///The data type that is generated by this retriever
      const std::string typeName;
  };
    
}

#endif
