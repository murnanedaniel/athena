/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef JIVEXML_VERTEXRETRIEVER_H
#define JIVEXML_VERTEXRETRIEVER_H

#include "JiveXML/IDataRetriever.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "TrkParameters/TrackParameters.h"

namespace JiveXML{

  /**
   * @class VertexRetriever
   * @brief Retrieves all @c Trk::VxCandidate objects
   *
   *  - @b Properties
   *    - <em>PrimaryVertexCollection</em><tt> =  'VxPrimaryCandidate'</tt>: @copydoc m_primaryVertexKey
   *    - <em>IgnoreTrackAssociations</em><tt> =  False</tt>: @copydoc m_ignoreTrackAssociations
   *    - <em>DoWritePrimAndSecVertexOnly</em><tt> =  False</tt>: @copydoc m_doWritePrimAndSecVertexOnly
   *    - <em>DoWriteHLT</em><tt> =  False</tt>: @copydoc m_doWriteHLT
   *
   *  - @b Retrieved @b Data
   *    - <em> x, y, z </em> : coordinates of vertex
   *    - @e primVxCand : boolean for vertex candidate
   *    - @e chi2 : @f$\chi^2@f$ of vertex fit
   *    - @e covMatrix : covariance matrix of vertex fit
   *    - @e numTracks : number of tracks associated with the vertex
   *    - @e tracks : indices of associated tracks in track collection
   *    .
   */
  class VertexRetriever : public IDataRetriever,
                          public AthAlgTool
  {
    public:
    
      /// Standard constructor
      VertexRetriever(const std::string& t,const std::string& n,const IInterface* p);
      
      /// Retrieve all the data
      virtual StatusCode retrieve(ToolHandle<IFormatTool> &FormatTool); 

      /// Return the name of the data type
      virtual std::string dataTypeName() const { return typeName; };

    private:

      //@name Property members
      //@{
      /// StoreGate key for primary vertex candidate collection
      std::string m_primaryVertexKey;
      /// StoreGate key for secondary vertex candidate collection
      std::string m_secondaryVertexKey;
      /// StoreGate key for conversion candidate collection
      std::string m_conversionVertexKey;
      /// wether to write HLTAutoKey objects
      bool m_doWriteHLT;
      /// write primary and secondary vertizes only - placeholder, to be removed
      bool m_doWritePrimAndSecVertexOnly;
      //@}
      /// StoreGate key for track collection for association
      std::string m_trackCollection;
      /// Chi^2 over NumberOfDegreesOfFreedom cut
      float m_chi2Cut;

      ///The data type that is generated by this retriever
      const std::string typeName;
      
      std::vector<const Trk::Perigee*> perigeeVector;
      virtual StatusCode fillPerigeeList(); 
  };

}
#endif
