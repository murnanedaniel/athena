///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// TileCosmicMuon_p2.h 
// Header file for class TileCosmicMuon_p2
// Author: Jose Maneira <maneira@lip.pt>
// Date:   July 2009
/////////////////////////////////////////////////////////////////// 
#ifndef TILETPCNV_TILECOSMICMUON_P2_H
#define TILETPCNV_TILECOSMICMUON_P2_H

#include <vector>

// forward declarations
class TileCosmicMuonCnv_p2;

class TileCosmicMuon_p2 {

  // Make the AthenaPoolCnv class our friend
  friend class TileCosmicMuonCnv_p2;

public:

  /** Default constructor: 
   */
  TileCosmicMuon_p2() : m_time(0.0), m_positionX(0.0),
                        m_positionY(0.0), m_positionZ(0.0),
                        m_directionPhi(0.0), m_directionTheta(0.0),
                        m_fitQuality(0.0), m_fitNCells(0),
                        m_pathTop(), m_pathBottom(),
                        m_energyTop(), m_energyBottom(),
                        m_trackCellHash(), m_segmentPath(),
                        m_segmentPartitionModuleSampling()   {}

private:

  float m_time;//!< Time of track at selected plane (y=0 for cosmics z=0 for beam)
  float m_positionX; //!< X coordinate of point in track at selected plane (y=0 for cosmics z=0 for beam)
  float m_positionY; //!< Y coordinate of point in track at selected plane (y=0 for cosmics z=0 for beam)
  float m_positionZ; //!< Z coordinate of point in track at selected plane (y=0 for cosmics z=0 for beam)
  float m_directionPhi; //!< Phi angle of track direction
  float m_directionTheta; //!< Theta angle of track direction
  float m_fitQuality; //!< Fit parameter: 0= no fit; (Hough) 1=fit ok; (Minuit) >0 chi-square 
  int   m_fitNCells;  //!< Number of cells used in fit
  
  /** Vector with length of track within Tile on top modules [0]:sampling A; [1]: BC; [2]: D */
  std::vector<float> m_pathTop; 
  /** Vector with length of track within Tile on bottom modules [0]:sampling A; [1]: BC; [2]: D */
  std::vector<float> m_pathBottom; 
  /** Vector with sum energy of cells close to track on top modules [0]:sampling A; [1]: BC; [2]: D */
  std::vector<float> m_energyTop; 
  /** Vector with sum energy of cells close to track on bottom modules [0]:sampling A; [1]: BC; [2]: D */
  std::vector<float> m_energyBottom; 
 
  /** Vector with list of Identifier Hash of cells close to track.*/
  std::vector<unsigned int> m_trackCellHash; 
  /** Vector with length of track within Tile on a given segment */
  std::vector<float> m_segmentPath; 
  /** Vector with segment partition/module/sampling - one byte for each */
  std::vector<unsigned int> m_segmentPartitionModuleSampling; 
};

#endif //> TILETPCNV_TILECOSMICMUON_P2_H
