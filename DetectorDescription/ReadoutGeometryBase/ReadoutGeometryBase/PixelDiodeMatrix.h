/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// PixelDiodeMatrix.h
///////////////////////////////////////////////////////////////////

#ifndef READOUTGEOMETRYBASE_PIXELDIODEMATRIX_H
#define READOUTGEOMETRYBASE_PIXELDIODEMATRIX_H

#include "ReadoutGeometryBase/SiCellId.h"
#include "GeoPrimitives/GeoPrimitives.h"
#include "TrkDetElementBase/TrkDetElementBase.h"

#include <memory>

namespace InDetDD {

   /** @class PixelDiodeMatrix
   
    Class  used to describe the segmentation of the pixel and allow for conversion between cell id and position.
     
    The class PixelDiodeMatrix represents what it calls a cell. In the context
    of this class a cell can either be just a simple cell as one would normally
    think of a pixel cell (eg a 50x400um or 50x600um cell) or it can be a
    container of multiple cells. To allow for different size cells and layout,
    it contains a lower cell, an upper cell and multiple middle
    cells. The middle cells are all the same (so all have the same size), but
    the upper and lower cells can be different and can have different size or be left out
    if they are not needed. The direction of replication is also specified (eta
    or phi direction). The size of the cells in the direction orthoganal to the
    replication direction must all be the same (there is no check for this though).
    
    To help understand the structure and its use, consider the current pixel
    geometry. The pixel module has two cell sizes short and long cells. These
    are refered to as normal and big bellow.  The long cells are in the region
    at the edge of the the FE chips in order to cover the gap between FE chips.
    
    The assembly of this structure is done when building the sensitive detector
    in GeoModel in the class GeoPixelSiCrystal.
    
     - normalCell (short pixel cell) 50x400um
     - bigCell    (long pixel cell): 50x600um 
    
     - singleChipRow (represents one row associated to a FE)
       - Size: phiWidth = 50um, etaWidth = 600*2+400*16) um = 7.6 CLHEP::mm 
       - upper cell: bigCell
       - middle cells: 16 x normalCell
       - lower cell: bigcell
       - direction: eta
     
     - singleRow (represents a row of cells covering the full length of the 
       module)
       - Size: phiWidth = 50um, etaWidth = 8 * (600*2+400*16) um = 60.8 CLHEP::mm
       - upper cell: none
       - middle cells: 8 x singleChipRow
       - lower cell: none
       - direction: eta
    
     - fullMatrix (represents the full module)
       - Size: phiWidth = 50 * 328 um = 16.4 CLHEP::mm, etaWidth = 60.8 CLHEP::mm
       - upper cell: none
       - middle cells: 328 x singleRow
       - lower cell: non
       - direction: phi
    
    Each cell knows its total size and total number of child cells in the phi and eta
    direction. This is determined by totaling up these quantities
    when the object is constructed.
    
    The structure allows for efficient navigation from cell number to 
    position and vice-versa. There are two main methods which allow for this navigation:
    cellIdOfPosition() (from cell id to position) and positionOfCell() (from position to cell id).
    
    As it is assumed that child cells fill up completely its parent, there is
    no range check in these methods. The methods are called from PixelDiodeMap
    where it first checks that is within range of the top level cell (refered to
    as the top level matrix).

    PixelDiodeMatrix might inherit std::shared_ptr but needs to return std::shared_ptr of this pointer.
    To do this, std::enable_shared_from_this and shared_from_this() have to be used.
    https://en.cppreference.com/w/cpp/memory/enable_shared_from_this
    
    See the description of these methods for more details.
    
    @image html pixeldiodematrix.gif

    @author Grant Gorfine
    - modified & maintained: Nick Styles, Andreas Salzburger
    **/
    
    class PixelDiodeMatrix {
      
    public:
      
      enum Direction {phiDir, etaDir};
      
      /// Destructor
      ~PixelDiodeMatrix() = default;
      
      /// Construct method for just a single cell.
      static std::shared_ptr<const PixelDiodeMatrix> construct(double phiWidth, double etaWidth);

      /** Construct method with multiple matrices.
          There may be zero or one lower matrix (pointer can be 0).
          Multiple (numCells) middle matrics (must pass a non zero pointer).
          There may be zero or one upper matrix (pointer can be 0).
          These must all have the same width and cells in the non replicating direction. */
      static std::shared_ptr<const PixelDiodeMatrix> construct(Direction direction,  // phi or eta
                                                               std::shared_ptr<const PixelDiodeMatrix> lowerCell,
                                                               std::shared_ptr<const PixelDiodeMatrix> middleCells,
                                                               int numCells,
                                                               std::shared_ptr<const PixelDiodeMatrix> upperCells);

      /** Return cell Id corresponding to a relative position within the matrix.
          The cellId is modified and returns the id relative to the passed cellId.
          That is, it adds the relative cellId to the cellId passed to the function. 
          A pointer to the correspond cell is returned. This can be used to get the
          size of the cell. */
      const PixelDiodeMatrix *cellIdOfPosition(const Amg::Vector2D & position, SiCellId & cellId) const;
      
      /** Return position correspong to cell with relative id withing the matrix.
          Returns the relative position added to the position passed in.
          A pointer to the correspond cell is returned. This can be used to get the
          size of the cell. */
      const PixelDiodeMatrix *positionOfCell(const SiCellId & cellId, Amg::Vector2D & position) const;
      
      /** Width in phi (x) direction.*/
      double phiWidth() const;

      /** Inverse of width in phi (x) direction.*/
      double phiWidthInverse() const;

      /** Width in eta (y) direction. */
      double etaWidth() const;

      /** Inverse of width in eta (y) direction. */
      double etaWidthInverse() const;

      /** Number of cells in phi (x) direction. */
      int phiCells() const;
      
      /** Number of cells in eta (y) direction. */
      int etaCells() const;
      
      /** Query wether the matrix is just a single cell. Mainly for internal use. */
      bool singleCell() const;
      
    private:
      /** Hidden constructor */
      PixelDiodeMatrix() {};

      /// Initialize for just a single cell.
      void initialize(double phiWidth, double etaWidth);
      
      /** Initialize for multiple matrices.
          There may be zero or one lower matrix (pointer can be 0).
          Multiple (numCells) middle matrics (must pass a non zero pointer).
          There may be zero or one upper matrix (pointer can be 0).
          These must all have the same width and cells in the non replicating direction. */
      void initialize(Direction direction,  // phi or eta
                      std::shared_ptr<const PixelDiodeMatrix> lowerCell,
                      std::shared_ptr<const PixelDiodeMatrix> middleCells,
                      int numCells,
                      std::shared_ptr<const PixelDiodeMatrix> upperCells);

      double        m_phiWidth;
      double        m_phiWidthInverse;
      double        m_etaWidth;
      double        m_etaWidthInverse;
      int           m_phiCells;    
      int           m_etaCells;
      Direction     m_direction;
      int m_numCells; // number of MiddleCells along m_direction direction.
      std::shared_ptr<const PixelDiodeMatrix> m_lowerCell;
      std::shared_ptr<const PixelDiodeMatrix> m_middleCells;
      std::shared_ptr<const PixelDiodeMatrix> m_upperCell;
      bool m_singleCell;
    
    };
    
    inline double PixelDiodeMatrix::phiWidth() const
    { 
      return m_phiWidth;
    }

    inline double PixelDiodeMatrix::phiWidthInverse() const
    { 
      return m_phiWidthInverse;
    }
    
    inline double PixelDiodeMatrix::etaWidth() const
    { 
      return m_etaWidth;
    }
    
    inline double PixelDiodeMatrix::etaWidthInverse() const
    { 
      return m_etaWidthInverse;
    }
    
    inline int PixelDiodeMatrix::phiCells() const
    { 
      return m_phiCells;
    }
    
    inline int PixelDiodeMatrix::etaCells() const
    { 
      return m_etaCells;
    }
    
    inline bool PixelDiodeMatrix::singleCell() const
    {
      return m_singleCell;
    }

} // end namespace

#endif //READOUTGEOMETRYBASE_PIXELDIODEMATRIX_H
