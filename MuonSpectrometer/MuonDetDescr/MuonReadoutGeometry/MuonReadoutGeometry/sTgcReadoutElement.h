/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUONREADOUTGEOMETRY_STGCREADOUTELEMENT_H
#define MUONREADOUTGEOMETRY_STGCREADOUTELEMENT_H

#include <string>
#include <utility>

#include "MuonIdHelpers/sTgcIdHelper.h"
#include "MuonReadoutGeometry/MuonChannelDesign.h"
#include "MuonReadoutGeometry/MuonClusterReadoutElement.h"
#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonReadoutGeometry/MuonPadDesign.h"

class BLinePar;
class GeoVFullPhysVol;

namespace MuonGM {
    /**
       An sTgcReadoutElement corresponds to a single STGC module; therefore
       typicaly a barrel muon station contains:
    */

    class sTgcReadoutElement final : public MuonClusterReadoutElement {
    public:
        /** constructor */

        sTgcReadoutElement(GeoVFullPhysVol* pv, const std::string& stName, int zi, int fi, int mL, bool is_mirrored, MuonDetectorManager* mgr);

        /** destructor */
        ~sTgcReadoutElement();

        /** function to be used to check whether a given Identifier is contained in the readout element */
        virtual bool containsId(const Identifier& id) const override final;

        /** distance to readout.
            If the local position is outside the active volume, the function first shift the position back into the active volume */
        virtual double distanceToReadout(const Amg::Vector2D& pos, const Identifier& id) const override final;

        /** strip number corresponding to local position.
            Should be renamed to channelNumber : the only public access for all hit types */
        virtual int stripNumber(const Amg::Vector2D& pos, const Identifier& id) const override final;

        /** Channel pitch. Gives full pitch for strips, width of a full wire group
        Gives the Height of a pad */
        double channelPitch(const Identifier& id) const;

        /** strip position - should be renamed to channel position
            If the strip number is outside the range of valid strips, the function will return false */
        virtual bool stripPosition(const Identifier& id, Amg::Vector2D& pos) const override final;

        bool stripGlobalPosition(const Identifier& id, Amg::Vector3D& gpos) const;

        /** pad number corresponding to local position */
        int padNumber(const Amg::Vector2D& pos, const Identifier& id) const;

        /** wire number corresponding to local position */
        int wireNumber(const Amg::Vector2D& pos, const Identifier& id) const;

        /** single wire pitch.
         *  sTGC wire pitch is the same for all chambers,
         *  so the default gas gap is set to the 1st gap */
        double wirePitch(int gas_gap = 1) const;

        /** Get the local position of the first wire of the chamber corresponding to the identifier */
        double positionFirstWire(const Identifier& id) const;

        /** Get the total number of wires (single wires) of a chamber **/
        int numberOfWires(const Identifier& id) const;

        /** pad position */
        bool padPosition(const Identifier& id, Amg::Vector2D& pos) const;

        /** pad corners */
        bool padCorners(const Identifier& id, std::vector<Amg::Vector2D>& corners) const;

        /** is eta=0 of QL1 or QS1? */
        bool isEtaZero(const Identifier& id, double posY) const;

        /** number of layers in phi/eta projection */
        virtual int numberOfLayers(bool) const override final;

        /** number of strips per layer */
        virtual int numberOfStrips(const Identifier& layerId) const override final;
        virtual int numberOfStrips(int, bool measuresPhi) const override final;

        /** space point position for a given pair of phi and eta identifiers
            The LocalPosition is expressed in the reference frame of the phi
           surface. If one of the identifiers is outside the valid range, the
           function will return false */
        virtual bool spacePointPosition(const Identifier& phiId, const Identifier& etaId, Amg::Vector2D& pos) const override final;

        /** Global space point position for a given pair of phi and eta identifiers
            If one of the identifiers is outside the valid range, the function will
           return false */
        virtual bool spacePointPosition(const Identifier& phiId, const Identifier& etaId, Amg::Vector3D& pos) const override final;

        /** TrkDetElementInterface */
        virtual Trk::DetectorElemType detectorType() const override final { return Trk::DetectorElemType::sTgc; }

        /** space point position for a pair of phi and eta local positions and a layer identifier
            The LocalPosition is expressed in the reference frame of the phi projection.
        */
        void spacePointPosition(const Amg::Vector2D& phiPos, const Amg::Vector2D& etaPos, Amg::Vector2D& pos) const;

        /** space point position, corrected for chamber deformations (b-lines), if b-lines are enabled.
            Accepts a precision (x) coordinate and a y-seed, in the local layer frame, and returns a 3D position, 
            in the same frame so that sTgcReadoutElement::transform() can be directly cast on it. Accounts for:
            a) PCB deformations (as-built), if as-built conditions are enabled
            b) Chamber deformations (b-lines), if b-lines are enabled
        */
        void spacePointPosition(const Identifier& layerId, double locXpos, double locYpos, Amg::Vector3D& pos) const;

        /** simHit local (SD) To Global position - to be used by MuonGeoAdaprors only      */
        Amg::Vector3D localToGlobalCoords(const Amg::Vector3D& locPos, Identifier id) const;

        /** @brief function to fill tracking cache */
        virtual void fillCache() override final;
        virtual void refreshCache() override final {
            clearCache();
            fillCache();
        }

        /** @brief returns the hash to be used to look up the surface and transform in the MuonClusterReadoutElement tracking cache */
        virtual int surfaceHash(const Identifier& id) const override final;

        /** @brief returns the hash to be used to look up the surface and transform in the MuonClusterReadoutElement tracking cache */
        int surfaceHash(int gasGap, int channelType) const;

        /** @brief returns the hash to be used to look up the normal and center in the MuonClusterReadoutElement tracking cache */
        virtual int layerHash(const Identifier& id) const override;

        /** @brief returns the hash to be used to look up the normal and center in the MuonClusterReadoutElement tracking cache */
        // int layerHash( int gasGap) const;        // does not fit in the scheme ( layer hash needs to follow surface hash )

        /** returns the hash function to be used to look up the surface boundary for a given identifier */
        virtual int boundaryHash(const Identifier& id) const override final;

        /** @brief returns whether the current identifier corresponds to a phi measurement */
        virtual bool measuresPhi(const Identifier& id) const override final;

        /** @brief initialize the design classes for this readout element */
        void initDesign(double largeX, double smallX, double lengthY, double stripPitch, double wirePitch, double stripWidth,
                        double wireWidth, double thickness);

        /** returns the MuonChannelDesign class for the given identifier */
        const MuonChannelDesign* getDesign(const Identifier& id) const;

        /** returns the MuonChannelDesign class  */
        const MuonChannelDesign* getDesign(int gasGap, int channelType) const;

        /** returns the MuonChannelDesign class for the given identifier */
        const MuonPadDesign* getPadDesign(const Identifier& id) const;
        MuonPadDesign* getPadDesign(const Identifier& id);

        /** returns the MuonChannelDesign */
        const MuonPadDesign* getPadDesign(int gasGap) const;

        /** set methods only to be used by MuonGeoModel */
        void setIdentifier(const Identifier& id);

        /** set methods only to be used by MuonGeoModel */
        void setChamberLayer(int ml) { m_ml = ml; }

        // double getSectorOpeningAngle(bool isLargeSector);

        /** read A-line parameters and include the chamber rotation/translation 
            in the local-to-global (ATLAS) reference frame transformaton */
        void setDelta(const ALinePar& aline);
        void setDelta(MuonDetectorManager* mgr);

        /** read B-line (chamber-deformation) parameters */
        void setBLinePar(const BLinePar& bLine);
        void setBLinePar(MuonDetectorManager* mgr);
    
        /** transform a position (in chamber-frame coordinates) to the deformed-chamber geometry */
        void posOnDefChamber(Amg::Vector3D& locPosML) const;

        bool  has_ALines() const { return (m_ALinePar != nullptr); }
        bool  has_BLines() const { return (m_BLinePar != nullptr); }
        const ALinePar* getALinePar() const { return m_ALinePar; }
        const BLinePar* getBLinePar() const { return m_BLinePar; }
        void  clearALinePar();
        void  clearBLinePar() { m_BLinePar = nullptr; }
        
    private:
        std::vector<MuonChannelDesign> m_phiDesign;
        std::vector<MuonChannelDesign> m_etaDesign;
        std::vector<MuonPadDesign> m_padDesign;

        std::vector<int> m_nStrips;
        std::vector<int> m_nWires;
        std::vector<int> m_nPads;
        int    m_nlayers{0};
        int    m_ml{0};
        double m_offset{0.};
        
        double m_sWidthChamber{0.}; // bottom base length (full chamber)
        double m_lWidthChamber{0.}; // top base length (full chamber)
        double m_lengthChamber{0.}; // radial size (full chamber)
        double m_tckChamber{0.};    // thickness (full chamber)
        bool   m_diamondShape{false};
        const ALinePar*  m_ALinePar{nullptr};
        const BLinePar*  m_BLinePar{nullptr};
        Amg::Transform3D m_delta{Amg::Transform3D::Identity()};

        // const double m_largeSectorOpeningAngle = 28.0;
        // const double m_smallSectorOpeningAngle = 17.0;

        // surface dimensions for strips
        std::vector<double> m_halfX;
        std::vector<double> m_minHalfY;
        std::vector<double> m_maxHalfY;
        // surface dimensions for pads and wires
        std::vector<double> m_PadhalfX;
        std::vector<double> m_PadminHalfY;
        std::vector<double> m_PadmaxHalfY;

        // transforms (RE->layer)
        Amg::Transform3D m_Xlg[4];
    };

    inline void sTgcReadoutElement::clearALinePar() {
        if (has_ALines()) {
            m_ALinePar = nullptr; 
            m_delta = Amg::Transform3D::Identity(); 
            refreshCache();
        }
    }

    inline int sTgcReadoutElement::surfaceHash(const Identifier& id) const {
        return surfaceHash(manager()->stgcIdHelper()->gasGap(id), manager()->stgcIdHelper()->channelType(id));
    }

    inline int sTgcReadoutElement::surfaceHash(int gasGap, int channelType) const {
        return (gasGap - 1) * 3 + (2 - channelType);  // assumes channelType=2 (wires), 1(strips), 0(pads)
    }

    inline int sTgcReadoutElement::layerHash(const Identifier& id) const {
        return surfaceHash(id);  // don't have a choice here : rewrite MuonClusterReadoutElement first
    }

    inline int sTgcReadoutElement::boundaryHash(const Identifier& id) const {
        int iphi = manager()->stgcIdHelper()->channelType(id) != 1 ? 1 : 0;  // wires and pads have locX oriented along phi
        if (std::abs(getStationEta()) < 3) iphi += 2 * (manager()->stgcIdHelper()->gasGap(id) - 1);
        return iphi;
    }

    inline bool sTgcReadoutElement::measuresPhi(const Identifier& id) const { return (manager()->stgcIdHelper()->channelType(id) != 1); }

    inline const MuonChannelDesign* sTgcReadoutElement::getDesign(const Identifier& id) const {
        if (manager()->stgcIdHelper()->channelType(id) == 1) return &(m_etaDesign[manager()->stgcIdHelper()->gasGap(id) - 1]);
        if (manager()->stgcIdHelper()->channelType(id) == 2) return &(m_phiDesign[manager()->stgcIdHelper()->gasGap(id) - 1]);
        return 0;
    }

    inline const MuonPadDesign* sTgcReadoutElement::getPadDesign(const Identifier& id) const {
        if (manager()->stgcIdHelper()->channelType(id) == 0) return &(m_padDesign[manager()->stgcIdHelper()->gasGap(id) - 1]);
        return 0;
    }

    inline MuonPadDesign* sTgcReadoutElement::getPadDesign(const Identifier& id) {
        if (manager()->stgcIdHelper()->channelType(id) == 0) return &(m_padDesign[manager()->stgcIdHelper()->gasGap(id) - 1]);
        return 0;
    }

    inline const MuonChannelDesign* sTgcReadoutElement::getDesign(int gasGap, int channelType) const {
        if (channelType == 1) return &(m_etaDesign[gasGap - 1]);
        if (channelType == 2) return &(m_phiDesign[gasGap - 1]);
        return 0;
    }

    inline const MuonPadDesign* sTgcReadoutElement::getPadDesign(int gasGap) const { return &(m_padDesign[gasGap - 1]); }

    inline double sTgcReadoutElement::distanceToReadout(const Amg::Vector2D& pos, const Identifier& id) const {
        const MuonChannelDesign* design = getDesign(id);
        if (!design) return -10000.;
        return design->distanceToReadout(pos);
    }

    inline int sTgcReadoutElement::stripNumber(const Amg::Vector2D& pos, const Identifier& id) const {
        if (manager()->stgcIdHelper()->channelType(id) == 0) return padNumber(pos, id);

        const MuonChannelDesign* design = getDesign(id);
        if (!design) return -1;
        return design->channelNumber(pos);
    }

    inline bool sTgcReadoutElement::stripPosition(const Identifier& id, Amg::Vector2D& pos) const {
        if (manager()->stgcIdHelper()->channelType(id) == 0) return padPosition(id, pos);

        const MuonChannelDesign* design = getDesign(id);
        if (!design) return 0;
        return design->channelPosition(manager()->stgcIdHelper()->channel(id), pos);
    }

    inline bool sTgcReadoutElement::stripGlobalPosition(const Identifier& id, Amg::Vector3D& gpos) const {
        Amg::Vector2D lpos(0., 0.);
        if (!stripPosition(id, lpos)) return false;
        surface(id).localToGlobal(lpos, Amg::Vector3D(0., 0., 0.), gpos);
        return true;
    }

    inline bool sTgcReadoutElement::padPosition(const Identifier& id, Amg::Vector2D& pos) const {
        const MuonPadDesign* design = getPadDesign(id);
        if (!design) return false;

        int padEta = manager()->stgcIdHelper()->padEta(id);
        int padPhi = manager()->stgcIdHelper()->padPhi(id);

        return design->channelPosition(std::pair<int, int>(padEta, padPhi), pos);
    }

    inline bool sTgcReadoutElement::padCorners(const Identifier& id, std::vector<Amg::Vector2D>& corners) const {
        const MuonPadDesign* design = getPadDesign(id);
        if (!design) return false;

        int padEta = manager()->stgcIdHelper()->padEta(id);
        int padPhi = manager()->stgcIdHelper()->padPhi(id);

        return design->channelCorners(std::pair<int, int>(padEta, padPhi), corners);
    }

    inline bool sTgcReadoutElement::isEtaZero(const Identifier& id, double posY) const {
        // This function returns true if we are in the eta 0 region of QL1/QS1

        const MuonChannelDesign* design = getDesign(id);
        if (!design) return false;

        if (design->wireCutout == 0.) return false;  // Not QL1 / QS1

        if (posY < 0.5 * design->xSize - design->wireCutout) return true;

        return false;
    }

    inline int sTgcReadoutElement::numberOfLayers(bool) const { return m_nlayers; }

    inline int sTgcReadoutElement::numberOfStrips(const Identifier& layerId) const {
        return numberOfStrips(manager()->stgcIdHelper()->gasGap(layerId) - 1, manager()->stgcIdHelper()->channelType(layerId) == 2);
    }

    inline int sTgcReadoutElement::numberOfStrips(int lay, bool measPhi) const {
        if (lay > -1 && lay < m_nlayers) { return measPhi ? m_nWires[lay] : m_nStrips[lay]; }
        return -1;
    }

    inline bool sTgcReadoutElement::spacePointPosition(const Identifier& phiId, const Identifier& etaId, Amg::Vector2D& pos) const {
        Amg::Vector2D phiPos;
        Amg::Vector2D etaPos;
        if (!stripPosition(phiId, phiPos) || !stripPosition(etaId, etaPos)) return false;
        spacePointPosition(phiPos, etaPos, pos);
        return true;
    }

    inline bool sTgcReadoutElement::spacePointPosition(const Identifier& phiId, const Identifier& etaId, Amg::Vector3D& pos) const {
        Amg::Vector2D lpos;
        spacePointPosition(phiId, etaId, lpos);
        surface(phiId).localToGlobal(lpos, pos, pos);
        return true;
    }

    inline void sTgcReadoutElement::spacePointPosition(const Amg::Vector2D& phiPos, const Amg::Vector2D& etaPos, Amg::Vector2D& pos) const {
        pos[0] = phiPos.x();
        pos[1] = etaPos.x();
    }

}  // namespace MuonGM

#endif  // MUONREADOUTGEOMETRY_STGCREADOUTELEMENT_H
