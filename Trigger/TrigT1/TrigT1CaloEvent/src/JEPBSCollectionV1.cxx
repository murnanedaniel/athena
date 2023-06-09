
#include "TrigT1CaloEvent/CMMEtSums.h"
#include "TrigT1CaloEvent/CMMJetHits.h"
#include "TrigT1CaloEvent/JEMEtSums.h"
#include "TrigT1CaloEvent/JEMHits.h"
#include "TrigT1CaloEvent/JetElement.h"
#include "TrigT1CaloEvent/JEPBSCollectionV1.h"

namespace LVL1 {

JEPBSCollectionV1::JEPBSCollectionV1(
                 const DataVector<JetElement>* jeCollection,
                 const DataVector<JEMHits>*    hitCollection,
     const DataVector<JEMEtSums>*  etCollection,
                 const DataVector<CMMJetHits>* cmmHitCollection,
     const DataVector<CMMEtSums>*  cmmEtCollection)
                 : m_jeCollection(jeCollection),
       m_hitCollection(hitCollection),
       m_etCollection(etCollection),
       m_cmmHitCollection(cmmHitCollection),
       m_cmmEtCollection(cmmEtCollection)
{
}

} // end namespace
