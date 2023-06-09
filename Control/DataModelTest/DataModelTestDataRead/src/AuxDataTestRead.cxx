/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id$
/**
 * @file DataModelTestDataRead/src/AuxDataTestRead.cxx
 * @author snyder@bnl.gov
 * @date May 2014
 * @brief Algorithm to test reading @c DataVector with auxiliary data.
 */


#include "AuxDataTestRead.h"
#include "DataModelTestDataCommon/B.h"
#include "DataModelTestDataCommon/BAux.h"
#include "DataModelTestDataCommon/BAuxStandalone.h"
#include "DataModelTestDataCommon/BAuxVec.h"
#include "AthContainers/AuxTypeRegistry.h"
#include "AthContainers/AuxStoreInternal.h"
//#include "AthContainers/PackedElement.h"
#include "AthLinks/ElementLink.h"
#include "AthenaKernel/errorcheck.h"
#include "CxxUtils/StrFormat.h"
#include <memory>


namespace DMTest {


/**
 * @brief Constructor.
 * @param name The algorithm name.
 * @param svc The service locator.
 */
AuxDataTestRead::AuxDataTestRead (const std::string &name,
                                  ISvcLocator *pSvcLocator)
  : AthAlgorithm (name, pSvcLocator),
    m_count(0)
{
  declareProperty ("ReadPrefix",  m_readPrefix);
  declareProperty ("WritePrefix", m_writePrefix);
}
  

/**
 * @brief Algorithm initialization; called at the beginning of the job.
 */
StatusCode AuxDataTestRead::initialize()
{
  return StatusCode::SUCCESS;
}


/**
 * @brief Algorithm event processing.
 */
StatusCode AuxDataTestRead::execute()
{
  ++m_count;
  std::cout << m_count << "\n";

  const SG::AuxTypeRegistry& r = SG::AuxTypeRegistry::instance();

  static BAux::Accessor<int> anInt1 ("anInt1");
  static BAux::Accessor<float> aFloat1 ("aFloat1");
  static BAux::Accessor<ElementLink<BAuxVec> > anEL ("anEL");
  static BAux::Accessor<DMTest::B> aB ("aB");
  static BAux::Accessor<float> dFloat1 ("dFloat1");
  static BAux::Accessor<int> dInt1 ("dInt1");
  static BAux::Accessor<int> dInt2 ("dInt2");
  //static BAux::Accessor<SG::PackedElement<unsigned int> > pInt ("pint");
  //static BAux::Accessor<SG::PackedElement<float> > pFloat ("pfloat");
  //static BAux::Accessor<SG::PackedElement<std::vector<int> > > pvint ("pvint");
  //static BAux::Accessor<SG::PackedElement<std::vector<float> > > pvfloat ("pvfloat");
  static BAux::Accessor<unsigned int> pInt ("pint");
  static BAux::Accessor<float> pFloat ("pfloat");
  static BAux::Accessor<std::vector<int> > pvint ("pvint");
  static BAux::Accessor<std::vector<float> > pvfloat ("pvfloat");

  const BAuxVec* vec = 0;
  CHECK( evtStore()->retrieve (vec, m_readPrefix + "bauxvec") );

  // Ordering of auxid is not reliable.  Sort by name.
  std::vector<std::string> names;
  for (SG::auxid_t auxid : vec->getAuxIDs())
    names.push_back (r.getName(auxid));
  std::sort (names.begin(), names.end());
  for (const std::string& n : names)
    std::cout << n << " ";
  std::cout << "\n";
  for (const BAux* belt : *vec) {
    std::cout << " anInt1: " << anInt1(*belt)
              << " aFloat1: " << aFloat1(*belt)
              << " pInt: " << pInt(*belt)
              << " pFloat: " << CxxUtils::strformat ("%.2f", pFloat(*belt))
              << " aB: " << aB(*belt).m_x
              << " dFloat1: " << dFloat1(*belt);
    if (dInt1.isAvailable(*belt))
      std::cout << " dInt1: " << dInt1(*belt);
    if (dInt2.isAvailable(*belt))
      std::cout << " dInt2: " << dInt2(*belt);
    std::cout << "\n";

    const std::vector<int>& pvi = pvint(*belt);
    std::cout << "  pvInt: [";
    for (auto ii : pvi)
      std::cout << ii << " ";
    std::cout << "]\n";

    const std::vector<float>& pvf = pvfloat(*belt);
    std::cout << "  pvFloat: [";
    for (auto ii : pvf)
      std::cout << CxxUtils::strformat ("%.3f", ii) << " ";
    std::cout << "]\n";
  }

  const BAux* b = 0;
  CHECK( evtStore()->retrieve (b, m_readPrefix + "b") );
  std::cout << "b anInt1: " << anInt1(*b) 
            << " aFloat1: " << CxxUtils::strformat ("%.1f", aFloat1(*b))
            << " anEL: " << anEL(*b).dataID() << "[" << anEL(*b).index() << "]"
            << " aB: " << aB(*b).m_x 
            << " dFloat1: " << dFloat1(*b);
  if (dInt1.isAvailable(*b))
    std::cout << " dInt1: " << dInt1(*b);
  if (dInt2.isAvailable(*b))
    std::cout << " dInt2: " << dInt2(*b);
  std::cout << "\n";
    
  if (!m_writePrefix.empty()) {
    // Passing this as the third arg of record will make the object const.
    bool LOCKED = false;

    std::unique_ptr<BAuxVec> vecnew (new BAuxVec);
    std::unique_ptr<SG::AuxStoreInternal> store (new SG::AuxStoreInternal);
    vecnew->setStore (store.get());
    for (size_t i = 0; i < vec->size(); i++) {
      vecnew->push_back (new BAux);
      *vecnew->back() = *(*vec)[i];
    }
    CHECK (evtStore()->record (std::move(vecnew), m_writePrefix + "bauxvec", LOCKED));
    CHECK (evtStore()->record (std::move(store), m_writePrefix + "bauxvecAux.", LOCKED));

    std::unique_ptr<BAuxStandalone> bnew (new BAuxStandalone);
    bnew->makePrivateStore (*b);
    CHECK (evtStore()->record (std::move(bnew), m_writePrefix + "b", LOCKED));
  }

  return StatusCode::SUCCESS;
}


/**
 * @brief Algorithm finalization; called at the end of the job.
 */
StatusCode AuxDataTestRead::finalize()
{
  return StatusCode::SUCCESS;
}


} // namespace DMTest

