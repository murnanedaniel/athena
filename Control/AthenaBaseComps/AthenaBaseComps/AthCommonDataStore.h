/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/** @class  AthCommonDataStore
 *  @file   AthenaBaseComps/AthCommonDataStore.h
 *  @author Charles Leggett
 *  @date   June 2018
 *  @brief  Templated class that provides access to Athena event stores
 *          and ability to set data dependencies via Properites.
 *          Implemented to reduce code duplication in AthAlgorithm, 
 *          AthAlgTool, and AthReEntrantAlgorithm
 */

#ifndef ATHENABASECOMPS_ATHCOMMONDATASTORE_H
#define ATHENABASECOMPS_ATHCOMMONDATASTORE_H

#include <string>
#include <type_traits>

// Need to do this very early so parser for VarHandleKey picked up
#include "GaudiKernel/StatusCode.h"
namespace SG {
  class VarHandleKey;
  class VarHandleKeyArray;
  class VarHandleBase;
}
namespace Gaudi {
  namespace Parsers {
    StatusCode parse(SG::VarHandleKey& v, const std::string& s);
    StatusCode parse(SG::VarHandleKeyArray& v, const std::string& s);
    StatusCode parse(SG::VarHandleBase& v, const std::string& s);
  }
}

#include "AthenaBaseComps/AthMsgStreamMacros.h"
#include "AthenaBaseComps/AthCheckMacros.h"
#include "AthenaBaseComps/HandleClassifier.h"

#include "GaudiKernel/ServiceHandle.h"
#include "StoreGate/StoreGateSvc.h"
#include "StoreGate/VarHandleProperty.h"
#include "StoreGate/VarHandleKeyProperty.h"
#include "StoreGate/VarHandleKey.h"
#include "StoreGate/VarHandleBase.h"
#include "StoreGate/VarHandleKeyArray.h"
#include "StoreGate/VarHandleKeyArrayProperty.h"



template <class PBASE>
class AthCommonDataStore : public PBASE {
public:
  template <typename... T>
  AthCommonDataStore(const std::string& name, T... args)
    : PBASE(name, args...),
      m_evtStore    ( "StoreGateSvc/StoreGateSvc",  name ),
      m_detStore    ( "StoreGateSvc/DetectorStore", name ),
      m_varHandleArraysDeclared (false)
  {

    this->declareProperty( "EvtStore",
                     m_evtStore = StoreGateSvc_t ("StoreGateSvc", name),
                     "Handle to a StoreGateSvc instance: it will be used to "
                     "retrieve data during the course of the job" );
    
    this->declareProperty( "DetStore",
                     m_detStore = StoreGateSvc_t ("StoreGateSvc/DetectorStore", name),
                     "Handle to a StoreGateSvc/DetectorStore instance: it will be used to "
                     "retrieve data during the course of the job" );

    auto props = this->getProperties();
    for( Gaudi::Details::PropertyBase* prop : props ) {
      if (prop->name() == "ExtraOutputs" || prop->name() == "ExtraInputs") {
        prop->declareUpdateHandler
          (&AthCommonDataStore<PBASE>::extraDeps_update_handler, this);
      }
    }
  }


  /** @brief The standard @c StoreGateSvc (event store)
   * Returns (kind of) a pointer to the @c StoreGateSvc
   */
  ServiceHandle<StoreGateSvc>& evtStore() { return m_evtStore; }

  /** @brief The standard @c StoreGateSvc (event store)
   * Returns (kind of) a pointer to the @c StoreGateSvc
   */
  const ServiceHandle<StoreGateSvc>& evtStore() const { return m_evtStore; }

  /** @brief The standard @c StoreGateSvc/DetectorStore
   * Returns (kind of) a pointer to the @c StoreGateSvc
   */
  const ServiceHandle<StoreGateSvc>& detStore() const { return m_detStore; }


  /**
   * @brief Perform system initialization for an algorithm.
   *
   * We override this to declare all the elements of handle key arrays
   * at the end of initialization.
   * See comments on updateVHKA.
   */
  virtual StatusCode sysInitialize() override;

  /**
   * @brief Handle START transition.
   *
   * We override this in order to make sure that conditions handle keys
   * can cache a pointer to the conditions container.
   */
  virtual StatusCode sysStart() override;
  
  /**
   * @brief Return this algorithm's input handles.
   *
   * We override this to include handle instances from key arrays
   * if they have not yet been declared.
   * See comments on updateVHKA.
   */
  virtual std::vector<Gaudi::DataHandle*> inputHandles() const override;


  /**
   * @brief Return this algorithm's output handles.
   *
   * We override this to include handle instances from key arrays
   * if they have not yet been declared.
   * See comments on updateVHKA.
   */
  virtual std::vector<Gaudi::DataHandle*> outputHandles() const override;



  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  public:
  /////////////////////////////////////////////////////////////////
  //
  //// Enable use of Gaudi::Property<Foo> m_foo {this,"NAME",init,"doc"};
  //   style properties in AthAlgorithms
  //

  template <class T>
  Gaudi::Details::PropertyBase& declareProperty(Gaudi::Property<T> &t) {
    typedef typename SG::HandleClassifier<T>::type htype;
    return AthCommonDataStore<PBASE>::declareGaudiProperty(t, htype());
  }

  private:
  /**
   * @brief specialization for handling Gaudi::Property<SG::VarHandleKey>
   *
   */
  template <class T>
  Gaudi::Details::PropertyBase& declareGaudiProperty(Gaudi::Property<T> &hndl,
                                 const SG::VarHandleKeyType&)
  {
    return *AthCommonDataStore<PBASE>::declareProperty(hndl.name(),
                                                       hndl.value(), 
                                                       hndl.documentation());

  }

  /**
   * @brief specialization for handling Gaudi::Property<SG::VarHandleKeyArray>
   *
   */
  template <class T>
  Gaudi::Details::PropertyBase& declareGaudiProperty(Gaudi::Property<T> &hndl, 
                                 const SG::VarHandleKeyArrayType&)
  {
    return *AthCommonDataStore<PBASE>::declareProperty(hndl.name(),
                                                       hndl.value(), 
                                                       hndl.documentation());

  }

  /**
   * @brief specialization for handling Gaudi::Property<SG::VarHandleBase>
   *
   */
  template <class T>
  Gaudi::Details::PropertyBase& declareGaudiProperty(Gaudi::Property<T> &hndl, 
                                 const SG::VarHandleType&)
  {
    return *AthCommonDataStore<PBASE>::declareProperty(hndl.name(),
                                                       hndl.value(), 
                                                       hndl.documentation());
  }


  /**
   * @brief specialization for handling everything that's not a
   * Gaudi::Property<SG::VarHandleKey> or a <SG::VarHandleKeyArray>
   *
   */
  template <class T>
  Gaudi::Details::PropertyBase& declareGaudiProperty(Gaudi::Property<T> &t, const SG::NotHandleType&)
  {
    return PBASE::declareProperty(t);
  }


  /////////////////////////////////////////////////////////////////
  //
  //// For automatic registration of Handle data products
  //

public:
  /**
   * @brief Declare a new Gaudi property.
   * @param name Name of the property.
   * @param hndl Object holding the property value.
   * @param doc Documentation string for the property.
   *
   * This is the version for types that derive from @c SG::VarHandleKey.
   * The property value object is put on the input and output lists as
   * appropriate; then we forward to the base class.
   */
  Gaudi::Details::PropertyBase* declareProperty(const std::string& name,
                            SG::VarHandleKey& hndl,
                            const std::string& doc,
                            const SG::VarHandleKeyType&)
  {
    this->declare(hndl);
    hndl.setOwner(this);

    return PBASE::declareProperty(name,hndl,doc);
  }


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  /**
   * @brief Declare a new Gaudi property.
   * @param name Name of the property.
   * @param hndl Object holding the property value.
   * @param doc Documentation string for the property.
   *
   * This is the version for types that derive from @c SG::VarHandleBase.
   * The property value object is put on the input and output lists as
   * appropriate; then we forward to the base class.
   */
  Gaudi::Details::PropertyBase* declareProperty(const std::string& name,
                            SG::VarHandleBase& hndl,
                            const std::string& doc,
                            const SG::VarHandleType&)
  {
    this->declare(hndl.vhKey());
    hndl.vhKey().setOwner(this);

    return PBASE::declareProperty(name,hndl,doc);
  }


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  Gaudi::Details::PropertyBase* declareProperty(const std::string& name,
                            SG::VarHandleKeyArray& hndArr,
                            const std::string& doc,
                            const SG::VarHandleKeyArrayType&)
  {

    // std::ostringstream ost;
    // ost << Algorithm::name() << " VHKA declareProp: " << name 
    //     << " size: " << hndArr.keys().size() 
    //     << " mode: " << hndArr.mode() 
    //     << "  vhka size: " << m_vhka.size()
    //     << "\n";
    // debug() << ost.str() << endmsg;

    hndArr.setOwner(this);
    m_vhka.push_back(&hndArr);

    Gaudi::Details::PropertyBase* p =  PBASE::declareProperty(name, hndArr, doc);
    if (p != 0) {
      p->declareUpdateHandler(&AthCommonDataStore<PBASE>::updateVHKA, this);
    } else {
      ATH_MSG_ERROR("unable to call declareProperty on VarHandleKeyArray " 
                    << name);
    }

    return p;

  }


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

  // Since the contents of the VarHandleKeyArrays have not been read 
  // in from the configurables by the time that declareProperty is
  // executed, we must cache them and loop through them later to
  // register the data dependencies.
  //
  // However, we cannot actually call declare() on the key instances
  // until we know that the vector cannot change size anymore --- otherwise,
  // the pointers given to declare() may become invalid.  That basically means
  // that we can't call declare() until the derived class's initialize()
  // completes.  So instead of doing it here (which would be too early),
  // we override sysInitialize() and do it at the end of that.  But,
  // Algorithm::sysInitialize() wants to have the handle lists after initialize()
  // completes in order to do dependency analysis.  It gets these lists
  // solely by calling inputHandles() and outputHandles(), so we can get this
  // to work by overriding those methods and adding in the current contents
  // of the arrays.

  void updateVHKA(Gaudi::Details::PropertyBase& /*p*/) {
    // debug() << "updateVHKA for property " << p.name() << " " << p.toString() 
    //         << "  size: " << m_vhka.size() << endmsg;
    for (auto &a : m_vhka) {
      std::vector<SG::VarHandleKey*> keys = a->keys();
      for (auto k : keys) {
        k->setOwner(this);
      }
    }
  }

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


  /**
   * @brief Declare a new Gaudi property.
   * @param name Name of the property.
   * @param property Object holding the property value.
   * @param doc Documentation string for the property.
   *
   * This is the generic version, for types that do not derive
   * from  @c SG::VarHandleKey.  It just forwards to the base class version
   * of @c declareProperty.
   */
  template <class T>
  Gaudi::Details::PropertyBase* declareProperty(const std::string& name,
                            T& property,
                            const std::string& doc,
                            const SG::NotHandleType&)
  {
    return PBASE::declareProperty(name, property, doc);
  }


  /**
   * @brief Declare a new Gaudi property.
   * @param name Name of the property.
   * @param property Object holding the property value.
   * @param doc Documentation string for the property.
   *
   * This dispatches to either the generic @c declareProperty or the one
   * for VarHandle/Key/KeyArray.
   */
  template <class T>
  Gaudi::Details::PropertyBase* declareProperty(const std::string& name,
                            T& property,
                            const std::string& doc="none") 
  {
    typedef typename SG::HandleClassifier<T>::type htype;
    return declareProperty (name, property, doc, htype());
  }

  

protected:
    /// remove all handles from I/O resolution
  void renounceArray( SG::VarHandleKeyArray& handlesArray ) {
    handlesArray.renounce();
  }

private:
  typedef ServiceHandle<StoreGateSvc> StoreGateSvc_t;
  /// Pointer to StoreGate (event store by default)
  StoreGateSvc_t m_evtStore;

  /// Pointer to StoreGate (detector store by default)
  StoreGateSvc_t m_detStore;


private:
  // to keep track of VarHandleKeyArrays for data dep registration
  std::vector<SG::VarHandleKeyArray*> m_vhka;
  bool m_varHandleArraysDeclared;


protected:
  /**
 * @brief Add StoreName to extra input/output deps as needed
 *
 * use the logic of the VarHandleKey to parse the DataObjID keys
 * supplied via the ExtraInputs and ExtraOuputs Properties to add
 * the StoreName if it's not explicitly given
 */
  void extraDeps_update_handler( Gaudi::Details::PropertyBase& ExtraDeps ); 


};

#include "AthCommonDataStore.icc"

#endif
