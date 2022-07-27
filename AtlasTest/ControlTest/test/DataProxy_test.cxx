/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#undef NDEBUG
#include "SGTools/DataProxy.h"

#include "TestTools/initGaudi.h"
#include "ToyConversion/FooBar.h"
#include "ToyConversion/ToyConversionSvc.h"
#include "GaudiKernel/GenericAddress.h"
#include "AthenaKernel/IProxyProviderSvc.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "AthenaKernel/StorableConversions.h"
#include "AthenaKernel/CLASS_DEF.h"
#include "SGTools/T2pMap.h"
#include "SGTools/TransientAddress.h"
#include <iostream>

#include <cassert>

const bool CREATEIF(true);

struct Bla{
  Bla(): i(0) {}
  Bla(int m): i(m) {}
  int i;
  void doNothing() const {};
};
CLASS_DEF( Bla, 7890, 0) 

using std::cerr;
using std::cout;
using std::endl;

using SG::DataProxy_cast;
using SG::DataProxy;
using SG::TransientAddress;

using namespace Athena_test;

int main() {
  cout << "*** DataProxy_test BEGINS ***" <<endl;

  ISvcLocator* pSvcLoc;
  if (!initGaudi("DataProxy_test.txt", pSvcLoc)) {
    cerr << "This test can not be run" << endl;
    return 0;
  }  

  assert(pSvcLoc);

  DataProxy emptyProxy;
  assert( !emptyProxy.isValid() );

  //cerr << "Now we expect to see a warning message:" << endl
  //     << "----Warning Message Starts--->>" << endl; 
  assert( !emptyProxy.accessData() );
  //cerr << "<<---Warning Message Ends-------" << endl;
  assert( emptyProxy.errNo() == DataProxy::NOCNVSVC );
  assert( !emptyProxy.object() );
  assert( !emptyProxy.address() );
  assert( emptyProxy.name().empty() );
  assert( emptyProxy.identifier().empty() );
  assert( 0 == emptyProxy.clID() );
  assert( 0 == emptyProxy.transientID().size() );
  assert( emptyProxy.isResetOnly() );
  assert( !emptyProxy.isConst() );
  

  Bla *pBla(new Bla(77));
  DataProxy transientProxy(SG::asStorable(pBla),
			   new TransientAddress(ClassID_traits<Bla>::ID(), "bla"));
  assert( transientProxy.isValid() );
  assert( transientProxy.accessData() );
  assert( transientProxy.object() );
  assert( !transientProxy.address() );
  assert( transientProxy.name() == "bla" );
  assert( transientProxy.identifier() == "bla" );
  assert( ClassID_traits<Bla>::ID() == transientProxy.clID() );
  assert( transientProxy.transientID().size() == 1);
  assert( transientProxy.isResetOnly() );
  assert( !transientProxy.isConst() );
  assert( (DataProxy_cast<Bla>(&transientProxy))->i == 77 );


  IConversionSvc* pIConvSvc(nullptr);
  StatusCode cs((pSvcLoc->service("EventPersistencySvc", pIConvSvc, CREATEIF)));
  assert(cs.isSuccess());
  assert(pIConvSvc);
  // create a transient address with IOA.
  TransientAddress* tGen = new TransientAddress(ClassID_traits<Foo>::ID(),
						"foo",
						new GenericAddress(ToyConversionSvc::storageType(), ClassID_traits<Foo>::ID(), "foo"));

  SG::T2pMap t2p;
  DataProxy addressProxy(tGen, pIConvSvc);
  addressProxy.setT2p(&t2p);

  assert( addressProxy.isValid() );
  assert( !addressProxy.object() );
  assert( addressProxy.address() );
  //the proxy provider mechanism is tested in ProxyProviderSvc_test
  assert( addressProxy.accessData() );
  assert( addressProxy.object() );
  assert( addressProxy.name() == "foo" );
  assert( addressProxy.identifier() == "foo" );
  assert( ClassID_traits<Foo>::ID() == addressProxy.clID() );
  assert( addressProxy.transientID().size() == 1);
  assert( addressProxy.isResetOnly() );
  assert( !addressProxy.isConst() );
  Foo* fptr(DataProxy_cast<Foo>(&addressProxy));
  fptr->doNothing(); //remove warning


  IProxyProviderSvc* pIPPSvc(nullptr);
  StatusCode psc(pSvcLoc->service("ProxyProviderSvc", pIPPSvc, CREATEIF));
  assert( psc.isSuccess() );
  assert( pIPPSvc );

  TransientAddress* tad = new TransientAddress(ClassID_traits<Bla>::ID(), "bla");

  DataProxy identifiedProxy(tad, pIConvSvc);

  assert( !identifiedProxy.isValid() );
  assert( !identifiedProxy.object() );
  assert( !identifiedProxy.address() );
  //the proxy provider mechanism is tested in ProxyProviderSvc_test
  cerr << "Now we expect to see a warning message:" << endl
       << "----Warning Message Starts--->>" << endl; 
  assert( !identifiedProxy.accessData() ); 
  assert( identifiedProxy.errNo() == DataProxy::NOIOA );
  cerr << "<<---Warning Message Ends-------" << endl;
  assert( identifiedProxy.name() == "bla" );
  assert( identifiedProxy.identifier() == "bla" );
  assert( ClassID_traits<Bla>::ID() == identifiedProxy.clID() );
  assert( identifiedProxy.transientID().size() == 1);
  assert( identifiedProxy.isResetOnly() );
  assert( !identifiedProxy.isConst() );


  cout << "*** DataProxy_test OK ***" <<endl;
  return 0;
}












