/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id$
/**
 * @file AthContainers/test/ViewVector_test.cxx
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2016
 * @brief Tests for ViewVector.
 */


#undef NDEBUG
#include <boost/preprocessor/stringize.hpp>
#include "AthContainers/ViewVector.h"
#include "AthContainers/ConstDataVector.h"
#include "AthContainers/DataVector.h"
#include "SGTools/BaseInfo.h"
#include "TestTools/expect_exception.h"
#include <iostream>
#include <cassert>


class C {};
typedef DataVector<C> DV;
typedef ViewVector<DV> VV;


void checkit (const VV& vv, const DV& dv, int n = -1)
{
  if (n < 0)
    n = dv.size();
  assert (vv.ownPolicy() == SG::VIEW_ELEMENTS);
  assert (static_cast<int>(vv.size()) == n);
  for (int i=0; i < n; i++)
    assert (vv[i] == dv[i]);
}

void test1()
{
  std::cout << "test1\n";

  DV dv;
  for (int i=0; i < 10; i++)
    dv.push_back (new C);

  VV vv1;
  assert (vv1.ownPolicy() == SG::VIEW_ELEMENTS);
  VV vv2 (SG::VIEW_ELEMENTS);
  assert (vv2.ownPolicy() == SG::VIEW_ELEMENTS);
  EXPECT_EXCEPTION (SG::ExcViewVectorNotView, VV vvx(SG::OWN_ELEMENTS));

  VV vv3 (10);
  assert (vv3.ownPolicy() == SG::VIEW_ELEMENTS);
  assert (vv3.size() == 10);
  VV vv4 (10, SG::VIEW_ELEMENTS);
  assert (vv4.ownPolicy() == SG::VIEW_ELEMENTS);
  assert (vv4.size() == 10);
  EXPECT_EXCEPTION (SG::ExcViewVectorNotView, VV(10, SG::OWN_ELEMENTS));

  VV vv5 (dv.begin(), dv.end());
  checkit (vv5, dv);

  VV vv6 (vv5);
  checkit (vv6, dv);

  VV vv7 (std::move (vv5));
  assert (vv5.size() == 0);
  checkit (vv7, dv);

  VV vv8 (dv);
  checkit (vv8, dv);

  DV dv2 (dv);
  VV vv9 (std::move(dv2));
  assert (dv2.size() == 0);
  checkit (vv9, dv);

  DV dv3;
  EXPECT_EXCEPTION (SG::ExcViewVectorNotView, VV(std::move(dv3)));

  VV vv10 {dv[0], dv[1]};
  checkit (vv10, dv, 2);

  vv2 = vv10;
  checkit (vv2, dv, 2);

  vv2 = dv;
  checkit (vv2, dv);

  vv2 = std::move(vv10);
  assert (vv10.size() == 0);
  checkit (vv2, dv, 2);

  dv2 = dv;
  vv2 = std::move(dv2);
  assert (dv2.size() == 0);
  checkit (vv2, dv);
  EXPECT_EXCEPTION (SG::ExcViewVectorNotView, vv2 = std::move(dv3));

  vv2 = {dv[0], dv[1]};
  checkit (vv2, dv, 2);

  vv2.clear();
  checkit (vv2, dv, 0);
  vv2.resize(10);
  vv2.clear (SG::VIEW_ELEMENTS);
  checkit (vv2, dv, 0);
  EXPECT_EXCEPTION (SG::ExcViewVectorNotView, vv2.clear (SG::OWN_ELEMENTS));

  vv2.resize(10);
  vv2.clear (SG::VIEW_ELEMENTS, SG::NEVER_TRACK_INDICES);
  checkit (vv2, dv, 0);
  EXPECT_EXCEPTION (SG::ExcViewVectorNotView, vv2.clear (SG::OWN_ELEMENTS,
                                                         SG::NEVER_TRACK_INDICES));
}


void test2()
{
  std::cout << "test2\n";

  const SG::BaseInfoBase* bib = SG::BaseInfoBase::find (typeid (VV));
  assert (bib != nullptr);
  assert (bib->is_base (typeid (DV)));
}


class D {};
typedef ViewVector<DataVector<D> > DView;
VIEWVECTOR_CLASS_DEF(DView, 9823475)


void test3()
{
  std::cout << "test3\n";

  EXPECT_EXCEPTION (SG::ExcMissingViewVectorCLID, ClassID_traits<VV>::ID());
  assert (ClassID_traits<DView>::ID() == 9823475);
  assert (ClassID_traits<DView>::typeName() == "DView");
  assert (ClassID_traits<DView>::typeInfo() == typeid(DView));
}


int main()
{
  test1();
  test2();
  test3();
  return 0;
}


class X {};
template class DataVector<X>;
template class ViewVector<DataVector<X> >;
template class ConstDataVector<ViewVector<DataVector<X> > >;
template class ViewVector<ConstDataVector<DataVector<X> > >;
