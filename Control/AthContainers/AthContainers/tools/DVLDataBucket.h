// This file's extension implies that it's C, but it's really -*- C++ -*-.

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// $Id: DVLDataBucket.h,v 1.2 2008-06-17 00:44:17 ssnyder Exp $

/**
 * @file  AthContainers/tools/DVLDataBucket.h
 * @author scott snyder
 * @date Mar 2008
 * @brief A @c DataBucket specialized for @c DataVector/@c DataList.
 */


#ifndef ATHCONTAINERS_TOOLS_DVLDATABUCKET_H
#define ATHCONTAINERS_TOOLS_DVLDATABUCKET_H

#include "SGTools/DataBucket.h"
#include "AthContainers/tools/DVLInfo.h"
#include <memory>


namespace SG {


/**
 * @brief A @c DataBucket specialized for @c DataVector/@c DataList.
 *
 * @c DataVector and @c DataList have some special conversion requirements.
 *
 * In particular, we may have a situation where we're asked to convert
 * from @c DataVector\<D> to @c DataVector\<B>, where these two classes
 * are not related by inheritance, but @c D derived from @B.
 * (In practice, this arises when one has multiple inheritance,
 * and can't mirror both arms with @c DataVector/@c DataList.)
 *
 * This @c DataBucket allows such conversions to happen by making
 * a copy of the container (as a view container) and converting the contained
 * pointers.  The new instance remains owned by the @c DataBucket.
 */
template <class T>
class DVLDataBucket
  : public SG::DataBucket<T>
{
public:
  /**
   * @brief Default constructor.
   */
  DVLDataBucket() {}


  /**
   * @brief Constructor from a payload object.
   * @param data Object to hold in the bucket.
   */
  DVLDataBucket (T* data);


#if __cplusplus > 201100
  /**
   * @brief Constructor from a payload object.
   * @param data Object to hold in the bucket.
   */
  DVLDataBucket(std::unique_ptr<T> data);
#endif


  /**
   * @brief Copy constructor.
   * @param other Bucket to copy.
   */
  DVLDataBucket (const DVLDataBucket& other);


  /**
   * @brief Destructor.
   */
  ~DVLDataBucket();


  /**
   * @brief Return the contents of the @c DataBucket,
   *        converted to type given by @a clid.  Note that only
   *        derived->base conversions are allowed here.
   * @param clid The class ID to which to convert.
   * @param irt To be called if we make a new instance.
   * @param isConst True if the object being converted is regarded as const.
   */
  virtual void* cast (CLID clid, IRegisterTransient* irt = 0,
                      bool isConst = true) const;
    

  /**
   * @brief Return the contents of the @c DataBucket,
   *        converted to type given by @a std::type_info.  Note that only
   *        derived->base conversions are allowed here.
   * @param clid The @a std::type_info of the type to which to convert.
   * @param irt To be called if we make a new instance.
   * @param isConst True if the object being converted is regarded as const.
   */
  virtual void* cast (const std::type_info& tinfo,
                      IRegisterTransient* irt = 0,
                      bool isConst = true) const;


  /**
   * @brief Return a new @c DataBucket whose payload has been cloned from the
   *        original one.
   */
  virtual DVLDataBucket* clone() const;


private:
  typedef std::pair<DataModel_detail::DVLInfoBase*, void*> ent_t;
  typedef std::vector<ent_t> vec_t;
  mutable vec_t m_copies;

  // Avoid coverity warning.
  DVLDataBucket& operator= (const DVLDataBucket&);
};


} // namespace SG


#include "AthContainers/tools/DVLDataBucket.icc"


#endif // not ATHCONTAINERS_TOOLS_DVLDATABUCKET_H
