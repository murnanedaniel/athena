//  ConfigException.h
//  TopoCore
//  Created by Joerg Stelzer on 11/18/12.
//  Copyright (c) 2012 Joerg Stelzer. All rights reserved.

#ifndef __TopoCore__Exception__
#define __TopoCore__Exception__

#include <iostream>
#include <sstream>

#define TCS_EXCEPTION(MSG) \
{ std::stringstream o;   \
o << MSG;\
throw TCS::Exception(o.str()); }

namespace TCS {
   
   class Exception : virtual public std::exception {
   public:
      enum type_t {
         CONFIG,
         RUNTIME
      };
      
      Exception(const std::string& msg) :
         m_msg(msg)
      {}
      
			virtual ~Exception() throw() {};

      virtual char const* what() const throw() { return m_msg.data(); }
      
   private:
      std::string m_msg;
      /*type_t m_type; // not used yet, commenting out to suppress compilation warning */
   };
   
}

#endif /* defined(__TopoCore__ConfigException__) */
