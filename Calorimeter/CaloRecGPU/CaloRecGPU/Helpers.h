// Dear emacs, this is -*- c++ -*-
//
// Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration

#ifndef CALORECGPU_HELPERS_H
#define CALORECGPU_HELPERS_H

#include <utility>
#include <type_traits>
#include <cstring>
//For memcpy, of all things...
#include <string>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <thread>
#include <mutex>
#include <memory>
#include <vector>

#if __cpp_lib_math_constants
  #include <numbers>
#endif
//This is the best way to have pi,
//but we provide a more manual alternative.
//Of course, there's also M_PI,
//but we wanted to ensure the type matched
//to prevent any GPU-based shenanigans.

namespace CaloRecGPU
{

#ifndef CUDA_AVAILABLE

  #ifdef __CUDA_ARCH__
    #define CUDA_AVAILABLE 1
  #elif __CUDA__
    #define CUDA_AVAILABLE 1
  #elif __CUDACC__
    #define CUDA_AVAILABLE 1
  #else
    #define CUDA_AVAILABLE 0
  #endif

#endif

#if CUDA_AVAILABLE



  /*!
    \brief Provides a way to report errors in \c CUDA kernels.

    \remark Standard CUDA definition that can be found almost everywhere...
  */
  inline void CUDA_gpu_assert(cudaError_t code, const char * file, int line, bool abort = true)
  {
    if (code != cudaSuccess)
      {
        printf("GPU Error: %s (%s %d)\n", cudaGetErrorString(code), file, line);
        if (abort)
          {
            exit(code);
          }
      }
  }

#define CUDA_HOS_DEV __host__ __device__

  /*!
    \brief Wraps up a `CUDA_gpu_assert` using `__FILE__` and `__LINE__`
    to improve error reporting (and debugging) capabilities.
  */
#define CUDA_ERRCHECK(ans) { CUDA_gpu_assert((ans), __FILE__, __LINE__); }

#else

#define CUDA_HOS_DEV

#endif

  namespace CUDA_Helpers
  {
    /*!
      \brief Allocates and returns the address of \p num bytes from GPU memory.
    */
    void * allocate(const size_t num);

    /*!
      \brief Deallocates \p address in GPU memory.
    */
    void deallocate(void * address);

    /*!
      \brief Copies \p num byte from \p source in GPU memory to \p dest in CPU memory.
    */
    void GPU_to_CPU(void * dest, const void * const source, const size_t num);

    /*!
      \brief Copies \p num byte from \p source in CPU memory to \p dest in GPU memory.
    */
    void CPU_to_GPU(void * dest, const void * const source, const size_t num);

    /*!
      \brief Copies \p num byte from \p source to \p dest, both in GPU memory.
    */
    void GPU_to_GPU(void * dest, const void * const source, const size_t num);
  }

  namespace Helpers
  {

    /// \brief Returns the ceiling of num/denom, with proper rounding.
    inline constexpr int int_ceil_div(const int num, const int denom)
    {
      return num / denom + (num % denom != 0);
    }

    /// \brief Returns 2 to the power of \p exp.
    template <class Base = float, class Exp = int>
    inline constexpr Base compile_time_pow2(const Exp exp)
    {
      Base ret = 1;
      if (exp < 0)
        {
          for (Exp i = 0; i < -exp; ++i)
            {
              ret /= Base(2);
            }
        }
      else
        {
          for (Exp i = 0; i < exp; ++i)
            {
              ret *= Base(2);
            }
        }
      return ret;
    }
    //Though we could possibly bit-hack stuff due to IEEE-754 reliance elsewhere,
    //it's not valid and type-safe C++...
    //Since it's compile-time, this being a trifle slower is meaningless.


    ///! \brief Just a wapper around constants with fallback.
    namespace Constants
    {
#ifdef __cpp_lib_math_constants
  template <class T>
  inline constexpr T pi = std::numbers::pi_v<T>;
#else
  template <class T>
  inline constexpr T pi = T(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624L);
#endif
    }

    ///! Holds dummy classes just to identify the place in which memory lives.
    namespace MemoryContext
    {
      struct CPU
      {
        constexpr static char const * name = "CPU";
      };
      struct CUDAGPU
      {
        constexpr static char const * name = "CUDA GPU";
      };
    }

    ///! Handles allocation of a type \p T, using \p indexer as the integer type to indicate sizes.
    template <class T, class indexer>
    class MemoryManagement
    {
     private:
      template <class C, class dummy = void> struct unary_helper;

      template <class dummy> struct unary_helper<MemoryContext::CPU, dummy>
      {
        static inline T * allocate(const indexer size)
        {
          return new T[size];
        }

        static inline void deallocate(T *& arr)
        {
          delete[] arr;
        }

      };

      template <class dummy> struct unary_helper<MemoryContext::CUDAGPU, dummy>
      {
        static inline T * allocate(const indexer size)
        {
          return static_cast<T *>(CUDA_Helpers::allocate(sizeof(T) * size));
        }

        static inline void deallocate(T *& arr)
        {
          CUDA_Helpers::deallocate(arr);
        }
      };

      template <class C1, class C2, class dummy = void> struct copy_helper;

      template <class dummy> struct copy_helper<MemoryContext::CPU, MemoryContext::CPU, dummy>
      {
        static inline void copy (T * dest, const T * const source, const indexer sz)
        {
          std::memcpy(dest, source, sizeof(T) * sz);
        }
      };

      template <class dummy> struct copy_helper<MemoryContext::CPU, MemoryContext::CUDAGPU, dummy>
      {
        static inline void copy (T * dest, const T * const source, const indexer sz)
        {
          CUDA_Helpers::GPU_to_CPU(dest, source, sizeof(T) * sz);
        }
      };

      template <class dummy> struct copy_helper<MemoryContext::CUDAGPU, MemoryContext::CUDAGPU, dummy>
      {
        static inline void copy (T * dest, const T * const source, const indexer sz)
        {
          CUDA_Helpers::GPU_to_GPU(dest, source, sizeof(T) * sz);
        }
      };

      template <class dummy> struct copy_helper<MemoryContext::CUDAGPU, MemoryContext::CPU, dummy>
      {
        static inline void copy (T * dest, const T * const source, const indexer sz)
        {
          CUDA_Helpers::CPU_to_GPU(dest, source, sizeof(T) * sz);
        }
      };


      template <class C1, class C2, class dummy = void> struct move_helper;

      template <class C1, class C2, class dummy> struct move_helper
      {
        inline static void move(T *& dest,  T *& source, const indexer sz)
        {
          dest = MemoryManagement<T, indexer>::template allocate<C1>(sz);
          MemoryManagement<T, indexer>::template copy<C1, C2>(dest, source, sz);
          MemoryManagement<T, indexer>::template deallocate<C2>(source);
        }
      };

      template <class C, class dummy> struct move_helper<C, C, dummy>
      {
        inline static void move(T *& dest,  T *& source, const indexer)
        {
          dest = source;
          source = nullptr;
        }
      };

     public:
      ///! \brief Allocates \p size elements from memory context \p Context
      template <class Context> static inline T * allocate(const indexer size)
      {
        T * ret = nullptr;
        if (size > 0)
          {
            ret = unary_helper<Context>::allocate(size);
          }
#if CALOREC_HELPERS_DEBUG
        std::cerr << "ALLOCATED " << size << " in " << Context::name << ": " << ret << std::endl;
#endif
        return ret;
      }

      ///! \brief Deallocates \p arr from memory context \p Context
      template <class Context> static inline void deallocate(T *& arr)
      {
        if (arr == nullptr)
        //This check is to ensure the code behaves on non-CUDA enabled platforms
        //where some destructors might still be called with nullptr.
        {
          return;
        }
        unary_helper<Context>::deallocate(arr);
#if CALOREC_HELPERS_DEBUG
        std::cerr << "DEALLOCATED in " << Context::name << ": " << arr << std::endl;
#endif
        arr = nullptr;
      }


      ///! \brief Copies \p sz bytes from \p source in \p SourceContext to \p dest in \p DestContext
      template <class DestContext, class SourceContext>
      static inline void copy(T * dest, const T * const source, const indexer sz)
      {
        if (sz > 0 && source != nullptr)
          {
            copy_helper<DestContext, SourceContext>::copy(dest, source, sz);
          }
#if CALOREC_HELPERS_DEBUG
        std::cerr << "COPIED " << sz << " from " << SourceContext::name << " to " << DestContext::name << ": " << source << " to " << dest << std::endl;
#endif
      }


      /**! \brief Moves \p sz bytes from \p source in \p SourceContext to \p dest in \p DestContext
            (performing the necessary data transfers and deallocations when the contexts are different).
            \p source is set to null and \p dest may be nulled if \p sz is 0 or \p source is null too.

      */
      template <class DestContext, class SourceContext>
      static inline void move(T *& dest,  T *& source, const indexer sz)
      {
#if CALOREC_HELPERS_DEBUG
        std::cerr << "MOVED " << sz << " from " << SourceContext::name << " to " << DestContext::name << ": " << source << " to " << dest;
#endif
        if (sz > 0 && source != nullptr)
          {
            move_helper<DestContext, SourceContext, std::is_same<DestContext, SourceContext>>::move(dest, source, sz);
          }
        else
          {
            dest = nullptr;
            deallocate<SourceContext>(source);
          }
#if CALOREC_HELPERS_DEBUG
        std::cerr << " | " << source << " to " << dest << std::endl;
#endif
      }

    };

    /*! Holds a run-time amount of objects of type \T, measuring sizes with \p indexer,
        in memory context \p Context. Automatically handles memory transfers to and from
        other contexts in the constructors and assignment operators.
        If \p hold_arrays is \p false, functions as a simple non-holding pointer,
        just keeping track of the \p Context for possible copies/transfers between memory contexts.
    */
    template <class T, class indexer, class Context, bool hold_arrays = true>
    class SimpleContainer;

    template <class T, class indexer, class Context>
    class SimpleContainer<T, indexer, Context, true>
    {
      static_assert(std::is_trivially_copyable<T>::value, "SimpleContainer only works with a trivially copyable type.");
      T * m_array;
      indexer m_size;

      template <class a, class b, class c, bool d> friend class SimpleContainer;

      using Manager = MemoryManagement<T, indexer>;

     public:

      CUDA_HOS_DEV inline indexer size() const
      {
        return m_size;
      }

      CUDA_HOS_DEV inline T & operator[] (const indexer i)
      {
        return m_array[i];
      }

      CUDA_HOS_DEV inline const T & operator[] (const indexer i) const
      {
        return m_array[i];
      }

      inline void clear()
      {
        Manager::deallocate(m_array);
        m_size = 0;
      }

      inline void resize(const indexer new_size)
      {
        if (new_size == 0)
          {
            clear();
          }
        else if (new_size != m_size)
          {
            T * temp = m_array;
            m_array = Manager::template allocate<Context>(new_size);
            Manager::template copy<Context, Context>(m_array, temp, std::min(m_size, new_size));
            Manager::template deallocate<Context>(temp);
            m_size = new_size;
          }
      }

      SimpleContainer() : m_array(nullptr), m_size(0)
      {
      }

      SimpleContainer(const indexer sz)
      {
        m_array = Manager::template allocate<Context>(sz);
        m_size = sz;
      }

      /*!
        \warning We assume the pointer is in a valid memory location!
      */
      SimpleContainer(T * other_array, const indexer sz)
      {
        m_array = Manager::template allocate<Context>(sz);
        Manager::template copy<Context, Context>(m_array, other_array, sz);
        m_size = sz;
      }

      SimpleContainer(const SimpleContainer & other)
      {
        m_size = other.m_size;
        m_array = Manager::template allocate<Context>(m_size);
        Manager::template copy<Context, Context>(m_array, other.m_array, m_size);
      }

      SimpleContainer(SimpleContainer && other)
      {
        m_size = other.m_size;
        m_array = nullptr;
        Manager::template move<Context, Context>(m_array, other.m_array, m_size);
        other.m_size = 0;
      }

      template <class other_indexer, class other_context, bool other_hold>
      SimpleContainer(const SimpleContainer<T, other_indexer, other_context, other_hold> & other)
      {
        m_size = other.m_size;
        m_array = Manager::template allocate<Context>(m_size);
        Manager::template copy<Context, other_context>(m_array, other.m_array, m_size);
      }

      template <class other_indexer, class other_context>
      SimpleContainer(SimpleContainer<T, other_indexer, other_context, true> && other)
      {
        m_size = other.m_size;
        m_array = nullptr;
        Manager::template move<Context, other_context>(m_array, other.m_array, m_size);
        other.m_size = 0;
      }

      SimpleContainer & operator= (const SimpleContainer & other)
      {
        if (this == &other)
          {
            return (*this);
          }
        else
          {
            resize(other.size());
            Manager::template copy<Context, Context>(m_array, other.m_array, m_size);
            return (*this);
          }
      }

      SimpleContainer & operator= (SimpleContainer && other)
      {
        if (this == &other)
          {
            return (*this);
          }
        else
          {
            clear();
            Manager::template move<Context, Context>(m_array, other.m_array, other.size());
            m_size = other.m_size;
            other.m_size = 0;
            return (*this);
          }
      }


      template <class other_indexer, class other_context, bool other_hold>
      SimpleContainer & operator= (const SimpleContainer<T, other_indexer, other_context, other_hold> & other)
      {
        resize(other.m_size);
        Manager::template copy<Context, other_context>(m_array, other.m_array, m_size);
        return (*this);
      }

      template <class other_indexer, class other_context>
      SimpleContainer & operator= (SimpleContainer<T, other_indexer, other_context, true> && other)
      {
        clear();
        Manager::template move<Context, other_context>(m_array, other.m_array, other.m_size);
        m_size = other.m_size;
        other.m_size = 0;
        return (*this);
      }

      ~SimpleContainer()
      {
        Manager::template deallocate<Context>(m_array);
        m_size = 0;
      }

      CUDA_HOS_DEV operator const T * () const
      {
        return m_array;
      }

      CUDA_HOS_DEV operator T * ()
      {
        return m_array;
      }

      template <class stream, class str = std::basic_string<typename stream::char_type> >
      void textual_output(stream & s, const str & separator = " ") const
      {
        if (std::is_same<Context, MemoryContext::CPU>::value)
          {
            s << m_size << separator;
            for (indexer i = 0; i < m_size - 1; ++i)
              {
                s << m_array[i] << separator;
              }
            s << m_array[m_size - 1];
          }
        else
          {
            SimpleContainer<T, indexer, MemoryContext::CPU, true> other(*this);
            other.textual_output(s, separator);
          }
      }

      template <class stream>
      void textual_input(stream & s)
      {
        if (std::is_same<Context, MemoryContext::CPU>::value)
          {
            indexer new_size;
            s >> new_size >> std::ws;
            if (s.fail())
              {
                //Throw errors, perhaps? Don't know if we can/should use exceptions...
                std::cerr << "FAILED READING " << this << "!" << std::endl;
                new_size = 0;
              }
            resize(new_size);
            for (indexer i = 0; i < m_size - 1; ++i)
              {
                s >> m_array[i];
                s >> std::ws;
              }
            s >> m_array[m_size - 1];
          }
        else
          {
            SimpleContainer<T, indexer, MemoryContext::CPU, true> other;
            other.textual_input(s);
            (*this) = other;
          }
      }

      template <class stream>
      void binary_output(stream & s) const
      {
        if (std::is_same<Context, MemoryContext::CPU>::value)
          {
            s.write((char *) &m_size, sizeof(indexer));
            for (indexer i = 0; i < m_size; ++i)
              {
                s.write((char *) (m_array + i), sizeof(T));
              }
          }
        else
          {
            SimpleContainer<T, indexer, MemoryContext::CPU, true> other(*this);
            other.binary_output(s);
          }
      }

      template <class stream>
      void binary_input(stream & s)
      {
        if (std::is_same<Context, MemoryContext::CPU>::value)
          {
            indexer new_size;
            s.read((char *) &new_size, sizeof(indexer));
            if (s.fail())
              {
                //Throw errors, perhaps? Don't know if we can/should use exceptions...
                std::cerr << "FAILED READING " << this << "!" << std::endl;
                new_size = 0;
              }
            resize(new_size);
            for (indexer i = 0; i < m_size; ++i)
              {
                s.read((char *) (m_array + i), sizeof(T));
              }
          }
        else
          {
            SimpleContainer<T, indexer, MemoryContext::CPU, true> other;
            other.binary_input(s);
            (*this) = other;
          }
      }

    };

    template <class T, class indexer, class Context>
    class SimpleContainer<T, indexer, Context, false>
    {
      static_assert(std::is_trivially_copyable<T>::value, "SimpleContainer only works with a trivially copyable type.");
      T * m_array;
      indexer m_size;

      using Manager = MemoryManagement<T, indexer>;

      template <class a, class b, class c, bool d> friend class SimpleContainer;

     public:

      CUDA_HOS_DEV inline indexer size() const
      {
        return m_size;
      }

      CUDA_HOS_DEV inline T & operator[] (const indexer i)
      {
        return m_array[i];
      }

      CUDA_HOS_DEV inline const T & operator[] (const indexer i) const
      {
        return m_array[i];
      }

      CUDA_HOS_DEV SimpleContainer() : m_array(nullptr), m_size(0)
      {
      }

      /*!
        \warning We assume the pointer is in a valid memory location!
      */
      CUDA_HOS_DEV SimpleContainer(T * other_array, const indexer sz)
      {
        m_array = other_array;
        m_size = sz;
      }

      template <class other_indexer, bool other_hold>
      CUDA_HOS_DEV SimpleContainer(const SimpleContainer<T, other_indexer, Context, other_hold> & other)
      {
        m_size = other.m_size;
        m_array = other.m_array;
      }

      CUDA_HOS_DEV SimpleContainer & operator= (const SimpleContainer & other)
      {
        if (this == &other)
          {
            return (*this);
          }
        else
          {
            m_array = other.m_array;
            m_size = other.m_size;
          }
      }

      template <class other_indexer, bool other_hold>
      CUDA_HOS_DEV SimpleContainer & operator= (const SimpleContainer<T, other_indexer, Context, other_hold> & other)
      {
        m_size = other.m_size;
        m_array = other.m_array;
        return (*this);
      }

      CUDA_HOS_DEV operator const T * () const
      {
        return m_array;
      }

      CUDA_HOS_DEV operator T * ()
      {
        return m_array;
      }
    };

    /// \brief Holds a run-time specified amount of objects of type \p T in CPU memory.
    template <class T, class indexer = unsigned int>
    using CPU_array = SimpleContainer<T, indexer, MemoryContext::CPU, true>;

    /// \brief Holds a run-time specified amount of objects of type \p T in CUDA GPU memory.
    template <class T, class indexer = unsigned int>
    using CUDA_array = SimpleContainer<T, indexer, MemoryContext::CUDAGPU, true>;

    /// \brief Non-owning pointer to an array of \p T in CUDA GPU memory.
    template <class T, class indexer = unsigned int>
    using CUDA_kernel_array = SimpleContainer<T, indexer, MemoryContext::CUDAGPU, false>;

    /*! Holds one objects of type \T in memory context \p Context.
        Automatically handles memory transfers to and from
        other contexts in the constructors and assignment operators.
        If \p hold_object is \p false, functions as a simple non-holding pointer,
        just keeping track of the \p Context for possible copies/transfers between memory contexts.
    */
    template <class T, class Context, bool hold_object = true>
    class SimpleHolder;

    template <class T, class Context>
    class SimpleHolder<T, Context, true>
    {
      static_assert(std::is_trivially_copyable<T>::value, "SimpleHolder only works with a trivially copyable type.");

      using indexer = unsigned int;

      T * m_object;

      using Manager = MemoryManagement<T, indexer>;

      template <class a, class b, bool c> friend class SimpleHolder;

     public:

      CUDA_HOS_DEV const T & operator *() const
      {
        return *m_object;
      }

      CUDA_HOS_DEV T & operator *()
      {
        return *m_object;
      }

      CUDA_HOS_DEV const T * operator ->() const
      {
        return m_object;
      }

      CUDA_HOS_DEV T * operator ->()
      {
        return m_object;
      }

      CUDA_HOS_DEV inline bool valid() const
      {
        return m_object != nullptr;
      }

      inline void clear()
      {
        Manager::template deallocate<Context>(m_object);
      }

      inline void allocate()
      {
        if (m_object == nullptr)
          {
            m_object = Manager::template allocate<Context>(1);
          }
      }

      SimpleHolder(): m_object(nullptr)
      {
      }

      SimpleHolder(const bool really_allocate)
      {
        if (really_allocate)
          {
            m_object = Manager::template allocate<Context>(1);
          }
        else
          {
            m_object = nullptr;
          }
      }

      /*!
        \warning We assume the pointer is in a valid memory location!
      */
      template < class X, class disabler = typename std::enable_if < std::is_base_of<T, X>::value || std::is_same<T, X>::value >::type >
      explicit SimpleHolder(X * other_p)
      {
        m_object = Manager::template allocate<Context>(1);
        Manager::template copy<Context, Context>(m_object, other_p, 1);
      }

      /*!
         \warning We assume the object is in a valid memory location!
       */
      template < class X, class disabler = typename std::enable_if < std::is_base_of<T, X>::value || std::is_same<T, X>::value >::type >
      SimpleHolder(const X & other_v) : SimpleHolder(&other_v)
      {


      }

      SimpleHolder(const SimpleHolder & other)
      {
        if (other.valid())
          {
            m_object = Manager::template allocate<Context>(1);
            Manager::template copy<Context, Context>(m_object, other.m_object, other.valid());
          }
        else
          {
            m_object = nullptr;
          }
      }

      template < class X, class other_context, bool other_hold,
                 class disabler = typename std::enable_if < std::is_base_of<T, X>::value || std::is_same<T, X>::value >::type >
      SimpleHolder(const SimpleHolder<X, other_context, other_hold> & other)
      {
        if (other.valid())
          {
            m_object = Manager::template allocate<Context>(1);
            Manager::template copy<Context, other_context>(m_object, other.m_object, other.valid());
          }
        else
          {
            m_object = nullptr;
          }
      }

      SimpleHolder(SimpleHolder && other)
      {
        m_object = nullptr;
        Manager::template move<Context, Context>(m_object, other.m_object, other.valid());
      }

      template < class X, class other_context,
                 class disabler = typename std::enable_if < std::is_base_of<T, X>::value || std::is_same<T, X>::value >::type >
      SimpleHolder(SimpleHolder<X, other_context, true> && other)
      {
        m_object = nullptr;
        Manager::template move<Context, other_context>(m_object, other.m_object, other.valid());
      }

      SimpleHolder & operator= (const SimpleHolder & other)
      {
        if (!valid() && other.valid())
          {
            allocate();
          }
        if (&other != this)
          {
            Manager::template copy<Context, Context>(m_object, other.m_object, other.valid());
          }
        return (*this);
      }

      template < class X, class other_context, bool other_hold,
                 class disabler = typename std::enable_if < std::is_base_of<T, X>::value || std::is_same<T, X>::value >::type >
      SimpleHolder & operator= (const SimpleHolder<X, other_context, other_hold> & other)
      {
        if (!valid() && other.valid())
          {
            allocate();
          }
        Manager::template copy<Context, other_context>(m_object, other.m_object, other.valid());
        return (*this);
      }

      SimpleHolder & operator= (SimpleHolder && other)
      {
        if (&other != this)
          {
            clear();
            Manager::template move<Context, Context>(m_object, other.m_object, other.valid());
          }
        return (*this);
      }

      template < class X, class other_context,
                 class disabler = typename std::enable_if < std::is_base_of<T, X>::value || std::is_same<T, X>::value >::type >
      SimpleHolder & operator= (SimpleHolder<X, other_context, true> && other)
      {
        clear();
        Manager::template move<Context, other_context>(m_object, other.m_object, other.valid());
        return (*this);
      }

      ~SimpleHolder()
      {
        Manager::template deallocate<Context>(m_object);
      }

      template < class X, class disabler = typename std::enable_if < std::is_base_of<X, T>::value || std::is_same<T, X>::value >::type >
      CUDA_HOS_DEV operator const X * () const
      {
        return m_object;
      }

      template < class X, class disabler = typename std::enable_if < std::is_base_of<X, T>::value || std::is_same<T, X>::value >::type >
      CUDA_HOS_DEV operator X * ()
      {
        return m_object;
      }

      template <class stream, class str = std::basic_string<typename stream::char_type> >
      void textual_output(stream & s, const str & separator = " ") const
      {
        if (std::is_same<Context, MemoryContext::CPU>::value)
          {
            if (m_object == nullptr)
              {
                s << 0;
              }
            else
              {
                s << 1 << separator << (*m_object);
              }
          }
        else
          {
            SimpleHolder<T, MemoryContext::CPU, true> other(*this);
            other.textual_output(s, separator);
          }
      }

      template <class stream>
      void textual_input(stream & s)
      {
        if (std::is_same<Context, MemoryContext::CPU>::value)
          {
            bool is_valid;
            s >> is_valid >> std::ws;
            if (s.fail())
              {
                //Throw errors, perhaps? Don't know if we can/should use exceptions...
                std::cerr << "FAILED READING " << this << "!" << std::endl;
                is_valid = false;
              }
            if (is_valid)
              {
                allocate();
                s >> (*m_object);
              }
            else
              {
                clear();
              }
          }
        else
          {
            SimpleHolder<T, MemoryContext::CPU, true> other;
            other.textual_input(s);
            (*this) = other;
          }
      }

      template <class stream>
      void binary_output(stream & s) const
      {
        if (m_object == nullptr)
          {
            return;
          }
        if (std::is_same<Context, MemoryContext::CPU>::value)
          {
            s.write((char *) m_object, sizeof(T));
          }
        else
          {
            SimpleHolder<T, MemoryContext::CPU, true> other(*this);
            other.binary_output(s);
          }
      }

      template <class stream>
      void binary_input(stream & s)
      {
        if (std::is_same<Context, MemoryContext::CPU>::value)
          {
            allocate();
            s.read((char *) m_object, sizeof(T));
          }
        else
          {
            SimpleHolder<T, MemoryContext::CPU, true> other;
            other.binary_input(s);
            (*this) = other;
          }
      }

    };

    template <class T, class Context>
    class SimpleHolder<T, Context, false>
    {
      static_assert(std::is_trivially_copyable<T>::value, "SimpleHolder only works with a trivially copyable type.");

      using indexer = unsigned int;

      T * m_object;

      using Manager = MemoryManagement<T, indexer>;

      template <class a, class b, bool c> friend class SimpleHolder;

     public:

      CUDA_HOS_DEV const T & operator *() const
      {
        return *m_object;
      }

      CUDA_HOS_DEV T & operator *()
      {
        return *m_object;
      }

      CUDA_HOS_DEV const T * operator ->() const
      {
        return m_object;
      }

      CUDA_HOS_DEV T * operator ->()
      {
        return m_object;
      }

      CUDA_HOS_DEV inline bool valid() const
      {
        return m_object != nullptr;
      }

      CUDA_HOS_DEV SimpleHolder() : m_object(nullptr)
      {
      }

      /*!
        \warning We assume the pointer is in a valid memory location!
      */
      template < class X, class disabler = typename std::enable_if < std::is_base_of<T, X>::value || std::is_same<T, X>::value >::type >
      CUDA_HOS_DEV SimpleHolder(X * other_p)
      {
        m_object = other_p;
      }

      template < class X, bool other_hold,
                 class disabler = typename std::enable_if < std::is_base_of<T, X>::value || std::is_same<T, X>::value >::type >
      CUDA_HOS_DEV SimpleHolder(const SimpleHolder<X, Context, other_hold> & other)
      {
        m_object = other.m_object;
      }

      template < class X, bool other_hold,
                 class disabler = typename std::enable_if < std::is_base_of<T, X>::value || std::is_same<T, X>::value >::type >
      CUDA_HOS_DEV SimpleHolder & operator= (const SimpleHolder<X, Context, other_hold> & other)
      {
        m_object = other.m_object;
        return (*this);
      }

      template < class X, class disabler = typename std::enable_if < std::is_base_of<X, T>::value || std::is_same<T, X>::value >::type >
      CUDA_HOS_DEV operator const X * () const
      {
        return m_object;
      }

      template < class X, class disabler = typename std::enable_if < std::is_base_of<X, T>::value || std::is_same<T, X>::value >::type >
      CUDA_HOS_DEV operator X * ()
      {
        return m_object;
      }
    };

    /// \brief Holds an object of type \p T in CPU memory.
    template <class T>
    using CPU_object = SimpleHolder<T, MemoryContext::CPU, true>;

    /// \brief Holds an object of type \p T in CUDA GPU memory.
    template <class T>
    using CUDA_object = SimpleHolder<T, MemoryContext::CUDAGPU, true>;

    /// \brief Non-owning pointer to an object of type \p T in CUDA GPU memory.
    template <class T>
    using CUDA_kernel_object = SimpleHolder<T, MemoryContext::CUDAGPU, false>;

    /*! \brief Manages objects of type \p T in a thread-safe way,
               ensuring that there's an object available for each separate thread
               while minimizing the number of allocations.

        Internally uses a \c std::vector of \c std::thread::id to manage the different threads,
        thus it is especially efficient for relatively small numbers of concurrent threads
        (which matches the use case, e. g. 64 threads )

    */
    template <class T>
    class separate_thread_holder
    {
     private:
      std::vector< std::unique_ptr<T> > m_held;
      std::vector< typename std::thread::id > m_thread_equivs;
      //For a sufficiently small number of threads
      //(not much more than 100 or so?)
      //it's faster to have linear search+insert
      //than any other addressing mode
      //(e. g. unordered_map)
      //We could still consider a more sophisticated solution...

      //Simple alternative: some sort of stack for non-assigned objects,
      //pushing and popping instead of linear searching.
      //(But with constant memory -> no (de)allocations.)

      std::mutex m_mutex;

     public:
      T & get_one()
      {
        std::lock_guard<std::mutex> lock_guard(m_mutex);
        std::thread::id this_id = std::this_thread::get_id();
        const std::thread::id invalid_id{};
        for (size_t i = 0; i < m_thread_equivs.size(); ++i)
          {
            if (m_thread_equivs[i] == invalid_id)
              {
                m_thread_equivs[i] = this_id;
                return *(m_held[i]);
              }
          }
        m_held.emplace_back(std::make_unique<T>());
        m_thread_equivs.emplace_back(this_id);
        return *(m_held.back());
      }

      void release_one()
      {
        std::lock_guard<std::mutex> lock_guard(m_mutex);
        std::thread::id this_id = std::this_thread::get_id();
        const std::thread::id invalid_id{};
        for (size_t i = 0; i < m_thread_equivs.size(); ++i)
          {
            if (m_thread_equivs[i] == this_id)
              {
                m_thread_equivs[i] = invalid_id;
              }
          }
      }

      void resize(const size_t new_size)
      {
        std::lock_guard<std::mutex> lock_guard(m_mutex);
        if (new_size < m_held.size())
          {
            m_held.resize(new_size);
            m_thread_equivs.resize(new_size);
          }
        else if (new_size > m_held.size())
          {
            const size_t to_add = new_size - m_held.size();
            const std::thread::id invalid_id{};
            for (size_t i = 0; i < to_add; ++i)
              {
                m_held.emplace_back(std::make_unique<T>());
                m_thread_equivs.emplace_back(invalid_id);
              }
          }
      }

      template <class F, class ... Args>
      void operate_on_all(F && f, Args && ... args)
      {
        std::lock_guard<std::mutex> lock_guard(m_mutex);
        for (std::unique_ptr<T> & obj : m_held)
          {
            f(*obj, std::forward<Args>(args)...);
          }
      }

      size_t held_size() const
      {
        return m_held.size();
      }

      size_t available_size() const
      {
        size_t count = 0;
        const std::thread::id invalid_id{};
        for (const auto & id : m_thread_equivs)
          {
            if (id == invalid_id)
              {
                ++count;
              }
          }
        return count;
      }

      size_t filled_size() const
      {
        return this->held_size() - this->available_size();
      }
    };

    /**! \brief Utility to access a \p separate_thread_holder to get a pointer and call \p release_one at scope exit.
    */
    template <class T>
    struct separate_thread_accessor
    {
     private:
      separate_thread_holder<T> & m_sth;
      T * m_held;
     public:
      separate_thread_accessor(separate_thread_holder<T> & s):
        m_sth(s), m_held(nullptr)
      {
      }
      T & get_one()
      {
        if (m_held == nullptr)
          {
            m_held = &(m_sth.get_one());
          }
        return *m_held;
      }
      void release_one()
      {
        if (m_held != nullptr)
          {
            m_sth.release_one();
            m_held = nullptr;
          }
      }
      ~separate_thread_accessor()
      {
        if (m_held != nullptr)
          {
            m_sth.release_one();
          }
      }
      separate_thread_accessor(separate_thread_holder<T> & s, T *& ptr):
        separate_thread_accessor(s)
      {
        get_one();
        ptr = m_held;
      }
    };
  }

}

#endif  // CALORECGPU_HELPERS_H