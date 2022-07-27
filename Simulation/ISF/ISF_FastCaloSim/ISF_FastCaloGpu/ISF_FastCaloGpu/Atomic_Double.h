/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ISF_FASTCALOGPU_ATOMIC_DOUBLE
#define ISF_FASTCALOGPU_ATOMIC_DOUBLE
#if !defined( __CUDA_ARCH__ ) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd( double* address, double val ) {
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int  old            = *address_as_ull, assumed;
  do {
    assumed = old;
    old     = atomicCAS( address_as_ull, assumed, __double_as_longlong( val + __longlong_as_double( assumed ) ) );
  } while ( assumed != old );
  return __longlong_as_double( old );
}
#endif
#endif
