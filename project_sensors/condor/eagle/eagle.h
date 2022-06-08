#pragma once

// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define EAGLE_HELPER_DLL_IMPORT __declspec(dllimport)
  #define EAGLE_HELPER_DLL_EXPORT __declspec(dllexport)
  #define EAGLE_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define EAGLE_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define EAGLE_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define EAGLE_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define EAGLE_HELPER_DLL_IMPORT
    #define EAGLE_HELPER_DLL_EXPORT
    #define EAGLE_HELPER_DLL_LOCAL
  #endif
#endif

// Now we use the generic helper definitions above to define EAGLE_API and EAGLE_LOCAL.
// EAGLE_API is used for the public API symbols. It either DLL imports or DLL exports (or does nothing for static build)
// EAGLE_LOCAL is used for non-api symbols.

#ifdef EAGLE_DLL // defined if EAGLE is compiled as a DLL
  #ifdef EAGLE_DLL_EXPORTS // defined if we are building the EAGLE DLL (instead of using it)
    #define EAGLE_API EAGLE_HELPER_DLL_EXPORT
  #else
    #define EAGLE_API EAGLE_HELPER_DLL_IMPORT
  #endif // EAGLE_DLL_EXPORTS
  #define EAGLE_LOCAL EAGLE_HELPER_DLL_LOCAL
#else // EAGLE_DLL is not defined: this means EAGLE is a static lib.
  #define EAGLE_API
  #define EAGLE_LOCAL
#endif // EAGLE_DLL

