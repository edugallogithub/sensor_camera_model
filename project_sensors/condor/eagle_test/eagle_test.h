#ifndef EAGLE_TEST_INCLUDED
#define EAGLE_TEST_INCLUDED

// 4251: some class or template needs to have dll-interface to be used by clients
// of this class. No problem when reffers to templates, clients can link succesfully
// because the class instantiation (the code) is in client side (so no link at all).
// Disabling 4251 can hide a problem when this class (which is dll exported): either
// (1) expose in non-private sections or (2) use in inline functions implementations,
// types witch are not dll exported
// 4996:  MSVC has gone into full security-paranoia mode. std::copy issues this warning
// when it is used with raw pointers, because when used incorrectly, it can result in
// buffer overflows. But if code is OK, no problem. Stupid Microsoft warning.
#if defined(_MSC_VER)
#	pragma warning(disable:4251)
#	pragma warning(disable:4996)
#endif

// Now we use the generic helper definitions above to define EAGLE_TEST_API and EAGLE_TEST_LOCAL.
// EAGLE_TEST_API is used for the public API symbols. It either DLL imports or DLL exports (or does nothing for static build)
// EAGLE_TEST_LOCAL is used for non-api symbols.

#ifdef EAGLE_TEST_DLL // defined if EAGLE_TEST is compiled as a DLL
    #ifdef EAGLE_TEST_DLL_EXPORTS // defined if we are building the EAGLE_TEST DLL (instead of using it)
        #define EAGLE_TEST_API EAGLE_HELPER_DLL_EXPORT
    #else
        #define EAGLE_TEST_API EAGLE_HELPER_DLL_IMPORT
    #endif // EAGLE_TEST_DLL_EXPORTS
    #define EAGLE_TEST_LOCAL EAGLE_HELPER_DLL_LOCAL
#else // EAGLE_TEST_DLL is not defined: this means EAGLE_TEST is a static lib.
    #define EAGLE_TEST_API
    #define EAGLE_TEST_LOCAL
#endif // EAGLE_TEST_DLL


#endif
