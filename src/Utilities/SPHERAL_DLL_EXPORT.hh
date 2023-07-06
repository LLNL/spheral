#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define SPHERAL_DLL_PUBLIC __attribute__ ((dllexport))
    #else
      #define SPHERAL_DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #else
    #ifdef __GNUC__
      #define SPHERAL_DLL_PUBLIC __attribute__ ((dllimport))
    #else
      #define SPHERAL_DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #endif
  #define SPHERAL_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define SPHERAL_DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define SPHERAL_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define SPHERAL_DLL_PUBLIC
    #define SPHERAL_DLL_LOCAL
  #endif
#endif
