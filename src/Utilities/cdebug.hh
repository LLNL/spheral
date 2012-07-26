//---------------------------------------------------------------------------
//
//                cdebug.hh -- A debugging cout/cerr.
//
// This ostream only sends messages to cerr if DEBUG is defined.
// Yeah, I know... this is a hack.  But it shouldn't impact performance when 
// you use a modern compiler, and it beats bracketing all your debugging 
// statements with #ifdef DEBUG ... #endif.
//
//---------------------------------------------------------------------------

#ifndef cdebug

#include <iostream>

#ifdef DEBUG
#define cdebug std::cerr 
#else
#define cdebug if (false) std::cerr
#endif

#endif
