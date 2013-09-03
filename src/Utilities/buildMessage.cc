#ifndef buildMessage_HH
#define buildMessage_HH

#include <iostream>
#include <sstream>
#include <string>

namespace DBC {
std::string buildMessage(const char* fileName, 
                         const int line,
                         const char* message) {
  std::stringstream str;
  str << "File " << fileName << ", line " << line << "." << std::endl
      << message << '\0';
  return str.str();
}
}
#endif
