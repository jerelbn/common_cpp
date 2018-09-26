// Name: Jerel Nielsen
// Date: 15 June 2017
// Desc: Container for common functions.

#include "common_cpp/common.h"

namespace common
{


// create string of red text
std::string redText(const std::string& input)
{
  return "\033[31m" + input + "\033[0m";
}


// create string of green text
std::string greenText(const std::string& input)
{
  return "\033[32m" + input + "\033[0m";
}


} // namespace common


