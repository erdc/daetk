#ifndef PLATFORM_H
#define PLATFORM_H

#include "Definitions.h"

namespace Daetk 
{
#ifdef NO_BOOL_TYPE
typedef int bool;
const bool true(1);
const bool false(0);
#endif
}
#endif
