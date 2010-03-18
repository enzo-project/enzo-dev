#ifdef DARWIN
// OS X has a uuid_t type defined in unistd.h,
// unless it's told otherwise
#define _POSIX_C_SOURCE 200112L
#endif

#include "kashmir/uuid.h"
#include "kashmir/devrand.h"

#include <iostream>
#include <fstream>

#include <unistd.h>
#include <sstream>

namespace{
  using kashmir::uuid_t;
  using kashmir::system::DevRand;
  using std::ostream;
  using std::ofstream;
}

void get_uuid(char *buffer){
  DevRand devrandom;
  kashmir::uuid_t uuid;
  devrandom >> uuid;
  std::stringstream s;
  s << uuid;
  sprintf(buffer, "%s", s.str().c_str());
}
