#include "uuid/uuid.h"

void get_uuid(char *buffer){
  uuid_t uuid_buf;
  uuid_generate(uuid_buf);
  uuid_unparse(uuid_buf, buffer);
}
