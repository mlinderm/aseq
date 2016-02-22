//
// Created by Michael Linderman on 2/22/16.
//

#pragma once

#include <memory>

// htslib
typedef struct __faidx_t faidx_t;

namespace boost {
namespace filesystem {
class path;
}
}

namespace aseq {
namespace io {

class ReferenceSource {
 public:
  ReferenceSource(const boost::filesystem::path&);

 protected:
  ReferenceSource();

 private:
  std::unique_ptr<faidx_t, void (*)(faidx_t*)> faidx_;
};
}
}