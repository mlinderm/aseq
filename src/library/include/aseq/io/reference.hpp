//
// Created by Michael Linderman on 2/22/16.
//

#pragma once

#include <memory>
#include "aseq/model/region.hpp"

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
  ReferenceSource(ReferenceSource&) = delete;
  ReferenceSource(const boost::filesystem::path&);
  ReferenceSource(const std::string&);
  virtual ~ReferenceSource() = default;

  virtual std::string Sequence(const model::Contig& ctg, model::Pos pos, model::Pos end);

 protected:
  ReferenceSource();

 private:
  std::unique_ptr<faidx_t, void (*)(faidx_t*)> faidx_;
};
}
}