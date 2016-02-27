//
// Created by Michael Linderman on 2/22/16.
//

#pragma once

#include <memory>
#include <aseq/model/contig.hpp>

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

  // TODO: Create shared region types
  virtual std::string Sequence(const model::Contig& ctg, int64_t pos, int64_t end);

 protected:
  ReferenceSource();

 private:
  std::unique_ptr<faidx_t, void (*)(faidx_t*)> faidx_;
};
}
}