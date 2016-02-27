//
// Created by Michael Linderman on 2/22/16.
//

#pragma once

#include "gmock/gmock.h"

namespace aseq {
namespace io {
namespace testing {

class MockReferenceSource : public ReferenceSource {
 public:
  MOCK_METHOD3(Sequence, std::string(const model::Contig&, int64_t, int64_t));
};
}
}
}