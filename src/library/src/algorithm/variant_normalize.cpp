//
// Created by Michael Linderman on 3/6/16.
//

#include <glog/logging.h>
#include <aseq/io/vcf.hpp>

#include "aseq/util/exception.hpp"
#include "aseq/model/variant_context.hpp"
#include "aseq/algorithm/variant.hpp"

namespace aseq {
namespace algorithm {

using model::Allele;
using model::VariantContext;

VariantContext UpdateREFAllele(io::ReferenceSource& ref, VariantContext&& cxt) {
  if (cxt.ref() == model::Allele::N) {
    // Update reference allele
    std::string seq = ref.Sequence(cxt.contig(), cxt.pos(), cxt.pos());
    return VariantContext(std::move(cxt), model::Allele(seq));
  } else
    return std::move(cxt);
}

namespace {

template <class T>
size_t MaxToTrimFromHead(const T& alleles) {
  size_t shortest = std::numeric_limits<size_t>::max();
  for (auto a : alleles) {
    if (a.IsSymbolic()) return 0;
    shortest = std::min(shortest, a.size());
  }
  return (shortest < 2) ? 0 : shortest - 1;
}

template <class T>
size_t HeadOfAllelesToTrim(const T& alleles) {
  size_t clip_to = MaxToTrimFromHead(alleles);
  auto ref_first = alleles[0].begin();
  for (size_t a = 1; clip_to > 0 && a < alleles.size(); a++) {
    auto r = std::mismatch(ref_first, ref_first + clip_to, alleles[a].begin());
    clip_to = std::min(clip_to, static_cast<size_t>(r.first - ref_first));
  }
  return clip_to;
}

template <class T>
size_t MaxToTrimFromTail(const T& alleles) {
  size_t shortest = std::numeric_limits<size_t>::max();
  bool any_symbolic = false;
  for (auto a : alleles) {
    if (a.IsSymbolic())
      any_symbolic = true;
    else
      shortest = std::min(shortest, a.size());
  }
  if (any_symbolic) {
    // Don't trim all the bases if there a symbolic allele
    return (shortest < 2) ? 0 : shortest - 1;
  } else
    return shortest;  // Can trim all the bases of shortest allele
}

template <class T>
size_t TailOfAllelesToTrim(const T& alleles) {
  size_t clip_to = MaxToTrimFromTail(alleles);
  auto ref_first = alleles[0].rbegin();
  for (size_t a = 1; clip_to > 0 && a < alleles.size(); a++) {
    if (alleles[a].IsSymbolic()) continue;
    auto r = std::mismatch(ref_first, ref_first + clip_to, alleles[a].rbegin());
    clip_to = std::min(clip_to, static_cast<size_t>(r.first - ref_first));
  }
  return clip_to;
}

}  // anonymous namespace

VariantContext LeftAlignAndTrimAlleles(io::ReferenceSource& ref, VariantContext&& cxt) {
  if (cxt.NumAltAlleles() == 0) return std::move(cxt);  // No need to trim if only a single allele

  VariantContext::Alleles alleles{cxt.ref()};
  alleles.insert(alleles.end(), std::begin(cxt.alts()), std::end(cxt.alts()));
  for (auto& a : alleles) {
    if (a.IsSymbolic()) {
      LOG_FIRST_N(WARNING, 1) << "Variants with symbolic alleles can't currently be left-aligned";
      return std::move(cxt);
    }
  }

  size_t tail_clip = TailOfAllelesToTrim(alleles), head_clip = HeadOfAllelesToTrim(alleles);
  if (tail_clip == 0 && head_clip == 0)
    return std::move(cxt);  // Common case with nothing to do, e.g. SNV or well-formed INS or DEL

  VariantContext::Alleles new_alleles(alleles);
  model::Pos new_pos = cxt.pos();
  while (tail_clip > 0) {
    bool any_empty = false;
    for (auto& a : new_alleles) {
      a = a.SubAllele(0, a.size() - tail_clip);
      any_empty |= a.size() == 0;
    }

    if (any_empty) {
      // Need to prepend bases to all alleles, query reference for prepend sequence
      model::Pos new_end = std::max(static_cast<model::Pos>(1), new_pos - 1);
      new_pos = std::max(static_cast<model::Pos>(1), new_pos - 10);
      std::string prepend = ref.Sequence(cxt.contig(), new_pos, new_end);
      for (auto& a : new_alleles) {
        a = prepend + a;
      }
    }

    tail_clip = TailOfAllelesToTrim(new_alleles);
  }

  // Update the clipping after tail trimming
  head_clip = HeadOfAllelesToTrim(new_alleles);
  if (head_clip > 0) {
    for (auto& a : new_alleles) {
      a = a.SubAllele(head_clip);
    }
  }

  model::Contig new_contig = cxt.contig();
  VariantContext new_cxt(std::move(cxt), new_contig, new_pos + head_clip, new_alleles.front(),
                         new_alleles.begin() + 1, new_alleles.end());
  return new_cxt;
};

VariantContext Normalize(io::ReferenceSource& ref, VariantContext&& cxt) {
  return LeftAlignAndTrimAlleles(ref, UpdateREFAllele(ref, std::move(cxt)));
}
}
}