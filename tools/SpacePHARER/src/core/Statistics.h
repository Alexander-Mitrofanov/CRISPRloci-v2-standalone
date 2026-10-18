// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>
namespace spacepharer {
struct HitEvidence {
    double log_p;
    double identity;
};
struct CombinedScore {
    double score;
    double mean_identity;
};
// Preserve the published engine's order-sensitive, below-threshold fallback.
// Changing this would be a separate scientific algorithm revision.
inline CombinedScore combine_score(const std::vector<HitEvidence> &hits, unsigned int orf_count) {
    const double threshold = std::log(1.0 / (orf_count + 1));
    double minimum = 0, minimum_identity = 0, sum = 0, identities = 0;
    std::size_t count = 0;
    for (const auto &hit : hits) {
        if (hit.log_p < minimum) {
            minimum = hit.log_p == 0 ? hit.log_p : -hit.log_p;
            minimum_identity = hit.identity;
        }
        if (hit.log_p < threshold) {
            sum -= hit.log_p - threshold;
            identities += hit.identity;
            ++count;
        }
    }
    return count == 0 ? CombinedScore{minimum, minimum_identity}
                      : CombinedScore{sum - threshold, identities / count};
}
inline double empirical_tail(const std::vector<double> &sorted_null, double score) {
    auto lower = std::lower_bound(sorted_null.begin(), sorted_null.end(), score);
    return (std::distance(lower, sorted_null.end()) + 0.5) / (sorted_null.size() + 1.0);
}
inline double combine_evalues(double protein, double nucleotide) {
    const double average = std::log(protein) * 0.5 + std::log(nucleotide) * 0.5;
    return average < std::log(nucleotide) ? std::exp(average) : nucleotide;
}
} // namespace spacepharer
