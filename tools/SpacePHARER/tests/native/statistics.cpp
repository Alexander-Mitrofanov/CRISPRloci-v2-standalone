#include "core/Statistics.h"
#include <cassert>
#include <cmath>
int main() {
    using namespace spacepharer;
    auto result = combine_score({{-10, .9}, {-8, .8}}, 4);
    assert(std::abs(result.score - (18 + std::log(.2))) < 1e-12);
    assert(std::abs(result.mean_identity - .85) < 1e-12);
    // Compatibility regression: legacy fallback is order-sensitive by design.
    assert(combine_score({{-.2, .7}, {-.1, .8}}, 4).score == .1);
    assert(combine_score({{-.1, .8}, {-.2, .7}}, 4).score == .2);
    assert(empirical_tail({}, 1) == .5);
    assert(empirical_tail({1, 2, 2, 3}, 2) == .7);
    assert(std::abs(combine_evalues(.01, .0001) - .0001) < 1e-15);
}
