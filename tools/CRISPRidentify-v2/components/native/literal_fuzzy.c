#include <limits.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PATTERN 63
#define MAX_ERRORS 6
#define MAX_PER_TYPE 3
#define SEED_BLOCK_COUNT (MAX_ERRORS + 1)
#define MIN_SEED_BLOCK_LENGTH 4
#define QUORUM_TARGET_BLOCK_LENGTH 3
#define MAX_QUORUM_BLOCKS (MAX_PATTERN / QUORUM_TARGET_BLOCK_LENGTH)
#define QUORUM_CODE_COUNT (64 + 256)

enum {
    FUZZY_SUB = 0,
    FUZZY_INS = 1,
    FUZZY_DEL = 2
};

typedef struct {
    int start;
    int end;
    int counts[3];
    int change_count;
    uint8_t change_type[MAX_ERRORS];
    int change_pos[MAX_ERRORS];
} InternalMatch;

typedef struct {
    int start;
    int end;
    int substitutions;
    int insertions;
    int deletions;
    int substitution_pos[MAX_PER_TYPE];
    int insertion_pos[MAX_PER_TYPE];
    int deletion_pos[MAX_PER_TYPE];
} LiteralMatch;

typedef struct {
    const unsigned char *pattern;
    int pattern_length;
    const unsigned char *text;
    int end;
    int search_anchor;
    int error_limit;
    int match_start;
    int counts[3];
    int change_count;
    uint8_t change_type[MAX_ERRORS];
    int change_pos[MAX_ERRORS];
} SearchState;

int literal_fuzzy_abi_version(void);
int literal_fuzzy_finditer(const unsigned char *pattern, int pattern_length,
                           const unsigned char *text, int text_length,
                           LiteralMatch *matches, int capacity);

static int total_errors(const SearchState *state) {
    return state->counts[FUZZY_SUB]
        + state->counts[FUZZY_INS]
        + state->counts[FUZZY_DEL];
}

static void push_change(SearchState *state, int type, int pos) {
    state->change_type[state->change_count] = (uint8_t)type;
    state->change_pos[state->change_count] = pos;
    state->change_count += 1;
    state->counts[type] += 1;
}

static void pop_change(SearchState *state, int type) {
    state->counts[type] -= 1;
    state->change_count -= 1;
}

/*
 * Literal-only equivalent of regex 2019.3.9's first successful fuzzy
 * backtracking path. Exact characters are consumed greedily. At a mismatch,
 * regex tries substitution, insertion, then deletion.
 */
static int anchored_dfs(SearchState *state, int pattern_pos, int text_pos,
                        InternalMatch *result) {
    int errors;

    if (pattern_pos == state->pattern_length) {
        result->start = state->match_start;
        result->end = text_pos;
        memcpy(result->counts, state->counts, sizeof(result->counts));
        result->change_count = state->change_count;
        memcpy(result->change_type, state->change_type,
               (size_t)state->change_count * sizeof(result->change_type[0]));
        memcpy(result->change_pos, state->change_pos,
               (size_t)state->change_count * sizeof(result->change_pos[0]));
        return 1;
    }

    if (text_pos < state->end
            && state->pattern[pattern_pos] == state->text[text_pos]) {
        return anchored_dfs(
            state, pattern_pos + 1, text_pos + 1, result
        );
    }

    errors = total_errors(state);
    if (errors >= state->error_limit || errors >= MAX_ERRORS)
        return 0;

    if (text_pos < state->end
            && state->counts[FUZZY_SUB] < MAX_PER_TYPE) {
        push_change(state, FUZZY_SUB, text_pos);
        if (anchored_dfs(
                state, pattern_pos + 1, text_pos + 1, result))
            return 1;
        pop_change(state, FUZZY_SUB);
    }

    /*
     * regex deliberately suppresses insertion exactly at search_anchor when
     * searching, because it can instead advance the candidate start.
     */
    if (text_pos < state->end
            && text_pos != state->search_anchor
            && state->counts[FUZZY_INS] < MAX_PER_TYPE) {
        push_change(state, FUZZY_INS, text_pos);
        if (anchored_dfs(state, pattern_pos, text_pos + 1, result))
            return 1;
        pop_change(state, FUZZY_INS);
    }

    if (state->counts[FUZZY_DEL] < MAX_PER_TYPE) {
        push_change(state, FUZZY_DEL, text_pos - 1);
        if (anchored_dfs(state, pattern_pos + 1, text_pos, result))
            return 1;
        pop_change(state, FUZZY_DEL);
    }

    return 0;
}

/* Exact Levenshtein prefilter for pattern versus  m-3 .. m+3 text bases. */
static int myers_possible(const uint64_t peq[256], uint64_t high_bit,
                          int pattern_length, const unsigned char *text,
                          int start, int end) {
    uint64_t pv = ~(uint64_t)0;
    uint64_t mv = 0;
    int minimum_length;
    int maximum_length;
    int available;
    int score;
    int index;

    if (pattern_length <= 0 || pattern_length > MAX_PATTERN)
        return 1;

    available = end - start;
    minimum_length = pattern_length - MAX_PER_TYPE;
    maximum_length = pattern_length + MAX_PER_TYPE;
    if (available < minimum_length)
        return 0;
    if (maximum_length > available)
        maximum_length = available;

    score = pattern_length;

    if (minimum_length <= 0 && pattern_length <= MAX_PER_TYPE)
        return 1;

    for (index = 0; index < maximum_length; index++) {
        uint64_t eq = peq[text[start + index]];
        uint64_t xv = eq | mv;
        uint64_t xh = ((((eq & pv) + pv) ^ pv) | eq);
        uint64_t ph = mv | ~(xh | pv);
        uint64_t mh = pv & xh;

        if (ph & high_bit)
            score += 1;
        if (mh & high_bit)
            score -= 1;

        ph = (ph << 1) | 1;
        mh <<= 1;
        pv = mh | ~(xv | ph);
        mv = ph & xv;

        if (index + 1 >= minimum_length && score <= MAX_ERRORS)
            return 1;
    }

    return 0;
}

static int seed_prefilter_supported(const unsigned char *pattern,
                                    int pattern_length) {
    int index;

    if (pattern_length < SEED_BLOCK_COUNT * MIN_SEED_BLOCK_LENGTH)
        return 0;

    for (index = 0; index < pattern_length; index++) {
        unsigned char base = pattern[index];
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T')
            return 0;
    }
    return 1;
}

static void build_kmp_prefix(const unsigned char *block, int block_length,
                             int prefix[MAX_PATTERN]) {
    int index;
    int matched = 0;

    prefix[0] = 0;
    for (index = 1; index < block_length; index++) {
        while (matched > 0 && block[index] != block[matched])
            matched = prefix[matched - 1];
        if (block[index] == block[matched])
            matched += 1;
        prefix[index] = matched;
    }
}

static int dna_block_code(const unsigned char *block, int block_length) {
    int code = 0;
    int index;

    for (index = 0; index < block_length; index++) {
        int base_code;
        if (block[index] == 'A')
            base_code = 0;
        else if (block[index] == 'C')
            base_code = 1;
        else if (block[index] == 'G')
            base_code = 2;
        else if (block[index] == 'T')
            base_code = 3;
        else
            return -1;
        code = (code << 2) | base_code;
    }

    if (block_length == 4)
        code += 64;
    return code;
}

static void mark_quorum_block_hit(unsigned char *votes,
                                  int *last_marked,
                                  int block_offset, int hit_start,
                                  int text_length) {
    int center = hit_start - block_offset;
    int first = center - MAX_PER_TYPE;
    int last = center + MAX_PER_TYPE;
    int candidate;

    if (first < 0)
        first = 0;
    if (last > text_length)
        last = text_length;
    if (first <= *last_marked)
        first = *last_marked + 1;
    if (first > last)
        return;

    for (candidate = first; candidate <= last; candidate++)
        votes[candidate] += 1;
    *last_marked = last;
}

static int candidate_meets_quorum(
        const unsigned char *pattern, int pattern_length,
        const unsigned char *text, int text_length, int candidate) {
    int block_count = pattern_length / QUORUM_TARGET_BLOCK_LENGTH;
    int base_block_length = pattern_length / block_count;
    int long_block_count = pattern_length % block_count;
    int quorum = block_count - MAX_ERRORS;
    int block_offset = 0;
    int votes = 0;
    int block_index;

    for (block_index = 0; block_index < block_count; block_index++) {
        int block_length = base_block_length
            + (block_index < long_block_count ? 1 : 0);
        int exact = 0;
        int displacement;

        for (displacement = -MAX_PER_TYPE;
                displacement <= MAX_PER_TYPE;
                displacement++) {
            int64_t hit_start = (int64_t)candidate
                + block_offset + displacement;
            if (hit_start >= 0
                    && hit_start + block_length <= text_length
                    && memcmp(
                        pattern + block_offset, text + hit_start,
                        (size_t)block_length) == 0) {
                exact = 1;
                break;
            }
        }

        votes += exact;
        block_offset += block_length;
        if (votes >= quorum)
            return 1;
        if (votes + block_count - block_index - 1 < quorum)
            return 0;
    }

    return 0;
}

/*
 * Intersect ``mask`` with a second, independently complete seed certificate.
 * Splitting a pattern into K disjoint blocks means that at least K-6 blocks
 * remain exact under any accepted alignment.  A candidate therefore needs a
 * vote from at least K-6 distinct blocks.  Each block is approximately three
 * bases long (and hence exactly three or four bases for the supported range).
 *
 * The compact code buckets turn all K exact-block searches into one pass over
 * the text.  Overlapping occurrences are retained, while merged +/-3 hit
 * intervals ensure that a single block votes at most once for any start.
 * Allocation failure conservatively leaves the original seven-block mask.
 */
static size_t intersect_quorum_seed_mask(
        unsigned char *mask, size_t mask_length, size_t marked_count,
        const unsigned char *pattern, int pattern_length,
        const unsigned char *text, int text_length) {
    unsigned char *votes;
    int block_heads[QUORUM_CODE_COUNT];
    int next_block[MAX_QUORUM_BLOCKS];
    int block_offsets[MAX_QUORUM_BLOCKS];
    int last_marked[MAX_QUORUM_BLOCKS];
    int block_count = pattern_length / QUORUM_TARGET_BLOCK_LENGTH;
    int base_block_length;
    int long_block_count;
    int block_offset = 0;
    int quorum;
    int index;
    int hit_start;
    size_t retained_count = 0;

    if (block_count <= MAX_ERRORS || block_count > MAX_QUORUM_BLOCKS)
        return marked_count;

    /* Sparse seven-block masks are cheaper to verify candidate by candidate. */
    if (marked_count < mask_length / 16) {
        size_t candidate;
        for (candidate = 0; candidate < mask_length; candidate++) {
            if (!mask[candidate])
                continue;
            if (candidate_meets_quorum(
                    pattern, pattern_length, text, text_length,
                    (int)candidate)) {
                retained_count += 1;
            }
            else {
                mask[candidate] = 0;
            }
        }
        return retained_count;
    }

    votes = (unsigned char *)calloc(mask_length, sizeof(*votes));
    if (!votes)
        return marked_count;

    for (index = 0; index < QUORUM_CODE_COUNT; index++)
        block_heads[index] = -1;

    base_block_length = pattern_length / block_count;
    long_block_count = pattern_length % block_count;
    for (index = 0; index < block_count; index++) {
        int block_length = base_block_length
            + (index < long_block_count ? 1 : 0);
        int code = dna_block_code(pattern + block_offset, block_length);

        block_offsets[index] = block_offset;
        last_marked[index] = -1;
        next_block[index] = block_heads[code];
        block_heads[code] = index;
        block_offset += block_length;
    }

    for (hit_start = 0; hit_start + 3 <= text_length; hit_start++) {
        int code4 = hit_start + 4 <= text_length
            ? dna_block_code(text + hit_start, 4) : -1;
        int code = code4 >= 0
            ? (code4 - 64) >> 2
            : dna_block_code(text + hit_start, 3);
        int block_index = code >= 0 ? block_heads[code] : -1;

        while (block_index >= 0) {
            mark_quorum_block_hit(
                votes, &last_marked[block_index],
                block_offsets[block_index], hit_start, text_length
            );
            block_index = next_block[block_index];
        }

        if (code4 >= 0) {
            block_index = block_heads[code4];
            while (block_index >= 0) {
                mark_quorum_block_hit(
                    votes, &last_marked[block_index],
                    block_offsets[block_index], hit_start, text_length
                );
                block_index = next_block[block_index];
            }
        }
    }

    quorum = block_count - MAX_ERRORS;
    for (index = 0; (size_t)index < mask_length; index++) {
        mask[index] = (unsigned char)(mask[index] && votes[index] >= quorum);
        retained_count += mask[index];
    }

    free(votes);
    return retained_count;
}

/*
 * Any match accepted by the matcher has at most six edits. Partitioning the
 * pattern into seven disjoint blocks therefore leaves at least one block
 * untouched and present exactly in the text. At most three insertions or
 * three deletions can precede that block, so its observed position identifies
 * the match start within +/-3 text bases.
 *
 * The seven-block union is intersected with the independently complete
 * three-base quorum above.  NULL means "do not filter": short/non-ACGT
 * patterns, allocation failure, or a final mask covering at least 75% of
 * starts all retain the original full scan.
 */
static unsigned char *build_seed_candidate_mask(
        const unsigned char *pattern, int pattern_length,
        const unsigned char *text, int text_length) {
    unsigned char *mask;
    size_t mask_length;
    size_t marked_count = 0;
    size_t dense_threshold;
    int base_block_length;
    int long_block_count;
    int block_offset = 0;
    int block_index;

    if (!seed_prefilter_supported(pattern, pattern_length))
        return NULL;

    mask_length = (size_t)text_length + 1;
    mask = (unsigned char *)calloc(mask_length, sizeof(*mask));
    if (!mask)
        return NULL;

    /* ceil(3 * mask_length / 4), without multiplication overflow. */
    dense_threshold = mask_length - mask_length / 4;
    base_block_length = pattern_length / SEED_BLOCK_COUNT;
    long_block_count = pattern_length % SEED_BLOCK_COUNT;

    for (block_index = 0;
            block_index < SEED_BLOCK_COUNT;
            block_index++) {
        const unsigned char *block = pattern + block_offset;
        int block_length = base_block_length
            + (block_index < long_block_count ? 1 : 0);
        int prefix[MAX_PATTERN];
        int matched = 0;
        int text_index;

        build_kmp_prefix(block, block_length, prefix);

        for (text_index = 0; text_index < text_length; text_index++) {
            while (matched > 0 && text[text_index] != block[matched])
                matched = prefix[matched - 1];
            if (text[text_index] == block[matched])
                matched += 1;

            if (matched == block_length) {
                int hit_start = text_index - block_length + 1;
                int center = hit_start - block_offset;
                int first = center - MAX_PER_TYPE;
                int last = center + MAX_PER_TYPE;
                int candidate;

                if (first < 0)
                    first = 0;
                if (last > text_length)
                    last = text_length;

                for (candidate = first; candidate <= last; candidate++) {
                    if (!mask[candidate]) {
                        mask[candidate] = 1;
                        marked_count += 1;
                        if (marked_count >= dense_threshold) {
                            free(mask);
                            return NULL;
                        }
                    }
                }

                /* Preserve overlapping exact block hits. */
                matched = prefix[matched - 1];
            }
        }

        block_offset += block_length;
    }

    if (marked_count > 0)
        marked_count = intersect_quorum_seed_mask(
            mask, mask_length, marked_count,
            pattern, pattern_length, text, text_length
        );

    if (marked_count >= dense_threshold) {
        free(mask);
        return NULL;
    }

    return mask;
}

static int first_search(const unsigned char *pattern, int pattern_length,
                        const uint64_t peq[256], uint64_t high_bit,
                        const unsigned char *candidate_mask,
                        const unsigned char *text, int pos, int end,
                        int error_limit, InternalMatch *result) {
    SearchState state;
    int start;

    memset(&state, 0, sizeof(state));
    state.pattern = pattern;
    state.pattern_length = pattern_length;
    state.text = text;
    state.end = end;
    /* Keep the incoming anchor; do not replace it with a selected start. */
    state.search_anchor = pos;
    state.error_limit = error_limit;

    for (start = pos;; start++) {
        if ((!candidate_mask || candidate_mask[start])
                && myers_possible(
                peq, high_bit, pattern_length, text, start, end)) {
            state.match_start = start;
            memset(state.counts, 0, sizeof(state.counts));
            state.change_count = 0;
            if (anchored_dfs(&state, 0, start, result))
                return 1;
        }
        if (start == end)
            break;
    }

    return 0;
}

static int enhanced_search(const unsigned char *pattern, int pattern_length,
                           const uint64_t peq[256], uint64_t high_bit,
                           const unsigned char *candidate_mask,
                           const unsigned char *text, int pos, int end,
                           InternalMatch *result) {
    InternalMatch candidate;
    InternalMatch best;
    int have_best = 0;
    int text_pos = pos;
    int slice_end = end;
    int error_limit = MAX_ERRORS;

    for (;;) {
        int errors;
        if (!first_search(
                pattern, pattern_length, peq, high_bit, candidate_mask,
                text, text_pos, slice_end, error_limit, &candidate))
            break;

        errors = candidate.counts[FUZZY_SUB]
            + candidate.counts[FUZZY_INS]
            + candidate.counts[FUZZY_DEL];
        if (have_best) {
            int best_errors = best.counts[FUZZY_SUB]
                + best.counts[FUZZY_INS]
                + best.counts[FUZZY_DEL];
            if (errors >= best_errors)
                break;
        }

        best = candidate;
        have_best = 1;
        if (errors == 0)
            break;

        text_pos = candidate.start;
        slice_end = candidate.end;
        error_limit = errors - 1;
    }

    if (have_best)
        *result = best;
    return have_best;
}

static void export_match(const InternalMatch *source, LiteralMatch *target) {
    int sub_index = 0;
    int ins_index = 0;
    int del_index = 0;
    int deletion_offset = 0;
    int index;

    memset(target, 0, sizeof(*target));
    target->start = source->start;
    target->end = source->end;
    target->substitutions = source->counts[FUZZY_SUB];
    target->insertions = source->counts[FUZZY_INS];
    target->deletions = source->counts[FUZZY_DEL];

    for (index = 0; index < source->change_count; index++) {
        int type = source->change_type[index];
        int pos = source->change_pos[index];
        if (type == FUZZY_SUB)
            target->substitution_pos[sub_index++] = pos;
        else if (type == FUZZY_INS)
            target->insertion_pos[ins_index++] = pos;
        else {
            target->deletion_pos[del_index++] = pos + deletion_offset;
            deletion_offset += 1;
        }
    }
}

int literal_fuzzy_abi_version(void) {
    return 1;
}

/*
 * Returns the number of matches written, or -(required count) when capacity
 * is insufficient. Patterns must be literal bytes of length <=63.
 */
int literal_fuzzy_finditer(const unsigned char *pattern, int pattern_length,
                           const unsigned char *text, int text_length,
                           LiteralMatch *matches, int capacity) {
    uint64_t peq[256] = {0};
    uint64_t high_bit;
    unsigned char *candidate_mask;
    int count = 0;
    int pos = 0;
    int index;

    if (!pattern || !text || pattern_length <= 0
            || pattern_length > MAX_PATTERN || text_length < 0
            || text_length > INT_MAX - MAX_PER_TYPE || capacity < 0
            || (capacity > 0 && !matches))
        return 0;

    for (index = 0; index < pattern_length; index++)
        peq[pattern[index]] |= (uint64_t)1 << index;
    high_bit = (uint64_t)1 << (pattern_length - 1);
    candidate_mask = build_seed_candidate_mask(
        pattern, pattern_length, text, text_length
    );

    for (;;) {
        InternalMatch match;
        if (!enhanced_search(
                pattern, pattern_length, peq, high_bit, candidate_mask,
                text, pos, text_length, &match))
            break;
        if (count < capacity)
            export_match(&match, &matches[count]);
        count += 1;
        if (match.end > pos)
            pos = match.end;
        else if (pos < text_length)
            pos += 1;
        else
            break;
    }

    free(candidate_mask);
    if (count > capacity)
        return -count;
    return count;
}
