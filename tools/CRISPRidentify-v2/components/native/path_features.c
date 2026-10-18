#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#define BITMASK_16 UINT32_C(65535)
#define BITMASK_32 UINT64_C(4294967295)
#define PY_HASH_MODULUS ((UINT64_C(1) << 61) - 1)
#define TUPLE_INITIAL_16 (UINT32_C(0x345678) & BITMASK_16)
#define TUPLE_MULTIPLIER_16 (UINT32_C(1000003) & BITMASK_16)
#define TUPLE_INCREMENT UINT32_C(82520)
#define TUPLE_FINAL_16 (UINT32_C(97531) & BITMASK_16)

static uint32_t feature_hash3(uint32_t first, uint32_t second,
                              uint32_t third) {
    uint32_t result;
    uint32_t multiplier;

    result = ((TUPLE_INITIAL_16 ^ (first & BITMASK_16))
              * TUPLE_MULTIPLIER_16) & BITMASK_16;
    multiplier = (TUPLE_MULTIPLIER_16 + TUPLE_INCREMENT + 4)
                 & BITMASK_16;
    result = ((result ^ (second & BITMASK_16)) * multiplier) & BITMASK_16;
    multiplier = (multiplier + TUPLE_INCREMENT + 2) & BITMASK_16;
    result = ((result ^ (third & BITMASK_16)) * multiplier) & BITMASK_16;
    return ((result + TUPLE_FINAL_16) & BITMASK_16) + 1;
}

static uint32_t feature_hash4(uint32_t first, uint32_t second,
                              uint32_t third, uint32_t fourth) {
    uint32_t result;
    uint32_t multiplier;

    result = ((TUPLE_INITIAL_16 ^ (first & BITMASK_16))
              * TUPLE_MULTIPLIER_16) & BITMASK_16;
    multiplier = (TUPLE_MULTIPLIER_16 + TUPLE_INCREMENT + 6)
                 & BITMASK_16;
    result = ((result ^ (second & BITMASK_16)) * multiplier) & BITMASK_16;
    multiplier = (multiplier + TUPLE_INCREMENT + 4) & BITMASK_16;
    result = ((result ^ (third & BITMASK_16)) * multiplier) & BITMASK_16;
    multiplier = (multiplier + TUPLE_INCREMENT + 2) & BITMASK_16;
    result = ((result ^ (fourth & BITMASK_16)) * multiplier) & BITMASK_16;
    return ((result + TUPLE_FINAL_16) & BITMASK_16) + 1;
}

int path_features_abi_version(void) {
    return 4;
}

/*
 * Emit packed (group, full feature, half feature) events in the exact order
 * consumed by sequence_path_vectorizer._direct_feature_dict. Neighborhoods
 * contain four hashes per sequence vertex, with zero marking absent radii.
 * Frontiers contain five packed pairs per vertex; each 16-bit coordinate is
 * stored as coordinate+1 so zero can mark an absent member.
 */
int path_feature_events(const uint32_t *neighborhoods,
                        const uint32_t *frontiers,
                        int length,
                        uint64_t *events,
                        int capacity) {
    int required = 0;
    int vertex_v;

    if (!neighborhoods || !frontiers || length <= 0 || capacity < 0
            || (capacity > 0 && !events))
        return 0;

    for (vertex_v = 0; vertex_v < length; vertex_v++) {
        const uint32_t *labels_v = neighborhoods + vertex_v * 4;
        int labels_v_length = 0;
        int distance;
        while (labels_v_length < 4 && labels_v[labels_v_length] != 0)
            labels_v_length += 1;

        for (distance = 0; distance < 5; distance++) {
            uint32_t packed_frontier = frontiers[vertex_v * 5 + distance];
            uint32_t encoded_vertices[2];
            int frontier_count = 0;
            int frontier_index;
            encoded_vertices[0] = packed_frontier & BITMASK_16;
            encoded_vertices[1] = packed_frontier >> 16;
            if (encoded_vertices[0] != 0)
                frontier_count += 1;
            if (encoded_vertices[1] != 0)
                frontier_count += 1;

            for (frontier_index = 0;
                    frontier_index < frontier_count;
                    frontier_index++) {
                int vertex_u = (int)encoded_vertices[frontier_index] - 1;
                const uint32_t *labels_u = neighborhoods + vertex_u * 4;
                int labels_u_length = 0;
                int limit;
                int radius;
                while (labels_u_length < 4 && labels_u[labels_u_length] != 0)
                    labels_u_length += 1;
                limit = labels_v_length < labels_u_length
                        ? labels_v_length : labels_u_length;

                for (radius = 0; radius < limit; radius++) {
                    uint32_t vertex_v_hash = labels_v[radius];
                    uint32_t vertex_u_hash = labels_u[radius];
                    uint32_t first_hash = vertex_v_hash;
                    uint32_t second_hash = vertex_u_hash;
                    uint32_t expanded_radius = (uint32_t)(radius * 2);
                    uint32_t expanded_distance = (uint32_t)(distance * 2);
                    uint32_t feature;
                    uint32_t half_feature;
                    uint32_t group = (uint32_t)(distance * 4 + radius);
                    if (first_hash > second_hash) {
                        uint32_t temporary = first_hash;
                        first_hash = second_hash;
                        second_hash = temporary;
                    }
                    feature = feature_hash4(
                        first_hash, second_hash,
                        expanded_radius, expanded_distance
                    );
                    half_feature = feature_hash3(
                        vertex_u_hash, expanded_radius, expanded_distance
                    );
                    if (required < capacity) {
                        events[required] = ((uint64_t)group << 34)
                                           | ((uint64_t)feature << 17)
                                           | (uint64_t)half_feature;
                    }
                    required += 1;
                }
            }
        }
    }

    return required <= capacity ? required : -required;
}

typedef struct {
    uint32_t key;
    uint32_t record_plus_one;
} FeatureCountSlot;

static int add_feature_count(uint64_t *counts,
                             int capacity,
                             int *count,
                             FeatureCountSlot *slots,
                             uint32_t slot_mask,
                             uint32_t group,
                             uint32_t feature) {
    uint32_t key = (group << 17) | feature;
    uint32_t slot_index = (key * UINT32_C(2654435761)) & slot_mask;

    for (;;) {
        FeatureCountSlot *slot = &slots[slot_index];
        if (slot->record_plus_one == 0) {
            int record_index;
            if (*count >= capacity)
                return 0;
            record_index = *count;
            *count += 1;
            counts[record_index] = ((uint64_t)group << 49)
                                   | ((uint64_t)feature << 32)
                                   | UINT64_C(1);
            slot->key = key;
            slot->record_plus_one = (uint32_t)record_index + 1;
            return 1;
        }
        if (slot->key == key) {
            counts[slot->record_plus_one - 1] += UINT64_C(1);
            return 1;
        }
        slot_index = (slot_index + 1) & slot_mask;
    }
}

/*
 * Aggregate the exact event stream into first-encounter-ordered
 * (group, feature, count) records. The order matches Python dict insertion
 * order, while integer counts remain exact until Python performs the frozen
 * floating-point normalization.
 */
int path_feature_counts(const uint32_t *neighborhoods,
                        const uint32_t *frontiers,
                        int length,
                        uint64_t *counts,
                        int capacity) {
    uint64_t *events;
    FeatureCountSlot *slots;
    uint32_t slot_capacity = 1;
    uint32_t slot_mask;
    int event_capacity;
    int event_count;
    int count = 0;
    int index;

    if (!neighborhoods || !frontiers || length <= 0 || capacity <= 0
            || !counts)
        return 0;
    event_capacity = length * 40;
    events = malloc((size_t)event_capacity * sizeof(*events));
    if (!events)
        return -1;
    event_count = path_feature_events(
        neighborhoods, frontiers, length, events, event_capacity
    );
    if (event_count < 0) {
        free(events);
        return event_count;
    }

    while (slot_capacity < (uint32_t)capacity * 2)
        slot_capacity <<= 1;
    slot_mask = slot_capacity - 1;
    slots = calloc((size_t)slot_capacity, sizeof(*slots));
    if (!slots) {
        free(events);
        return -1;
    }

    for (index = 0; index < event_count; index++) {
        uint64_t packed = events[index];
        uint32_t group = (uint32_t)((packed >> 34) & UINT64_C(31));
        uint32_t feature = (uint32_t)((packed >> 17) & UINT64_C(0x1ffff));
        uint32_t half_feature = (uint32_t)(packed & UINT64_C(0x1ffff));
        if (!add_feature_count(
                counts, capacity, &count, slots, slot_mask,
                group, feature)
                || !add_feature_count(
                    counts, capacity, &count, slots, slot_mask,
                    group, half_feature)) {
            free(slots);
            free(events);
            return -capacity;
        }
    }

    free(slots);
    free(events);
    return count;
}


typedef struct {
    uint32_t feature;
    double value;
} PackedFeatureValue;

typedef struct {
    uint32_t key;
    uint32_t record_plus_one;
} FinalFeatureSlot;

static int compare_packed_feature_value(const void *first,
                                        const void *second) {
    const PackedFeatureValue *left = first;
    const PackedFeatureValue *right = second;
    if (left->feature < right->feature)
        return -1;
    if (left->feature > right->feature)
        return 1;
    return 0;
}

/*
 * Return the final normalized sparse row rather than the intermediate
 * per-(distance,radius) counts. The loops reproduce Python insertion and
 * floating-point order exactly. Build without fast-math or contraction.
 */
int path_feature_row(const uint32_t *neighborhoods,
                     const uint32_t *frontiers,
                     int length,
                     uint32_t *feature_ids,
                     double *feature_values,
                     int capacity) {
    uint64_t *counts = NULL;
    uint32_t *group_order = NULL;
    uint8_t group_seen[20] = {0};
    PackedFeatureValue *final_values = NULL;
    FinalFeatureSlot *slots = NULL;
    uint32_t slot_capacity = 1;
    uint32_t slot_mask;
    int count_capacity;
    int count;
    int group_count = 0;
    int final_count = 0;
    int index;
    int group_index;
    double total_squared_norm = 0.0;
    double total_norm;

    if (!neighborhoods || !frontiers || length <= 0 || capacity < 0
            || (capacity > 0 && (!feature_ids || !feature_values)))
        return 0;

    count_capacity = length * 80;
    if (count_capacity < 1)
        count_capacity = 1;
    counts = malloc((size_t)count_capacity * sizeof(*counts));
    group_order = malloc(20 * sizeof(*group_order));
    final_values = malloc((size_t)count_capacity * sizeof(*final_values));
    if (!counts || !group_order || !final_values)
        goto allocation_failure;

    count = path_feature_counts(
        neighborhoods, frontiers, length, counts, count_capacity
    );
    if (count < 0)
        goto allocation_failure;

    while (slot_capacity < (uint32_t)count_capacity * 2)
        slot_capacity <<= 1;
    slot_mask = slot_capacity - 1;
    slots = calloc((size_t)slot_capacity, sizeof(*slots));
    if (!slots)
        goto allocation_failure;

    for (index = 0; index < count; index++) {
        uint32_t group = (uint32_t)((counts[index] >> 49) & UINT64_C(31));
        if (group >= 20)
            goto invalid_result;
        if (!group_seen[group]) {
            group_seen[group] = 1;
            group_order[group_count++] = group;
        }
    }

    for (group_index = 0; group_index < group_count; group_index++) {
        uint32_t group = group_order[group_index];
        double squared_norm = 0.0;
        double norm;

        for (index = 0; index < count; index++) {
            uint32_t record_group =
                (uint32_t)((counts[index] >> 49) & UINT64_C(31));
            if (record_group == group) {
                double value = (double)(counts[index] & BITMASK_32);
                squared_norm += value * value;
            }
        }
        norm = sqrt(squared_norm);
        if (!(norm > 0.0))
            goto invalid_result;

        for (index = 0; index < count; index++) {
            uint32_t record_group =
                (uint32_t)((counts[index] >> 49) & UINT64_C(31));
            uint32_t feature;
            uint32_t slot_index;
            double value;
            if (record_group != group)
                continue;
            feature = (uint32_t)((counts[index] >> 32) & UINT64_C(0x1ffff));
            value = (double)(counts[index] & BITMASK_32) / norm;
            slot_index = (feature * UINT32_C(2654435761)) & slot_mask;
            for (;;) {
                FinalFeatureSlot *slot = &slots[slot_index];
                if (slot->record_plus_one == 0) {
                    slot->key = feature;
                    slot->record_plus_one = (uint32_t)final_count + 1;
                    final_values[final_count].feature = feature;
                    final_values[final_count].value = value;
                    final_count += 1;
                    break;
                }
                if (slot->key == feature) {
                    final_values[slot->record_plus_one - 1].value = value;
                    break;
                }
                slot_index = (slot_index + 1) & slot_mask;
            }
        }
    }

    for (index = 0; index < final_count; index++) {
        double value = final_values[index].value;
        total_squared_norm += value * value;
    }
    total_norm = sqrt(total_squared_norm);
    if (!(total_norm > 0.0))
        goto invalid_result;
    for (index = 0; index < final_count; index++)
        final_values[index].value /= total_norm;

    if (final_count > capacity) {
        free(slots);
        free(final_values);
        free(group_order);
        free(counts);
        return -final_count;
    }
    qsort(
        final_values, (size_t)final_count, sizeof(*final_values),
        compare_packed_feature_value
    );
    for (index = 0; index < final_count; index++) {
        feature_ids[index] = final_values[index].feature;
        feature_values[index] = final_values[index].value;
    }

    free(slots);
    free(final_values);
    free(group_order);
    free(counts);
    return final_count;

invalid_result:
allocation_failure:
    free(slots);
    free(final_values);
    free(group_order);
    free(counts);
    return -1;
}

/* Accumulate in the exact sorted-feature order used by _sparse_dot. */
double path_feature_dot(const uint32_t *first_ids,
                        const double *first_values,
                        int first_count,
                        const uint32_t *second_ids,
                        const double *second_values,
                        int second_count) {
    int first_index = 0;
    int second_index = 0;
    double result = 0.0;

    if (!first_ids || !first_values || first_count < 0
            || !second_ids || !second_values || second_count < 0)
        return 0.0;
    while (first_index < first_count && second_index < second_count) {
        uint32_t first_feature = first_ids[first_index];
        uint32_t second_feature = second_ids[second_index];
        if (first_feature < second_feature) {
            first_index += 1;
        } else if (second_feature < first_feature) {
            second_index += 1;
        } else {
            result += first_values[first_index] * second_values[second_index];
            first_index += 1;
            second_index += 1;
        }
    }
    return result;
}

static int64_t unsigned_as_signed_64(uint64_t value) {
    uint64_t magnitude;

    if (value <= (uint64_t)INT64_MAX)
        return (int64_t)value;
    magnitude = (~value) + UINT64_C(1);
    if (magnitude == (UINT64_C(1) << 63))
        return INT64_MIN;
    return -(int64_t)magnitude;
}

static int64_t legacy_hash_integer(int64_t value) {
    uint64_t remainder;
    int64_t result;

    if (value >= 0) {
        remainder = (uint64_t)value % PY_HASH_MODULUS;
        result = (int64_t)remainder;
    } else {
        uint64_t magnitude = (~(uint64_t)value) + UINT64_C(1);
        remainder = magnitude % PY_HASH_MODULUS;
        result = -(int64_t)remainder;
    }
    return result == -1 ? -2 : result;
}

static int64_t legacy_hash_int_tuple(const int64_t *values, int count) {
    uint64_t result = UINT64_C(0x345678);
    uint64_t multiplier = UINT64_C(1000003);
    int remaining = count;
    int index;

    for (index = 0; index < count; index++) {
        int64_t item_hash;
        remaining -= 1;
        item_hash = legacy_hash_integer(values[index]);
        result = (result ^ (uint64_t)item_hash) * multiplier;
        multiplier += UINT64_C(82520) + (uint64_t)(remaining * 2);
    }
    result += UINT64_C(97531);
    if (result == UINT64_MAX)
        return -2;
    return unsigned_as_signed_64(result);
}

static int hash_low32_plus_one(int64_t hash_value, uint32_t *output) {
    uint32_t low_bits = (uint32_t)(uint64_t)hash_value;

    /* 2**32 cannot be represented in the ABI-v3 uint32 neighborhood row. */
    if (low_bits == UINT32_MAX)
        return 0;
    *output = low_bits + UINT32_C(1);
    return 1;
}

/*
 * Build Python-3.7-compatible path neighborhoods directly from ASCII bytes,
 * then reuse the exact ABI-v3 packed-row implementation. ``label_hashes`` is
 * the frozen 128-entry one-byte ASCII table supplied by legacy_hash.py.
 * Return -2 when an exact intermediate cannot fit the uint32 packed contract;
 * callers must use the Python fallback in that case.
 */
int path_feature_ascii_row(const uint8_t *sequence,
                           int length,
                           const uint32_t *label_hashes,
                           const uint32_t *frontiers,
                           uint32_t *feature_ids,
                           double *feature_values,
                           int capacity) {
    uint32_t *expanded_vertex_hashes = NULL;
    uint32_t *neighborhoods = NULL;
    int last_coordinate;
    int position;
    int root_coordinate;
    int result;

    if (!sequence || !label_hashes || !frontiers || length <= 0
            || length > (int)BITMASK_16 || capacity < 0
            || (capacity > 0 && (!feature_ids || !feature_values)))
        return 0;
    for (position = 0; position < length; position++) {
        if (sequence[position] >= 128)
            return 0;
    }

    last_coordinate = 2 * length - 2;
    expanded_vertex_hashes = malloc(
        (size_t)(last_coordinate + 1) * sizeof(*expanded_vertex_hashes)
    );
    neighborhoods = calloc(
        (size_t)length * 4, sizeof(*neighborhoods)
    );
    if (!expanded_vertex_hashes || !neighborhoods) {
        result = -1;
        goto cleanup;
    }

    if (length == 1) {
        int64_t values[2];
        values[0] = (int64_t)label_hashes[sequence[0]];
        values[1] = 0;
        if (!hash_low32_plus_one(
                legacy_hash_int_tuple(values, 2),
                &expanded_vertex_hashes[0])) {
            result = -2;
            goto cleanup;
        }
    } else {
        int64_t edge_values[2];
        uint32_t edge_vertex_hash;
        int last_position = length - 1;

        edge_values[0] = (int64_t)label_hashes[(uint8_t)'-'];
        edge_values[1] = 2;
        if (!hash_low32_plus_one(
                legacy_hash_int_tuple(edge_values, 2),
                &edge_vertex_hash)) {
            result = -2;
            goto cleanup;
        }
        for (position = 0; position < length; position++) {
            int64_t values[2];
            int coordinate = position * 2;
            values[0] = (int64_t)label_hashes[sequence[position]];
            values[1] = (
                position == 0 || position == last_position ? 1 : 2
            );
            if (!hash_low32_plus_one(
                    legacy_hash_int_tuple(values, 2),
                    &expanded_vertex_hashes[coordinate])) {
                result = -2;
                goto cleanup;
            }
            if (coordinate < last_coordinate)
                expanded_vertex_hashes[coordinate + 1] = edge_vertex_hash;
        }
    }

    for (root_coordinate = 0;
            root_coordinate <= last_coordinate;
            root_coordinate += 2) {
        int eccentricity = root_coordinate;
        int other_side = last_coordinate - root_coordinate;
        int depth_limit;
        int distance;
        int radius_index = 0;
        int64_t running_hash = (int64_t)UINT32_C(0xAAAAAAAA);
        uint32_t *root_neighborhood =
            neighborhoods + (root_coordinate / 2) * 4;

        if (other_side > eccentricity)
            eccentricity = other_side;
        depth_limit = eccentricity < 6 ? eccentricity : 6;
        for (distance = 0; distance <= depth_limit; distance++) {
            int left = root_coordinate - distance;
            int right = root_coordinate + distance;
            int64_t distance_hash_signed;
            uint32_t distance_hash;
            int64_t running_values[3];

            if (left >= 0 && right <= last_coordinate && right != left) {
                int64_t values[2];
                uint32_t first = expanded_vertex_hashes[left];
                uint32_t second = expanded_vertex_hashes[right];
                if (first > second) {
                    uint32_t temporary = first;
                    first = second;
                    second = temporary;
                }
                values[0] = (int64_t)first;
                values[1] = (int64_t)second;
                distance_hash_signed = legacy_hash_int_tuple(values, 2);
            } else {
                int64_t value;
                int coordinate = left >= 0 ? left : right;
                value = (int64_t)expanded_vertex_hashes[coordinate];
                distance_hash_signed = legacy_hash_int_tuple(&value, 1);
            }
            if (!hash_low32_plus_one(
                    distance_hash_signed, &distance_hash)) {
                result = -2;
                goto cleanup;
            }

            running_values[0] = running_hash;
            running_values[1] = (int64_t)distance_hash;
            running_values[2] = (int64_t)distance;
            running_hash = unsigned_as_signed_64(
                (uint64_t)running_hash
                ^ (uint64_t)legacy_hash_int_tuple(running_values, 3)
            );
            if (distance % 2 == 0) {
                if (!hash_low32_plus_one(
                        running_hash,
                        &root_neighborhood[radius_index])) {
                    result = -2;
                    goto cleanup;
                }
                radius_index += 1;
            }
        }
    }

    result = path_feature_row(
        neighborhoods, frontiers, length,
        feature_ids, feature_values, capacity
    );

cleanup:
    free(neighborhoods);
    free(expanded_vertex_hashes);
    return result;
}
