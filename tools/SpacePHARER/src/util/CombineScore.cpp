#include "../core/CheckedInput.h"
#include "../core/Statistics.h"
#include "Aggregation.h"
#include "Debug.h"
#include "LocalParameters.h"
#include "Util.h"
#include "itoa.h"
#include <memory>

#ifdef OPENMP
#include <omp.h>
#endif

class PvalueAggregator : public Aggregation {
  public:
    PvalueAggregator(std::string queryDbName, std::string targetDbName, const std::string &resultDbName,
                     const std::string &outputDbName, unsigned int threads, unsigned int compressed)
        : Aggregation(targetDbName, resultDbName, outputDbName, threads, compressed) {

        std::string sizeDBName = queryDbName + "_set_size";
        std::string sizeDBIndex = queryDbName + "_set_size.index";
        querySizeReader = std::make_unique<DBReader<unsigned int>>(
            sizeDBName.c_str(), sizeDBIndex.c_str(), threads,
            DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        querySizeReader->open(DBReader<unsigned int>::NOSORT);

        sizeDBName = targetDbName + "_set_size";
        sizeDBIndex = targetDbName + "_set_size.index";
        targetSizeReader = std::make_unique<DBReader<unsigned int>>(
            sizeDBName.c_str(), sizeDBIndex.c_str(), threads,
            DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
        targetSizeReader->open(DBReader<unsigned int>::NOSORT);
    }

    ~PvalueAggregator() {
        targetSizeReader->close();

        querySizeReader->close();
    }

    void prepareInput(unsigned int MAYBE_UNUSED(querySetKey), unsigned int MAYBE_UNUSED(thread_idx)) {}

    // Get all result of a single Query Set VS a Single Target Set and return the multiple-match p-value for
    // it
    std::string aggregateEntry(std::vector<std::vector<std::string>> &dataToAggregate,
                               unsigned int querySetKey, unsigned int targetSetKey, unsigned int thread_idx) {

        std::string buffer;
        char keyBuffer[255];
        char *tmpBuff = Itoa::u32toa_sse2(targetSetKey, keyBuffer);
        buffer.append(keyBuffer, tmpBuff - keyBuffer - 1);
        buffer.append("\t");

        const unsigned int orfCount =
            Util::fast_atoi<unsigned int>(querySizeReader->getDataByDBKey(querySetKey, thread_idx));
        std::vector<spacepharer::HitEvidence> hits;
        hits.reserve(dataToAggregate.size());
        for (const auto &row : dataToAggregate) {
            spacepharer::require_columns(row.size(), 3);
            hits.push_back(
                {spacepharer::finite_number(row[1].c_str()), spacepharer::finite_number(row[2].c_str())});
        }
        const auto result = spacepharer::combine_score(hits, orfCount);
        buffer.append(SSTR(result.score));
        buffer.append("\t");
        buffer.append(SSTR(result.mean_identity));
        return buffer;
    }

  private:
    std::unique_ptr<DBReader<unsigned int>> querySizeReader;
    std::unique_ptr<DBReader<unsigned int>> targetSizeReader;
};

int combinescore(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    PvalueAggregator aggregation(par.db1, par.db2, par.db3, par.db4, (unsigned int)par.threads,
                                 par.compressed);
    return aggregation.run();
}
