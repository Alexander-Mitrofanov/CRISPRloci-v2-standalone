#include "../core/CheckedInput.h"
#include "../core/Statistics.h"
#include "../core/RecordIndex.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int combineprotnuclaln(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> protAlnReader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                         DBReader<unsigned int>::USE_DATA |
                                             DBReader<unsigned int>::USE_INDEX);
    protAlnReader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    spacepharer::ReaderGuard<decltype(protAlnReader)> protAlnReaderGuard(protAlnReader);

    DBReader<unsigned int> nuclAlnReader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                         DBReader<unsigned int>::USE_DATA |
                                             DBReader<unsigned int>::USE_INDEX);
    nuclAlnReader.open(DBReader<unsigned int>::NOSORT);
    spacepharer::ReaderGuard<decltype(nuclAlnReader)> nuclAlnReaderGuard(nuclAlnReader);

    DBWriter dbw(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed,
                 protAlnReader.getDbtype());
    dbw.open();

    Debug::Progress progress(protAlnReader.getSize());

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        const char *protEntry[255];
        const char *nuclEntry[255];

        std::string buffer;
        buffer.reserve(1024 * 1024);

#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < protAlnReader.getSize(); ++id) {
            progress.updateProgress();

            unsigned int dbKey = protAlnReader.getDbKey(id);

            size_t nuclId = nuclAlnReader.getId(dbKey);
            if (nuclId == UINT_MAX) {
                Debug(Debug::WARNING) << "Missing nucleotide alignment result for key " << dbKey << " \n";
                continue;
            }

            char *nuclData = nuclAlnReader.getData(nuclId, thread_idx);
            const auto nucleotideRows = spacepharer::index_records(
                nuclData, Util::fast_atoi<unsigned int>, Util::skipLine);
            char *protData = protAlnReader.getData(id, thread_idx);
            while (*protData != '\0') {
                size_t columns = Util::getWordsOfLine(protData, protEntry, 255);
                if (columns < Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                    Debug(Debug::ERROR) << "Invalid alignment result record\n";
                    EXIT(EXIT_FAILURE);
                }
                protData = Util::skipLine(protData);

                size_t targetSetId = Util::fast_atoi<size_t>(protEntry[0]);
                double protEval = spacepharer::finite_number(protEntry[3]);

                double updatedEval = FLT_MAX;
                double nuclSeqId = 0;
                const auto matching = targetSetId <= UINT_MAX
                                          ? nucleotideRows.find(static_cast<unsigned int>(targetSetId))
                                          : nucleotideRows.end();
                if (matching != nucleotideRows.end()) {
                    // The historical final matching nucleotide row wins. Keep
                    // every matching row in input order and retain its arithmetic.
                    for (char *nuclCurrent : matching->second) {
                        columns = Util::getWordsOfLine(nuclCurrent, nuclEntry, 255);
                        if (columns < Matcher::ALN_RES_WITHOUT_BT_COL_CNT) {
                            Debug(Debug::ERROR) << "Invalid alignment result record\n";
                            EXIT(EXIT_FAILURE);
                        }

                        nuclSeqId = spacepharer::finite_number(nuclEntry[2]);
                        double nuclEval = spacepharer::finite_number(nuclEntry[3]);
                        updatedEval = spacepharer::combine_evalues(protEval, nuclEval);
                    }
                }
                buffer.append(protEntry[0], protEntry[2] - protEntry[0]);
                // substitute protein sequence identity with nucl sequence identity to later calculated
                // average seqid per set
                buffer.append(SSTR(nuclSeqId));
                buffer.append("\t");
                buffer.append(SSTR(updatedEval));
                buffer.append("\t");
                buffer.append(protEntry[4], (protData - 1) - protEntry[4]);
                buffer.append("\n");
            }
            dbw.writeData(buffer.c_str(), buffer.length(), dbKey, thread_idx);
            buffer.clear();
        }
    }
    dbw.close();

    return EXIT_SUCCESS;
}