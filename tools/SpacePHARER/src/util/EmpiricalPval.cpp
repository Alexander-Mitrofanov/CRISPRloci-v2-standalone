#include "../core/CheckedInput.h"
#include "../core/Statistics.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "FastSort.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "Util.h"

#ifdef OPENMP
#include <omp.h>
#endif

int empiricalpval(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> posEvalDb(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                     DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    posEvalDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    spacepharer::ReaderGuard<decltype(posEvalDb)> posEvalDbGuard(posEvalDb);

    Debug::Progress progress(posEvalDb.getSize());

    // sort cEval of neg evalDB
    DBReader<unsigned int> negEvalDb(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                     DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA);
    negEvalDb.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    spacepharer::ReaderGuard<decltype(negEvalDb)> negEvalDbGuard(negEvalDb);

    std::vector<double> negToSort;

#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif

        std::vector<double> threadNegToSort;
        const char *entry[255];
#pragma omp for schedule(dynamic, 10)
        // TODO: check if control cEval list is empty
        for (size_t id = 0; id < negEvalDb.getSize(); id++) {
            progress.updateProgress();
            char *data = negEvalDb.getData(id, thread_idx);
            while (*data != '\0') {
                size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid result record.\n";
                    EXIT(EXIT_FAILURE);
                }
                double eval = spacepharer::finite_number(entry[1]);
                threadNegToSort.emplace_back(eval);
                data = Util::skipLine(data);
            }
        }
#pragma omp critical
        negToSort.insert(negToSort.end(), threadNegToSort.begin(), threadNegToSort.end());
    }
    SORT_PARALLEL(negToSort.begin(), negToSort.end());

    // compute empirical P-val for each pos entry based on neg evalDB
    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed,
                    Parameters::DBTYPE_GENERIC_DB);
    writer.open();
    progress.reset(posEvalDb.getSize());
#pragma omp parallel
    {
        int thread_idx = 0;
#ifdef OPENMP
        thread_idx = omp_get_thread_num();
#endif
        std::string buffer;
        buffer.reserve(1000000);
        const char *entry[255];
        // const double evalCoef = 1.0 ;// * numTargetSets / numControlSets;

        double empPval;
#pragma omp for schedule(dynamic, 10)
        for (size_t id = 0; id < posEvalDb.getSize(); id++) {
            progress.updateProgress();
            char *data = posEvalDb.getData(id, thread_idx);
            while (*data != '\0') {
                size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns < 2) {
                    Debug(Debug::ERROR) << "Invalid result record.\n";
                    EXIT(EXIT_FAILURE);
                }
                double eval = spacepharer::finite_number(entry[1]);
                empPval = spacepharer::empirical_tail(negToSort, eval);
                const char *current = data;
                data = Util::skipLine(data);
                size_t length = data - current;
                std::string line(current, length - 1);
                buffer.append(line);
                buffer.append("\t");
                buffer.append(SSTR(empPval));
                buffer.append("\n");
            }
            writer.writeData(buffer.c_str(), buffer.length(), posEvalDb.getDbKey(id), thread_idx);
            buffer.clear();
        }
    }
    writer.close();

    return EXIT_SUCCESS;
}
