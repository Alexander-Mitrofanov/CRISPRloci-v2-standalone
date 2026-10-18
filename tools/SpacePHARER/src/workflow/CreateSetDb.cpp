#include "CommandCaller.h"
#include "DBReader.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "Util.h"
#include "Workflow.h"
#include "createsetdb.sh.h"

void setcreatesetdbDefaults(Parameters *p) { p->orfMinLength = 30; }

int createsetdb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setcreatesetdbDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, 0);

    // check if temp dir exists and if not, try to create it:
    std::string tmpDir = par.filenames.back();
    par.filenames.pop_back();
    std::string hash = SSTR(par.hashParameter(command.databases, par.filenames, par.createsetdbworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);

    std::string outDb = par.filenames.back();
    par.filenames.pop_back();

    CommandCaller cmd;
    cmd.addVariable("OUTDB", outDb.c_str());
    cmd.addVariable("TMP_PATH", tmpDir.c_str());
    cmd.addVariable("REMOVE_TMP", spacepharer::enabled(par.removeTmpFiles));
    cmd.addVariable("REVERSE_FRAGMENTS", spacepharer::enabled(par.reverseFragments == 1));
    cmd.addVariable("EXTRACTORFS_SPACER", spacepharer::enabled(par.extractorfsSpacer == 1));
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("EXTRACTORFS_PAR", par.createParameterString(par.extractorfs).c_str());
    cmd.addVariable("TRANSLATENUCS_PAR", par.createParameterString(par.translatenucs).c_str());
    cmd.addVariable("SWAPDB_PAR", par.createParameterString(par.swapdb).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("VERBOSITY_PAR", par.createParameterString(par.onlyverbosity).c_str());
    std::vector<MMseqsParameter *> createtaxdb;
    createtaxdb.push_back(&(par.PARAM_TAX_MAPPING_FILE));
    createtaxdb.push_back(&(par.PARAM_NCBI_TAX_DUMP));
    cmd.addVariable("TAXMAPPING", par.taxMappingFile.empty() == false ? par.taxMappingFile.c_str() : nullptr);
    cmd.addVariable("CREATETAXDB_PAR", par.createParameterString(createtaxdb).c_str());

    std::string program(tmpDir + "/createsetdb.sh");
    return spacepharer::execute_workflow(cmd, program, createsetdb_sh, createsetdb_sh_len, par.filenames);
}
