#include "macojo.h"


int main(int argc, char** argv) 
{    
    auto start = steady_clock::now();

    LOGGER << setprecision(12);

    int code = set_read_process_output_options(argc, argv);
    if (code != 0) return code;

    HyperParams& params = get_params();

    for (int chr : params.chr_list) {
        params.curr_chr = chr;
        LOGGER << "\nProcessing chromosome " << chr << endl;

        MACOJO macojo;
        if (!macojo.read_input_files()) {
            LOGGER.w("No valid SNPs found on chromosome " + to_string(chr) + ", skipping chromosome\n");
            continue;
        }

        macojo.entry_function();
        LOGGER.i("===========================================\n");
    }
    
    if (params.if_output_all && bad_SNP_reason.size() > 0) {
        string badSNPfile = params.output_name + ".badsnps";
        ofstream badSNPout(badSNPfile);
        if (!badSNPout) LOGGER.e("Cannot open the file [" + badSNPfile + "] to write");

        badSNPout << "SNP\tReason\n";
        for (const auto& kv : bad_SNP_reason)
            badSNPout << kv.first << "\t" << kv.second << "\n";

        badSNPout.close();
        LOGGER.i("List of bad SNPs saved into [" + badSNPfile + "]");
    }

    auto end = steady_clock::now();

    #ifdef __linux__
        float vmem = roundf(1000.0 * getVMPeakKB() / 1024 / 1024) / 1000.0; 
        float mem = roundf(1000.0 * getMemPeakKB() / 1024 / 1024) / 1000.0;     
        LOGGER << setprecision(3) << "Peak memory: " << mem << " GB; Virtual memory: " << vmem << " GB." << endl;
    #endif

    double secs = duration<double>(end-start).count();
    int hours = int(secs / 3600);
    int mins = int((secs - hours * 3600) / 60);

    LOGGER << "Overall computational time: ";
    if (hours) LOGGER << hours << " hour" << (hours == 1 ? " " : "s ");
    if (mins)  LOGGER << mins  << " minute" << (mins == 1 ? " " : "s ");
    LOGGER << fixed << setprecision(3) << secs - hours * 3600 - mins * 60 << " seconds" << endl;

    LOGGER.close();
    return 0;
}
