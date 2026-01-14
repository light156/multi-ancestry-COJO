#include "include/macojo.h"


int main(int argc, char** argv) 
{    
    auto start = steady_clock::now();

    // PLINK score mode, hidden from user
    bool score_mode = false;

    for (int i = 1; i < argc; i++) {
        string arg_str = argv[i];
        if (arg_str == "--score") {
            score_mode = true;
            break;
        }
    }

    if (score_mode) {
        SharedData score_shared;
        Cohort c(score_shared, 0);
        c.calc_polygenic_score(argc, argv);
    } else {
        LOGGER << setprecision(10);

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
