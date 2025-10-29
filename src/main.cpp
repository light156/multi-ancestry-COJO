#include "macojo.h"


int main(int argc, char** argv) 
{    
    LOGGER.ts("main");
    LOGGER << setprecision(12);

    MACOJO macojo;

    int code = macojo.set_read_process_output_options(argc, argv);
    if (code != 0) return code;

    macojo.read_cojo_PLINK_files();
    macojo.entry_function();

    #ifdef __linux__
        float vmem = roundf(1000.0 * getVMPeakKB() / 1024/1024) / 1000.0; 
        float mem = roundf(1000.0 * getMemPeakKB() / 1024/1024) / 1000.0;     
        LOGGER << setprecision(3) << "Peak memory: " << mem << " GB; Virtual memory: " << vmem << " GB." << endl;
    #endif

    float duration = LOGGER.tp("main");
    int hours = (int) duration / 3600;
    string time_str = (hours == 0) ? "" : (to_string(hours) + " hour" + ((hours == 1) ? " ": "s "));
    int mins = (int) (duration - 3600 * hours) / 60;
    time_str += (mins == 0) ? "" : (to_string(mins) + " minute" + ((mins == 1) ? " ": "s "));
    float seconds = duration - 3600 * hours - 60 * mins;
    time_str = time_str + to_string(seconds) + " seconds";
    LOGGER << "Overall computational time: " << time_str << endl;
    
    LOGGER.close();
    return 0;
}
