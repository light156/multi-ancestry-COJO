#include "macojo.h"


int main(int argc, char** argv) 
{    
    LOGGER << setprecision(12);

    auto start = chrono::steady_clock::now();

    MACOJO macojo;

    int code = macojo.set_read_process_output_options(argc, argv);
    if (code != 0) return code;

    macojo.read_input_files();
    macojo.entry_function();

    auto end = chrono::steady_clock::now();

    #ifdef __linux__
        float vmem = roundf(1000.0 * getVMPeakKB() / 1024/1024) / 1000.0; 
        float mem = roundf(1000.0 * getMemPeakKB() / 1024/1024) / 1000.0;     
        LOGGER << setprecision(3) << "Peak memory: " << mem << " GB; Virtual memory: " << vmem << " GB." << endl;
    #endif

    double secs = chrono::duration<double>(end-start).count();
    int hours = int(secs / 3600);
    int mins = int((secs - hours * 3600) / 60);

    LOGGER << "Overall computational time: ";
    if (hours) LOGGER << hours << " hour" << (hours == 1 ? " " : "s ");
    if (mins)  LOGGER << mins  << " minute" << (mins == 1 ? " " : "s ");
    LOGGER << fixed << setprecision(3) << secs - hours * 3600 - mins * 60 << " seconds" << endl;

    LOGGER.close();
    return 0;
}
