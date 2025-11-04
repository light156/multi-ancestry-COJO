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

    LOGGER.tp("main");    
    LOGGER.close();
    return 0;
}
