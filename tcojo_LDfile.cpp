/*
#include "tcojo.h"


void Cohort::count_LD_number(string PLINKfile) 
{   
    string bimFile = PLINKfile+".bim";
    ifstream Bim(bimFile.c_str());
    if (!Bim) LOGGER.e(0, "cannot open BIM file [" + bimFile + "] to read");

    // Read bim file
    int SNP_num = 0;
    string SNP_buf, str_buf;
    int i_buf;

    map<int, string> LD_map;

    while (Bim) {
        Bim >> str_buf;
        if (Bim.eof()) break;
        Bim >> SNP_buf >> str_buf >> i_buf >> str_buf >> str_buf;
        LD_map.insert(LD_map.end(), make_pair(i_buf, SNP_buf));
        SNP_num++;
    }

    Bim.clear();
    Bim.close();
    
    auto iter1 = LD_map.begin(), iter2 = LD_map.begin();
    int temp_num = 0;
    double temp_num_sum = 0;

    while (iter1 != LD_map.end()) {
        while (iter2 != LD_map.end() && iter2->first < iter1->first + 2000000) {
            temp_num++; iter2++;
        }

        temp_num_sum += temp_num;
        iter1++;
        temp_num--;
        
    }

    cout << LD_map.size() << " " << temp_num_sum / LD_map.size() << endl;
}
*/