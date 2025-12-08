#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <ios>
#include <iostream>
#include <iomanip>

#ifdef _WIN32
#include <io.h>
#include <stdio.h>
#define ISATTY(x) _isatty(_fileno(x))
#else 
#include <unistd.h>
#define ISATTY(x) isatty(fileno(x))
#endif

#define LOGGER Logger::instance()

using std::string;
using std::cout;
using std::endl;


class Logger {
public:
    enum Level{INFO, PROMPT, WARN, ERROR};

    Logger(){
        #ifndef _WIN32
            use_color = ISATTY(stdout);
        #else
            use_color = false;
        #endif
    }

    // Singleton accessor
    static Logger& instance() {
        static Logger inst;
        return inst;
    }

    void e(const string& message) {
        *this << ERROR << "Error: " << INFO << message << endl;
        exit(EXIT_FAILURE);
    }    
    
    void i(const string& message) {
        *this << INFO << message << endl;
    }
    
    void w(const string &message) {
        *this << WARN << "Warning: " << INFO << message << endl;
    }

    template <typename T>
    void e(const string& message, const T& title) {
        *this << ERROR << title << " " << INFO << message << endl;
        exit(EXIT_FAILURE);
    }

    template <typename T>
    void i(const string& message, const T& title) {
        *this << PROMPT << title << " " << INFO << message << endl;
    }

    template <typename T>
    void w(const string &message, const T& title) {
        *this << WARN << title << " " << INFO << message << endl;
    }

    Logger& operator<<(Level level) {
        if (use_color) {
            switch (level) {
                case INFO:   cout << "\033[0m";     break; // reset
                case PROMPT: cout << "\033[0;32m";  break; // green
                case WARN:   cout << "\033[0;33m";  break; // yellow
                case ERROR:  cout << "\033[0;31m";  break; // red
                default:     cout << "\033[0m";     break;
            }
        }
        return *this;
    }

    template<typename T>
    Logger& operator<<(const T& t) {
        cout << t;
        logFile << t;
        return *this;
    }

    Logger& operator<<(std::ostream& (*op)(std::ostream&)) {
        (*op)(cout);
        (*op)(logFile);
        return *this;
    }

    Logger& operator<<(std::ios& (*pf)(std::ios&)) {
        cout << pf;
        logFile << pf;
        return *this;
    }

    Logger& operator<<(std::ios_base& (*pf)(std::ios_base&)) {
        cout << pf;
        logFile << pf;
        return *this;
    }

    void open(string ofile){
        if (ofile.empty()) e("Please provide valid log file name");
        
        logFile.open(ofile, std::ios::out);
        if(!logFile.is_open()) e("Cannot write to log file [" + ofile + "], please check file/folder permission or disk quota");
    }

    void close() {
        logFile.close();
    }

private:
    bool use_color;
    std::ofstream logFile;
};


inline int getMemPeakKB() {
    std::ifstream file("/proc/self/status");
    if (!file.is_open()) return -1;

    std::string line;
    while (std::getline(file, line)) {
        if (line.rfind("VmHWM:", 0) == 0) { // starts with "VmHWM:"
            auto pos = line.find_first_of("0123456789");
            if (pos == std::string::npos) return -1;

            std::istringstream iss(line.substr(pos));
            int value = 0;
            iss >> value;
            return value;
        }
    }
    return -1;
}


inline int getVMPeakKB() {
    std::ifstream file("/proc/self/status");
    if (!file.is_open()) return -1;

    std::string line;
    while (std::getline(file, line)) {
        if (line.rfind("VmPeak:", 0) == 0) { // starts with "VmPeak:"
            auto pos = line.find_first_of("0123456789");
            if (pos == std::string::npos) return -1;

            std::istringstream iss(line.substr(pos));
            int value = 0;
            iss >> value;
            return value;
        }
    }
    return -1;
}
