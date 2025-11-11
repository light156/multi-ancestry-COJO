#pragma once

#include <string>
#include <fstream>
#include <ios>
#include <iostream>

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
using std::endl;
using std::cout;


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

    const char* color_code(Level level) const {
        if (!use_color) return "";
        switch (level) {
            case INFO:   return "\033[0m";     // reset
            case PROMPT: return "\033[0;32m";  // green
            case WARN:   return "\033[0;33m";  // yellow
            case ERROR:  return "\033[0;31m";  // red
            default:     return "\033[0m";
        }
    }

    void e(const string& message, const string& title = "") {
        string head = title.empty() ? "Error: " : (title+" ");
        *this << ERROR << head << INFO << message << endl;
        *this << INFO << "" << INFO << "An error occurs, please check the options or data" << endl;
        exit(EXIT_FAILURE);
    }

    void i(const string& message, const string& title = "") {
        string head = title.empty() ? "" : (title + " ");
        *this << PROMPT << head << INFO << message << endl;
    }

    void w(const string &message, const string& title = "") {
        string head = title.empty() ? "Warning: " : (title + " ");
        *this << WARN << head << INFO << message << endl;
    }

    Logger& operator<<(Level level) {
        cout << color_code(level);
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

    int precision(int p) {
        cout.precision(p);
        logFile.precision(p);
        return cout.precision();
    }

    string setprecision(int p) {
        cout.precision(p);
        logFile.precision(p);
        return "";
    }

    void open(string ofile){
        if(logFile.is_open() && !filename.empty()){
            cout << "Logger has been set, not support to set another time" << endl;
        } else {
            logFile.open(ofile, std::ios::out);
            filename = ofile;
            if(!(logFile.is_open() && !filename.empty())){
                e("can't write to log file [" + ofile + "].\nPlease check file/folder permission or disk quota.");
            }
        }
    }

    void close() {
        logFile.close();
    }

private:
    bool use_color;
    string filename;
    std::ofstream logFile;
};
