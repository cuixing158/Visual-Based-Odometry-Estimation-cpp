/**
* @file        :CommandParser.h
* @brief       :在on inux,Windows,MacOS平台上解析命令行参数，类似OpenCV的CommandLineParser
* @details     :This is the detail description.
* @date        :2022/03/26 10:32:01
* @author      :cuixingxing(cuixingxing150@gmail.com)
* @version     :1.0
*
* @copyright Copyright (c) 2022
*
*/

// Example, C++11
//
// int main(int argc, const char **argv) {
//     auto show_help = [&]() {
//         std::cout << argv[0] << " [-h|--help|-?] [-f=file|--file=file] [-o=dstpath|--outpath=dstpath] [-p=srcpath|--path=srcpath]" << std::endl;
//         exit(0);
//     };

//     // Simple functional api. No initialization required.
//     bool help = CommandParser(false, "-h", "--help", "-?");
//     std::string outpath = CommandParser("results/", "-o", "--outpath");// 第一个为默认项
//     std::string imgPath = CommandParser("", "-p", "--path");
//     if (imgPath.empty()) {
//         show_help();
//     }

//     // OOP map-based api. Explicit (argc, argv) initialization required.
//     struct CommandParser args(argc, argv);
//     if (args.has("-h") || args.has("--help") || args.has("-?") || args.size() == 1) {
//         show_help();
//     }
//     if (args.has("-d") || args.has("--depth") || args.has("--max-depth")) {
//         std::string arg = args["-d"];
//         if (arg.empty()) arg = args["--depth"];
//         if (arg.empty()) arg = args["--max-depth"];
//         int depth = atoi(arg.c_str());
//         std::cout << "provided depth: " << depth << std::endl;
//     }
//     if (args.has("-f") || args.has("--file")) {
//         std::string arg = args["-f"];
//         if (arg.empty()) arg = args["--file"];
//         std::string fname = arg;
//         std::cout << "provided file: " << fname << std::endl;
//     }
//     std::cout << "---" << std::endl;
//     std::cout << args.cmdline() << std::endl;
//     // std::cout << args.size() << " provided args: " << args.str() << std::endl;
// }

#pragma once
#include <map>
#include <sstream>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <io.h>
// #include <shellapi.h>
// #include <winsock2.h>
#include <windows.h>
#pragma comment(lib, "Shell32.lib")
#else
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <fstream>
#include <sstream>
#endif

namespace CommandParser_utils {
// string conversion

template <typename T>
inline T as(const std::string &self) {
    T t;
    return (std::istringstream(self) >> t) ? t : (T)(self.size() && (self != "0") && (self != "false"));
}

template <>
inline char as(const std::string &self) {
    return self.size() == 1 ? (char)(self[0]) : (char)(as<int>(self));
}
template <>
inline signed char as(const std::string &self) {
    return self.size() == 1 ? (signed char)(self[0]) : (signed char)(as<int>(self));
}
template <>
inline unsigned char as(const std::string &self) {
    return self.size() == 1 ? (unsigned char)(self[0]) : (unsigned char)(as<int>(self));
}

template <>
inline const char *as(const std::string &self) {
    return self.c_str();
}
template <>
inline std::string as(const std::string &self) {
    return self;
}

// token split

inline size_t split(std::vector<std::string> &tokens, const std::string &self, const std::string &delimiters) {
    std::string str;
    tokens.clear();
    for (auto &ch : self) {
        if (delimiters.find_first_of(ch) != std::string::npos) {
            if (str.size()) tokens.push_back(str), str = "";
            tokens.push_back(std::string() + ch);
        } else
            str += ch;
    }
    return str.empty() ? tokens.size() : (tokens.push_back(str), tokens.size());
};

// portable cmdline

inline std::vector<std::string> cmdline() {
    std::vector<std::string> args;
    std::string arg;
#ifdef _WIN32
    int argv;
    auto *list = CommandLineToArgvW(GetCommandLineW(), &argv);
    if (list) {
        for (int i = 0; i < argv; ++i) {
            std::wstring ws(list[i]);
            args.push_back(std::string(ws.begin(), ws.end()));
        }
        LocalFree(list);
    }
#else
    pid_t pid = getpid();

    char fname[32] = {};
    sprintf(fname, "/proc/%d/cmdline", pid);
    std::ifstream ifs(fname);
    if (ifs.good()) {
        std::stringstream ss;
        ifs >> ss.rdbuf();
        arg = ss.str();
    }
    for (auto end = arg.size(), i = end - end; i < end; ++i) {
        auto st = i;
        while (i < arg.size() && arg[i] != '\0') ++i;
        args.push_back(arg.substr(st, i - st));
    }
#endif
    return args;
}
}  // namespace CommandParser_utils

// main map class; explicit initialization

struct CommandParser : public std::map<std::string, std::string> {
    using super = std::map<std::string, std::string>;

    CommandParser(int argc, const char **argv) : super() {
        // reconstruct vector
        std::vector<std::string> args(argc, std::string());
        for (int i = 0; i < argc; ++i) {
            args[i] = argv[i];
        }
        // create key=value and key= args as well
        for (auto &it : args) {
            std::vector<std::string> tokens;
            auto size = CommandParser_utils::split(tokens, it, "=");

            if (size == 3 && tokens[1] == "=")
                (*this)[tokens[0]] = tokens[2];
            else if (size == 2 && tokens[1] == "=")
                (*this)[tokens[0]] = true;
            else if (size == 1 && tokens[0] != argv[0])
                (*this)[tokens[0]] = true;
        }
        // recreate args
        while (argc--) {
            (*this)[std::to_string(argc)] = std::string(argv[argc]);
        }
    }

    CommandParser(const std::vector<std::string> &args) : super() {
        std::vector<const char *> argv;
        for (auto &it : args) {
            argv.push_back(it.c_str());
        }
        *this = CommandParser(argv.size(), argv.data());
    }

    size_t size() const {
        unsigned i = 0;
        while (has(std::to_string(i))) ++i;
        return i;
    }

    bool has(const std::string &op) const {
        return this->find(op) != this->end();
    }

    std::string str() const {
        std::stringstream ss;
        std::string sep;
        for (auto &it : *this) {
            ss << sep << it.first << "=" << it.second;
            sep = ',';
        }
        return ss.str();
    }

    std::string cmdline() const {
        std::stringstream cmd;
        std::string sep;
        // concatenate args
        for (auto end = size(), arg = end - end; arg < end; ++arg) {
            cmd << sep << this->find(std::to_string(arg))->second;
            sep = ' ';
        }
        return cmd.str();
    }
};

// variadic syntax sugars {

template <typename T>
inline T CommandParser(const T &defaults, const char *argv) {
    static struct CommandParser map(CommandParser_utils::cmdline());
    return map.has(argv) ? CommandParser_utils::as<T>(map[argv]) : defaults;
}

template <typename T, typename... Args>
inline T CommandParser(const T &defaults, const char *arg0, Args... argv) {
    T t = CommandParser<T>(defaults, arg0);
    return t == defaults ? CommandParser<T>(defaults, argv...) : t;
}

inline const char *CommandParser(const char *defaults, const char *argv) {
    static struct CommandParser map(CommandParser_utils::cmdline());
    return map.has(argv) ? CommandParser_utils::as<const char *>(map[argv]) : defaults;
}

template <typename... Args>
inline const char *CommandParser(const char *defaults, const char *arg0, Args... argv) {
    const char *t = CommandParser(defaults, arg0);
    return t == defaults ? CommandParser(defaults, argv...) : t;
}
