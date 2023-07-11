/**
 * @file           : path.h
 * @target         : windows or linux or mac os
 * @details        : 在on inux,Windows,MacOS平台上操作文件路径
 * @author         : cuixingxing
 * @email          : cuixingxing150@gmail.com
 * @date           : 17-Mar-2022 12:02:52
 * @version        : V1.0.0
 * @copyright      : Copyright (C) 2022 cuixingxing.All rights reserved.
 */

//////////////////////////////////Example////////////////////////////////////////////////////////////
// #include <iostream>
// #include "c_cpp_utils/path.h"

// using namespace std;
// using namespace filesystem;

// int main(int argc, char **argv) {
// #if !defined(WIN32)
//     path path1("/dir 1/dir 2/");
// #else
//     path path1("C:\\dir 1\\dir 2\\");
// #endif
//     path path2("dir 3");

//     cout << path1.exists() << endl;
//     cout << path1 << endl;
//     cout << (path1 / path2) << endl;
//     cout << (path1 / path2).parent_path() << endl;
//     cout << (path1 / path2).parent_path().parent_path() << endl;
//     cout << (path1 / path2).parent_path().parent_path().parent_path() << endl;
//     cout << (path1 / path2).parent_path().parent_path().parent_path().parent_path() << endl;
//     cout << path().parent_path() << endl;
//     cout << "some/path.ext:operator==() = " << (path("some/path.ext") == path("some/path.ext")) << endl;
//     cout << "some/path.ext:operator==()
//         (unequal) = " << (path(" some / path.ext ") == path(" another / path.ext ")) << endl;

//                                                                         cout
//                     << "nonexistant:exists = " << path("nonexistant").exists() << endl;
//     cout << "nonexistant:is_file = " << path("nonexistant").is_file() << endl;
//     cout << "nonexistant:is_directory = " << path("nonexistant").is_directory() << endl;
//     cout << "nonexistant:filename
//         = " << path(" nonexistant
//           ").filename() << endl; cout <<
//           "nonexistant:extension = "
//           << path("nonexistant").extension() << endl;
//     cout << "filesystem/path.h:exists = " << path("filesystem/path.h").exists() << endl;
//     cout << "filesystem/path.h:is_file = " << path("filesystem/path.h").is_file() << endl;
//     cout << "filesystem/path.h:is_directory = " << path("filesystem/path.h").is_directory() << endl;
//     cout << "filesystem/path.h:filename = " << path("filesystem/path.h").filename() << endl;
//     cout << "filesystem/path.h:extension = " << path("filesystem/path.h").extension() << endl;
//     cout << "filesystem/path.h:make_absolute = " << path("filesystem/path.h").make_absolute() << endl;
//     cout << "../filesystem:exists = " << path("../filesystem").exists() << endl;
//     cout
//         << "../filesystem:is_file = " << path("../filesystem").is_file() << endl;
//     cout << "../filesystem:is_directory = " << path("../filesystem").is_directory() << endl;
//     cout << "../filesystem:extension = " << path("../filesystem").extension() << endl;
//     cout << "../filesystem:filename = " << path("../filesystem").filename() << endl;
//     cout << "../filesystem:make_absolute = " << path("../filesystem").make_absolute() << endl;

//     cout << "resolve(filesystem/path.h) = " << resolver().resolve("filesystem/path.h") << endl;
//     cout << "resolve(nonexistant) = " << resolver().resolve("nonexistant") << endl;
//     return 0;
// }
//
////////////////////////////////////////////////////////////////////////////////////////////
#pragma once
#include <dirent.h>

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

//#include "fwd.h"

#if defined(_WIN32)
#include <ShlObj.h>
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <sys/stat.h>

#if defined(__linux)
#include <linux/limits.h>
#endif

namespace filesystem {

/**
 * @brief Simple class for manipulating paths on Linux/Windows/Mac OS
 *
 * This class is just a temporary workaround to avoid the heavy boost
 * dependency until boost::filesystem is integrated into the standard template
 * library at some point in the future.
 */
class path {
   public:
    enum path_type {
        windows_path = 0,
        posix_path = 1,
#if defined(_WIN32)
        native_path = windows_path
#else
        native_path = posix_path
#endif
    };

    path() : m_type(native_path), m_absolute(false), m_smb(false) {}

    path(const path &path)
        : m_type(path.m_type),
          m_path(path.m_path),
          m_absolute(path.m_absolute),
          m_smb(path.m_smb) {}

    path(path &&path)
        : m_type(path.m_type),
          m_path(std::move(path.m_path)),
          m_absolute(path.m_absolute),
          m_smb(path.m_smb) {}

    path(const char *string) : m_smb(false) { set(string); }

    path(const std::string &string) : m_smb(false) { set(string); }

#if defined(_WIN32)
    path(const std::wstring &wstring) { set(wstring); }
    path(const wchar_t *wstring) { set(wstring); }
#endif

    // 获取path路径有多少个
    size_t length() const { return m_path.size(); }

    // 判断path是否有路径
    bool empty() const { return m_path.empty(); }

    //判断是否为绝对路径
    bool is_absolute() const { return m_absolute; }

    // 获取绝对路径
    path make_absolute() const {
#if !defined(_WIN32)
        char temp[PATH_MAX];
        if (realpath(str().c_str(), temp) == NULL)
            throw std::runtime_error("Internal error in realpath(): " +
                                     std::string(strerror(errno)));
        return path(temp);
#else
        std::wstring value = wstr(), out(MAX_PATH_WINDOWS, '\0');
        DWORD length =
            GetFullPathNameW(value.c_str(), MAX_PATH_WINDOWS, &out[0], NULL);
        if (length == 0)
            throw std::runtime_error("Internal error in realpath(): " +
                                     std::to_string(GetLastError()));
        return path(out.substr(0, length));
#endif
    }

    // 判断路径或者文件是否存在
    bool exists() const {
#if defined(_WIN32)
        return GetFileAttributesW(wstr().c_str()) != INVALID_FILE_ATTRIBUTES;
#else
        struct stat sb;
        return stat(str().c_str(), &sb) == 0;
#endif
    }

    // 获取文件大小
    size_t file_size() const {
#if defined(_WIN32)
        struct _stati64 sb;
        if (_wstati64(wstr().c_str(), &sb) != 0)
            throw std::runtime_error("path::file_size(): cannot stat file \"" +
                                     str() + "\"!");
#else
        struct stat sb;
        if (stat(str().c_str(), &sb) != 0)
            throw std::runtime_error("path::file_size(): cannot stat file \"" +
                                     str() + "\"!");
#endif
        return (size_t)sb.st_size;
    }
    // 判断path是否为路径
    bool is_directory() const {
#if defined(_WIN32)
        DWORD result = GetFileAttributesW(wstr().c_str());
        if (result == INVALID_FILE_ATTRIBUTES) return false;
        return (result & FILE_ATTRIBUTE_DIRECTORY) != 0;
#else
        struct stat sb;
        if (stat(str().c_str(), &sb)) return false;
        return S_ISDIR(sb.st_mode);
#endif
    }

    // 判断path是否为文件
    bool is_file() const {
#if defined(_WIN32)
        DWORD attr = GetFileAttributesW(wstr().c_str());
        return (attr != INVALID_FILE_ATTRIBUTES &&
                (attr & FILE_ATTRIBUTE_DIRECTORY) == 0);
#else
        struct stat sb;
        if (stat(str().c_str(), &sb)) return false;
        return S_ISREG(sb.st_mode);
#endif
    }

    // 获取扩展名
    std::string extension() const {
        const std::string &name = filename();
        size_t pos = name.find_last_of(".");
        if (pos == std::string::npos) return "";
        return name.substr(pos + 1);
    }

    // 获取文件名，包含扩展名
    std::string filename() const {
        if (empty()) return "";
        const std::string &last = m_path[m_path.size() - 1];
        return last;
    }

    // 获取文件名，不包含扩展名
    std::string filenameNoExt() const {
        if (empty()) return "";
        const std::string &last = m_path[m_path.size() - 1];
        size_t pos = last.find_last_of(".");
        if (pos == std::string::npos) return "";
        return last.substr(0, pos);
    }

    // 获取父路径
    path parent_path() const {
        path result;
        result.m_absolute = m_absolute;
        result.m_smb = m_smb;

        if (m_path.empty()) {
            if (!m_absolute) result.m_path.push_back("..");
        } else {
            size_t until = m_path.size() - 1;
            for (size_t i = 0; i < until; ++i) result.m_path.push_back(m_path[i]);
        }
        return result;
    }

    // 重载文件分隔符
    path operator/(const path &other) const {
        if (other.m_absolute)
            throw std::runtime_error("path::operator/(): expected a relative path!");
        if (m_type != other.m_type)
            throw std::runtime_error(
                "path::operator/(): expected a path of the same type!");

        path result(*this);

        for (size_t i = 0; i < other.m_path.size(); ++i)
            result.m_path.push_back(other.m_path[i]);

        return result;
    }

    // path对象转换为字符串路径
    std::string str(path_type type = native_path) const {
        std::ostringstream oss;

        if (m_absolute) {
            if (m_type == posix_path)
                oss << "/";
            else {
                size_t length = 0;
                for (size_t i = 0; i < m_path.size(); ++i)
                    // No special case for the last segment to count the NULL character
                    length += m_path[i].length() + 1;
                if (m_smb) length += 2;

                // Windows requires a \\?\ prefix to handle paths longer than MAX_PATH
                // (including their null character). NOTE: relative paths >MAX_PATH are
                // not supported at all in Windows.
                if (length > MAX_PATH_WINDOWS_LEGACY) {
                    if (m_smb)
                        oss << "\\\\?\\UNC\\";
                    else
                        oss << "\\\\?\\";
                } else if (m_smb)
                    oss << "\\\\";
            }
        }

        for (size_t i = 0; i < m_path.size(); ++i) {
            oss << m_path[i];
            if (i + 1 < m_path.size()) {
                if (type == posix_path)
                    oss << '/';
                else
                    oss << '\\';
            }
        }

        return oss.str();
    }

    void set(const std::string &str, path_type type = native_path) {
        m_type = type;
        if (type == windows_path) {
            std::string tmp = str;

            // Long windows paths (sometimes) begin with the prefix \\?\. It should
            // only be used when the path is >MAX_PATH characters long, so we remove
            // it for convenience and add it back (if necessary) in str()/wstr().
            static const std::string LONG_PATH_PREFIX = "\\\\?\\";
            if (tmp.length() >= LONG_PATH_PREFIX.length() &&
                std::mismatch(std::begin(LONG_PATH_PREFIX),
                              std::end(LONG_PATH_PREFIX), std::begin(tmp))
                        .first == std::end(LONG_PATH_PREFIX)) {
                tmp.erase(0, LONG_PATH_PREFIX.length());
            }

            // Special-case handling of absolute SMB paths, which start with the
            // prefix "\\".
            if (tmp.length() >= 2 && tmp[0] == '\\' && tmp[1] == '\\') {
                m_path = {};
                tmp.erase(0, 2);

                // Interestingly, there is a special-special case where relative paths
                // may be specified as beginning with a "\\" when a non-SMB file with a
                // more-than-260-characters-long absolute _local_ path is
                // double-clicked. This seems to only happen with single-segment
                // relative paths, so we can check for this condition by making sure no
                // further path separators are present.
                if (tmp.find_first_of("/\\") != std::string::npos)
                    m_absolute = m_smb = true;
                else
                    m_absolute = m_smb = false;

                // Special-case handling of absolute SMB paths, which start with the
                // prefix "UNC\"
            } else if (tmp.length() >= 4 && tmp[0] == 'U' && tmp[1] == 'N' &&
                       tmp[2] == 'C' && tmp[3] == '\\') {
                m_path = {};
                tmp.erase(0, 4);
                m_absolute = true;
                m_smb = true;
                // Special-case handling of absolute local paths, which start with the
                // drive letter and a colon "X:\" So that UTF-8 works, do not call
                // std::isalpha if the high bit is set, as that causes an assert on
                // Windows.
            } else if (tmp.length() >= 3 && ((unsigned char)tmp[0] < 0x80) &&
                       std::isalpha(tmp[0]) && tmp[1] == ':' &&
                       (tmp[2] == '\\' || tmp[2] == '/')) {
                m_path = {tmp.substr(0, 2)};
                tmp.erase(0, 3);
                m_absolute = true;
                m_smb = false;
                // Relative path
            } else {
                m_path = {};
                m_absolute = false;
                m_smb = false;
            }

            std::vector<std::string> tokenized = tokenize(tmp, "/\\");
            m_path.insert(std::end(m_path), std::begin(tokenized),
                          std::end(tokenized));
        } else {
            m_path = tokenize(str, "/");
            m_absolute = !str.empty() && str[0] == '/';
        }
    }

    path &operator=(const path &path) {
        m_type = path.m_type;
        m_path = path.m_path;
        m_absolute = path.m_absolute;
        m_smb = path.m_smb;
        return *this;
    }

    path &operator=(path &&path) {
        if (this != &path) {
            m_type = path.m_type;
            m_path = std::move(path.m_path);
            m_absolute = path.m_absolute;
            m_smb = path.m_smb;
        }
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const path &path) {
        os << path.str();
        return os;
    }

    bool remove_file() {
#if !defined(_WIN32)
        return std::remove(str().c_str()) == 0;
#else
        return DeleteFileW(wstr().c_str()) != 0;
#endif
    }

    bool resize_file(size_t target_length) {
#if !defined(_WIN32)
        return ::truncate(str().c_str(), (off_t)target_length) == 0;
#else
        HANDLE handle = CreateFileW(wstr().c_str(), GENERIC_WRITE, 0, nullptr, 0,
                                    FILE_ATTRIBUTE_NORMAL, nullptr);
        if (handle == INVALID_HANDLE_VALUE) return false;
        LARGE_INTEGER size;
        size.QuadPart = (LONGLONG)target_length;
        if (SetFilePointerEx(handle, size, NULL, FILE_BEGIN) == 0) {
            CloseHandle(handle);
            return false;
        }
        if (SetEndOfFile(handle) == 0) {
            CloseHandle(handle);
            return false;
        }
        CloseHandle(handle);
        return true;
#endif
    }

    static path getcwd() {
#if !defined(_WIN32)
        char temp[PATH_MAX];
        if (::getcwd(temp, PATH_MAX) == NULL)
            throw std::runtime_error("Internal error in getcwd(): " +
                                     std::string(strerror(errno)));
        return path(temp);
#else
        std::wstring temp(MAX_PATH_WINDOWS, '\0');
        if (!_wgetcwd(&temp[0], MAX_PATH_WINDOWS))
            throw std::runtime_error("Internal error in getcwd(): " +
                                     std::to_string(GetLastError()));
        return path(temp.c_str());
#endif
    }

#if defined(_WIN32)
    std::wstring wstr(path_type type = native_path) const {
        std::string temp = str(type);
        int size =
            MultiByteToWideChar(CP_UTF8, 0, &temp[0], (int)temp.size(), NULL, 0);
        std::wstring result(size, 0);
        MultiByteToWideChar(CP_UTF8, 0, &temp[0], (int)temp.size(), &result[0],
                            size);
        return result;
    }

    void set(const std::wstring &wstring, path_type type = native_path) {
        std::string string;
        if (!wstring.empty()) {
            int size = WideCharToMultiByte(CP_UTF8, 0, &wstring[0],
                                           (int)wstring.size(), NULL, 0, NULL, NULL);
            string.resize(size, 0);
            WideCharToMultiByte(CP_UTF8, 0, &wstring[0], (int)wstring.size(),
                                &string[0], size, NULL, NULL);
        }
        set(string, type);
    }

    path &operator=(const std::wstring &str) {
        set(str);
        return *this;
    }
#endif

    bool operator==(const path &p) const { return p.m_path == m_path; }
    bool operator!=(const path &p) const { return p.m_path != m_path; }

   protected:
    static std::vector<std::string> tokenize(const std::string &string,
                                             const std::string &delim) {
        std::string::size_type lastPos = 0,
                               pos = string.find_first_of(delim, lastPos);
        std::vector<std::string> tokens;

        while (lastPos != std::string::npos) {
            if (pos != lastPos)
                tokens.push_back(string.substr(lastPos, pos - lastPos));
            lastPos = pos;
            if (lastPos == std::string::npos || lastPos + 1 == string.length()) break;
            pos = string.find_first_of(delim, ++lastPos);
        }

        return tokens;
    }

   protected:
#if defined(_WIN32)
    static const size_t MAX_PATH_WINDOWS = 32767;
#endif
    static const size_t MAX_PATH_WINDOWS_LEGACY =
        260;                          // windows平台路径字符串最大长度
    path_type m_type;                 // 特定平台的路径类型，枚举类型，0为windows，1为Linux
    std::vector<std::string> m_path;  // 按照文件分隔符分隔的字符串向量
    bool m_absolute;                  // 是否为绝对路径
    bool m_smb;                       // Unused, except for on Windows
};

// 创建单个目录
inline bool create_directory(const path &p) {
#if defined(_WIN32)
    return CreateDirectoryW(p.wstr().c_str(), NULL) != 0;
#else
    return mkdir(p.str().c_str(), S_IRWXU) == 0;
#endif
}

//创建多层目录级别
inline bool create_directories(const path &p) {
#if defined(_WIN32)
    return SHCreateDirectory(nullptr, p.make_absolute().wstr().c_str()) ==
           ERROR_SUCCESS;
#else
    if (create_directory(p.str().c_str())) return true;

    if (p.empty()) return false;

    if (errno == ENOENT) {
        if (create_directory(p.parent_path()))
            return mkdir(p.str().c_str(), S_IRWXU) == 0;
        else
            return false;
    }
    return false;
#endif
}

// 获取path路径下（包含子路径、子孙路径）所有以wildchard为后缀的全路径图片名字,返回图片数量
inline size_t getFullNames(const path &p,
                           std::vector<std::string> &saveFileNames,
                           const char *wildchard = ".jpg") {
    bool gIgnoreHidden = true;  // 忽略隐藏文件
    size_t numsFiles = 0;
    saveFileNames.clear();
    path fullpath = p.make_absolute();
    std::string curr_path = p.str();

    // std::stringstream stream;
    // stream << fullpath;
    // stream >> curr_path;
    const char *srcPath = curr_path.c_str();

    std::regex reg_obj(wildchard, std::regex::icase);
    DIR *dirFile = opendir(srcPath);
    if (dirFile) {
        struct dirent *hFile;
        while ((hFile = readdir(dirFile)) != NULL) {
            if (!strcmp(hFile->d_name, ".")) continue;
            if (!strcmp(hFile->d_name, "..")) continue;
            if (gIgnoreHidden && (hFile->d_name[0] == '.'))
                continue;  // in linux hidden files all start with '.'

            path tmp_path = path(curr_path + "/" + std::string(hFile->d_name));
            if (!tmp_path.is_file()) {
                std::vector<std::string> tempFileNames;
                size_t tempNums = getFullNames(tmp_path, tempFileNames, wildchard);
                numsFiles += tempNums;
                saveFileNames.insert(saveFileNames.end(), tempFileNames.begin(),
                                     tempFileNames.end());
            } else {
                if (std::regex_search(hFile->d_name, reg_obj)) {
                    std::string full_path = curr_path + "/" + hFile->d_name;
                    saveFileNames.push_back(full_path);
                }
            }
        }
        closedir(dirFile);
    }
    return saveFileNames.size();
}

/**
 * This convenience class looks for a file or directory given its name
 * and a set of search paths. The implementation walks through the
 * search paths in order and stops once the file is found.
 */
class resolver {
   public:
    typedef std::vector<path>::iterator iterator;
    typedef std::vector<path>::const_iterator const_iterator;

    resolver() { m_paths.push_back(path::getcwd()); }

    size_t size() const { return m_paths.size(); }

    iterator begin() { return m_paths.begin(); }
    iterator end() { return m_paths.end(); }

    const_iterator begin() const { return m_paths.begin(); }
    const_iterator end() const { return m_paths.end(); }

    void erase(iterator it) { m_paths.erase(it); }

    void prepend(const path &path) { m_paths.insert(m_paths.begin(), path); }
    void append(const path &path) { m_paths.push_back(path); }
    const path &operator[](size_t index) const { return m_paths[index]; }
    path &operator[](size_t index) { return m_paths[index]; }

    path resolve(const path &value) const {
        for (const_iterator it = m_paths.begin(); it != m_paths.end(); ++it) {
            path combined = *it / value;
            if (combined.exists()) return combined;
        }
        return value;
    }

    friend std::ostream &operator<<(std::ostream &os, const resolver &r) {
        os << "resolver[" << std::endl;
        for (size_t i = 0; i < r.m_paths.size(); ++i) {
            os << "  \"" << r.m_paths[i] << "\"";
            if (i + 1 < r.m_paths.size()) os << ",";
            os << std::endl;
        }
        os << "]";
        return os;
    }

   private:
    std::vector<path> m_paths;
};

}  // namespace filesystem
