#pragma once

#include <cstdio>
#include <format>
#include <ostream>
#include <string>
#include <utility>

template<class... Args>
void myprint(std::format_string<Args...> fmt, Args&&... args)
{
    const auto s = format(fmt, std::forward<Args>(args)...);
    fwrite(s.data(), 1, s.size(), stdout);
}

template<class... Args>
void myprintln(std::format_string<Args...> fmt, Args&&... args)
{
    myprint("{}\n", format(fmt, std::forward<Args>(args)...));
}

inline void myprintln()
{
    fputc('\n', stdout);
}

template<class... Args>
void myprint(std::ostream& out, std::format_string<Args...> fmt, Args&&... args)
{
    const auto s = format(out.getloc(), fmt, std::forward<Args>(args)...);
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
}

template<class... Args>
void myprintln(std::ostream& out, std::format_string<Args...> fmt, Args&&... args)
{
    const auto s = format(out.getloc(), fmt, std::forward<Args>(args)...);
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
    out.put('\n');
}

inline void myprintln(std::ostream& out)
{
    out.put('\n');
}


