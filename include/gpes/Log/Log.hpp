#pragma once 


#include <string>
#include <string_view>
#include <memory>
#include <source_location>
#include <format>
#include <concepts>

#include "Core/utility.hpp"

namespace gpes::log {

    enum class LogVerbosity : uint8_t {
        NoLogging = 0,
        Display,
        Warning,
        Error,
        Log,
        Fatal
    };

    struct LogCategory {
        explicit LogCategory(const std::string& name) : name_(name) {}
        std::string name() const {return name_;}
    private:
        const std::string name_;
    };

    class Log {
    public:
        static Log& getInstance(){
            static Log instance;
            return instance;
        }

    // !! Do not use the direct call; it's unsafe. Use macro version with static checks !!
    void log(const LogCategory& category,  //
        LogVerbosity verbosity,            //
        const std::string& message,        //
        bool showLocation = false,         //
        const std::source_location location = std::source_location::current()) const;

    private: 
        Log();
        ~Log();

        class Impl;
        std::unique_ptr<Impl> pImpl_;
    };






    constexpr LogVerbosity c_minVerbosity = LogVerbosity::Display;
    constexpr LogVerbosity c_maxVerbosity = LogVerbosity::Fatal;

    // concepts
    template <typename T>
    concept ValidLogCategory = std::constructible_from<LogCategory, T>;

    template <typename T>
    concept LoggableMessage = std::convertible_to<T, std::string> || std::convertible_to<T, std::string_view>;

    template <LogVerbosity V>
    concept ValidVerbosityLevel = V == LogVerbosity::NoLogging   //
                                || V == LogVerbosity::Display  //
                                || V == LogVerbosity::Warning  //
                                || V == LogVerbosity::Error    //
                                || V == LogVerbosity::Log      //
                                || V == LogVerbosity::Fatal;

}


#define DEFINE_LOG_CATEGORY_STATIC(logName)       \
    namespace                                     \
    {                                             \
        const gpes::log::LogCategory logName(#logName); \
    }

#define GPES_LOG_IMPL(categoryName, verbosity, showLocation, formatStr, ...)                                                        \
    do                                                                                                                            \
    {                                                                                                                             \
        if constexpr (gpes::log::LogVerbosity::verbosity >= gpes::log::c_minVerbosity &&                                              \
                      gpes::log::LogVerbosity::verbosity <= gpes::log::c_maxVerbosity)                                                \
        {                                                                                                                         \
            static_assert(gpes::log::ValidVerbosityLevel<gpes::log::LogVerbosity::verbosity>,                                         \
                "Verbosity must be one of: NoLogging, Display, Warning, Error, Log, Fatal");                                      \
            static_assert(gpes::log::ValidLogCategory<decltype(categoryName)>, "Category must be of type LogCategory");             \
            static_assert(                                                                                                        \
                gpes::log::LoggableMessage<decltype(formatStr)>, "Message must be convertible to std::string or std::string_view"); \
            gpes::log::Log::getInstance().log(                                                                                      \
                categoryName, gpes::log::LogVerbosity::verbosity, std::format(formatStr __VA_OPT__(, ) __VA_ARGS__), showLocation); \
        }                                                                                                                         \
    } while (0)

#define GPES_LOG(categoryName, verbosity, formatStr, ...) GPES_LOG_IMPL(categoryName, verbosity, false, formatStr, __VA_ARGS__)
#define GPES_LOG_DEBUG(categoryName, verbosity, formatStr, ...) GPES_LOG_IMPL(categoryName, verbosity, true, formatStr, __VA_ARGS__)