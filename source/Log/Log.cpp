
#include "gpes/Log/Log.hpp"

#include <memory>
#include <unordered_map>
#include <filesystem>
#include <chrono>


#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"


using namespace gpes::log;
namespace fs = std::filesystem;

namespace {
    const std::unordered_map<LogVerbosity, spdlog::level::level_enum> c_verbosityMap = {
        {LogVerbosity::NoLogging, spdlog::level::off},
        {LogVerbosity::Display, spdlog::level::info},
        {LogVerbosity::Warning, spdlog::level::warn},
        {LogVerbosity::Error, spdlog::level::err},
        {LogVerbosity::Log, spdlog::level::info},
        {LogVerbosity::Fatal, spdlog::level::critical}
    };
    constexpr const char* c_logPattern = "[%H:%M:%S.%e] [%^%l%$] %v";
    const fs::path c_logDirectory = fs::path("logs");
    constexpr const char* c_logFilePrefix = "GPES";
    constexpr const char* c_logFileExtension = "txt";
    constexpr const char* c_timestampFormat = "{:%Y.%d.%m-%H.%M.%S}";
};

class Log::Impl
{
public:
    Impl()
    {
        using namespace spdlog;

        const auto consoleSink = std::make_shared<sinks::stdout_color_sink_mt>();
        m_consoleLogger = std::make_unique<logger>("console", consoleSink);
        m_consoleLogger->set_pattern(c_logPattern);

        setLogFile(makeLogFile());
    }

    void log(LogVerbosity verbosity, const std::string& message) const
    {
        const auto spdLevel = c_verbosityMap.at(verbosity);
        if (verbosity != LogVerbosity::Log && m_consoleLogger->should_log(spdLevel))
        {
            m_consoleLogger->log(spdLevel, message);
        }

        if (m_fileLogger->should_log(spdLevel))
        {
            m_fileLogger->log(spdLevel, message);
        }

        if (verbosity == LogVerbosity::Fatal)
        {
            // PLATFORM_BREAK();
        }
    }

    void setLogFile(const fs::path& logFile)
    {
        if (logFile.empty()) {
            return;
        }
        if (logFile.has_parent_path()) {
            fs::create_directories(logFile.parent_path());
        }
        const auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logFile.string(), true);
        m_fileLogger = std::make_unique<spdlog::logger>("file", fileSink);
        m_fileLogger->set_pattern(c_logPattern);
    }

private:
    std::unique_ptr<spdlog::logger> m_consoleLogger;
    std::unique_ptr<spdlog::logger> m_fileLogger;

    fs::path makeLogFile() const
    {
        fs::create_directory(c_logDirectory);
        const auto now = std::chrono::system_clock::now();
        const auto nowSeconds = std::chrono::floor<std::chrono::seconds>(now);
        const std::string timestamp = std::format(c_timestampFormat, nowSeconds);
        const std::string logName = std::format("{}-{}.{}", c_logFilePrefix, timestamp, c_logFileExtension);
        return c_logDirectory / logName;
    }
};

// interface

Log::Log() : pImpl_(std::make_unique<Impl>()) {}
Log::~Log() = default;

void Log::setLogFile(const fs::path& path)
{
    pImpl_->setLogFile(path);
}

void Log::log(const LogCategory& category, LogVerbosity verbosity, const std::string& message, bool showLocation,
    const std::source_location location) const
{
    const std::string fmtMsg = showLocation
                                   ? std::format("[{}] [{}:{}] {}", category.name(), location.function_name(), location.line(), message)
                                   : std::format("[{}] {}", category.name(), message);

    pImpl_->log(verbosity, fmtMsg);
}
