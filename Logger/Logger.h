#pragma once

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE 

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <iostream>

#define LOG_INFO(...)  SPDLOG_INFO(__VA_ARGS__)
#define LOG_DEBUG(...) SPDLOG_DEBUG(__VA_ARGS__)
#define LOG_ERROR(...) SPDLOG_ERROR(__VA_ARGS__)

class Logger 
{
public:
    static void init(const std::string& log_file_path = "log/app.log") 
    {
        try 
        {
            auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
            console_sink->set_level(spdlog::level::info);

            auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file_path, true);
            file_sink->set_level(spdlog::level::trace);

            std::vector<spdlog::sink_ptr> sinks {console_sink, file_sink};
            auto logger = std::make_shared<spdlog::logger>("multi_sink", sinks.begin(), sinks.end());

            logger->set_pattern("[%Y-%m-%d %H:%M:%S:%e.%f UTC%z] [%^%l%$] [%s:%#] %v");
            logger->set_level(spdlog::level::trace);
            spdlog::flush_every(std::chrono::seconds(1));
            spdlog::set_default_logger(logger);
            spdlog::flush_on(spdlog::level::err); 
        }
        catch (const spdlog::spdlog_ex& ex) 
        {
            std::cerr << "Log initialization failed: " << ex.what() << std::endl;
        }
    }
};