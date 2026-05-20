#include "TestRunner.h"
#include "../../../Logger/Logger.h"
#include <spdlog/spdlog.h>
#include <clocale>

int main() 
{
    setlocale(LC_ALL, "Russian"); 
    Logger::init("../../../log/app.log");
    
    LOG_INFO("AGO Application - Starting");
    
    RunConsoleBenchmarks();

    
    LOG_INFO("AGO Application - Finished");
    spdlog::shutdown();
    return 0;
}