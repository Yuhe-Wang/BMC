#include "zeusConfig.h"
#include <cstdarg>
#include <cstdio>
#include <string>


using namespace std;

namespace Zeus {
    namespace LogProviders {

        void InformationLog(const char *fmt, ...)
        {
            char logBuffer[4096];
            va_list params;
            va_start(params, fmt);
            int end = vsprintf(logBuffer, fmt, params);
            if (end > 0)
            {
                end--;
                if(logBuffer[end] == '\n')
                {
                    logBuffer[end] = '\0';
                }
            }
            else
            {
                sprintf(logBuffer,  "InformationLog: Invalid Information log parameters");
            }
            printf( "%s", logBuffer );
            printf( "\n" );
            fflush(stdout);
            va_end(params);
        }

        void WarningLog(const char *fmt, ...)
        {
            char logBuffer[4096];
            va_list params;
            va_start(params, fmt);
            int end = vsprintf(logBuffer, fmt, params);
            if (end > 0)
            {
                end--;
                if(logBuffer[end] == '\n')
                {
                    logBuffer[end] = '\0';
                }
            }
            else
            {
                sprintf(logBuffer, "WarningLog: Invalid Warning log parameters");
            }
            printf( "%s", logBuffer );
            printf( "\n" );
            va_end(params);
        }

        void ErrorLog(const char *fmt, ...)
        {
            char logBuffer[4096];
            va_list params;
            va_start(params, fmt);
            int end = vsprintf(logBuffer, fmt, params);
            if (end > 0)
            {
                end--;
                if(logBuffer[end] == '\n')
                {
                    logBuffer[end] = '\0';
                }
            }
            else
            {
                sprintf(logBuffer, "ErrorLog: Invalid Error log parameters");
            }
            printf( "%s", logBuffer );
            printf( "\n" );
            va_end(params);
        }

        void FatalLog(const char *fmt, ...)
        {
            char logBuffer[4096];
            va_list params;
            va_start(params, fmt);
            int end = vsprintf(logBuffer, fmt, params);
            if (end > 0)
            {
                end--;		    
                if(logBuffer[end] == '\n')
                {
                    logBuffer[end] = '\0';
                }
            }
            else
            {
                sprintf(logBuffer, "FatalLog: Invalid Fatal log parameters");
            }
            printf( "%s", logBuffer );
            printf( "\n" );
            va_end(params);
        }

    } // Namespace LogProviders
} // Namespace Zeus

#ifdef AE_SIT_UNITTEST
#include "WinUnit.h"

#include "zeusLogProviders.h"
using namespace Zeus::LogProviders;
BEGIN_TEST( Logging )
{
    InformationLog( "Information logging test" );
    WarningLog( "Warning logging test" );
    FatalLog( "Fatal logging test" );
}
END_TEST
#endif
