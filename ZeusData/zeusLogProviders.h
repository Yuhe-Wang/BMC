#pragma once

// IK comment: The implementation of these logging methods does not use any serialization mechanism. 
//             They can be used from different threads and this implies relying on the clib implementation 
//             of the various printf methods to function properly. I think, we should provide a locking 
//             mechanism, so output from different threads is not mixed or even destroyed. 
//
namespace Zeus {
    /**
    * \brief Namespace containing the log providers.
    *
    * These methods currently just print to stdout.
    */
    namespace LogProviders {
        /**
        * Print informational message.
        */
        void InformationLog(const char *fmt, ...);
        /**
        * Print warning message.
        */
        void WarningLog(const char *fmt, ...);
        /**
        * Print error message.
        */
        void ErrorLog(const char *fmt, ...);
        /**
        * Print fatal error message and exit (or throw exception)
        */
        void FatalLog(const char *fmt, ...);

    } // Namespace LogProviders
}
