#pragma once

#include <iostream>

/**
 * The following message logger macros are available:
 *
 * VERBOSE_MSG(message)  prints a verbose message
 * DEBUG_MSG(message)    prints a debug message
 * INFO_MSG(message)     prints information
 * WARNING_MSG(message)  prints a warning
 * ERROR_MSG(message)    prints an error
 * FATAL_MSG(message)    prints an error and throws an exception
 *
 * Verbose messages are only printed in ENABLE_VERBOSE is defined.
 * Debug messages are only printed in ENABLE_DEBUG is defined.
 * Colored output can be enabled if ENABLE_COLORS is defined.
 * All output can be disabled if ENABLE_SILENT is defined.
 */

/**
 * The following definitions may be put into a config.h file,
 * generated by cmake.
 */
#define ENABLE_COLORS
#define ENABLE_DEBUG
#define ENABLE_VERBOSE
#undef  ENABLE_SILENT

namespace himalaya {

/**
 * Enum representing the kind of output.
 */
enum ELogLevel { kVerbose, kDebug, kInfo, kWarning, kError, kFatal };

} // namespace himalaya

// All macros print to std::cerr by default in order to not interfere
// with normal Himalaya output
#define LOG_OUTPUT_STREAM std::cerr

#ifdef ENABLE_VERBOSE
   #define VERBOSE_MSG(message) LOG(himalaya::kVerbose, message)
#else
   #define VERBOSE_MSG(message)
#endif

#ifdef ENABLE_DEBUG
   #define DEBUG_MSG(message) LOG(himalaya::kDebug,   message)
#else
   #define DEBUG_MSG(message)
#endif

#define INFO_MSG(message)    LOG(himalaya::kInfo,    message)
#define WARNING_MSG(message) LOG(himalaya::kWarning, message)
#define ERROR_MSG(message)   LOG(himalaya::kError,   message)

#ifdef ENABLE_SILENT
   #define FATAL_MSG(message)                         \
      do {                                            \
         exit(EXIT_FAILURE);                          \
         throw std::runtime_error(message);           \
      } while (false)
#else
   #define FATAL_MSG(message)                                         \
      do {                                                            \
         LOG(himalaya::kFatal, message);                              \
         throw std::runtime_error(message);                           \
      } while (false)
#endif

#ifdef ENABLE_SILENT
   #define PRINT_PREFIX(level)
#else
   #define PRINT_PREFIX(level)                                        \
      do {                                                            \
         switch (level) {                                             \
         case himalaya::kVerbose: LOG_OUTPUT_STREAM << "Himalaya verbose: "; break; \
         case himalaya::kDebug:   LOG_OUTPUT_STREAM << "Himalaya debug: "; break; \
         case himalaya::kInfo:    LOG_OUTPUT_STREAM << "Himalaya info: "; break; \
         case himalaya::kWarning: LOG_OUTPUT_STREAM << "Himalaya warning: "; break; \
         case himalaya::kError:   LOG_OUTPUT_STREAM << "Himalaya error: "; break; \
         case himalaya::kFatal:   LOG_OUTPUT_STREAM << "Himalaya fatal error: "; break; \
         default:                                                     \
            break;                                                    \
         }                                                            \
      } while (false)
#endif

#ifdef ENABLE_SILENT
   #define PRINT_FILE_LINE(level)
#else
   #define PRINT_FILE_LINE(level)                                     \
      do {                                                            \
         switch (level) {                                             \
         case himalaya::kFatal:                                       \
            LOG_OUTPUT_STREAM << "(file: " << __FILE__                \
                              << ", line: " << __LINE__ << ") ";      \
            break;                                                    \
         default:                                                     \
            break;                                                    \
         }                                                            \
      } while (false)
#endif

#ifdef ENABLE_SILENT
   #define PRINT_COLOR_CODE(level)
#else
   #define PRINT_COLOR_CODE(level)                                  \
      do {                                                          \
         switch (level) {                                           \
         case himalaya::kVerbose: LOG_OUTPUT_STREAM << "\033[0;36m"; break; \
         case himalaya::kDebug:   LOG_OUTPUT_STREAM << "\033[0;34m"; break; \
         case himalaya::kInfo:    LOG_OUTPUT_STREAM << "\033[1;34m"; break; \
         case himalaya::kWarning: LOG_OUTPUT_STREAM << "\033[0;31m"; break; \
         case himalaya::kError:   LOG_OUTPUT_STREAM << "\033[1;31m"; break; \
         case himalaya::kFatal:   LOG_OUTPUT_STREAM << "\033[41;1;37m"; break; \
         default:                                                   \
            break;                                                  \
         }                                                          \
      } while (false)
#endif

#ifdef ENABLE_SILENT
   #define RESET_COLOR(level)
#else
   #define RESET_COLOR(level)                                    \
      do {                                                       \
         LOG_OUTPUT_STREAM << "\033[0m";                         \
      } while (false)
#endif

#ifdef ENABLE_SILENT
   #define PRINT_MESSAGE(level, message)
#else
   #define PRINT_MESSAGE(level, message)                         \
      do {                                                       \
         LOG_OUTPUT_STREAM << message;                           \
      } while (false)
#endif

#ifdef ENABLE_SILENT
   #define PRINT_ENDL(level)
#else
   #define PRINT_ENDL(level)                                     \
      do {                                                       \
         LOG_OUTPUT_STREAM << std::endl;                         \
      } while (false)
#endif

#ifdef ENABLE_SILENT
   #define LOG(level, message)
#else
   #ifdef ENABLE_COLORS
      #define LOG(level, message)                                \
      do {                                                       \
         PRINT_COLOR_CODE(level);                                \
         PRINT_PREFIX(level);                                    \
         PRINT_FILE_LINE(level);                                 \
         RESET_COLOR(level);                                     \
         PRINT_MESSAGE(level, message);                          \
         PRINT_ENDL(level);                                      \
      } while (false)
   #else
      #define LOG(level, message)                  \
      do {                                         \
         PRINT_PREFIX(level);                      \
         PRINT_FILE_LINE(level);                   \
         PRINT_MESSAGE(level, message);            \
         PRINT_ENDL(level);                        \
      } while (false)
   #endif
#endif
