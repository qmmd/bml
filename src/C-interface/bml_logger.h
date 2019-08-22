/** \file */

#ifndef __BML_LOGGER_H
#define __BML_LOGGER_H

#include "bml_types.h"

#include <stdlib.h>

/** The log-levels. */
typedef enum
{
    /** Debugging messages. */
    BML_LOG_DEBUG,
    /** Info messages. */
    BML_LOG_INFO,
    /** Warning messages. */
    BML_LOG_WARNING,
    /** Error messages. */
    BML_LOG_ERROR
} bml_log_level_t;

void bml_log(
    bml_log_level_t log_level,
    char *format,
    ...);

void bml_log_location(
    bml_log_level_t log_level,
    char *filename,
    int linenumber,
    char *format,
    ...);

/** Convenience macro to write a BML_LOG_DEBUG level message. */
#define LOG_DEBUG(format, ...) \
    bml_log_location(BML_LOG_DEBUG, __FILE__, __LINE__, format, ##__VA_ARGS__)
/** Convenience macro to write a BML_LOG_INFO level message. */
#define LOG_INFO(format, ...) \
    bml_log(BML_LOG_INFO, format, ##__VA_ARGS__)
/** Convenience macro to write a BML_LOG_WARNING level message. */
#define LOG_WARN(format, ...) \
    bml_log_location(BML_LOG_WARNING, __FILE__, __LINE__, format, ##__VA_ARGS__)
/** Convenience macro to write a BML_LOG_ERROR level message. */
#define LOG_ERROR(format, ...) \
    bml_log_location(BML_LOG_ERROR, __FILE__, __LINE__, format, ##__VA_ARGS__)

#endif
