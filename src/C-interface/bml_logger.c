#define _POSIX_C_SOURCE 200809L

#include "bml_logger.h"

#include <execinfo.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

static bml_log_level_t global_log_level = GLOBAL_DEBUG;

#define FMT_LENGTH 10000
#define TRACE_DEPTH 20
#define BUFFER_SIZE 512

#ifndef ADDR2LINE
#define ADDR2LINE "false"
#endif

/** Resolve line number in backtrace.
 */
static void
resolve_linenumber(
    void *trace,
    char *string)
{
    char buffer[BUFFER_SIZE];
    char *line;

    /* Find first occurence of '(' or ' ' in strings[i] and assume
     * everything before that is the file name. (Don't go beyond 0 though
     * (string terminator)
     */
    int p = 0;
    while (string[p] != '(' && string[p] != ' ' && string[p] != 0)
    {
        ++p;
    }

    char syscom[BUFFER_SIZE];
    snprintf(syscom, BUFFER_SIZE - 1, ADDR2LINE " %p -e %.*s", trace, p,
             string);

    /* last parameter is the file name of the symbol */
    FILE *output = popen(syscom, "r");
    if (!output)
    {
        fprintf(stderr, "error executing %s\n", syscom);
        exit(EXIT_FAILURE);
    }

    while ((line = fgets(buffer, BUFFER_SIZE, output)) != NULL)
    {
        printf("%s", line);
    }
    pclose(output);
}

/** Print out backtrace.
 */
static void
print_backtrace(
    void)
{
    void *buffer[TRACE_DEPTH];
    size_t size;
    char **strings;
    size_t i;

    size = backtrace(buffer, TRACE_DEPTH);
    strings = backtrace_symbols(buffer, size);

    printf("Obtained %zd stack frames.\n", size);

    for (i = 3; i < size; i++)
    {
#ifdef HAVE_ADDR2LINE
        resolve_linenumber(buffer[i], strings[i]);
#else
        printf("%s\n", strings[i]);
#endif
    }

    free(strings);
}

/** Log a message.
 *
 * \param log_level The log level.
 * \param format The format (as in printf()).
 * \param ap The variadic argument list.
 */
void
bml_log_real(
    bml_log_level_t log_level,
    char *format,
    va_list ap)
{
    char new_format[FMT_LENGTH];

    if (log_level >= global_log_level)
    {
        if (log_level == BML_LOG_DEBUG)
        {
            snprintf(new_format, FMT_LENGTH - 1, "[DEBUG] %s", format);
        }
        else
        {
            strncpy(new_format, format, FMT_LENGTH - 1);
        }
        vprintf(new_format, ap);

        if (log_level == BML_LOG_ERROR)
        {
            print_backtrace();
            exit(EXIT_FAILURE);
        }
    }
}

/** Log a message.
 *
 * \param log_level The log level.
 * \param format The format (as in printf()).
 */
void
bml_log(
    bml_log_level_t log_level,
    char *format,
    ...)
{
    va_list ap;

    va_start(ap, format);
    bml_log_real(log_level, format, ap);
    va_end(ap);
}

/** Log a message with location, i.e. filename and linenumber..
 *
 * \param log_level The log level.
 * \param filename The filename to log.
 * \param linenumber The linenumber.
 * \param format The format (as in printf()).
 */
void
bml_log_location(
    bml_log_level_t log_level,
    char *filename,
    int linenumber,
    char *format,
    ...)
{
    va_list ap;
    char new_format[FMT_LENGTH];

    snprintf(new_format, FMT_LENGTH - 1, "[%s:%d] %s", filename, linenumber,
             format);
    va_start(ap, format);
    bml_log_real(log_level, new_format, ap);
    va_end(ap);
}
