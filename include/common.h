#ifndef DISCRETE_ELASTIC_RODS_COMMON_H
#define DISCRETE_ELASTIC_RODS_COMMON_H

#ifndef PRINT_DEBUG
// #define PRINT_DEBUG
#endif

#ifdef PRINT_DEBUG
#define print_debug(str) std::cout << "[" << __func__ << "] " << (str) << std::endl
#else
#define print_debug(str) \
    do                   \
    {                    \
    } while (0)
#endif

#endif // DISCRETE_ELASTIC_RODS_COMMON_H

#ifndef KEEP_TURNING
    #define KEEP_TURNING
#endif

