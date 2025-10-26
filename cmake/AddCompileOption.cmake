# - Custom function for checking compiler flag
#
# This module provides the function enable_supported_compiler_flags(target).
# It checks if the compiler supports a list of given flags and sets them as
# compile options for the given target.
#    
#    enable_supported_compiler_flags(<target>
#                                    [CHECK_C_FLAGS] [PUBLIC]
#                                    [FLAGS <flag_1> [<flag_2> ...]]
#    )
#
# The default is to check CXX flags, if the `CHECK_C_FLAGS` option
# is given C flags are checked instead.
# If the option `PUBLIC` is given, the flags are added publicly to
# `<target>`. The default is private.
# All arguments listed after `FLAGS` will be checked and possibly
# added as a compile option.

include(CheckCompilerFlag)

function(enable_supported_compiler_flag target)
    cmake_parse_arguments(PARSE_ARGV 1 arg "CHECK_C_FLAGS;PUBLIC" "" "FLAGS")

    if(NOT DEFINED arg_FLAGS)
        message(STATUS "No flags given")
        return()
    endif()

    if(NOT TARGET "${target}")
        message(WARNING "\"${target}\" is not a target")
        return()
    endif()

    if(arg_CHECK_C_FLAGS)
        set(LANG C)
    else()
        set(LANG CXX)
    endif()

    if(arg_PUBLIC)
        set(VISIBILITY PUBLIC)
    else()
        set(VISIBILITY PRIVATE)
    endif()

    foreach(flag IN LISTS arg_FLAGS)
        get_target_property(old_flags "${target}" COMPILE_OPTIONS)
        string(FIND "${old_flags}" "${flag}" flag_already_set)

        if (NOT (flag_already_set EQUAL -1))
            message(STATUS "Flag \"${flag}\" already set for target \"${target}")
            continue()
        endif()

        check_compiler_flag(${LANG} "${flag}" flag_supported)

        if(flag_supported)
            target_compile_options("${target}" ${VISIBILITY} "${flag}")
            message(STATUS "Added flag \"${flag}\" to target \"${target}\"")
        else()
            message(STATUS "Flag \"${flag}\" not supported")
        endif()

        unset(flag_supported CACHE)
    endforeach()
endfunction()
