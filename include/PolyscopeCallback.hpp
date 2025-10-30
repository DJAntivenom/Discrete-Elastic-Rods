/**
 * \file PolyscopeCallback.hpp
 * \brief Defines the callback passed to polyscope.
 */

#ifndef __POLYSCOPE_CALLBACK_HPP__
#define __POLYSCOPE_CALLBACK_HPP__

/**
 * @brief Override the cout rdbuf to see log in GUI
 * 
 * Once this function is called, std::cout prints messages not only to stdout,
 * but also to the GUI in the Analysis window.
 */
void overrideCoutBuffer();

/**
 * @brief Called each frame by polyscope.
 *
 * Responsible for generating the GUI elements and rendering
 * the simulated rods.
 */
void polyscopeCallback();

#endif
