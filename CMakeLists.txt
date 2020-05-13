# Minimum CMake version required
cmake_minimum_required(VERSION 3.7 FATAL_ERROR)

set (ProjectName MeshlessAnimationFramework)

# Project name
project(${ProjectName} LANGUAGES CXX)

#Specify compiler features
set(CMAKE_CXX_STANDARD 17)

# Compiler options
#if(MSVC)
#    message("Compiler: ${CMAKE_CXX_COMPILER_ID}")
#elseif (CMAKE_CXX_COMPILER_ID MATCHES GNU)
#    message("Compiler: ${CMAKE_CXX_COMPILER_ID}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Werror")
#elseif (CMAKE_CXX_COMPILER_ID MATCHES Clang)
#  message("Compiler: ${CMAKE_CXX_COMPILER_ID}")

#  # Default compile flags (read: turn on all checks)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-documentation -pedantic -Werror")

#  # Extra compile flags to turn of stuff we don't want to handle at the moment
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-exit-time-destructors")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-conversion")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-padded")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-weak-vtables")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-global-constructors")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-double-promotion")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-inconsistent-missing-override")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-float-equal")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-qualifiers")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-non-virtual-dtor")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-covered-switch-default")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-gnu-statement-expression")
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-conditional-uninitialized")
#else()
#  message(WARNING "Compiler <${CMAKE_CXX_COMPILER_ID}> is not a supported C++ compiler.")
#endif()

if(MSVC AND NOT VCPKG_TOOLCHAIN)
  message(FATAL_ERROR
    "Enable VCPKG toolchain on MSVC platform\n"
    "  Star over; CMAKE_TOOLCHAIN_FILE definition missing!\n"
    "  i.e.: cmake -DCMAKE_TOOLCHAIN_FILE=<path_to>\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake <path_to>\\source\\"
    )
endif()

# Options

#files
file(GLOB_RECURSE HEADER_FILES
    ${CMAKE_SOURCE_DIR}/src/*.h
    ${CMAKE_SOURCE_DIR}/src/*.hpp)

file(GLOB_RECURSE SOURCE_FILES
    ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Must be changed to VK_USE_PLATFORM_WAYLAND_KHR for linux I think.
add_definitions(-DVK_USE_PLATFORM_WIN32_KHR)
add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX)

# Required packages
find_package(Vulkan REQUIRED)
find_package(glm REQUIRED)
find_package(OpenMesh REQUIRED)

add_executable(${ProjectName} ${SOURCE_FILES})# ${HEADER_FILES})

target_include_directories(${ProjectName}
    PRIVATE Vulkan::Vulkan glm OpenMeshCore ${ImGui_DIR})

set(LIBS
    Vulkan::Vulkan
    glm
    OpenMeshCore
)

target_link_libraries(${ProjectName} ${LIBS})
