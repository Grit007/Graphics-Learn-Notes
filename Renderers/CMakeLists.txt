cmake_minimum_required(VERSION 3.24)
project(Renderers)

set(CMAKE_CXX_STANDARD 11)

add_executable(Renderers main.cpp core/camera.cpp platform/win32.cpp utils/maths.cpp utils/model.cpp utils/tgaimage.cpp)
