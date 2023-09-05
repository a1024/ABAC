# CMakeList.txt : CMake project for lossy, include source and define project specific logic here.

include_directories(
	"C:/Projects/a/libjpeg-turbo"
)
link_directories(
	"C:/Projects/a/libjpeg-turbo/build/Release"
)
add_executable(lossy
	"lodepng.c"
	"lossy.c"
	"tests.c"
	"transforms.c"
	"util.c"
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET lossy PROPERTY CXX_STANDARD 20)
endif()
