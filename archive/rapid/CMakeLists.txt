
include_directories(
	"C:/Projects/a/OpenCL-Headers"
#	"C:/Projects/Vulkan-Headers/include"

	"E:/C/OpenCL-Headers"
#	"E:/C/a/Vulkan-Headers/include"
)
link_directories(
	"C:/Projects/a/OpenCL-Headers"
#	"C:/Projects/Vulkan-Headers"

	"E:/C/OpenCL-Headers"
#	"E:/C/a/Vulkan-Headers"
)
add_executable(rapid
	r01.c
	rapid.c
	util.c
	
	htscodecs/arith_dynamic.c
	htscodecs/fqzcomp_qual.c
	htscodecs/htscodecs.c
	htscodecs/pack.c
	htscodecs/rANS_static.c
	htscodecs/rANS_static4x16pr.c
	htscodecs/rANS_static32x16pr.c
	htscodecs/rANS_static32x16pr_avx2.c
	htscodecs/rANS_static32x16pr_avx512.c
	htscodecs/rANS_static32x16pr_neon.c
	htscodecs/rANS_static32x16pr_sse4.c
	htscodecs/rle.c
	htscodecs/tokenise_name3.c
	htscodecs/utils.c
)
