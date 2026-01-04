include_directories(
	"."
	"zydist_amalgamated410"
)
add_definitions(-DZYDIS_DEFINES_H=)
add_definitions(-DZYDIS_EXPORT=)
add_definitions(-DZYDIS_NO_EXPORT=)

#	add_definitions(-DPROFILER=)
#	add_definitions(-DENABLE_PROFILER_DISASSEMBLY=)

#	add_definitions(-DRELEASE=)
#	add_definitions(-DCODEC_EXT=c34)

add_executable(codec
	c01.c
	c02.c
	c03.c
	c04.c
	c05.c
	c06.c
	c07.c
	c08.c
	c09.c
	c10.c
	c11.c
	c12.c
	c13.c
	c14.c
	c15.c
	c16.c
	c17.c
	c18.c
	c19.c
	c20.c
	c21.c
	c22.c
	c23.c
	c24.c
	c25.c
	c26.c
	c27.c
	c28.c
	c29.c
	c30.c
	c31.c
	c32.c
	c33.c
	c34.c
	c35.c
	c36.c
	c37.c
	c38.c
	c39.c
	c40.c
	c41.c
	c42.c
	c43.c
	c44.c
	c45.c
	codec.c
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

	zydist_amalgamated410/Zydis.c
)
