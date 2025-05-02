add_executable(alic
	a01.c
	a02.c
	a03.c
	a04.c
	a05.c

	fpc/fpc.c

	fse/debug.c
	fse/entropy_common.c
	fse/fse_compress.c
	fse/fse_decompress.c
	fse/fseU16.c
	fse/hist.c
	fse/huf_compress.c
	fse/huf_decompress.c

	sli.c
)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-DRELWITHDEBINFO")
