SRC = c03.c codec.c util.c
#SRC = c07.c codec.c util.c
#SRC = c05.c codec.c util.c
#SRC = c01.c c02.c c03.c c04.c c05.c c06.c c07.c codec.c util.c

ifeq ($(OS),Windows_NT)
LIBS += -static
else
LIBS += -lm
endif

.PHONY: c01 c02 c03 c03-clang c05 c07 c09 c10 c11 c12 c13 c13-debug c13-profile

build:
	gcc -O3 -Wall -mavx2 -mbmi -mlzcnt $(SRC) -o codec -static

build-clang:
	clang -O3 -Wall -mavx2 -mbmi -mlzcnt $(SRC) -o codec -static

c01:
	gcc -O3 -Wall -DCODEC_FUNC=c01_codec -mavx2 -mbmi -mlzcnt  c01.c codec.c util.c  -o c01 $(LIBS)

c02:
	gcc -O3 -Wall -DCODEC_FUNC=c02_codec -mavx2 -mbmi -mlzcnt  c02.c codec.c util.c  -o c02 $(LIBS)

c03:
	gcc -O3 -Wall -DCODEC_FUNC=c03_codec -mavx2 -mbmi -mlzcnt  c03.c codec.c util.c  -o c03 $(LIBS)

c03-clang:
	clang -O3 -Wall -DCODEC_FUNC=c03_codec -mavx2 -mbmi -mlzcnt  c03.c codec.c util.c  -o c03-clang $(LIBS)

c05:
	gcc -O3 -Wall -DCODEC_FUNC=c05_codec -mavx2 -mbmi -mlzcnt  c05.c codec.c util.c  -o c05 $(LIBS)

c07:
	gcc -O3 -Wall -DCODEC_FUNC=c07_codec -mavx2 -mbmi -mlzcnt  c07.c codec.c util.c  -o c07 $(LIBS)

c09:
	gcc -O3 -Wall -DCODEC_FUNC=c09_codec -mavx2 -mbmi -mlzcnt  c09.c codec.c util.c  -o c09 $(LIBS)

c10:
	gcc -O3 -Wall -DCODEC_FUNC=c10_codec -mavx2 -mbmi -mlzcnt  c10.c codec.c util.c  -o c10 $(LIBS)

c11:
	gcc -O3 -Wall -DCODEC_FUNC=c11_codec -mavx2 -mbmi -mlzcnt  c11.c codec.c util.c  -o c11 $(LIBS)

c12:
	gcc -O3 -Wall -DCODEC_FUNC=c12_codec -mavx2 -mbmi -mlzcnt  c12.c codec.c util.c  -o c12 $(LIBS)

c13:
	gcc -O3 -Wall -DCODEC_FUNC=c13_codec -mavx2 -mbmi -mlzcnt  c13.c codec.c util.c  -o c13 $(LIBS)

c13-debug:
	gcc -g -Wall -DCODEC_FUNC=c13_codec -mavx2 -mbmi -mlzcnt  c13.c codec.c util.c  -o c13-debug $(LIBS)

c13-profile:
	gcc -O -g -Wall -DCODEC_FUNC=c13_codec -mavx2 -mbmi -mlzcnt  c13.c codec.c util.c  -o c13-profile $(LIBS)

#c13-gprofile:
#	gcc -g -Wall -DCODEC_FUNC=c13_codec -mavx2 -mbmi -mlzcnt  c13.c codec.c util.c  -o c13-gprofile $(LIBS) -Wl,--no-as-needed,-lprofiler,--as-needed

c14:
	gcc -O3 -Wall -DCODEC_FUNC=c14_codec -mavx2 -mbmi -mlzcnt  c14.c codec.c util.c  -o c14 $(LIBS)

c15:
	gcc -O3 -Wall -DCODEC_FUNC=c15_codec -mavx2 -mbmi -mlzcnt  c15.c codec.c util.c  -o c15 $(LIBS)

c16:
	gcc -O3 -Wall -DCODEC_FUNC=c16_codec -mavx2 -mbmi -mlzcnt  c16.c codec.c util.c  -o c16 $(LIBS)

c17:
	gcc -O3 -Wall -DCODEC_FUNC=c17_codec -mavx2 -mbmi -mlzcnt  c17.c codec.c util.c  -o c17 $(LIBS)

c18:
	gcc -O3 -Wall -DCODEC_FUNC=c18_codec -mavx2 -mbmi -mlzcnt  c18.c codec.c util.c  -o c18 $(LIBS)
