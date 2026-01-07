
ifeq ($(OS), Windows_NT)
LIBS += -static
else
LIBS += -lm
endif

hist:
	gcc -Wall -O3 -s hist.c -o hist $(LIBS)

psnr:
	gcc -Wall -O3 -s psnr.c -o psnr $(LIBS)

ssim:
	gcc -Wall -O3 -s ssim.c -o ssim $(LIBS)

wavpsnr:
	gcc -Wall -O3 -s wavpsnr.c -o wavpsnr $(LIBS)

diff:
	gcc -Wall -O3 -s diff.c -o diff $(LIBS)

chism:
	gcc -Wall -O3 -mavx2 -mbmi -mlzcnt -s chism.c -o chism $(LIBS)
