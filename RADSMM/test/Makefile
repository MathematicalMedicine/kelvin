CFLAGS = -Wall -L.. -I..
all:	tst1 tst2 tst3 tst4 tst1f tst2f

tst1:	tst1.c ../libRADSMM.a
	cc $(CFLAGS) tst1.c  -o tst1 -lRADSMM -lm

tst2:	tst2.c ../libRADSMM.a
	cc $(CFLAGS) tst2.c -o tst2 -lRADSMM -lm

tst2f:		tst2f.c ../libRADSMM.a
	cc tst2f.c -O -o tst2f -I.. -L.. -lRADSMM -lm

tst3:	tst3.o ../libRADSMM.a
	cc tst3.o -o tst3 -L.. -lRADSMM -lm

tst3.o:	tst3.c ../RADSMM.h
	cc -O4 -c -I.. -otst3.o tst3.c


tst4:	tst4.o ../libRADSMM.a
	cc tst4.o -o tst4 -I.. -L.. -lRADSMM -lm

tst4.o:	tst4.c ../RADSMM.h
	cc -O4 -c -I.. -otst4.o tst4.c


tst1f:	tst1f.c ../libRADSMM.a
	cc tst1f.c -O -o tst1f -I.. -L.. -lRADSMM -lm



