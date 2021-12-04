
DEBUG=-g -pg

CC=gcc

MEMCHECK=-fsanitize=address,undefined -fno-omit-frame-pointer

CFLAGS=-O3 -Wall -I./ $(DEBUG) $(MEMCHECK)

LIBS=-lm -lpthread

all : libpso.a test-vector test-fill_initialrnd test-mfunc test-workunit test-workitems # test-robotarm

pso.o : pso.h

pso_util.o : pso_util.h

test-mfunc.o : pso.h pso_util.h

libpso.a : CFLAGS+=-fPIC

libpso.a : pso_util.o update_range.o workitems.o pso.o
	ar r $@ $^

libpso.so : CFLAGS+=-fPIC

libpso.so : pso_util.o update_range.o workitems.o pso.o
	$(CC) -shared -o $@ $^

test-fill_initialrnd : LIBS+=-lpso -lm

test-fill_initialrnd : LDFLAGS+=-L./

test-fill_initialrnd : libpso.a test-fill_initialrnd.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS) $(LIBS)

test-mfunc : LDFLAGS+=-L./

test-mfunc : LIBS+=-lpso -lpthread

test-mfunc : libpso.a sincfunc.o mfuncs.o test-mfunc.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS) $(LIBS)

mfunc.asc : NUM_PARTICLES=2500000
mfunc.asc : NUM_GENERATIONS=15000
mfunc.asc : NUM_THREADS=8

mfunc.asc : test-mfunc
	@./$^ $(NUM_PARTICLES) $(NUM_GENERATIONS) $(NUM_THREADS) $@

pyramid.asc : NUM_PARTICLES=16000
pyramid.asc : NUM_GENERATIONS=7500
pyramid.asc : NUM_THREADS=8

pyramid.asc : test-mfunc
	@./$^ $(NUM_PARTICLES) $(NUM_GENERATIONS) $(NUM_THREADS) $@ "pyramid"

bumps.asc : NUM_PARTICLES=128000
bumps.asc : NUM_GENERATIONS=2000
bumps.asc : NUM_THREADS=8

bumps.asc : test-mfunc
	@./$^ $(NUM_PARTICLES) $(NUM_GENERATIONS) $(NUM_THREADS) $@ "bumps"

sinc.asc : NUM_PARTICLES=16000
sinc.asc : NUM_GENERATIONS=7500
sinc.asc : NUM_THREADS=8

sinc.asc : test-mfunc
	@./$^ $(NUM_PARTICLES) $(NUM_GENERATIONS) $(NUM_THREADS) $@ "sinc"

quicktest.asc : NUM_PARTICLES=8000
quicktest.asc : NUM_GENERATIONS=5
quicktest.asc : NUM_THREADS=3

quicktest.asc : test-mfunc
	@DEBUG=1 ./$^ $(NUM_PARTICLES) $(NUM_GENERATIONS) $(NUM_THREADS) $@ "bumps"

torus.asc : NUM_PARTICLES=50000
torus.asc : NUM_GENERATIONS=1750
torus.asc : NUM_THREADS=8

torus.asc : test-mfunc
	@./$^ $(NUM_PARTICLES) $(NUM_GENERATIONS) $(NUM_THREADS) $@ "torus"

tube.asc : NUM_PARTICLES=64000
tube.asc : NUM_GENERATIONS=1500
tube.asc : NUM_THREADS=8

tube.asc : test-mfunc
	@./$^ $(NUM_PARTICLES) $(NUM_GENERATIONS) $(NUM_THREADS) $@ "tube"

stairs.asc : NUM_PARTICLES=128000
stairs.asc : NUM_GENERATIONS=1500
stairs.asc : NUM_THREADS=8

stairs.asc : test-mfunc
	@./$^ $(NUM_PARTICLES) $(NUM_GENERATIONS) $(NUM_THREADS) $@ "stairs"

test-workitems.o : workitems.h pso.h

workitems.o : workitems.h

test-workitems : workitems.o test-workitems.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS) $(LIBS)

test-robotarm : dot.o ga.o sincfunc.o particlepack.o pso.o pswarm.o
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS) $(LIBS)

clean :
	rm libpso.a libpso.so *.o test-mfunc test-robotarm
