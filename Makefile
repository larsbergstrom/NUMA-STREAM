default:
	cc -std=c99 -O3 -DNUM_THREADS=64 -DN=5000000 -o stream-pthreads stream-pthreads.c -lnuma -pthread
