sources = ldpc.cpp decode.cpp encode.cpp itime.cpp mathfun.cpp main.cpp mtrand/mtrand.cpp mtrand/MTRand_gaussian.cpp

basic = g++ -o ldpc -O3 $(sources)

radio:
	$(basic) -D NTHREADS=2 -pthread -march=core2

sprout:
	$(basic) -D NTHREADS=64 -mcpu=ultrasparc -mtune=ultrasparc

