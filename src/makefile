sources = ldpc.cpp decode.cpp encode.cpp itime.cpp mathfun.cpp main.cpp mtrand/mtrand.cpp mtrand/MTRand_gaussian.cpp

basic = g++ -o ldpc $(sources)

radio:
	$(basic) -D NTHREADS=2 -pthread -O3 -march=native
radio-debug:
	$(basic) -D NTHREADS=2 -pthread -g
radio-debugopt:
	$(basic) -D NTHREADS=2 -pthread -O3 -march=native -g

sprout:
	$(basic) -D NTHREADS=64 -pthreads -O3 -mcpu=ultrasparc -mtune=ultrasparc

