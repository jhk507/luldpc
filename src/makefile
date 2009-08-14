sources = ldpc.cpp decode.cpp encode.cpp itime.cpp mathfun.cpp main.cpp mtrand/mtrand.cpp mtrand/MTRand_gaussian.cpp

basic = g++ -O3 -pthread $(sources)

radio:
	$(basic) -march=core2

sprout:
	$(basic) -mcpu=ultrasparc -mtune=ultrasparc
