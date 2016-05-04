CFLAGS += -I include --std=c++14 -Wall -Wextra -Werror -O3

test: include/kdbush.hpp test.cpp Makefile
	$(CXX) test.cpp $(CFLAGS) -o test
	./test

bench: include/kdbush.hpp bench.cpp Makefile
	$(CXX) bench.cpp $(CFLAGS) -o bench
	./bench

format:
	clang-format include/*.hpp *.cpp -i

clean:
	rm test

default: test
