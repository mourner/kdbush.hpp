CFLAGS += -I include --std=c++14 -Wall -Wextra -Werror -O3
DEPS = Makefile include/kdbush.hpp

build/test: test.cpp $(DEPS)
	mkdir -p build
	$(CXX) test.cpp $(CFLAGS) -o build/test

build/bench: bench.cpp $(DEPS)
	mkdir -p build
	$(CXX) bench.cpp $(CFLAGS) -o build/bench

run-bench: build/bench
	./build/bench

run-test: build/test
	./build/test

format:
	clang-format include/*.hpp *.cpp -i

clean:
	rm -rf build

default: run-test
