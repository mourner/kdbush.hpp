CFLAGS += -I include --std=c++14 -Wall -Wextra -Werror

test: include/kdbush.hpp test.cpp Makefile
	$(CXX) test.cpp $(CFLAGS) -o test
	./test

format: include/kdbush.hpp test.cpp
	clang-format include/kdbush.hpp test.cpp -i

clean:
	rm test

default: test
