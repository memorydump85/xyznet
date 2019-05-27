CFLAGS = -std=c++17 -O3 -I ../inc/      \
         -Wall -Wextra -Werror          \
         -DEIGEN_NO_DEBUG

INCLUDES =                  \
    -I ../deps/cppkit/inc   \
    `pkg-config eigen3 --cflags`

.PHONY: all clean run_tests

all: run_simulation.x simulator_test.x run_tests

simulator.o: simulator.hh simulator.cc
	c++ $(CFLAGS) $(INCLUDES) -c -o $@ simulator.cc

run_simulation.x: simulator.hh simulator.o run_simulation.cc
	c++ $(CFLAGS) $(INCLUDES) -o $@ run_simulation.cc simulator.o

simulator_test.x: simulator.hh simulator.o simulator_test.cc
	c++ $(CFLAGS) $(INCLUDES) -UEIGEN_NO_DEBUG -o $@ simulator_test.cc simulator.o

run_tests: simulator_test.x
	@echo
	@./simulator_test.x
	@echo

clean:
	rm -f *.x *.o
