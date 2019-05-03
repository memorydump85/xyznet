CFLAGS = -std=c++17 -O3 -I ../inc/ -Wall -Wextra -Werror

INCLUDES = \
	-I ../deps/cppkit/inc \
	`pkg-config eigen3 --cflags`

.PHONY: all clean

all: run_simulation.x simulator_test.x

run_simulation.x: simulator.hh run_simulation.cc
	c++ $(CFLAGS) $(INCLUDES) -o $@ run_simulation.cc

simulator_test.x: simulator.hh simulator_test.cc
	c++ $(CFLAGS) $(INCLUDES) -o $@ simulator_test.cc

clean: *
	rm -f simulator_test.x run_simulation.x
