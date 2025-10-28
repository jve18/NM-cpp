CPPC = g++
CPPCFLAGS = -c -Wall

all: test_fzero_bisect

test_fzero_bisect: main.o root.o
	$(CPPC) main.o root.o -o test_fzero_bisect

main.o: main.cpp
	$(CPPC) $(CPPCFLAGS) main.cpp

root.o: ../../../source/root.cpp
	$(CPPC) $(CPPCFLAGS) ../../../source/root.cpp

clean: rm -rf *o test_fzero_bisect

