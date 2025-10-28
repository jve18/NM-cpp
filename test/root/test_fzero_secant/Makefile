CPPC = g++
CPPCFLAGS = -c -Wall

all: test_fzero_secant

test_fzero_secant: main.o root.o
	$(CPPC) main.o root.o -o test_fzero_secant

main.o: main.cpp
	$(CPPC) $(CPPCFLAGS) main.cpp

root.o: ../../../source/root.cpp
	$(CPPC) $(CPPCFLAGS) ../../../source/root.cpp

clean: rm -rf *o test_fzero_secant

