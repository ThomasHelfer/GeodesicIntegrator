CC = g++
CFLAGS = -O2 -Wall -std=c++14 

output: integration.cpp
	$(CC) $(CFLAGS) integration.cpp -o output

clean: 
	rm output
