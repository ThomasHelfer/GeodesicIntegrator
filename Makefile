CC = g++

CFLAGS = -O2 -Wall 

output: integration.cpp
	$(CC) $(CFLAGS) integration.cpp 


