FLAGS = -c -g3 -Wall -O3
OBJS = *.o 
TARGETS = search 
CC = g++

search: search_main.o Functionality.o 
	g++ -Wall -o search search_main.o Functionality.o  -O3

cluster: clustering_main.o Functionality.o 
	g++ -Wall -o cluster clustering_main.o Functionality.o  -O3

search_main.o: ./search_main.cpp 
	$(CC) $(FLAGS) ./search_main.cpp

clustering_main.o: ./clustering_main.cpp 
	$(CC) $(FLAGS) ./clustering_main.cpp

Functionality.o: ./src/Functionality.cpp ./headers/Functionality.hpp
	$(CC) $(FLAGS) ./src/Functionality.cpp

clean:
	rm -rf $(TARGETS) $(OBJS)