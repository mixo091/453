FLAGS = -c -g3 -Wall -O3
OBJS = *.o 
TARGETS = search 
CC = g++


search: a3_main.o Functionality.o curve.o point.o frechet.o config.o simplification.o interval.o
	g++ -o search a3_main.o Functionality.o curve.o point.o frechet.o config.o simplification.o interval.o -O3 -g3 -Wall

a3_main.o: ./a3_main.cpp ./headers/frechet/frechet.hpp ./headers/frechet/types.hpp ./headers/frechet/point.hpp ./headers/frechet/curve.hpp
	g++ -g3 -c -Wall ./a3_main.cpp ./src/frechet.cpp ./src/point.cpp ./src/curve.cpp -O3 

Functionality.o: ./src/Functionality.cpp ./headers/Functionality.hpp ./headers/frechet/frechet.hpp ./headers/frechet/types.hpp ./headers/frechet/point.hpp ./headers/frechet/curve.hpp
	$(CC) $(FLAGS) ./src/Functionality.cpp ./src/frechet.cpp ./src/point.cpp ./src/curve.cpp

frechet.o: ./src/frechet.cpp ./headers/frechet/frechet.hpp ./headers/frechet/curve.hpp ./headers/frechet/point.hpp ./headers/frechet/types.hpp ./headers/frechet/interval.hpp
	$(CC) $(FLAGS) ./src/frechet.cpp ./src/point.cpp ./src/interval.cpp ./src/curve.cpp

point.o: ./src/point.cpp ./headers/frechet/point.hpp ./headers/frechet/types.hpp ./headers/frechet/interval.hpp
	$(CC) $(FLAGS) ./src/point.cpp 

config.o: ./src/config.cpp ./headers/frechet/config.hpp
	$(CC) $(FLAGS) ./src/config.cpp 

interval.o: ./src/interval.cpp ./headers/frechet/interval.hpp
	$(CC) $(FLAGS) ./src/interval.cpp 

simplification.o: ./src/simplification.cpp ./headers/frechet/simplification.hpp
	$(CC) $(FLAGS) ./src/simplification.cpp 

curve.o: ./src/curve.cpp ./headers/frechet/curve.hpp ./headers/frechet/simplification.hpp
	$(CC) $(FLAGS) ./src/curve.cpp ./src/simplification.cpp

clean:
	rm -rf $(TARGETS) $(OBJS)

