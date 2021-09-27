run: main.cpp
	g++ -o run main.cpp -pthread -lbcm2835 -lfftw3

clean:
	rm *.o run -rf
