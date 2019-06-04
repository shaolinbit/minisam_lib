# minisam

因子图库
gcc compiler version:5.4.0 20160609 (Ubuntu 5.4.0-6ubuntu1~16.04.11)
Usage:
1.install eigen3 

2.Copy libminisam.so to library path,eg:/usr/lib

3.build the example: 
	g++ -fpic -std=c++11 -o main main.cpp -lminisam -I /usr/include/eigen3/
4.run the example:
	./main

