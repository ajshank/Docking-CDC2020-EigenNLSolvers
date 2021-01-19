EIGENDIR=$(shell pkg-config --cflags eigen3)

cdctest: cdc-test.cpp
	g++ cdc-test.cpp $(EIGENDIR) -O3 -ffast-math -o cdctest
