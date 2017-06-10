// if the macro is not defined, the tests will not be compiled

#ifndef DEF_TEST_HPP
#define DEF_TEST_HPP

#include "LargeScaleLinearProgramming.hpp"


#ifdef COMPILE_TEST

#include <cstdlib>



class Timer {
public:
	Timer(){}
	void begin_timing() {
		started = true;
		start = chrono::system_clock::now();
	}
	void stop_timing() {
		end = chrono::system_clock::now();
		if (!started) {
			cerr << "warning: stop_timing called before begin_timing";
		}
		started = false;
	}
	float get_duration() {
		auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
		return (float)duration.count() / 1000.0f;
	}
private:
	bool started = false;
	chrono::time_point<chrono::system_clock> start{};
	chrono::time_point<chrono::system_clock> end{};
};


// test basic functions of OnDiskMatrix class
void test_OnDiskMatrix();

void test_GenerateRandomMatrix();

// test writing an 3000x3000 matrix and reading each row of it
void test_OnDiskMatrix_ReadingTime();

void test_SimplexMethod();


#endif

#endif // !DEF_TEST_HPP
