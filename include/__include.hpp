#ifndef DEF___INCLUDE_HPP
#define DEF___INCLUDE_HPP

// Dependent on the Eigen library
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <string>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <memory>
#include <limits>
#include <map>

using namespace std;


inline void expr_check(bool expr, const char* info) {
	if (!expr) {
		cerr << info << "\n";
		throw runtime_error{ info };
	}
}

// error tolerance
const double mach_eps = 1.0e-8;


template<typename V>
inline bool fpeq(V left, V right) {
	if (abs(right - left) < (V)mach_eps)return true;
	else return false;
}

#define COMPILE_TEST


#endif // !DEF___INCLUDE_HPP


