#ifndef TOOLS_H_
#define TOOLS_H_

#include "global.h"
#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
using namespace std;

int gettimeofday(struct timeval *tv, void * /*tzv*/);

double get_current_time()
{
	timeval t;
	gettimeofday(&t, 0);
	return (double)t.tv_sec + (double)t.tv_usec / 1000000;
	//	return 0;
}

template <typename T>
void print_array(T* array, T len)
{
	cout << "{";
	for (int i = 0; i < len; i++)
	{
		cout << array[i];
		if (i < len - 1)
		{
			cout << ", ";
		}
	}
	cout << "}" << endl;
}

template <typename T>
void print_vertor(vector<T> vec)
{
	size_t len = vec.size();
	cout << "{";
	for (size_t i = 0; i < len; i++)
	{
		cout << vec[i] << ", ";
	}
	cout << "}" << endl;
}

int gettimeofday(struct timeval *tv, void * /*tzv*/) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  tv->tv_sec  = ts.tv_sec;
  tv->tv_usec = ts.tv_nsec / 1000;
  return 0;
}

void print_time_elapsed(std::string desc, struct timeval* start, struct
												timeval* end)
{
	struct timeval elapsed;
	if (start->tv_usec > end->tv_usec) {
		end->tv_usec += 1000000;
		end->tv_sec--;
	}
	elapsed.tv_usec = end->tv_usec - start->tv_usec;
	elapsed.tv_sec = end->tv_sec - start->tv_sec;
	float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
	std::cout << desc << "Total Time Elapsed: " << std::to_string(time_elapsed) << 
		"seconds" << std::endl;
}

/*
 * Returns the 1-based index of the last (i.e., most significant) bit set in x.
 */
inline uint64_t last_bit_set (uint64_t x) {
  assert (x > 0);
  return (sizeof (uint64_t) * 8 - __builtin_clzll (x));
}

inline uint64_t floor_lg (uint64_t x) {
  return (last_bit_set (x) - 1);
}

inline uint64_t ceil_lg (uint64_t x) {
  return (last_bit_set (x - 1));
}

/* Returns the largest power of 2 not greater than x
 * (i.e., $2^{\lfloor \lg x \rfloor}$).
 */
inline uint64_t hyperfloor (uint64_t x) {
  return (1 << floor_lg (x));
}

/* Returns the smallest power of 2 not less than x
 * (i.e., $2^{\lceil \lg x \rceil}$).
 */
inline uint64_t hyperceil (uint64_t x) {
  return (1 << ceil_lg (x));
}

inline uint64_t ceil_div (uint64_t x, uint64_t y) {
  assert (x > 0);
  return (1 + ((x - 1) / y));
}

// size_t binary_search(size_t *array, size_t start, size_t end, size_t target)
// {
// 	return 1;
// }

#endif

