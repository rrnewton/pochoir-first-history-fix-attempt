/*
 **********************************************************************************
 *  Copyright (C) 2010  Massachusetts Institute of Technology
 *  Copyright (C) 2010  Yuan Tang <yuantang@csail.mit.edu>
 * 		                Charles E. Leiserson <cel@mit.edu>
 * 	 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

/* Test bench - 2D heat equation, Periodic version */
#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>

#include <pochoir.hpp>

using namespace std;
#define SIMPLE 0
/* N_RANK includes both time and space dimensions */
#define N_RANK 2
// #define N_SIZE 555 
// #define T_SIZE 555
#define TOLERANCE (1e-6)

void check_result(int t, int j, int i, double a, double b)
{
	if (abs(a - b) < TOLERANCE) {
//		printf("a(%d, %d, %d) == b(%d, %d, %d) == %f : passed!\n", t, j, i, t, j, i, a);
	} else {
		printf("a(%d, %d, %d) = %f, b(%d, %d, %d) = %f : FAILED!\n", t, j, i, a, t, j, i, b);
	}

}

    Pochoir_Boundary_2D(heat_bv_2D, arr, t, i, j)
        /* this is non-periodic boundary value */
        /* we already shrinked by using range I, J, K,
         * so the following code to set boundary index and
         * boundary rvalue is not necessary!!! 
         */
        int new_i = i, new_j = j;
        if (new_i < 0) 
            new_i += arr.size(1);
        else if (new_i >= arr.size(1))
            new_i -= arr.size(1);

        if (new_j < 0) 
            new_j += arr.size(0);
        else if (new_j >= arr.size(0))
            new_j -= arr.size(0);
        return arr.get(t, new_i, new_j);
    Pochoir_Boundary_end

int main(int argc, char * argv[])
{
	const int BASE = 1024;
	int t;
	struct timeval start, end;
    int N_SIZE = 0, T_SIZE = 0;

    if (argc < 3) {
        printf("argc < 3, quit! \n");
        exit(1);
    }
    N_SIZE = StrToInt(argv[1]);
    T_SIZE = StrToInt(argv[2]);
    printf("N_SIZE = %d, T_SIZE = %d\n", N_SIZE, T_SIZE);
	/* data structure of Pochoir - row major */
	Pochoir_Array<double, N_RANK> a(N_SIZE, N_SIZE), b(N_SIZE, N_SIZE);
    Pochoir_Stencil<double, N_RANK> heat_2D;
    Pochoir_uRange I(0, N_SIZE-1), J(0, N_SIZE-1);
    Pochoir_Shape_info<2> heat_shape_2D[5] = {{1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, -1}, {0, 0, 1}};

	for (int i = 0; i < N_SIZE; ++i) {
	for (int j = 0; j < N_SIZE; ++j) {
        if (i == 0 || i == N_SIZE-1
            || j == 0 || j == N_SIZE-1) {
            a(0, i, j) = a(1, i, j) = 0;
        } else {
#if DEBUG 
		    a(0, i, j) = i * N_SIZE + j;
		    a(1, i, j) = 0;
#else
            a(0, i, j) = 1.0 * (rand() % BASE); 
            a(1, i, j) = 0; 
#endif
        }
        b(0, i, j) = a(0, i, j);
        b(1, i, j) = 0;
	} }

	cout << "a(T+1, J, I) = 0.125 * (a(T, J+1, I) - 2.0 * a(T, J, I) + a(T, J-1, I)) + 0.125 * (a(T, J, I+1) - 2.0 * a(T, J, I) + a(T, J, I-1)) + a(T, J, I)" << endl;
    Pochoir_kernel_2D(heat_2D_fn, t, i, j)
	    a(t+1, i, j) = 0.125 * (a(t, i+1, j) - 2.0 * a(t, i, j) + a(t, i-1, j)) + 0.125 * (a(t, i, j+1) - 2.0 * a(t, i, j) + a(t, i, j-1)) + a(t, i, j);
    Pochoir_kernel_end

    /* we have to bind arrayInUse and Shape together 
     * => One arrayInUse, one shape[] => One slope[]
     * because each arrayInUse needs to know the slope to determine
     * the boundary region and when to call the user supplied boundary
     * value function
     */
    heat_2D.registerBoundaryFn(a, heat_bv_2D);
    heat_2D.registerArrayInUse(a);
    heat_2D.registerShape(heat_shape_2D);
    heat_2D.registerDomain(I, J);

	gettimeofday(&start, 0);
    heat_2D.run(T_SIZE, heat_2D_fn);
	gettimeofday(&end, 0);
	std::cout << "Pochoir ET: consumed time :" << 1.0e3 * tdiff(&end, &start) << "ms" << std::endl;

    b.registerShape(heat_shape_2D);
    b.registerBV(heat_bv_2D);

	gettimeofday(&start, 0);
    /* cilk_for + zero-padding */
	for (int t = 0; t < T_SIZE; ++t) {
    cilk_for (int i = 0; i <= N_SIZE-1; ++i) {
	for (int j = 0; j <= N_SIZE-1; ++j) {
        b(t+1, i, j) = 0.125 * (b(t, i+1, j) - 2.0 * b(t, i, j) + b(t, i-1, j)) + 0.125 * (b(t, i, j+1) - 2.0 * b(t, i, j) + b(t, i, j-1)) + b(t, i, j); } } }
	gettimeofday(&end, 0);
	std::cout << "Naive Loop: consumed time :" << 1.0e3 * tdiff(&end, &start) << "ms" << std::endl;

	t = T_SIZE;
	for (int i = 0; i <= N_SIZE-1; ++i) {
	for (int j = 0; j <= N_SIZE-1; ++j) {
		check_result(t, i, j, a.interior(t, i, j), b.interior(t, i, j));
	} } 

	return 0;
}
