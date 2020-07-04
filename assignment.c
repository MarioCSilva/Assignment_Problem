////////////////////////////////////////////////////////////////////////////////////////////////////
//
// AED, 2019/2020
//
// Leandro Emanuel Soares Alves da Silva 93446
// Mário Francisco Silva 93430
//
// Brute-force solution of the assignment problem (https://en.wikipedia.org/wiki/Assignment_problem)
//
// Compile with "cc -Wall -O2 assignment.c -lm" or equivalent
//
// In the assignment problem we will solve here we have n agents and n tasks; assigning agent
//   a
// to task
//   t
// costs
//   cost[a][t]
// The goal of the problem is to assign one agent to each task such that the total cost is minimized
// The total cost is the sum of the costs
//
// Things to do:
//   0. (mandatory)
//      Place the student numbers and names at the top of this file
//   1. (highly recommended)
//      Read and understand this code
//   2. (mandatory)
//      Modify the function generate_all_permutations to solve the assignment problem
//      Compute the best and worst solutions for all problems with sizes n=2,...,14 and for each
//      student number of the group
//   3. (mandatory)
//      Calculate and display an histogram of the number of occurrences of each cost
//      Does it follow approximately a normal distribution?
//      Note that the maximum possible cost is n * t_range
//   4. (optional)
//      For each problem size, and each student number of the group, generate one million (or more!)
//      random permutations and compute the best and worst solutions found in this way; compare
//      these solutions with the ones found in item 2
//      Compare the histogram computed in item 3 with the histogram computed using the random
//      permutations
//   5. (optional)
//      Try to improve the execution time of the program (use the branch-and-bound technique)
//   6. (optional)
//      Surprise us, by doing something more!
//   7. (mandatory)
//      Write a report explaining what you did and presenting your results
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#define NDEBUG  // uncomment to skip disable asserts (makes the code slightly faster)
#include <assert.h>


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// problem data
//
// max_n ........ maximum problem size
// cost[a][t] ... cost of assigning agent a to task t
//

//
// if your compiler complains about srandom() and random(), replace #if 0 by #if 1
//
#if 0
# define srandom srand
# define random  rand
#endif

#define max_n    32           // do not change this (maximum number of agents, and tasks)
#define range    20           // do not change this (for the pseudo-random generation of costs)
#define t_range  (3 * range)  // do not change this (maximum cost of an assignment)

static int cost[max_n][max_n];
static int seed; // place a student number here!


static void init_costs(int n)
{
	if(n == -3)
	{ // special case (example for n=3)
		cost[0][0] = 3; cost[0][1] = 8; cost[0][2] = 6;
		cost[1][0] = 4; cost[1][1] = 7; cost[1][2] = 5;
		cost[2][0] = 5; cost[2][1] = 7; cost[2][2] = 5;
		return;
	}
	if(n == -5)
	{ // special case (example for n=5)
		cost[0][0] = 27; cost[0][1] = 27; cost[0][2] = 25; cost[0][3] = 41; cost[0][4] = 24;
		cost[1][0] = 28; cost[1][1] = 26; cost[1][2] = 47; cost[1][3] = 38; cost[1][4] = 21;
		cost[2][0] = 22; cost[2][1] = 48; cost[2][2] = 26; cost[2][3] = 14; cost[2][4] = 24;
		cost[3][0] = 32; cost[3][1] = 31; cost[3][2] =  9; cost[3][3] = 41; cost[3][4] = 36;
		cost[4][0] = 24; cost[4][1] = 34; cost[4][2] = 30; cost[4][3] = 35; cost[4][4] = 45;
		return;
	}
	assert(n >= 1 && n <= max_n);
	srandom((unsigned int)seed * (unsigned int)max_n + (unsigned int)n);
	for(int a = 0;a < n;a++)
		for(int t = 0;t < n;t++)
			cost[a][t] = 3 + (random() % range) + (random() % range) + (random() % range); // [3,3*range]
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// code to measure the elapsed time used by a program fragment (an almost copy of elapsed_time.h)
//
// use as follows:
//
//   (void)elapsed_time();
//   // put your code to be time measured here
//   dt = elapsed_time();
//   // put morecode to be time measured here
//   dt = elapsed_time();
//
// elapsed_time() measures the CPU time between consecutive calls
//

#if defined(__linux__) || defined(__APPLE__)

//
// GNU/Linux and MacOS code to measure elapsed time
//

#include <time.h>

static double elapsed_time(void)
{
	static struct timespec last_time,current_time;

	last_time = current_time;
	if(clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&current_time) != 0)
		return -1.0; // clock_gettime() failed!!!
	return            ((double)current_time.tv_sec - (double)last_time.tv_sec)
				 + 1.0e-9 * ((double)current_time.tv_nsec - (double)last_time.tv_nsec);
}

#endif

#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)

//
// Microsoft Windows code to measure elapsed time
//

#include <windows.h>

static double elapsed_time(void)
{
	static LARGE_INTEGER frequency,last_time,current_time;
	static int first_time = 1;

	if(first_time != 0)
	{
		QueryPerformanceFrequency(&frequency);
		first_time = 0;
	}
	last_time = current_time;
	QueryPerformanceCounter(&current_time);
	return (double)(current_time.QuadPart - last_time.QuadPart) / (double)frequency.QuadPart;
}

#endif


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// function to generate a pseudo-random permutation
//

void random_permutation(int n,int t[n])
{
	assert(n >= 1 && n <= 1000000);
	for(int i = 0;i < n;i++)
		t[i] = i;
	for(int i = n - 1;i > 0;i--)
	{
		int j = (int)floor((double)(i + 1) * (double)random() / (1.0 + (double)RAND_MAX)); // range 0..i
		assert(j >= 0 && j <= i);
		int k = t[i];
		t[i] = t[j];
		t[j] = k;

	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// place to store best and worst solutions (also code to print them)
//

static int min_cost,min_cost_assignment[max_n], min_cost_row[max_n]; // smallest cost information
static int max_cost,max_cost_assignment[max_n], max_cost_row[max_n]; // largest cost information
static long n_visited; // number of permutations visited (examined)
// place your histogram global variable here
static double cpu_time;
static int histogram[max_n*t_range] = { 0 };


#define minus_inf  -1000000000  // a very small integer
#define plus_inf   +1000000000  // a very large integer

static void reset_solutions(void)
{
	min_cost = plus_inf;
	max_cost = minus_inf;
	n_visited = 0l;
	// place your histogram initialization code here
	for( int i = 0; i < 60*max_n; i++)
		histogram[i] = 0;
	

	cpu_time = 0.0;
}

#define show_info_1        (1 << 0)
#define show_info_2        (1 << 1)
#define show_costs         (1 << 2)
#define show_min_solution  (1 << 3)
#define show_max_solution  (1 << 4)
#define show_histogram     (1 << 5)
#define show_all           (0xFFFF)


double factorial(int n)
{
  double c;
  double result = 1;
 
  for (c = 1; c <= n; c++)
    result = result * c;
 
  return result;
}

static void show_solutions(int n,char *header,int what_to_show, int type)
{
	printf("%s\n",header);
	if((what_to_show & show_info_1) != 0)
	{
		printf("  seed .......... %d\n",seed);
		printf("  n ............. %d\n",n);
	}
	if((what_to_show & show_info_2) != 0)
	{
		printf("  visited ....... %ld\n",n_visited);
		printf("  cpu time ...... %.3fs\n",cpu_time);
	}
	if((what_to_show & show_costs) != 0)
	{
		printf("  costs .........");
		for(int a = 0;a < n;a++)
		{
			for(int t = 0;t < n;t++)
				printf(" %2d",cost[a][t]);
			printf("\n%s",(a < n - 1) ? "                 " : "");
		}
	}
	if((what_to_show & show_min_solution) != 0)
	{
		printf("  min cost ...... %d\n",min_cost);
		if(min_cost != plus_inf)
		{
			printf("  assignement ...");
			for(int i = 0;i < n;i++)
				printf(" %d",min_cost_assignment[i]);
			printf("\n");
		}
	}
	if((what_to_show & show_max_solution) != 0)
	{
		printf("  max cost ...... %d\n",max_cost);
		if(max_cost != minus_inf)
		{
			printf("  assignement ...");
			for(int i = 0;i < n;i++)
				printf(" %d",max_cost_assignment[i]);
			printf("\n");
		}
	}
	if((what_to_show & show_histogram) != 0)
	{
		// place your code to print the histogram here
		FILE *fp;	
		char name[60];
		char number[20];

		if (type==0){ // general permutations histogram
			strcpy(name, "GenN");
			snprintf(number, sizeof(number), "%d", n);
			strcat(name, number);
		} else if (type==1){ // random permutations histogram
			strcpy(name, "RanN");
			snprintf(number, sizeof(number), "%d", n);
			strcat(name, number);
		}

		fp=fopen(name, "w+");

		for(int i=0; i<=n*t_range;i++){
			fprintf(fp, "%d %d\n", i, histogram[i]);
		}
		if (type==1 || type == 0) {
			fclose(fp);
		}
		
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// code used to generate all permutations of n objects
//
// n ........ number of objects
// m ........ index where changes occur (a[0], ..., a[m-1] will not be changed)
// a[idx] ... the number of the object placed in position idx
//
// TODO: modify the following function to solve the assignment problem
//

static void generate_all_permutations(int n,int m,int a[n])
{
	if(m < n - 1)
	{
		for(int i = m;i < n;i++)
		{
#define swap(i,j)  do { int t = a[i]; a[i] = a[j]; a[j] = t; } while(0)
			swap(i,m);                            // exchange a[i] with a[m]
			generate_all_permutations(n,m + 1,a); // recurse
			swap(i,m);                            // undo the exchange of a[i] with a[m]
#undef swap
		}
	}
	else
	{
		n_visited++;

		int total_cost = 0;
		for( int i=0; i<n;i++){
			total_cost += cost[i][a[i]];
		}
		if (total_cost<min_cost){
			min_cost=total_cost;
			for(int j = 0; j < n; j++)
				min_cost_assignment[j] = a[j];
		}
		if (total_cost>max_cost){
			max_cost=total_cost;
			for(int j = 0; j < n; j++)
				max_cost_assignment[j] = a[j];
		}

		histogram[total_cost] += 1;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// code to generate min cost using branch-and-bound method of n objects
//
// n ........ number of objects
// m ........ index where changes occur (a[0], ..., a[m-1] will not be changed)
// a[idx] ... the number of the object placed in position idx
// parcial_cost ....... parcial cost

static void costs_row(int n)
{
	int min, max;
	for(int a = 0;a < n; a++){
		min = plus_inf; max = minus_inf;
		for(int t = 0; t < n; t++){
			if(min > cost[a][t])
				min = cost[a][t];
			if(max < cost[a][t])
				max = cost[a][t];
		}
		min_cost_row[a] = min;
		max_cost_row[a] = max;
	}
}

static void min_branch_and_bound(int n, int m,int a[n], int partial_cost)
{
	int min = 0;
	for(int i = m; i < n; i++){	min += min_cost_row[i]; }

	if(min_cost < partial_cost + min)
		return;

	if(m < n - 1)
	{
		for(int i = m;i < n;i++)
		{
#define swap(i,j)  do { int t = a[i]; a[i] = a[j]; a[j] = t; } while(0)
			swap(i,m);                            // exchange a[i] with a[m]
			min_branch_and_bound(n,m + 1,a, partial_cost + cost[m][a[m]]); // recurse
			swap(i,m);                            // undo the exchange of a[i] with a[m]
#undef swap
		}
	}
	else
	{
		int total_cost = partial_cost + cost[m][a[m]];

		if(min_cost > total_cost){
			min_cost = total_cost;
			for(int j = 0; j < n; j++)
				min_cost_assignment[j] = a[j];
		}

		n_visited++;
	}
}

static void max_branch_and_bound(int n, int m,int a[n], int partial_cost)
{
	int max = 0;
	for(int i = m; i < n; i++){	max += max_cost_row[i]; }

	if(max_cost > partial_cost + max)
		return;

	if(m < n - 1)
	{
		for(int i = m;i < n;i++)
		{
#define swap(i,j)  do { int t = a[i]; a[i] = a[j]; a[j] = t; } while(0)
			swap(i,m);                            // exchange a[i] with a[m]
			max_branch_and_bound(n,m + 1,a, partial_cost + cost[m][a[m]]); // recurse
			swap(i,m);                            // undo the exchange of a[i] with a[m]
#undef swap
		}
	}
	else
	{
		int total_cost = partial_cost + cost[m][a[m]];

		if(max_cost < total_cost){
			max_cost = total_cost;
			for(int j = 0; j < n; j++)
				max_cost_assignment[j] = a[j];
		}

		n_visited++;
	}
}

static void generate_all_permutations_branch_and_bound(int n, int m,int a[n], int partial_cost)
{
	costs_row(n);
	min_branch_and_bound(n, m, a, partial_cost);
	max_branch_and_bound(n, m, a, partial_cost);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// main program
//

int main(int argc,char **argv)
{
	if(argc == 2 && argv[1][0] == '-' && argv[1][1] == 'e')
	{
		seed = 0;
		{
			int n = 3;
			init_costs(-3); // costs for the example with n = 3
			int a[n];
			for(int i = 0;i < n;i++)
				a[i] = i;
			reset_solutions();
			(void)elapsed_time();
			generate_all_permutations(n,0,a);
			cpu_time = elapsed_time();
			show_solutions(n,"Example for n=3",show_all, 0);
			printf("\n");
		}
		{
			int n = 5;
			init_costs(-5); // costs for the example with n = 5
			int a[n];
			for(int i = 0;i < n;i++)
				a[i] = i;
			reset_solutions();
			(void)elapsed_time();
			generate_all_permutations(n,0,a);
			cpu_time = elapsed_time();
			show_solutions(n,"Example for n=5",show_all, 0);
			return 0;
		}
	}
	if(argc == 2)
	{
		seed = atoi(argv[1]); // seed = student number
		if(seed >= 0 && seed <= 1000000)
		{
			for(int n = 1;n <= max_n;n++)
			{
				init_costs(n);
				show_solutions(n,"-----------------------------------------------------------------------\nProblem statement",show_info_1 | show_costs, 0);
				//
				// 2.
				//
				if(n <= 14) // use a smaller limit here while developing your code
				{
					FILE *fp;
					FILE *map;
					FILE *mip;
					int a[n];
					fp=fopen("TimeGen", "a");
					mip=fopen("MinCostGen", "a");
					map=fopen("MaxCostGen", "a");
					for(int i = 0;i < n;i++)
						a[i] = i; // initial permutation
					reset_solutions();
					(void)elapsed_time();
					generate_all_permutations(n,0,a);
					cpu_time = elapsed_time();
					fprintf(fp, "%d %f\n", n, cpu_time);
					fclose(fp);
					fprintf(mip, "%d %d\n", n, min_cost);
					fclose(mip);
					fprintf(map, "%d %d\n", n, max_cost);
					fclose(map);
					show_solutions(n,"Brute force",show_info_2 | show_min_solution | show_max_solution | show_histogram,0);
				}

				if(n <= 16)
				{
					FILE *fp;
					FILE *map;
					FILE *mip;
					int a[n];
					fp=fopen("TimeBB", "a");
					mip=fopen("MinCostBB", "a");
					map=fopen("MaxCostBB", "a");
					for(int i = 0;i < n;i++)
						a[i] = i; // initial permutation
					reset_solutions();
					(void)elapsed_time();
					generate_all_permutations_branch_and_bound(n,0,a,0);
					cpu_time = elapsed_time();
					fprintf(fp, "%d %f\n", n, cpu_time);
					fclose(fp);
					fprintf(mip, "%d %d\n", n, min_cost);
					fclose(mip);
					fprintf(map, "%d %d\n", n, max_cost);
					fclose(map);
					show_solutions(n,"\nBrute force with branch-and-bound",show_info_2 | show_min_solution | show_max_solution, 2);	
				}
				FILE *fp;
				FILE *map;
				FILE *mip;
				int t[n];
				fp=fopen("TimeRan", "a");
				mip=fopen("MinCostRan", "a");
				map=fopen("MaxCostRan", "a");
				reset_solutions();
				(void)elapsed_time();
				for(int i=0;i<1000000;i++)
				{
					random_permutation(n, t);
					int total_cost = 0;
					for( int i=0; i<n;i++){
						total_cost += cost[i][t[i]];
					}
					if (total_cost<min_cost){
						min_cost=total_cost;
						for(int j = 0; j < n; j++)
							min_cost_assignment[j] = t[j];
					}
					if (total_cost>max_cost){
						max_cost=total_cost;
						for(int j = 0; j < n; j++)
							max_cost_assignment[j] = t[j];
					}
					histogram[total_cost] += 1;
				}
				cpu_time = elapsed_time();
				fprintf(fp, "%d %f\n", n, cpu_time);
				fprintf(mip, "%d %d\n", n, min_cost);
				fclose(mip);
				fprintf(map, "%d %d\n", n, max_cost);
				fclose(map);		
				n_visited=1000000;
				show_solutions(n,"\nRandom permutations",show_info_2 | show_min_solution | show_max_solution | show_histogram,1);
				fclose(fp);

				//
				// done
				//
				printf("\n");
			}
			return 0;
		}
	}
	fprintf(stderr,"usage: %s -e              # for the examples\n",argv[0]);
	fprintf(stderr,"usage: %s student_number\n",argv[0]);
	return 1;
}