/***************************************************************************
	changepoint.h  -  Header file for change point detection
	-----------------
	copyright            : (C) 2004 by Haw Yang
	email                : hawyang@princeton.edu
 ***************************************************************************/

#define CHANGEPOINT_VERSION "1.2-gaussian-"
#define COMPILE_DATE "20120715"

#ifdef __MINGW32__
	#ifndef PLATFORM
		#define PLATFORM "MinGW"
	#endif
#endif
#ifdef __GNUC__
	#ifndef PLATFORM
		#define PLATFORM "Linux-gnu"
	#endif
#endif

/*******************************
*         Data Structures      *
*******************************/
enum bool {false=0, true=1};
struct data {
	double        time;       // chronological time stamp
	double        dt;         // inter-photon duration
	double        value;      // numeric value of the data
	};

struct changepoint {
	size_t              i;    // change point index to the trajectory
	size_t              il;   // change point confidence left bound
	size_t              ir;   // change point confidence right bound
	struct changepoint  *left;
	struct changepoint  *right;
	int                 height;
	};
			
/*******************************
*      Function Declaration    *
*******************************/
int FindCP();
void CheckCP();
void MakeCPArray();
double zbrent();
double Vostrikova();
double C_ac();

