// C++ Gillespie acode for simulating NFkB gene regulatory model after Sneppen (2006), Hoffmann et al (2002) 

/*////////////////////// SUMMARY OF THE MODEL ///////////////////////////

Species={Nn:0 | In:1 | N:2 | I:3 | Cn:4 | C:5 | Im:6 | ON:7 | OFF:8 | Dc:9 | 10: DcNn} 

R0:
        ON > OFF + Nn
        koff*ON
R1:
        OFF + Nn > ON
        kon*Nn*OFF
R2:
        N > Nn
        k_Nin*N
R3:
        Nn + In > Cn
        k_fn*Nn*In
R4:
        Cn > Nn + In
        k_bn*Cn
R5:
        ON > Im + ON
        k_t*ON
R6:
        Im > $pool
        gm*Im
R7:
        Im > Im + I
        k_tl*Im
R8:
        N + I > C
        k_f*N*I
R9:
        C > N + I
        k_b*C
R10:
        I > In
        k_Iin*I
R11:
        In > I
        k_Iout*In
R12:
        C  > N
	a*C
R13:
        Cn > C
        k_Cout*Cn
        
R14:    Nn+Dc -> DcNn
            kd_on*Nn*Dc

R15:    DcNn -> Dc + Nn
            kd_off*DcNn

R16: In+DcNn -> Cn+Dc
            ks*In*Dc*Nn

R17: In+ON -> Cn+OFF
            ks*In*ON

*////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <limits.h>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <sys/time.h>

#define TRAJ      1          	// # trajectories
#define ENDTIME   20000.0	// end of time
#define TIMESTEP  1		// interval of output
#define N         18		// number of reactions
#define M         11		// number of chemical species

using namespace std;

int x[M];				// population of chemical species
double k[N];				// reaction rates
double p[N];				// propencities of reaction
int S[M][N];				// data structure for updating x[]


//////////////// Initial conditions /////////////////

void init(int x[], double k[], int S[][N], char *FF){

	ifstream inFile;
	inFile.open(FF);
	if (!inFile) {
  		cerr << "Can't open input file " << FF << endl;
  		exit(1);
	}
	char name[9];		
	double value;
	double params[M+N];
	int pp=0;
		
	while (inFile >> name >> value) {
		params[pp]=value;
		pp++;	
	}
	inFile.close();
	
	// population of chemical species
	
	for(int i = 0; i < M; i ++) {
		x[i] = static_cast<int>(params[i]);
	} 

	// reaction rates

	for(int i = M; i < M+N; i ++) {
		k[i-M]=params[i];
	} 

	// Stochiometric matrix		
	for(int i = 0; i < M; i ++) {
  	  for(int j = 0; j < N; j ++) S[i][j] = 0;
	}

	S[0][0] =   1; 	  			// Change of Species 1 by R 0
	S[0][1] =  -1;   			// Change of Species 1 by R 1
	S[0][2] =   1;   			// Change of Species 1 by R 2
	S[0][3] =  -1;   			// ...
	S[0][4] =   1;
	S[0][14] = -1; 
	S[0][15] = 1;
	S[1][3] =  -1;
	S[1][4] =   1;
	S[1][10] =  1;
	S[1][11] = -1;
	S[1][16] = -1;
	S[1][17] = -1;
	S[2][2] =  -1;
	S[2][8] =  -1;  
	S[2][9] =   1;
	S[2][12] =  1;
	S[3][7]  =  1;
	S[3][8]  = -1;
	S[3][9]  =  1;
	S[3][10] = -1;
	S[3][11] =  1; 
	S[4][3]  =  1;
	S[4][4]  = -1;
	S[4][13] = -1;
	S[4][16] = 1;
	S[4][17] = 1;
	S[5][8]  =  1;
	S[5][9]  = -1;
	S[5][12] = -1; 
	S[5][13] =  1;
	S[6][5]  =  1; 
	S[6][6]  = -1;
	S[7][0]  = -1; 
	S[7][1]  =  1;
	S[7][17] = -1;
	S[8][0]  =  1;
	S[8][1]  = -1;
	S[8][17] = 1;
	S[9][14] = -1;
	S[9][15] = 1;
	S[9][16] = 1;
	S[10][14] = 1;
	S[10][15] = -1;
	S[10][16] = -1;
}

//////////////// Prototypes /////////////////

void init(int x[], double k[], int S[][N]);
void update_p(double p[], double k[], int x[]);
void update_x(int x[], double k[], int S[][N], int reaction);
double sum(double a[], int n);
int select_reaction(double p[], int pn, double sum_propencity, double r);
string make_filename( const string& basename, const string& ext );

////////////////////////// Function definitions ///////////////////////////////////////

void update_p(double p[], double k[], int x[]){

	p[0] = k[0]*x[7];	// R1
	p[1] = k[1]*x[0]*x[8];	// R2
	p[2] = k[2]*x[2];	// R3
	p[3] = k[3]*x[0]*x[1];  // ...
	p[4] = k[4]*x[4]; 
	p[5] = k[5]*x[7];
	p[6] = k[6]*x[6];
	p[7] = k[7]*x[6];
	p[8] = k[8]*x[2]*x[3];
	p[9] = k[9]*x[5];
	p[10] = k[10]*x[3];
	p[11] = k[11]*x[1];
	p[12] = k[12]*x[5];
	p[13] = k[13]*x[4];
	p[14] = k[14]*x[0]*x[9];
	p[15] = k[15]*x[10];
	p[16] = k[16]*x[1]*x[10];
	p[17] = k[17]*x[1]*x[7];

}

void update_x(int x[], double k[], int S[][N], int reaction){

	for (int i=0; i < M; i++){
		x[i] = x[i] + S[i][reaction]; 	
	}

}

double sum(double a[], int n){
	int i;
	double s=0.0;
	for(i=0; i<n; i++) 
		s += a[i];
	return(s);
}

int select_reaction(double p[], int pn, double sum_propencity, double r){
	int reaction = -1;
	double sp = 0.0;
	int i;
	r = r * sum_propencity;
	for(i=0; i<pn; i++){
		sp += p[i];
		if(r < sp){
			reaction = i;
			break;
		}
	}
	return reaction;
}

string make_filename( const string& base, const string& ext ){
  ostringstream result;
  result << base << ext;
  return result.str();
  }

//////////////////////////////////////////////// Gillespie engine ///////////////////////////////////////////////////


int main(int argc, char *argv[]){

//	srand(time(NULL));
        struct timeval time;
        gettimeofday(&time,NULL);

     // microsecond has 1 000 000
     // Assuming you did not need quite that accuracy
     // Also do not assume the system clock has that accuracy.
        srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
     // 
     // The trouble here is that the seed will repeat every
     // 24 days or so.
     //
     // If you use 100 (rather than 1000) the seed repeats every 248 days.
     //
     // Do not make the MISTAKE of using just the tv_usec
     // This will mean your seed repeats every second.
     

	if (argc < 2 ) {
        cout << "ERROR: At least 1 file argument is required, 0 given!"<< endl;
	cout<< "PROPER USAGE:  1. Input file name, 2. Output file name(OPTIONAL)"<<endl;
	return 1;
    	}
	if (argv[2] == NULL){ 
		argv[2] = "gill_out_";
	}
	
	for( int traj=0; traj < TRAJ; traj++){
		
  		

		// initialize sim params
		double sum_propencity = 0.0;	// sum of propencities
		double tau=0.0;			// step of time
		int reaction=0;			// reaction number selected
		double t=0.0;			// time
		int tn = 0; 			// time interval counter
		double r1;			// random number-1
		double r2;			// random number-2

		// read init data
		init(x, k, S, argv[1]);
		
		FILE *out =  fopen(make_filename(argv[2], ".dat" ).c_str(), "w");

		// main loop 
		while (t < ENDTIME){
			
			// output			//output(out, t, x, M);
			
			if (t >= tn*TIMESTEP){	
				fprintf(out, "%f", t);
				for(int i=0; i<M; i++){
					fprintf(out, "\t%d", x[i]); 
					}
				fprintf(out, "\n");
				tn +=1;		
			}

			// Update propencity
			update_p(p, k, x);
			sum_propencity = sum(p, N);

			// Sample tau
			double r1 = rand()/(double)RAND_MAX;
			if(sum_propencity > 0){
				tau = -log(r1)/sum_propencity;
			}else{
				break;
			}

			// Sample reaction
			double r2 = rand()/(double)RAND_MAX;
			reaction = select_reaction(p, N, sum_propencity, r2);

			// Update species
			update_x(x, k, S, reaction);
			// Update time
			t += tau;
			
		}

		fclose(out);
		
	}
	
	return 0;
}
