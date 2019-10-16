#include <stdio.h>
#include <math.h>


//define a complex number struct
typedef struct complex {
	double a; //real part
	double b; //complex part
} complex;


//complex multiplication
complex compMult(complex A, complex B){
	complex output;
	
	output.a = A.a*B.a - A.b*B.b; //real part
	output.b = A.a*B.b + B.a*A.b; // complex part

	return output;
}


//complex addition
complex compAdd(complex A, complex B) {
	complex output;
	
	output.a = A.a + B.a; //real part
	output.b = A.b + B.b; // complex part
	
	return output;
}


//define function h_1(t)
complex h_1 (double t){
	complex output; 
	
	output.a = cos(t) + cos(5*t);
	output.b = sin(t) + sin(5*t);

	return output;
}


//define function h_2(t)
complex h_2 (double t){
	complex output;
	
	output.a = exp(-1*pow(t-M_PI,2)/2);
	output.b = 0;
	
	return output;
}


//function to print the imaginary and real parts of h_1 and h_2 to text files
//sampled 100 times between 0 and 2PI
void printToFile(){

	FILE *fp;
	
	//reset files upon each run of the program
	fp = fopen("h_1.txt","w");
	fp = fopen("h_2.txt","w");
	
	//print values of h_1 and h_2 sampled 100 times between 0 and 2PI to two separate files
	for(double i = 0; i < (99*M_PI/50); i += ((2*M_PI)/100)){
		
		//h_1 - store sampled values in h_1.txt
		fp = fopen("h_1.txt","a");
		fprintf(fp,"%f + %fi\n", h_1(i).a, h_1(i).b);
		fclose(fp);	
		
		//h_2 - store sampled values in h_2.txt
		fp = fopen("h_2.txt","a");
		fprintf(fp,"%f + %fi\n", h_2(i).a, h_2(i).b);
		fclose(fp);		
	}
}


//FF: Fourier Factor - from equation 7
//The complex number multiplied by the sampled function to produce the discrete Fourier values
complex FF (int n, int k, int N){
	complex output;
	
	output.a = cos(2*M_PI*n*k/N); //real part
	output.b = -1*sin(2*M_PI*n*k/N); //complex part
	
	return output;
}


//IFF: Inverse Fourier Factor - from equation 8
//The complex number multiplied by the discrete Fourier values to produce the inverse Fourier values
complex IFF (int n, int k, int N){
	complex output;
	
	output.a = cos(2*M_PI*n*k/N); //real part
	output.b = sin(2*M_PI*n*k/N); //imaginary part
	
	return output;
}


//Fourier Transform on an arbitrary function 
//Values are stored in an array 
void FourierTransform (complex (*p_func)(double), double Delta, int N, complex store[], int check){

	complex H_f_n; //f is an arbitrary function and n is the discrete Fourier transformed number
	
	for(int n = 0; n < N; n++){
	
		//reset variables after each H_f,n value is calculated and stored
		H_f_n.a = 0;
		H_f_n.b = 0;
		
		for(int k = 0; k < N; k++){
	
			complex temp1 = (*p_func)(k*Delta); //calculate the value of the function at k*Sampling interval
			complex temp2 = FF(n, k, N); //calculate the Fourier Factor for this n, k and N
			complex temp3 = compMult(temp1,temp2); //multiply the function value by the Fourier Factor
			H_f_n = compAdd(H_f_n, temp3); //Sum all of these values up to N-1 
	
		}
		
		store[n] = H_f_n; //store the calculated H_f,n value to the corresponding value of n
			
		switch(check){ //the check input is used to correctly output to the user which function is being used
		
			case 1: //h_1
				if(n == 0){
					printf("Discrete Fourier Transform values for h_1:\n");
				}
				
				printf("H_1,%d %lf + %lfi\n",n, H_f_n.a, H_f_n.b);
				break;
				
			case 2: //h_2
				if(n == 0){
					printf("Discrete Fourier Transform values for h_2:\n");
				}
				
				printf("H_2,%d %lf + %lfi\n",n, H_f_n.a, H_f_n.b);
				break;
					
			default: //any other function
				if(n == 0){
					printf("Discrete Fourier Transform values for h_f:\n");
				}
				
				printf("H_f,%d %lf + %lfi\n",n, H_f_n.a, H_f_n.b);
				break;	
		} 	
	}
}


//Inverse Fourier Transform 
//can be used on any discrete Fourier values stored inside an array
//for the three functions used the values are stored in text files
//for an unknown function values are output to the user
void IFT (int N, complex FTValues[], int check){
	
	complex h_f_k; //f is arbitrary function and f(Delta*k) is the recovered value of the function
	FILE *fp;

	for(int k = 0; k < N; k++){
	
		//reset variables after each h_f,k value is calculated and stored
		h_f_k.a = 0;
		h_f_k.b = 0;
	
		for(int n = 0; n < N; n++){
		
			complex temp1 = FTValues[n]; //retrieve the stored discrete Fourier value at n
			complex temp2 = IFF(n, k, N); //calculate the Inverse Fourier Factor for this n, k and N
			complex temp3 = compMult(temp1,temp2); //multiply the discrete Fourier value by the Inverse Fourier Factor
			h_f_k = compAdd(h_f_k, temp3); //Sum all of these values up to N-1 
	
		}
	
		h_f_k.a = h_f_k.a/N; //divide the real part by N
		h_f_k.b = h_f_k.b/N; //divide the imaginary part by N
		
		switch(check){ //the check input is used to correctly store the Inverse Fourier values in the right file
		
			case 1: //h_1 - store Inverse Fourier values in h_1'.txt
				if(k == 0){
					//reset file upon each run of the program
					fp = fopen("h_1'.txt","w");
				}
					
				fp = fopen("h_1'.txt","a");
				fprintf(fp,"%f + %fi\n", h_f_k.a, h_f_k.b);
				fclose(fp);	
				break;
			
			case 2: //h_2 - store Inverse Fourier values in h_2'.txt
				if(k == 0){
					//reset file upon each run of the program
					fp = fopen("h_2'.txt","w");
				}
					
				fp = fopen("h_2'.txt","a");
				fprintf(fp,"%f + %fi\n", h_f_k.a, h_f_k.b);
				fclose(fp);	
				break;
				
			case 3: //h_3 - store Inverse Fourier values in h_3'.txt
				if(k == 0){
					//reset file upon each run of the program
					fp = fopen("h_3'.txt","w");
				}
				
				fp = fopen("h_3'.txt","a");
				fprintf(fp,"%f + %fi\n", h_f_k.a, h_f_k.b);
				fclose(fp);	
				break;
			
			default: //any other Inverse Fourier Transfrom
				printf("h_f',%d: %lf + %lfi\n", k, h_f_k.a, h_f_k.b);
				break;
		}
	}  
}


//read the values from h3.txt into a 2d array
void ReadFromFile(double Store[200][4]){

	FILE *h3;
	h3 = fopen("h3.txt", "r"); //open the file for reading
	
	double *p; //declare a double pointer
	p = &Store[0][0]; //set the pointer to the first position of the 2d array
	
	for(int i = 0; i < 800; i+=4){ //cycles through each row of the file; skips 4 to allow each value to be read into separate array positions
	
		//scan the four elements in each row of the file and place them in consecutive positions in the 2d array
		fscanf(h3, "%lf%*c %lf%*c %lf%*c %lf", i+p, (i+p+1), (i+p+2), (i+p+3));
	}
}


//convert the real and imaginary parts from the h3.txt file to one complex array
void ConvertToComplex (complex Convert[], double Store[200][4]){
	
	//cycles through each row of the 2d array and combines the real and imaginary parts 
	//to produce one complex number and store it in another array
	for(int i = 0; i < 200; i++){ 
	
		Convert[i].a = Store[i][2];
		Convert[i].b = Store[i][3];
	}
}


//Fourier Transform for the data read from the h3.txt file
//Values are stored in an array 
void h3FourierTransform (int N, complex store[], complex h3Input[]){

	complex H_3_n; //n is the discrete Fourier Transformed number
	
		for(int n = 0; n < N; n++){
	
			//reset variables after each H_3,n value is calculated and stored
			H_3_n.a = 0;
			H_3_n.b = 0;
		
			for(int k = 0; k < N; k++){
	
				complex temp1 = h3Input[k]; //retrieve the stored h_3 value at k*Delta
				complex temp2 = FF(n, k, N); //calculate the Fourier Factor for this n, k and N
				complex temp3 = compMult(temp1,temp2); //multiply the function value by the Fourier Factor
				H_3_n = compAdd(H_3_n, temp3); //Sum all of these values up to N-1 
	
			}
		
			store[n] = H_3_n; //store the calculated H_3,n value to the corresponding value of n
			//printf("H_3,%d %.2f + %.2fi\n",n, H_3_n.a, H_3_n.b);
		}
}


//Picks the top four values of mod(H_3) and populates a new array with the values
void h3FourierTopFour (complex FourierStoreh3[], complex FourierStoreh3T4[]){

	//array to store the top four absolute H_3,n values and their index in the original Fourier Transfrom array
	double Top4[4][2] = {{0,0}, {0,0}, {0,0}, {0,0}}; 
	double tempValue1, tempIndex1, tempValue2, tempIndex2, tempValue3, tempIndex3;
	
	//Bubble Sort
	/* This for loop cycles through every Fourier transformed value of the h_3 function and compares
	the absolute value with the number stored in the last (4th) position of the Top4 array, if it is larger
	then it replaces the last value in the Top4 array with the current value and stores the index (value of i).
	It then checks to see if the last value is larger than the above value (3rd) and if it is then it will 
	swap those two values and their associated indexes in the array. This continues until the value is smaller
	than the value above or it has reached the top of the array. */
	
	for(int i = 0; i < 200; i++){
	
		if(fabs(FourierStoreh3[i].a) > fabs(Top4[3][0])){
		
			Top4[3][0] = FourierStoreh3[i].a;
			Top4[3][1] = i;
			
			if(fabs(Top4[3][0]) > fabs(Top4[2][0])){
			
				tempValue1 = Top4[2][0];
				tempIndex1 = Top4[2][1];
				
				Top4[2][0] = Top4[3][0];
				Top4[2][1] = Top4[3][1];
				
				Top4[3][0] = tempValue1;
				Top4[3][1] = tempIndex1;
				
				if(fabs(Top4[2][0]) > fabs(Top4[1][0])){
				
					tempValue2 = Top4[1][0];
					tempIndex2 = Top4[1][1];
				
					Top4[1][0] = Top4[2][0];
					Top4[1][1] = Top4[2][1];
				
					Top4[2][0] = tempValue2;
					Top4[2][1] = tempIndex2;
				
					if(fabs(Top4[1][0]) > fabs(Top4[0][0]) ){
				
						tempValue3 = Top4[0][0];
						tempIndex3 = Top4[0][1];
				
						Top4[0][0] = Top4[1][0];
						Top4[0][1] = Top4[1][1];
				
						Top4[1][0] = tempValue3;
						Top4[1][1] = tempIndex3;
					
					}
				}	
			}
		}
	}
	
	
	/* This for loop cycles through every value in the new 'Fourier Store h_3 Top 4' array and assigns 
	every value to be 0 + 0i unless the index (i) is equal to an index of one of the top 4 values stored 
	in the Top4 array where it subsequently inserts the value stored. */
	
	for(int i = 0; i < 200; i++){
	
		FourierStoreh3T4[i].a = 0.0;
		FourierStoreh3T4[i].b = 0.0;
		
		if(i == Top4[0][1]){
		
			FourierStoreh3T4[i].a = Top4[0][0];
			FourierStoreh3T4[i].b = 0.0;
		}
		
		if(i == Top4[1][1]){
		
			FourierStoreh3T4[i].a = Top4[1][0];
			FourierStoreh3T4[i].b = 0.0;
		}
		
		if(i == Top4[2][1]){
		
			FourierStoreh3T4[i].a = Top4[2][0];
			FourierStoreh3T4[i].b = 0.0;
		}
		
		if(i == Top4[3][1]){
		
			FourierStoreh3T4[i].a = Top4[3][0];
			FourierStoreh3T4[i].b = 0.0;
		}
	}
}


int main() {
	
	double Delta = 2*M_PI/100; //value of Delta for h_1 and h_2
	
	complex FourierStoreh_1[100]; //array used to store the Fourier Transformed values of h_1
	complex FourierStoreh_1_2[100]; //array used to store the Fourier Transformed values of h_1 with H_1,1 set to 0+0i
	
	complex FourierStoreh_2[100]; //array used to store the Fourier Transformed values of h_2
	complex FourierStoreh_2_2[100]; //array used to store the Fourier Transformed values of h_2 with H_2,0 set to 0+0i
	
	double Store[200][4]; //array used to store the values imported from the h3.txt file
	complex h_3Complex[200]; //array used to store the combined real and imaginary parts from the h3.txt file as a complex number
	complex FourierStoreh_3[200]; //array used to store the Fourier Transformed values of h_3
	complex FourierStoreh_3_T4[200]; //array used to store only the top 4 absolute Fourier Transformed values of h_3, the rest set to 0+0i
	
	complex (*p_h_1)(double); //pointer for the h_1 function
	complex (*p_h_2)(double); //pointer for the h_2 function
	p_h_1 = &h_1; //assign the pointer value to point to the h_1 function
	p_h_2 = &h_2; //assign the pointer value to point to the h_2 function
	
	
	printToFile(); //print the sampled values of h_1 and h_2 to h_1.txt and h_2.txt files
	
	
	//Discrete Fourier Transform on the h_1 function with the values stored in the FourierStoreh_1 array
	FourierTransform(p_h_1, Delta, 100, FourierStoreh_1, 1); 
	//Discrete Fourier Transform on the h_2 function with the values stored in the FourierStoreh_2 array
	FourierTransform(p_h_2, Delta, 100, FourierStoreh_2, 2);
	
	
	/* This for loop assigns all of the h_1 Fourier Transformed values from the FourierStoreh_1 array to the 
	FourierStoreh_1_2 array. */
	for(int i = 0; i < 100; i++){
		FourierStoreh_1_2[i] = FourierStoreh_1[i];
	}
	
	//Assigns the H_1,1 value to 0+0i
	FourierStoreh_1_2[1].a = 0;
	FourierStoreh_1_2[1].b = 0;
	
	
	/* This for loop assigns all of the h_2 Fourier Transformed values from the FourierStoreh_2 array to the 
	FourierStoreh_2_2 array. */
	for(int i = 0; i < 100; i++){
		FourierStoreh_2_2[i] = FourierStoreh_2[i];
	}
	
	//Assigns the H_2,0 value to 0+0i
	FourierStoreh_2_2[0].a = 0;
	FourierStoreh_2_2[0].b = 0;
	
	
	//Inverse Fourier Transform on the FourierStoreh_1_2 array; results printed to the h_1'.txt file 
	IFT(100, FourierStoreh_1_2, 1);
	//Inverse Fourier Transform on the FourierStoreh_2_2 array; results printed to the h_2'.txt file
	IFT(100, FourierStoreh_2_2, 2); 
	
	
	//Read the values from the h3.txt file and storing them in the 'Store' array
	ReadFromFile(Store);
	//Convert the real and imaginary parts from the h3.txt file to one complex number and storing them in the h_3Complex array
	ConvertToComplex(h_3Complex, Store);
	
	
	//Discrete Fourier Transform on the h_3Complex array with the values stored in the FourierStoreh_3 array
	h3FourierTransform(200, FourierStoreh_3, h_3Complex);
	/* Pick the top four absolute Fourier Transformed values from FourierStoreh_3 array and populate the FourierStoreh_3_T4 array
	with all 0+0i except the top four values */
	h3FourierTopFour(FourierStoreh_3, FourierStoreh_3_T4);
	
	
	//Inverse Fourier Transform on the FourierStoreh_3_T4 array; results printed to the h_3'.txt file 
	IFT(200, FourierStoreh_3_T4, 3); 
	
	return 0;
}