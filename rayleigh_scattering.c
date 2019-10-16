#include <stdio.h>
#include <math.h>

//Linear Congruential Sequence
//produces a random number between 0 and m based on the last produced value X
//a, c and m are the seed values for the random sequence
long lcs (long X, int a, int c, long m){

	long Xn;
	Xn = (a*X + c)%m;
	
	return Xn;
}


//Cumulative Method for producing random numbers with distribution 3/8(1+cos^2x)sinx
//argument is a random number between 0 and 1
double thetaHat (double rand){
	
	//method produces a cubic equation which is solved using the cubic formula
	double q,delta,x0;
	double p = 3.0;
	
	q = 8.0*rand - 4.0;	
	delta = pow(q,2.0)/4.0 + pow(p,3.0)/27.0;
	x0 = cbrt(-q/2.0 + sqrt(delta)) + cbrt(-q/2.0 - sqrt(delta));
	x0 = acos(x0);
	
	return x0;
}


//Cumulative Method for producing random numbers with distribution sinx
//argument is a random number between 0 and 1
double theta (double rand){
	
	double x0;
	
	x0 = acos(1-2.0*rand);
	
	return x0;
}


//Cumulative Method for producing random numbers with distribution e^-x
//argument is a random number between 0 and 1
double Tau (double rand){
	
	double x0;
	
	x0 = -log(1-rand);
	
	return x0;
}


//Function to find the angle between two 3x1 column vectors
double angle(double a1, double a2, double a3, double b1, double b2, double b3){
	
	double dot, magA, magB, theta;
	
	dot = a1*b1 + a2*b2 + a3*b3; //find the dot product of the two vectors
	magA = sqrt(pow(a1,2) + pow(a2,2) + pow(a3,2)); //find the magnitude of vector A
	magB = sqrt(pow(b1,2) + pow(b2,2) + pow(b3,2)); //find the magnitude of vector B
	
	theta = acos(dot/(magA*magB)); //calculate the angle between A and B
	
	return theta;
}


//Method to simulate isotropic scattering of a particle released directly up from (0,0,0)
//arguments: Zmin: lower bound of z axis, Zmax: upper bound of z axis (used for detection)
//TauMax: optical depth, N: number of photons to be simulated, noBins: the number of detection bins
//scatterProb: probability of a particle being scattered (not absorbed) 
void isotropic(double Zmin, double Zmax, double TauMax, int N, int noBins, double scatterProb){
	
	//declare the seed values used for random number generation 
	int a = 16807;
	int c = 0;
	long m = 2147483647;
	long X0 = 5;

	int count = 0; //photon detection counter (will equal N if particles cannot be absorbed)
	int bins[noBins]; //detection array
	double binMax[noBins]; //store the maximum value of each bin

	for(int i=0; i<noBins; i++){
		bins[i] = 0; //clear storage 
		binMax[i] = 0;
		
		//calculate and store the maximum value of each bin
		if(i == 0){
			binMax[i] = acos(1 - 1.0/20.0);
		}
		
		binMax[i] = acos(cos(binMax[i-1]) - 1.0/20.0);
	}
	
	//declare variables to be used
	double Z, X, Y, TauTemp, ThetaTemp, PhiTemp, L, randX, Xtot, Ytot;
	int binCheck = 0;
	int binCounter = 0;
	int scatterCheck = 0;
	
	//simulate N photons 
	for(int i = 0; i<N; i++){
		
		scatterCheck = 0; //reset scatter check
		
		X0 = lcs(X0,a,c,m);
		randX = (double)X0/m; //generate a random number between 0 and 1
		
		//move the particle directly up to a random distance Z (based on e^x distribution)
		X = 0.0;
		Y = 0.0;
		Z = (Tau(randX)/TauMax) * Zmax;
		
		X0 = lcs(X0,a,c,m);
		randX = (double)X0/m; //generate a random number between 0 and 1
		
		if(randX<=scatterProb){ //detect if the photon has been scattered or absorbed
		
			if(Z<Zmax){ //check the photon has not already escaped
				//repeat while the photon remains in the relevant space and has not been absorbed 
				while(Z<Zmax && Z>Zmin && scatterCheck==0){ 
				
					//generate a random theta, phi and distance based on the relevant probability distributions
					X0 = lcs(X0,a,c,m);
					randX = (double)X0/m;
					ThetaTemp = theta(randX);
				
					X0 = lcs(X0,a,c,m);
					randX = (double)X0/m;
					PhiTemp = randX*2*M_PI;
					
					X0 = lcs(X0,a,c,m);
					randX = (double)X0/m;
					TauTemp = Tau(randX);
					L = (TauTemp/TauMax) * Zmax;
			
					//update the postion of the photon
					X = X + L*sin(ThetaTemp)*cos(PhiTemp);
					Y = Y + L*sin(ThetaTemp)*sin(PhiTemp);
					Z = Z + L*cos(ThetaTemp);	
					
					//detect if the photon has been scattered or absorbed
					X0 = lcs(X0,a,c,m);
					randX = (double)X0/m;
					if(randX>scatterProb){
						scatterCheck = 1;
					}
				}
			} 
			
			if(scatterCheck == 0){ //if the photon has not been absorbed
				if(Z>Zmax){ //if the photon has escaped out of the positive detection side
			
					//check for anomalous photons that escape with theta>pi/2
					if(ThetaTemp > (M_PI/2.0)){
					
						i--; //set the photon counter back by 1
					} else{
					
						//record the magnitude of the X and Y distance away from the start point
						Xtot = Xtot + fabs(X);
						Ytot = Ytot + fabs(Y);
						count++; //increment the photon counter
					}
					
					//method to bin the detected angles
					//checks if the detected angle is less than each maximum bin value
					//breaks once the angle has been binned
					while(binCheck == 0){
						if(binCounter>(noBins-1)){
							binCheck++;
						} else if(ThetaTemp<=binMax[binCounter]){
							bins[binCounter]++;
							binCheck++;
						}	
						binCounter++;
					}
					
					//reset counters 
					binCheck = 0; 
					binCounter = 0;
			
				} else {
					i--; //if photon returns out of the bottom of the space set the photon counter back by 1
				}
			}
		} 
	}
	
	
	double centralBinAngle, binNEnergy, binError;
	
	printf("Isotropic Scattering\n");
	printf("Tau: %.2f Zmin: %.2f Zmax: %.2f\n",TauMax,Zmin,Zmax);
	printf("Bin No.   Central Bin Angle   Normalised Energy      Error\n");
	
	//method to output the information 
	for(int i=0; i<noBins; i++){
		
		//calculate the central values of each bin
		if(i == 0){
			centralBinAngle = binMax[i]/2.0;
		} else {
			centralBinAngle = (binMax[i-1]+binMax[i])/2.0;
		}
		
		//calculate and output the bin number, central bin angle, normalised energy and error
		binNEnergy = (double)bins[i]/count;
		binError = binNEnergy/sqrt((double)bins[i]);
		
		printf("%d              %.6f            %.6f         %.6f\n",i+1,centralBinAngle,binNEnergy,binError);
		
		binCounter = binCounter + bins[i];
	} 
	
	printf("Percentage of photons absorbed: %.2f%%\n", 100.0 - (double)count/N * 100.0);
	printf("Average Magnitude X: %.6f Y: %.6f\n", (double)Xtot/count, (double)Ytot/count);
}


//Method to simulate Rayleigh scattering of a particle released directly up from (0,0,0)
//arguments: Zmin: lower bound of z axis, Zmax: upper bound of z axis (used for detection)
//TauMax: optical depth, N: number of photons to be simulated, noBins: the number of detection bins
//scatterProb: probability of a particle being scattered (not absorbed) 
void Rayleigh(double Zmin, double Zmax, double TauMax, int N, int noBins, double scatterProb){

	//declare the seed values used for random number generation 
	int a = 16807;
	int c = 0;
	long m = 2147483647;
	long X0 = 5;

	int count = 0; //photon detection counter (will equal N if particles cannot be absorbed)
	int bins[noBins]; //detection array
	double binMax[noBins]; //store the maximum value of each bin

	for(int i=0; i<noBins; i++){
		bins[i] = 0; //clear storage 
		binMax[i] = 0;
		
		//calculate and store the maximum value of each bin
		if(i == 0){
			binMax[i] = acos(1 - 1.0/20.0);
		}
		
		binMax[i] = acos(cos(binMax[i-1]) - 1.0/20.0);
	}
	
	//declare variables to be used
	double Z, X, Y, TauTemp, ThetaTemp, PhiTemp, L, randX, Xtot, Ytot, Xtemp, Ytemp, Ztemp, alpha, beta, det;
	int binCheck = 0;
	int binCounter = 0;
	int scatterCheck = 0;
	
	//[x][y][z]
	double direction[3][2]; //array to store coordinates of previous two points
	double ZaxisDirection[3]; //the Z axis in the photons frame
	double NaxisDirection[3]; //vector normal to the new Z axis
	double tFC[3]; //tempFrameCoords
	int posCounter;
	
	//simulate N photons 
	for(int i = 0; i<N; i++){
		
		scatterCheck = 0; //reset scatter check 
		posCounter = 0; //reset the position counter
		
		X0 = lcs(X0,a,c,m);
		randX = (double)X0/m;
		
		//move the particle directly up to a random distance Z (based on e^x distribution)
		X = 0.0;
		Y = 0.0;
		Z = (Tau(randX)/TauMax) * Zmax;
		
		//store the first point in the direction array 
		direction[0][0] = X;
		direction[1][0] = Y;
		direction[2][0] = Z;
		posCounter++;
		
		X0 = lcs(X0,a,c,m);
		randX = (double)X0/m;
		
		if(randX<=scatterProb){ //detect if the photon has been scattered or absorbed
		
			if(Z<Zmax){//check the photon has not already escaped
				//repeat while the photon remains in the relevant space and has not been absorbed 
				while(Z<Zmax && Z>Zmin && scatterCheck == 0){
		
					//generate a random theta, phi and distance based on the relevant probability distributions
					X0 = lcs(X0,a,c,m);
					randX = (double)X0/m;
					ThetaTemp = thetaHat(randX);
					
					X0 = lcs(X0,a,c,m);
					randX = (double)X0/m;
					PhiTemp = randX*2*M_PI;
			
					X0 = lcs(X0,a,c,m);
					randX = (double)X0/m;
					TauTemp = Tau(randX);
					L = (TauTemp/TauMax) * Zmax;
					
					//if calculating the first point then the new Z axis direction is the same as the old z axis direction
					//therefore the positions can already be calculated
					if(posCounter == 1){
					
						//update the postion of the photon
						X = X + L*sin(ThetaTemp)*cos(PhiTemp);
						Y = Y + L*sin(ThetaTemp)*sin(PhiTemp);
						Z = Z + L*cos(ThetaTemp);
						
						//detect if the photon has been scattered or absorbed
						X0 = lcs(X0,a,c,m);
						randX = (double)X0/m;
						if(randX>scatterProb){
							scatterCheck = 1;
						}	
						
						//store the second point in the direction array
						direction[0][1] = X;
						direction[1][1] = Y;
						direction[2][1] = Z;
						posCounter++;
						
					}else {
						
						//calculate the direction vector of the Z axis by substituting the position vectors of the two previous points
						ZaxisDirection[0] = direction[0][1] - direction[0][0];
						ZaxisDirection[1] = direction[1][1] - direction[1][0];
						ZaxisDirection[2] = direction[2][1] - direction[2][0];
						
						//calculate the cross product of z and Z to produce the normal vector
						NaxisDirection[0] = -ZaxisDirection[1];
						NaxisDirection[1] = ZaxisDirection[0];
						NaxisDirection[2] = 0.0;
						
						//calculate the angle between z and Z
						beta = angle(0.0,0.0,1.0,ZaxisDirection[0],ZaxisDirection[1],ZaxisDirection[2]);
						//calculate the angle between x and N
						alpha = angle(1.0,0.0,0.0,NaxisDirection[0],NaxisDirection[1],NaxisDirection[2]);
						
						//generate a new point in the frame of the moving particle (Z)
						tFC[0] = L*sin(ThetaTemp)*cos(PhiTemp);
						tFC[1] = L*sin(ThetaTemp)*sin(PhiTemp);
						tFC[2] = L*cos(ThetaTemp);	
						
						//calculate the determinant of the rotation matrix
						det = pow(cos(alpha),2)*pow(cos(beta),2)+pow(sin(alpha),2)+sin(alpha)*cos(alpha)*sin(beta);
						
						//use the inverse rotation matrix to transform the coordinates into the original frame 
						Xtemp = (cos(alpha)*cos(beta)*tFC[0] - sin(alpha)*tFC[1] + cos(alpha)*sin(beta)*tFC[2])*(1.0/det);
						Ytemp = (sin(alpha)*cos(beta)*tFC[0] + cos(alpha)*tFC[1] + sin(alpha)*sin(beta)*tFC[2])*(1.0/det);
						Ztemp = (-sin(beta)*tFC[0] + cos(beta)*tFC[2])*(1.0/det);
						
						//calculate the value of theta in the original frame
						ThetaTemp = atan(Ytemp/Xtemp);
						
						//update the postion of the photon
						X = X + Xtemp;
						Y = Y + Ytemp;
						Z = Z + Ztemp;
						
						//detect if the photon has been scattered or absorbed
						X0 = lcs(X0,a,c,m);
						randX = (double)X0/m;
						if(randX>scatterProb){
							scatterCheck = 1;
						}
						
						//update the direction array to store the two most recent points
						direction[0][0] = direction[0][1];
						direction[1][0] = direction[1][1];
						direction[2][0] = direction[2][1];
						
						direction[0][1] = X;
						direction[1][1] = Y;
						direction[2][1] = Z;
						
					}
				}
			} 
			
			if(scatterCheck == 0){ //if the photon has not been absorbed
				if(Z>Zmax){ //if the photon has escaped out of the positive detection side
			
					//check for anomalous photons that escape with theta>pi/2
					if(ThetaTemp > (M_PI/2.0)){
					
						i--; //set the photon counter back by 1
					} else{
					
						//record the magnitude of the X and Y distance away from the start point
						Xtot = Xtot + fabs(X);
						Ytot = Ytot + fabs(Y);
						count++; //increment the photon counter
					}
					
					//method to bin the detected angles
					//checks if the detected angle is less than each maximum bin value
					//breaks once the angle has been binned
					while(binCheck == 0){
						if(binCounter>(noBins-1)){
							binCheck++;
						} else if(ThetaTemp<=binMax[binCounter]){
							bins[binCounter]++;
							binCheck++;
						}	
						binCounter++;
					}
					
					//reset counters 
					binCheck = 0; 
					binCounter = 0;
			
				} else {
					i--; //if photon returns out of the bottom of the space set the photon counter back by 1
				}
			}
		} 
	}
	
	
	double centralBinAngle, binNEnergy, binError;
	
	printf("Rayleigh Scattering\n");
	printf("Tau: %.2f Zmin: %.2f Zmax: %.2f\n",TauMax,Zmin,Zmax);
	printf("Bin No.   Central Bin Angle   Normalised Energy      Error\n");
	
	//method to output the information 
	for(int i=0; i<noBins; i++){
		
		//calculate the central values of each bin
		if(i == 0){
			centralBinAngle = binMax[i]/2.0;
		} else {
			centralBinAngle = (binMax[i-1]+binMax[i])/2.0;
		}
		
		//calculate and output the bin number, central bin angle, normalised energy and error
		binNEnergy = (double)bins[i]/count;
		binError = binNEnergy/sqrt((double)bins[i]);
		
		printf("%d              %.6f            %.6f         %.6f\n",i+1,centralBinAngle,binNEnergy,binError);
		
		binCounter = binCounter + bins[i];
	} 
	
	printf("Percentage of photons absorbed: %.2f%%\n", 100.0 - (double)count/N * 100.0);
	printf("Average Magnitude X: %.6f Y: %.6f\n", (double)Xtot/count, (double)Ytot/count);

}


int main(){
	
	//declare the seed values used for random number generation 
	int a = 16807;
	int c = 0;
	long m = 2147483647;
	long X = 2;	
	
	int N = 1000000; //number of points used to simulate the probability distributions
	int noBins = 20; //number of bins 
	
	//declare variables to be used
	double out,x0,y0;
	int countCum, countDisc, index;
	countCum = 0;
	countDisc = 0;
	
	int binsCum[noBins];
	int binsDisc[noBins];
	double dx = (double)1/noBins;
	
	for(int i=0; i<noBins; i++){
		binsCum[i] = 0; //clear storage 
		binsDisc[i] = 0;
	}
	
	//compare the discard and cumulative methods for producing random numbers with a distribution of 3/8(1+cos^2x)sinx
	for(int i=0; i<N; i++){
		
		//produce random number between 0 and pi (range of the probability distribution)
		X = lcs(X,a,c,m);
		out = (double)X/m;
		x0 = out*M_PI;
		
		//produce random number between 0 and 6^-0.5 (the top of the probability distribution)
		X = lcs(X,a,c,m);
		out = (double)X/m;
		y0 = out*1.0/sqrt(6);
	
		//if the point (x0,y0) lies below the distribution at x0 then bin the point
		if(y0<((0.375)*(1+pow(cos(x0),2))*sin(x0))){
			
			index = (int)((x0/(M_PI))/dx);
			binsDisc[index] = binsDisc[index]+1;
			countDisc++;
		}
		
		//produce a random number between 0 and pi weighted using the cumulative method
		X = lcs(X,a,c,m);
		out = (double)X/m;
		x0 = thetaHat(out);
		
		//bin the resulting output
		index = (int)((x0/(M_PI))/dx);
		binsCum[index] = binsCum[index]+1;
		countCum++;
	}
	
	//print the comparison of the results 
	printf("Bin No.   Bin Prob. Discard   Bin Prob. Cumulative   \n");
  	
  	for(int i=0; i<noBins; i++){

		printf("%d             %.6f             %.6f\n",i+1,(double)binsDisc[i]/countDisc,(double)binsCum[i]/countCum);
		
	} 
	
	printf("Number of points used in Discard: %d\n", countDisc);
	printf("Discard Efficiency: %.2f%%\n", (double)countDisc/N *100.0);
	printf("Number of points used in Cumulative: %d\n", countCum);
	printf("Cumulative Efficiency: %.2f%%\n", (double)countCum/N *100.0);
	printf("\n");
	
	//double Zmin, double Zmax, double TauMax, int N, int noBins, double scatterProb
	isotropic(0.0,2.0,10.0,1000000,20,1.0);
	printf("\n");
	Rayleigh(0.0,2.0,10.0,1000000,20,1.0);
	printf("\n");
	Rayleigh(0.0,2.0,0.1,1000000,20,1.0);
	
	return 0;
}