#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

//Function to solve the Riemann problem using the HLL scheme with a 2nd order time step and the Courant-Friedrichs-Lewy (CFL) stability condition.
//arguments(file name, number of grid spaces, polytropic index,  time to integrate to, CFL factor, left density, left velocity, left pressure, right desnsity, right velocity, right pressure)
void HLL(char fileName[], int N, double gamma, double t_max, double Cfl, double rhoL, double VL, double PL, double rhoR, double VR, double PR){
	
	//allocate the memory to be used throughout the program 	
	double **q = malloc(3 * sizeof(double *));
	double **qnew = malloc(3 * sizeof(double *));
	double **qtemp = malloc(3 * sizeof(double *));
	
	double **f = malloc(3 * sizeof(double *));
	double **fnew = malloc(3 * sizeof(double *));
	double **ftemp = malloc(3 * sizeof(double *));
	
	double **K1 = malloc(3 * sizeof(double *));
	double **K2 = malloc(3 * sizeof(double *));
	
	for(int i = 0; i < 3; i++){
		q[i] = malloc(N * sizeof(double));
		qnew[i] = malloc(N * sizeof(double));
		qtemp[i] = malloc(N * sizeof(double));
		
		f[i] = malloc(2*N * sizeof(double));
		fnew[i] = malloc(2*N * sizeof(double));
		ftemp[i] = malloc(2*N * sizeof(double));
		
		K1[i] = malloc(N * sizeof(double));
		K2[i] = malloc(N * sizeof(double));
	}
	
	double *P = malloc(N * sizeof(double));
	double *V = malloc(N * sizeof(double));
	double *xcoord = malloc(N * sizeof(double));
	
	//assign the variables to be used throughout the program 
	double Sn, SnMax, Pstar, qL, qR, SL, SR, PrL, PrR, CsL, CsR, deltaT;
	double deltaX = (double)1/N;
	double t = 0.0;
	int i = 0; 

	//set up the initial conditions of the problem 
	for(int i = 0;i<N;i++){
		
		xcoord[i] = (double)i/N; //store the x coordinates 
		
		//assign left and right initial conditions 
		if(i<(int)N/2){
		
			q[0][i] = rhoL;
			V[i] = VL;
			P[i] = PL;
		
		} else {
			
			q[0][i] = rhoR;
			V[i] = VR;
			P[i] = PR;
		
		}
		
		q[1][i] = q[0][i]*V[i]; //calculate initial momentum 
		q[2][i] = P[i]/(gamma-1.0) + 0.5*q[0][i]*pow(V[i],2.0); //calculate initial energy 
		
		//calculate initial fluxes
		f[0][2*i] = q[1][i];
		f[1][2*i] = pow(q[1][i],2.0)/q[0][i] + (gamma-1.0)*(q[2][i]-0.5*pow(q[1][i],2.0)/q[0][i]);
		f[2][2*i] = q[1][i]*(q[2][i]+(gamma-1.0)*(q[2][i]-0.5*pow(q[1][i],2.0)/q[0][i]))/q[0][i];
		
		//calculate the approximate wave speed
		Sn = sqrt((gamma/q[0][i])*(gamma-1.0)*(q[2][i]-0.5*pow(q[1][i],2.0)/q[0][i])) + fabs(V[i]);
		
		if(Sn>SnMax){
			SnMax = Sn; //store the maximum wave speed
		}
	}
	
	//use the maximum wave speed to calculate the initial time step as per the CFL condition 
	deltaT = Cfl*deltaX/SnMax;
	
	
	while(t<t_max){ //integrate until the maximum time 
		
		SnMax = 0.0; //reset the max wave speed 
		
		//This section calculates the left and right velocities used in the HLL method for the initial position (i = 0) for a first time evaluation (K1)
		PrL = (gamma-1.0)*(q[2][0]-0.5*pow(q[1][0],2.0)/q[0][0]); //calculate pressure at i=0
		PrR = (gamma-1.0)*(q[2][1]-0.5*pow(q[1][1],2.0)/q[0][1]); //calculate pressure at i=1
		CsL = sqrt(gamma*PrL/q[0][0]); //calculate the speed of sound at i=0
		CsR = sqrt(gamma*PrR/q[0][1]); //calculate the speed of sound at i=1
		
		//use these to calculate p*
		Pstar = 0.5*(PrL+PrR)-0.125*(q[1][1]/q[0][1]-q[1][0]/q[0][0])*(q[0][0]+q[0][1])*(CsL+CsR);
		
		//according to p* and left and right pressures calculate the values qL and qR
		if(Pstar <= PrL){
			qL = 1.0;
		}else{
			qL = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrL - 1.0));
		}
		
		if(Pstar <= PrR){
			qR = 1.0;
		}else{
			qR = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrR - 1.0));
		}
		
		//use pressures and speeds of sound to calculate the left and right velocities 
		SL = q[1][0]/q[0][0] - CsL*qL; 
		SR = q[1][1]/q[0][1] + CsR*qR;
		
		//dependent on the right and left velocities assign the correct fluxes at i=0.5
		if(0.0<=SL){
		
			for(int k=0; k<3; k++){
				f[k][1] = f[k][0];
			}
		
		} else if(SL<=0.0 && 0.0<=SR){
		
			for(int k=0; k<3; k++){
				f[k][1] = (SR*f[k][0] - SL*f[k][2] + SL*SR*(q[k][1] - q[k][0]))/(SR-SL);
			}
	
		} else if(0.0>=SR){
			
			for(int k=0; k<3; k++){
				f[k][1] = f[k][2];
			}
		}
		
		
		//loop through all the positions excluding the start and end point
		for(int j = 1;j<N-1;j++){
		
			//This section calculates the left and right velocities used in the HLL method for the position (i = j) for a first time evaluation (K1)
			PrL = (gamma-1.0)*(q[2][j]-0.5*pow(q[1][j],2.0)/q[0][j]); //calculate pressure at i=j
			PrR = (gamma-1.0)*(q[2][j+1]-0.5*pow(q[1][j+1],2.0)/q[0][j+1]); //calculate pressure at i=j+1
			CsL = sqrt(gamma*PrL/q[0][j]); //calculate the speed of sound at i=j
			CsR = sqrt(gamma*PrR/q[0][j+1]); //calculate the speed of sound at i=j+1
			
			//use these to calculate p*
			Pstar = 0.5*(PrL+PrR)-0.125*(q[1][j+1]/q[0][j+1]-q[1][j]/q[0][j])*(q[0][j]+q[0][j+1])*(CsL+CsR);
		
			//according to p* and left and right pressures calculate the values qL and qR
			if(Pstar <= PrL){
				qL = 1.0;
			}else{
				qL = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrL - 1.0));
			}
		
			if(Pstar <= PrR){
				qR = 1.0;
			}else{
				qR = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrR - 1.0));
			}
			
			//use pressures and speeds of sound to calculate the left and right velocities 
			SL = q[1][j]/q[0][j] - CsL*qL;
			SR = q[1][j+1]/q[0][j+1] + CsR*qR;
			
			//dependent on the right and left velocities assign the correct fluxes at i=j+0.5
			if(0.0<=SL){
				
				for(int k=0; k<3; k++){
					f[k][2*j+1] = f[k][2*j];
				}
				
			} else if(SL<=0.0 && 0.0<=SR){
			
				for(int k=0; k<3; k++){
					f[k][2*j+1] = (SR*f[k][2*j] - SL*f[k][2*j+2] + SL*SR*(q[k][j+1] - q[k][j]))/(SR-SL);
				}
				
			} else if(0.0>=SR){
				
				for(int k=0; k<3; k++){
					f[k][2*j+1] = f[k][2*j+2];
				}
			}
			
			//use the fluxes at i=j-0.5 and i=j+0.5 to calculate a temporary q and a K1 value as per the 2nd order Runge-Kutta time step
			for(int k=0; k<3; k++){
				qtemp[k][j] = q[k][j] - (deltaT/deltaX)*(f[k][2*j+1]-f[k][2*j-1]);
				K1[k][j] = -(deltaT/deltaX)*(f[k][2*j+1]-f[k][2*j-1]);
			}
			
			//use the temporary q values to calculate temporary fluxes at integar values of i
			ftemp[0][2*j] = qtemp[1][j];
			ftemp[1][2*j] = pow(qtemp[1][j],2.0)/qtemp[0][j] + (gamma-1.0)*(qtemp[2][j]-0.5*pow(qtemp[1][j],2.0)/qtemp[0][j]);
			ftemp[2][2*j] = qtemp[1][j]*(qtemp[2][j]+(gamma-1.0)*(qtemp[2][j]-(0.5*pow(qtemp[1][j],2.0)/qtemp[0][j])))/qtemp[0][j];
			
			//calculate the approximate wave speed
			Sn = sqrt((gamma/qtemp[0][j])*(gamma-1.0)*(qtemp[2][j]-0.5*pow(qtemp[1][j],2.0)/qtemp[0][j])) + fabs(qtemp[1][j]/qtemp[0][j]);
			
			//store the maximum wave speed 
			if(Sn>SnMax){
				SnMax = Sn;
			}
		}
		
		//apply outflow boundary conditions 
		for(int k=0; k<3; k++){
			qtemp[k][0] = q[k][1];
			qtemp[k][N-1] = q[k][N-2];
			
			K1[k][0] = qtemp[k][0] - q[k][0];
			K1[k][N-1] = qtemp[k][N-1] - q[k][N-1];
			
			ftemp[k][0] = f[k][2];
			ftemp[k][2*N-2] = f[k][2*N-4];
		}
		
		
		//This section calculates the left and right velocities used in the HLL method for the initial position (i = 0) for a second time evaluation (K2)
		PrL = (gamma-1.0)*(qtemp[2][0]-0.5*pow(qtemp[1][0],2.0)/qtemp[0][0]); //calculate pressure at i=0
		PrR = (gamma-1.0)*(qtemp[2][1]-0.5*pow(qtemp[1][1],2.0)/qtemp[0][1]); //calculate pressure at i=1
		CsL = sqrt(gamma*PrL/qtemp[0][0]); //calculate the speed of sound at i=0
		CsR = sqrt(gamma*PrR/qtemp[0][1]); //calculate the speed of sound at i=1
		
		//use these to calculate p*
		Pstar = 0.5*(PrL+PrR)-0.125*(qtemp[1][1]/qtemp[0][1]-qtemp[1][0]/qtemp[0][0])*(qtemp[0][0]+qtemp[0][1])*(CsL+CsR);
		
		//according to p* and left and right pressures calculate the values qL and qR
		if(Pstar <= PrL){
			qL = 1.0;
		}else{
			qL = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrL - 1.0));
		}
		
		if(Pstar <= PrR){
			qR = 1.0;
		}else{
			qR = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrR - 1.0));
		}
		
		//use pressures and speeds of sound to calculate the left and right velocities 
		SL = qtemp[1][0]/qtemp[0][0] - CsL*qL;
		SR = qtemp[1][1]/qtemp[0][1] + CsR*qR;
		
		//dependent on the right and left velocities assign the correct fluxes at i=0.5
		if(0.0<=SL){
		
			for(int k=0; k<3; k++){
				ftemp[k][1] = ftemp[k][0];
			}
		
		} else if(SL<=0.0 && 0.0<=SR){
		
			for(int k=0; k<3; k++){
				ftemp[k][1] = (SR*ftemp[k][0] - SL*ftemp[k][2] + SL*SR*(qtemp[k][1] - qtemp[k][0]))/(SR-SL);
			}
	
		} else if(0.0>=SR){
			
			for(int k=0; k<3; k++){
				ftemp[k][1] = ftemp[k][2];
			}
		}
		
		//loop through all the positions excluding the start and end point
		for(int j = 1;j<N-1;j++){
			
			//This section calculates the left and right velocities used in the HLL method for the position (i = j) for a second time evaluation (K2)
			PrL = (gamma-1.0)*(qtemp[2][j]-0.5*pow(qtemp[1][j],2.0)/qtemp[0][j]); //calculate pressure at i=j
			PrR = (gamma-1.0)*(qtemp[2][j+1]-0.5*pow(qtemp[1][j+1],2.0)/qtemp[0][j+1]); //calculate pressure at i=j+1
			CsL = sqrt(gamma*PrL/qtemp[0][j]); //calculate the speed of sound at i=j
			CsR = sqrt(gamma*PrR/qtemp[0][j+1]); //calculate the speed of sound at i=j+1
			
			//use these to calculate p*
			Pstar = 0.5*(PrL+PrR)-0.125*(qtemp[1][j+1]/qtemp[0][j+1]-qtemp[1][j]/qtemp[0][j])*(qtemp[0][j]+qtemp[0][j+1])*(CsL+CsR);
			
			//according to p* and left and right pressures calculate the values qL and qR
			if(Pstar <= PrL){
				qL = 1.0;
			}else{
				qL = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrL - 1.0));
			}
		
			if(Pstar <= PrR){
				qR = 1.0;
			}else{
				qR = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrR - 1.0));
			}
			
			//use pressures and speeds of sound to calculate the left and right velocities 
			SL = qtemp[1][j]/qtemp[0][j] - CsL*qL;
			SR = qtemp[1][j+1]/qtemp[0][j+1] + CsR*qR;
			
			//dependent on the right and left velocities assign the correct fluxes at i=j+0.5
			if(0.0<=SL){
				
				for(int k=0; k<3; k++){
					ftemp[k][2*j+1] = ftemp[k][2*j];
				}
				
			} else if(SL<=0.0 && 0.0<=SR){
			
				for(int k=0; k<3; k++){
					ftemp[k][2*j+1] = (SR*ftemp[k][2*j] - SL*ftemp[k][2*j+2] + SL*SR*(qtemp[k][j+1] - qtemp[k][j]))/(SR-SL);
				}
				
			} else if(0.0>=SR){
				
				for(int k=0; k<3; k++){
					ftemp[k][2*j+1] = ftemp[k][2*j+2];
				}
			}
			
			
			//use the fluxes at i=j-0.5 and i=j+0.5 to calculate a K2 value 
			//use the K1 and K2 values to calculate q at a time t+deltaT
			for(int k=0; k<3; k++){
				K2[k][j] = -(deltaT/deltaX)*(ftemp[k][2*j+1]-ftemp[k][2*j-1]);
				qnew[k][j] = q[k][j] + 0.5*(K1[k][j] + K2[k][j]);
			}
			
			//use the new q values to calculate the new fluxes at integar values of i
			fnew[0][2*j] = qnew[1][j];
			fnew[1][2*j] = pow(qnew[1][j],2.0)/qnew[0][j] + (gamma-1.0)*(qnew[2][j]-0.5*pow(qnew[1][j],2.0)/qnew[0][j]);
			fnew[2][2*j] = qnew[1][j]*(qnew[2][j]+(gamma-1.0)*(qnew[2][j]-(0.5*pow(qnew[1][j],2.0)/qnew[0][j])))/qnew[0][j];
			
			//calculate the approximate wave speed
			Sn = sqrt((gamma/qnew[0][j])*(gamma-1.0)*(qnew[2][j]-0.5*pow(qnew[1][j],2.0)/qnew[0][j])) + fabs(qnew[1][j]/qnew[0][j]);
			
			//store the maximum wave speed
			if(Sn>SnMax){
				SnMax = Sn;
			}
		}
		
		//apply outflow boundary conditions and copy the new q values into the original storage for the next iteration
		for(int k=0; k<3; k++){
			qnew[k][0] = q[k][1];
			qnew[k][N-1] = q[k][N-2];
			
			fnew[k][0] = f[k][2];
			fnew[k][2*N-2] = f[k][2*N-4];
			
			memcpy(q[k],qnew[k],N*sizeof(double));
			memcpy(f[k],fnew[k],2*N*sizeof(double));
		}
		
		
		t = t + deltaT; //move forward in time by delta T
		deltaT = Cfl*deltaX/SnMax; //calculate the new delta T according to the CFL condition
		i++;
		
		//printf("%.6f %d\n",t,i); //optional time and iteration print out, useful for long run-times
	} 
	

	FILE *fp;
	
	//reset files upon each run of the program
	fp = fopen(fileName,"w");
	
	//print out values to the console and to a text file for plotting
	printf("X coord   Density   Velocity   Pressure   Energy\n");
	
	for(int i = 0;i<N;i++){
		printf("%.4f    %.6f  %.6f   %.6f   %.6f\n",xcoord[i],q[0][i],q[1][i]/q[0][i],(gamma-1.0)*(q[2][i]-0.5*pow(q[1][i],2)/q[0][i]),(q[2][i]-0.5*pow(q[1][i],2)/q[0][i])/q[0][i]);
		fp = fopen(fileName,"a");
		fprintf(fp,"%.6f %.6f %.6f %.6f %.6f\n",xcoord[i],q[0][i],q[1][i]/q[0][i],(gamma-1.0)*(q[2][i]-0.5*pow(q[1][i],2)/q[0][i]),(q[2][i]-0.5*pow(q[1][i],2)/q[0][i])/q[0][i]);
		fclose(fp);	
	} 
	
	//free memory used throughout the program
	for(int i = 0; i < 3; i++){
		free(q[i]); 
		free(qnew[i]); 
		free(qtemp[i]); 
		
		free(f[i]); 
		free(fnew[i]);
		free(ftemp[i]); 
		
		free(K1[i]); 
		free(K2[i]); 
	}
	free(q);
	free(qnew);
	free(qtemp);
	free(f);
	free(fnew);
	free(ftemp);
	free(K1);
	free(K2);
	
	free(P);
	free(V);
	free(xcoord);
}

//Function to solve the Riemann problem in 2D using the HLL scheme with the Courant-Friedrichs-Lewy (CFL) stability condition.
//For the problem we are solving the initial conditions are split into top half and bottom half and hence the arguments for the function are likewise 
//The initial condition arguments could easily be changed to left and right or into four quadrants dependent on the problem 
//arguments(file name, number of gird spaces per axis, ploytropic index, time to integrate to, CFL factor, bottom density, upper density, bottom x velocity, upper x velocity, bottom pressure, upper pressure)
void HLL2D(char fileName[], int N, double gamma, double t_max, double Cfl, double rhoD, double rhoU, double VxiD, double VxiU, double PressureD, double PressureU){
	
	//allocate the memory to be used throughout the program 
	double ***q = malloc(4 * sizeof(double **));
	double ***qnew = malloc(4 * sizeof(double **));
	
	double ***fx = malloc(4 * sizeof(double **));
	double ***fxnew = malloc(4 * sizeof(double **));
	
	double ***fy = malloc(4 * sizeof(double **));
	double ***fynew = malloc(4 * sizeof(double **)); 
	
	double **P = malloc(N * sizeof(double *));
	double **Vx = malloc(N * sizeof(double *));
	double **Vy = malloc(N * sizeof(double *));
	double **xcoord = malloc(N * sizeof(double *));
	double **ycoord = malloc(N * sizeof(double *));
	
	for(int i = 0; i < N; i++){
		P[i] = malloc(N * sizeof(double));
		Vx[i] = malloc(N * sizeof(double));
		Vy[i] = malloc(N * sizeof(double));
		xcoord[i] = malloc(N * sizeof(double));
		ycoord[i] = malloc(N * sizeof(double));
	}
	
	for(int i = 0; i < 4; i++){
		q[i] = malloc(N * sizeof(double *));
		qnew[i] = malloc(N * sizeof(double *));
		
		fy[i] = malloc(N * sizeof(double *));
		fynew[i] = malloc(N * sizeof(double *));
		
		for(int j = 0; j < N; j++){
			q[i][j] = malloc(N * sizeof(double));
			qnew[i][j] = malloc(N * sizeof(double));
			
			fy[i][j] = malloc(2 * N * sizeof(double));
			fynew[i][j] = malloc(2 * N * sizeof(double));
		}
	}
	
	for(int i = 0; i < 4; i++){
		fx[i] = malloc(2 * N * sizeof(double *));
		fxnew[i] = malloc(2 * N * sizeof(double *));
		
		for(int j = 0; j < 2 * N; j++){
			fx[i][j] = malloc(N * sizeof(double));
			fxnew[i][j] = malloc(N * sizeof(double));
		}
	}
	
	//assign variables to be used throughout the program  
	double Sn, SnMax, Vtot, Pstar, qL, qR, SL, SR, PrL, PrR, CsL, CsR, deltaT;
	double deltaX = (double)1/N;
	double t = 0.0;
	int i = 0; 

	//set up the initial conditions of the problem 
	for(int i = 0;i<N;i++){
	
		for(int j = 0; j<N; j++){
		
			xcoord[i][j] = (double)i/N; //store the x coordinates
			ycoord[i][j] = (double)j/N; //store the y coordinates
			
			//assign the up and down initial conditions
			if(j<(int)N/2){
		
				q[0][i][j] = rhoD;
				Vx[i][j] = VxiD;
				P[i][j] = PressureD;
		
			} else {
			
				q[0][i][j] = rhoU;
				Vx[i][j] = VxiU;
				P[i][j] = PressureU;
			}
			
			
			Vy[i][j] = 0.1*sin((2*M_PI*i)/N); //calculate the initial y velocities
			
			Vtot = sqrt(pow(Vy[i][j],2.0)+pow(Vx[i][j],2.0)); //calculate the total velocity
		
			
			q[1][i][j] = q[0][i][j]*Vx[i][j]; //calculate initial x momentum
			q[2][i][j] = q[0][i][j]*Vy[i][j]; //calculate initial y momentum 
			q[3][i][j] = P[i][j]/(gamma-1.0) + 0.5*q[0][i][j]*pow(Vtot,2.0);//calculate initial energy
		
			//calculate initial x fluxes
			fx[0][2*i][j] = q[0][i][j]*Vx[i][j];
			fx[1][2*i][j] = pow(q[1][i][j],2.0)/q[0][i][j] + P[i][j];
			fx[2][2*i][j] = q[0][i][j]*Vx[i][j]*Vy[i][j];
			fx[3][2*i][j] = (P[i][j]/(gamma-1.0) + 0.5*q[0][i][j]*pow(Vtot,2.0) + P[i][j])*Vx[i][j];
			
			//calculate initial y fluxes 
			fy[0][i][2*j] = q[0][i][j]*Vy[i][j];
			fy[1][i][2*j] = q[0][i][j]*Vx[i][j]*Vy[i][j];
			fy[2][i][2*j] = pow(q[2][i][j],2.0)/q[0][i][j] + P[i][j];
			fy[3][i][2*j] = (P[i][j]/(gamma-1.0) + 0.5*q[0][i][j]*pow(Vtot,2.0) + P[i][j])*Vy[i][j];
			
			//calculate the approximate wave speen 
			Sn = sqrt(gamma*P[i][j]/q[0][i][j]) + fabs(Vtot);
		
			if(Sn>SnMax){
				SnMax = Sn; //store the maximum wave speed
			}
		}
	}
	
	//use the maximum wave speed to calculate the initial time step as per the CFL condition 
	deltaT = Cfl*deltaX/SnMax;


	while(t<t_max){ //integrate until the maximum time 
		
		
		SnMax = 0.0; //reset the max wave speed
		
		//first iteration across the grid considers the fluxes in the x direction
		for(int j = 0; j<N; j++){
		
			//This section calculates the left and right velocities used in the HLL method for the initial position (x = 0)
			PrL = (gamma-1.0)*(q[3][0][j]-0.5*pow(q[1][0][j],2.0)/q[0][0][j]); //calculate the pressure at x=0, y=j 
			PrR = (gamma-1.0)*(q[3][1][j]-0.5*pow(q[1][1][j],2.0)/q[0][1][j]); //calculate the pressure at x=1, y=j 
			CsL = sqrt(gamma*PrL/q[0][0][j]); //calculate the speed of sound at x=0, y=j
			CsR = sqrt(gamma*PrR/q[0][1][j]); //calculate the speed of sound at x=1, y=j
			
			//use these to calculate p*
			Pstar = 0.5*(PrL+PrR)-0.125*(q[1][1][j]/q[0][1][j]-q[1][0][j]/q[0][0][j])*(q[0][0][j]+q[0][1][j])*(CsL+CsR);
			
			//according to p* and left and right pressures calculate the values qL and qR
			if(Pstar <= PrL){
				qL = 1.0;
			}else{
				qL = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrL - 1.0));
			}
		
			if(Pstar <= PrR){
				qR = 1.0;
			}else{
				qR = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrR - 1.0));
			}	
			
			//use pressures and speeds of sound to calculate the left and right velocities 
			SL = q[1][0][j]/q[0][0][j] - CsL*qL;
			SR = q[1][1][j]/q[0][1][j] + CsR*qR;
			
			//dependent on the right and left velocities assign the correct fluxes at x=0.5, y=j
			if(0.0<=SL){
		
				for(int k=0; k<4; k++){
					fx[k][1][j] = fx[k][0][j];
				}
		
			} else if(SL<=0.0 && 0.0<=SR){
		
				for(int k=0; k<4; k++){
					fx[k][1][j] = (SR*fx[k][0][j] - SL*fx[k][2][j] + SL*SR*(q[k][1][j] - q[k][0][j]))/(SR-SL);
				}
	
			} else if(0.0>=SR){
			
				for(int k=0; k<4; k++){
					fx[k][1][j] = fx[k][2][j];
				}
			}
		
		
			//loop through all the x positions excluding the start and end point
			for(int i = 1;i<N-1;i++){
				
				//This section calculates the left and right velocities used in the HLL method for the position x=i, y=j
				PrL = (gamma-1.0)*(q[3][i][j]-0.5*pow(q[1][i][j],2.0)/q[0][i][j]); //calculate pressure at x=i, y=j
				PrR = (gamma-1.0)*(q[3][i+1][j]-0.5*pow(q[1][i+1][j],2.0)/q[0][i+1][j]); //calculate the pressure at x=i+1, y=j
				CsL = sqrt(gamma*PrL/q[0][i][j]); //calculate the speed of sound at x=i, y=j
				CsR = sqrt(gamma*PrR/q[0][i+1][j]); //calculate the speed of sound at x=i+1, y=j
		
				//use these to calculate p*
				Pstar = 0.5*(PrL+PrR)-0.125*(q[1][i+1][j]/q[0][i+1][j]-q[1][i][j]/q[0][i][j])*(q[0][i][j]+q[0][i+1][j])*(CsL+CsR);
				
				//according to p* and left and right pressures calculate the values qL and qR
				if(Pstar <= PrL){
					qL = 1.0;
				}else{
					qL = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrL - 1.0));
				}
		
				if(Pstar <= PrR){
					qR = 1.0;
				}else{
					qR = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrR - 1.0));
				}
				
				//use pressures and speeds of sound to calculate the left and right velocities 
				SL = q[1][i][j]/q[0][i][j] - CsL*qL;
				SR = q[1][i+1][j]/q[0][i+1][j] + CsR*qR;
				
				//dependent on the right and left velocities assign the correct fluxes at x=i+0.5, y=j
				if(0.0<=SL){
				
					for(int k=0; k<4; k++){
						fx[k][2*i+1][j] = fx[k][2*i][j];
					}
				
				} else if(SL<=0.0 && 0.0<=SR){
			
					for(int k=0; k<4; k++){
						fx[k][2*i+1][j] = (SR*fx[k][2*i][j] - SL*fx[k][2*i+2][j] + SL*SR*(q[k][i+1][j] - q[k][i][j]))/(SR-SL);
					}
				
				} else if(0.0>=SR){
				
					for(int k=0; k<4; k++){
						fx[k][2*i+1][j] = fx[k][2*i+2][j];
					}
				}
			
				//use the fluxes at x=i-0.5, y=j and x=i+0.5, y=j to calculate a new q state 
				for(int k=0; k<4; k++){
					qnew[k][i][j] = q[k][i][j] - (deltaT/deltaX)*(fx[k][2*i+1][j]-fx[k][2*i-1][j]);
				}
				
				//calculate total velocity
				Vtot = sqrt((pow(qnew[1][i][j],2.0)+pow(qnew[2][i][j],2.0))/pow(qnew[0][i][j],2.0));
				
				
				//calculate the approximate wave speed
				Sn = sqrt(gamma*(gamma-1.0)*(qnew[3][i][j]-0.5*qnew[0][i][j]*pow(Vtot,2.0))/qnew[0][i][j]) + fabs(Vtot);
			
				if(Sn>SnMax){
					SnMax = Sn; //store the maximum wave speed 
				}
			}
		
			//apply periodic boundary conditions in the x direction
			for(int k=0; k<4; k++){
				qnew[k][0][j] = q[k][N-2][j];
				qnew[k][N-1][j] = q[k][1][j];
			
				fxnew[k][0][j] = fx[k][2*N-4][j];
				fxnew[k][2*N-2][j] = fx[k][2][j];
			}
		}
		
		
		//second iteration across the grid considers the fluxes in the y direction using x-flux state as an input
		for(int i = 0; i<N; i++){
			
			//This section calculates the up and down velocities used in the HLL method for the initial position (y = 0)
			PrL = (gamma-1.0)*(qnew[3][i][0]-0.5*pow(qnew[2][i][0],2.0)/qnew[0][i][0]); //calculate the pressure at x=i, y=0
			PrR = (gamma-1.0)*(qnew[3][i][1]-0.5*pow(qnew[2][i][1],2.0)/qnew[0][i][1]); //calculate the pressure at x=i, y=1
			CsL = sqrt(gamma*PrL/qnew[0][i][0]); //calculate the speed of sound at x=i, y=0
			CsR = sqrt(gamma*PrR/qnew[0][i][1]); //calculate the speed of sound at x=i, y=1
			
			//use these to calculate p*
			Pstar = 0.5*(PrL+PrR)-0.125*(qnew[2][i][1]/qnew[0][i][1]-qnew[2][i][0]/qnew[0][i][0])*(qnew[0][i][0]+qnew[0][i][1])*(CsL+CsR);
		
			//according to p* and up and down pressures calculate the values qL and qR
			if(Pstar <= PrL){
				qL = 1.0;
			}else{
				qL = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrL - 1.0));
			}
		
			if(Pstar <= PrR){
				qR = 1.0;
			}else{
				qR = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrR - 1.0));
			}	
			
			//use pressures and speeds of sound to calculate the up and down velocities
			SL = qnew[2][i][0]/qnew[0][i][0] - CsL*qL;
			SR = qnew[2][i][1]/qnew[0][i][1] + CsR*qR;
			
			//dependent on the up and down velocities assign the correct fluxes at x=i, y=0.5
			if(0.0<=SL){
		
				for(int k=0; k<4; k++){
					fy[k][i][1] = fy[k][i][0];
				}
		
			} else if(SL<=0.0 && 0.0<=SR){
		
				for(int k=0; k<4; k++){
					fy[k][i][1] = (SR*fy[k][i][0] - SL*fy[k][i][2] + SL*SR*(qnew[k][i][1] - qnew[k][i][0]))/(SR-SL);
				}
	
			} else if(0.0>=SR){
			
				for(int k=0; k<4; k++){
					fy[k][i][1] = fy[k][i][2];
				}
			}
		
			//loop through all the y positions excluding the start and end point
			for(int j = 1;j<N-1;j++){
		
				//This section calculates the up and down velocities used in the HLL method for the position x=i, y=j
				PrL = (gamma-1.0)*(qnew[3][i][j]-0.5*pow(qnew[2][i][j],2.0)/qnew[0][i][j]); //calculate pressure at x=i, y=j
				PrR = (gamma-1.0)*(qnew[3][i][j+1]-0.5*pow(qnew[2][i][j+1],2.0)/qnew[0][i][j+1]); //calculate pressure at x=i, y=j+1
				CsL = sqrt(gamma*PrL/qnew[0][i][j]); //calculate the speed of sound at x=i, y=j
				CsR = sqrt(gamma*PrR/qnew[0][i][j+1]); //calculate the speed of sound at x=i, y=j+1
		
				//use these to calculate p*
				Pstar = 0.5*(PrL+PrR)-0.125*(qnew[2][i][j+1]/qnew[0][i][j+1]-qnew[2][i][j]/qnew[0][i][j])*(qnew[0][i][j]+qnew[0][i][j+1])*(CsL+CsR);
			
			//according to p* and up and down pressures calculate the values qL and qR
				if(Pstar <= PrL){
					qL = 1.0;
				}else{
					qL = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrL - 1.0));
				}
		
				if(Pstar <= PrR){
					qR = 1.0;
				}else{
					qR = sqrt(1.0+(gamma+1.0)/(2.0*gamma)*(Pstar/PrR - 1.0));
				}
				
				//use pressures and speeds of sound to calculate the up and down velocities 
				SL = qnew[2][i][j]/qnew[0][i][j] - CsL*qL;
				SR = qnew[2][i][j+1]/qnew[0][i][j+1] + CsR*qR;
				
				//dependent on the up and down velocities assign the correct fluxes at x=i, y=j+0.5
				if(0.0<=SL){
				
					for(int k=0; k<4; k++){
						fy[k][i][2*j+1] = fy[k][i][2*j];
					}
				
				} else if(SL<=0.0 && 0.0<=SR){
			
					for(int k=0; k<4; k++){
						fy[k][i][2*j+1] = (SR*fy[k][i][2*j] - SL*fy[k][i][2*j+2] + SL*SR*(qnew[k][i][j+1] - qnew[k][i][j]))/(SR-SL);
					}
				
				} else if(0.0>=SR){
				
					for(int k=0; k<4; k++){
						fy[k][i][2*j+1] = fy[k][i][2*j+2];
					}
				}
			
				//use the fluxes at x=i, y=j-0.5 and x=i, y=j+0.5 to calculate the new q state 
				for(int k=0; k<4; k++){
					q[k][i][j] = qnew[k][i][j] - (deltaT/deltaX)*(fy[k][i][2*j+1]-fy[k][i][2*j-1]);
				}
				
				//calculate total velocity
				Vtot = sqrt((pow(q[1][i][j],2.0)+pow(q[2][i][j],2.0))/pow(q[0][i][j],2.0));
				
				//use the updated state to calculate the new y fluxes 
				fxnew[0][2*i][j] = qnew[1][i][j];
				fxnew[1][2*i][j] = pow(qnew[1][i][j],2.0)/qnew[0][i][j] + (gamma-1.0)*(qnew[3][i][j]-0.5*qnew[0][i][j]*pow(Vtot,2.0));
				fxnew[2][2*i][j] = qnew[1][i][j]*qnew[2][i][j]/qnew[0][i][j];
				fxnew[3][2*i][j] = (qnew[3][i][j] + (gamma-1.0)*(qnew[3][i][j]-0.5*qnew[0][i][j]*pow(Vtot,2.0)))*(qnew[1][i][j]/qnew[0][i][j]);
				
				//use the updated state to calculate the new y fluxes 
				fynew[0][i][2*j] = q[2][i][j];
				fynew[1][i][2*j] = q[1][i][j]*q[2][i][j]/q[0][i][j];
				fynew[2][i][2*j] = pow(q[2][i][j],2.0)/q[0][i][j] + (gamma-1.0)*(q[3][i][j]-0.5*q[0][i][j]*pow(Vtot,2.0));
				fynew[3][i][2*j] = (q[3][i][j] + (gamma-1.0)*(q[3][i][j]-0.5*q[0][i][j]*pow(Vtot,2.0)))*(q[2][i][j]/q[0][i][j]);
				
				//calculate the approximate wave speed
				Sn = sqrt(gamma*(gamma-1.0)*(q[3][i][j]-0.5*q[0][i][j]*pow(Vtot,2.0))/q[0][i][j]) + fabs(Vtot);
			
				if(Sn>SnMax){
					SnMax = Sn; //store the maximum wave speed 
				}
			}
		
			//apply outflow boundary conditions in the y direction
			for(int k=0; k<4; k++){
				q[k][i][0] = qnew[k][i][1];
				q[k][i][N-1] = qnew[k][i][N-2];
			
				fynew[k][i][0] = fy[k][i][2];
				fynew[k][i][2*N-2] = fy[k][i][2*N-4];
			}
			
			//copy the new x fluxes over to the original memory storage for the next iteration 
			for(int i = 0; i < 4; i++){
				for(int j = 0; j < 2 * N; j++){
					memcpy(fx[i][j],fxnew[i][j], N*sizeof(double));
				}
			}
			
			//copy the new y fluxes over to the original memory storage for the next iteration 
			for(int i = 0; i < 4; i++){
				for(int j = 0; j < N; j++){
					memcpy(fy[i][j],fynew[i][j], 2*N*sizeof(double));
				}
			}
		}
		
		t = t + deltaT; //move forward in time by delta T
		deltaT = Cfl*deltaX/SnMax; //calculate the new delta T according to the CFL condition
		i++;
		
		printf("%.6f %d\n",t,i); //optional time and iteration print out, useful for long run-times
	} 
	
	FILE *fp;
	
	//reset files upon each run of the program
	fp = fopen(fileName,"w");
	
	//output to console and to a text file for plotting 
	printf("X coord   Y coord   Density\n");
	
	for(int i = 0;i<N;i++){
		for(int j = 0;j<N;j++){
		
			printf("%.4f    %.4f  %.6f   \n",xcoord[i][j],ycoord[i][j],q[0][i][j]);
			fp = fopen(fileName,"a");
			fprintf(fp,"%.8f    %.8f  %.8f   \n",xcoord[i][j],ycoord[i][j],q[0][i][j]);
			fclose(fp);	
		}
	}  
	
	//free all memory used throughout the program 
	for(int i = 0; i < N; i++){
		free(P[i]);
		free(Vx[i]); 
		free(Vy[i]);
		free(xcoord[i]);
		free(ycoord[i]);
	}
	free(P);
	free(Vx);
	free(Vy);
	free(ycoord);
	free(xcoord);
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < N; j++){
			free(q[i][j]); 
			free(qnew[i][j]); 
			free(fy[i][j]); 
			free(fynew[i][j]); 
		}
		free(q[i]); 
		free(qnew[i]); 
		free(fy[i]); 
		free(fynew[i]); 
	}
	free(q);
	free(qnew);
	free(fy);
	free(fynew);
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 2 * N; j++){
			free(fx[i][j]); 
			free(fxnew[i][j]); 
		}
		free(fx[i]); 
		free(fxnew[i]);
	}
	free(fx);
	free(fxnew);
}

int main(){
	
	//Produce the data for question 1
	HLL("HLLSetupA.txt",100,1.4,0.25,0.3,1.0,0.0,1.0,0.125,0.0,0.1);
	HLL("HLLSetupB.txt",100,1.4,0.15,0.3,1.0,-2.0,0.4,1.0,2.0,0.4);
	
	//Produce the data for question 2
	//The figure in the answers used N=300 grid points which takes approximately 40mins
	HLL2D("HLL2D300.txt",300, 1.4, 3.0, 0.3, 2.0,1.0,0.5,-0.5,2.5,2.5);
	return 0;
}
