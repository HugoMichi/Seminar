/****************************************************************************
 *                   FAU Erlangen SS14
 *                   Siwir2 Uebung 1 - Elliptic PDE with Multigrid
 *                   written by Michael Hildebrand, Tim Brendel
 *                   Mai 2014
 *****************************************************************************/

#include "header.h"
using namespace std;


int main(int argc, char *argv[]) {

    if (argc != 2) {
        cerr<<"error: wrong number of arguments"<<endl;
        cout<<"call ./mgsolve l"<<endl;
        cout<<"l: number of levels"<<endl;
        exit(EXIT_FAILURE);
    }
    l = atoi(argv[1]);
    if(l<3 || l>17){
        cerr<<"error: arguments out of range"<<endl;
        exit(EXIT_FAILURE);
    }
    NX = NY = pow(2,l)+1;
    H = 2./(NX-1);
    string your_alias= "broetchen_kinder";
    
    
    
    
    
//     int id;
    
     omp_set_num_threads(4);
//       #pragma omp parallel //private (id)
//  	 {
     //for(int i = 0; i<4; ++i){
//        id = omp_get_thread_num();
//        cout<<"ich bin thread no. "<<id<<endl;
//        if(id ==0)
//  	cout<<"es sind "<<omp_get_num_threads()<<" threads aktiv"<<endl;
//      }
//      cout<<endl;
// 	 }
    
	 

    
    
    
    
    

    double* u = new double[NX*NY];
    memset(u,0,sizeof(double)*NY*NX);

    init_polar(u, NX, NY);
//     initSemBD(u);
    save_in_file("init.dat", u, NX, NY);

    double* f = new double[NX*NY];
    memset(f,0,sizeof(double)*NY*NX);

    double* res = new double[NY*NX];
    memset(res,0,sizeof(double)*NY*NX);
    

    std::cout<<"Your Alias: "<<your_alias<<std::endl;
    struct timeval t0, t;
    gettimeofday(&t0, NULL);
    solveMG(u, f, res);
		 //cool smoothing; smooth cooling
// 	 do_gauss_seidel(u,f,NX,NY,1);
    gettimeofday(&t, NULL);
    std::cout << "Wall clock time of MG execution: " <<
    ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t)1000000 +
    (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3
    << " ms" << std::endl;

    save_in_file("solution.dat", u, NX, NY);
    delete[] res;
    delete[] f;
    delete[] u;
    exit(EXIT_SUCCESS);
}


void save_in_file(const char *str, double *matrix, const int n_x, const int n_y){

    ofstream file;
    file.open(str, ios::out);
    if(!(file.is_open())){
        printf("%p konnte nicht gespeichert werden\n", str);
        exit(1);
    }
    double hx_local = 2./(n_x-1);
    double hy_local = 2./(n_y-1);

    //file << setprecision(12);
      for(int yi =0 ; yi <n_y ; ++yi){
         for(int xj =0 ; xj <n_x ; ++xj){

            file << xj*hx_local - 1.0 << ' ';
            file << yi*hy_local - 1.0 << ' ';

            file << matrix[IDX(xj,yi)] << endl;
        }
        //file << endl;
    }
    file.close();
}


//prolongation von grob/coarse nach fein/fine
void prolongation(double *u_co, double *u_fi, const int n_x, const int n_y){

   int Nx_co=(n_x/2)+1;
   int Ny_co=(n_y/2)+1;

   /*
      for(int j = 0; j < Ny_co-1; ++j)
      {
      for(int i = 0; i < Nx_co-1; ++i)
      {
      if (i!=0 && j!=0)
      u_fi[IDX(2*i,2*j)]  	+= u_co[j * Nx_co+ i]; // centre 
      u_fi[IDX(2*i+1,2*j)]	+= 1./2. * (u_co[j * Nx_co+ i] + u_co[j * Nx_co+ i+1]);
      u_fi[IDX(2*i,2*j+1)]    	+= 1./2. * (u_co[j * Nx_co+ i] + u_co[(j+1) * Nx_co+i]);
      u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * (u_co[j * Nx_co+ i] + u_co[(j+1) * Nx_co+i]+ u_co[j * Nx_co+ i+1]+ u_co[(j+1) * Nx_co+ i+1]);
      } 
      }
      */

   //set four courners
   int i=0;
   int j=0;
   double center = u_co[j * Nx_co+ i];
   u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * center;

   i=Nx_co-1;

   center = u_co[j * Nx_co+ i];
   u_fi[IDX(2*i-1,2*j+1)] 	+= 1./4. * center;

   j=Ny_co-1;
   center = u_co[j * Nx_co+ i];
   u_fi[IDX(2*i-1,2*j-1)] 	+= 1./4. * center;

   i=0;
   center = u_co[j * Nx_co+ i];
   u_fi[IDX(2*i+1,2*j-1)] 	+= 1./4. * center;

   //y-border-columns
   for(j = 1; j < Ny_co-1; ++j)
   {
      i=0;
      center = u_co[j * Nx_co+ i];

      u_fi[IDX(2*i+1,2*j)]	+= 1./2. * center;
      u_fi[IDX(2*i+1,2*j-1)] 	+= 1./4. * center;
      u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * center;

      i=Nx_co-1;

      center = u_co[j * Nx_co+ i];

      u_fi[IDX(2*i-1,2*j)]	+= 1./2. * center;
      u_fi[IDX(2*i-1,2*j+1)] 	+= 1./4. * center;
      u_fi[IDX(2*i-1,2*j-1)] 	+= 1./4. * center;
   }

   //x-border-rows
   for(i = 1; i < Nx_co-1; ++i)
   {
      j=0;
      center = u_co[j * Nx_co+ i];


      u_fi[IDX(2*i,2*j+1)]    	+= 1./2. * center;
      u_fi[IDX(2*i-1,2*j+1)] 	+= 1./4. * center;
      u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * center;

      j=Ny_co-1;

      center = u_co[j * Nx_co+ i];

      u_fi[IDX(2*i , 2*j-1)]    	+= 1./2. * center;
      u_fi[IDX(2*i+1,2*j-1)] 	+= 1./4. * center;
      u_fi[IDX(2*i-1,2*j-1)] 	+= 1./4. * center;
   }

   //inner grid
   for(j = 1; j < Ny_co-1; ++j)
   {
      for(i = 1; i < Nx_co-1; ++i)
      {

	 center = u_co[j * Nx_co+ i];
	 u_fi[IDX(2*i,2*j)]  	+= center;

	 u_fi[IDX(2*i+1,2*j)]	+= 1./2. * center;
	 u_fi[IDX(2*i-1,2*j)]	+= 1./2. * center;
	 u_fi[IDX(2*i,2*j+1)]    	+= 1./2. * center;
	 u_fi[IDX(2*i,2*j-1)]    	+= 1./2. * center;

	 u_fi[IDX(2*i-1,2*j+1)] 	+= 1./4. * center;
	 u_fi[IDX(2*i+1,2*j-1)] 	+= 1./4. * center;
	 u_fi[IDX(2*i-1,2*j-1)] 	+= 1./4. * center;
	 u_fi[IDX(2*i+1,2*j+1)] 	+= 1./4. * center;
      } 
   }
}


void do_gauss_seidel(double *u, double *f, const int n_x, const int n_y, const int c){
	
// 	if( n_x <= 257 || n_y <= 257)	omp_set_num_threads(1);

	double h = 1./(n_x-1);
   for(int it=0; it<c; ++it){
/*     //red
       #pragma omp parallel for
      for (int y=1; y<n_y-1; y++)
      {
	 for (int x=(y%2)+1; x<n_x-1; x+=2)
	 {
		if(y==(n_y-1)/2 && x>=(n_x-1)/2){u[IDX(x,y)]=0.0;}
		else	
	        u[IDX(x,y)] = 1.0/4.0 * (h*h*f[IDX(x,y)] + (u[IDX(x,y-1)] + u[IDX(x,y+1)] + u[IDX(x-1,y)] + u[IDX(x+1,y)]));
	 }
      }
      //black
       #pragma omp parallel for
      for (int y=1; y<n_y-1; y++)
      {
        for (int x=((y+1)%2)+1; x<n_x-1; x+=2)
        {
            if(y==(n_y-1)/2 && x>=(n_x-1)/2){u[IDX(x,y)]=0.0;}
		else	
		u[IDX(x,y)] = 1.0/4.0 * (h*h*f[IDX(x,y)] + (u[IDX(x,y-1)] + u[IDX(x,y+1)] + u[IDX(x-1,y)] + u[IDX(x+1,y)]));
        }
      }
*/   
      //all red points of 1st row
      for (int x=2; x<n_x-1; x+=2)
      {
        u[IDX(x,1)] = 1.0/4.0 * (h*h*f[IDX(x,1)] + (u[IDX(x,0)] + u[IDX(x,2)] + u[IDX(x-1,1)] + u[IDX(x+1,1)]));

      }
      //dann
      //bla--red--bla 
      //red--bla--red
      //bla--red--bla
 
      for (int y=2; y<n_y-1; y++)
      {
	 for (int x=(y%2)+1; x<n_x-1; x+=2)
	 {
           if(y==(n_y-1)/2 && x>=(n_x-1)/2){u[IDX(x,y)]=0.0;}
           else{
           //red
           u[IDX(x,y)] = 1.0/4.0 * (h*h*f[IDX(x,y)] + (u[IDX(x,y-1)] + u[IDX(x,y+1)] + u[IDX(x-1,y)] + u[IDX(x+1,y)]));
           }
           if(y==((n_y-1)/2-1) && x>=(n_x-1)/2){u[IDX(x,y-1)]=0.0;}
           else{

           //black
           u[IDX(x,y-1)] = 1.0/4.0 * (h*h*f[IDX(x,y-1)] + (u[IDX(x,y-2)] + u[IDX(x,y)] + u[IDX(x-1,y-1)] + u[IDX(x+1,y-1)]));
          }
        }
      }
      //all black points of last row
      for(int x=1; x<n_x-1; x+=2) //wegen ungeradem gridsize x=1  <=>  x=((y+1)%2)+1
      {
        u[IDX(x,n_y-2)] = 1.0/4.0 * (h*h*f[IDX(x,n_y-2)] + (u[IDX(x,n_y-3)] + u[IDX(x,n_y-1)] + u[IDX(x-1,n_y-2)] + u[IDX(x+1,n_y-2)]));

      }
   }
}

void initSemBD(double* u){
	for(int i = 0; i < NX; ++i){//g(X,Y)= (X^2 + Y^2)^2/3 * sin(1/2atan2(x,y))
//x-richtung
		 u[(NY-1) * NX + i] =pow((i*H-1)*(i*H-1)+(1)*(1),3./2.)* sin(0.5*fabs(atan2(1.,(i*H-1.))) );
		 u[(0) * NX + i] = pow((i*H-1)*(i*H-1)+(-1)*(-1),3./2.)* sin(0.5*fabs(atan2(-1.,(i*H-1.))) );
		
//y-richtung		
		 u[(i) * NX + 0] = pow((-1)*(-1)+(i*H-1)*(i*H-1),3./2.)*sin(0.5*fabs(atan2((i*H-1),-1.)) );
		 u[(i) * NX + NX-1] = pow((1)*(1)+(i*H-1)*(i*H-1),3./2.)* sin(0.5*fabs(atan2((i*H-1),1.)) );
 
    }
	for(int i = NX/2; i < NX; ++i){//g(X,Y)= (X^2 + Y^2)^2/3 * sin(1/2atan2(x,y))
//x-richtung
		 u[(NY/2) * NX + i] = pow((i*H-1)*(i*H-1)+0,3./2.)* sin(0.5*fabs(atan2(0.,(i*H-1.))) );
	


 
    }
}


//do restriction from residual to f_coarse
void restriction(double* f_co, double* res, const int n_x, const int n_y){

    int Nx_co=(n_x/2)+1;
    int Ny_co=(n_y/2)+1;

    //x=0(left border) und x=1(right border)
//     #pragma omp parallel for
    for(int j=0; j<Ny_co;j++){
        f_co[j*Nx_co+0] = res[IDX(0,2*j)];
        f_co[j*Nx_co+Nx_co-1] = res[IDX(n_x-1,2*j)];
    }
    //y=0(lower border) und y=1(upper border)
//     #pragma omp parallel for
    for(int i=0; i<Nx_co;i++){
        f_co[0*Nx_co+i] = res[IDX(2*i,0)];
        f_co[(Ny_co-1)*Nx_co+i] = res[IDX(2*i,n_y-1)];
    }


//     #pragma omp parallel for
    //wird langsamer fuer 4, 8, 16  threads
    for(int j=1;j<Ny_co-1;j++){
        for(int i=1;i<Nx_co-1;i++){
            f_co[j*Nx_co+i] =
                RES_CENTER*res[IDX(2*i,2*j)]+ //restriction stencil
                RES_HORIZONTAL*(res[IDX(2*i-1,2*j)]+ res[IDX(2*i+1,2*j)])+
                RES_VERTICAL*(res[IDX(2*i,2*j-1)]+ res[IDX(2*i,2*j+1)])+
                RES_CORNER*(res[IDX(2*i-1,2*j-1)]+ res[IDX(2*i+1,2*j-1)]+
                        res[IDX(2*i-1,2*j+1)]+ res[IDX(2*i+1,2*j+1)]);
        }
    }
   /*for(int i=(Nx_co-1)/2; i<Nx_co;i++){
        f_co[((Ny_co-1)/2)*Nx_co+i] = res[IDX(2*i,(Ny_co-1)/2)];
    }*/

}


//recursive multigrid function
void mgm(double* u,double* f,int v1,int v2,int n_x, int n_y, int gamma){

    //Pre-smoothing
    do_gauss_seidel(u,f,n_x,n_y,v1);

    double* res = new double[n_y*n_x];
    memset(res,0,sizeof(double)*n_y*n_x);

    //residual calculation
    residuum(res, f, u, n_x, n_y);

    //Coarse grid size calculation
    int Nx_co=(n_x/2)+1;
    int Ny_co=(n_y/2)+1;

    //coarse rhs f
    double* f_co = new double[Ny_co*Nx_co];

    //full weighted restriction
    restriction(f_co,res, n_x, n_y);
    delete[] res;

    //coarse rhs c
    double* c_co = new double[Ny_co*Nx_co];
    memset(c_co,0,sizeof(double)*Ny_co*Nx_co);

    //level==coarsest
    if(Nx_co==3||Ny_co==3)
    {
        double h_co = 1./2.;

        c_co[1 * Nx_co + 1] = ( (h_co*h_co) * f[1 * Nx_co + 1]
                +  u[0 * Nx_co + 1]
                +  u[2 * Nx_co + 1]
                +  u[1 * Nx_co + 2]
                +  u[1 * Nx_co + 0]
                ) / 4.;
    }
    else
    {
        //recursive call
		for(int i = 0; i<gamma; ++i){
			mgm(c_co, f_co, v1, v2, Nx_co ,Ny_co, gamma);
		}
		delete[] f_co;
    }

    prolongation(c_co, u, n_x, n_y);
    delete[] c_co;

    //post-smoothing
    do_gauss_seidel(u, f, n_x, n_y, v2);

}

//calculation of  residual
void residuum(double* res,double* f, double* u, const int n_x,const int n_y){

    double hx_local = 1./(n_x-1);
    double hy_local = 1./(n_y-1);
    
    #pragma omp parallel for
    for(int j=1; j<n_y-1; j++){
        for(int i=1;i<n_x-1;i++){
// 			if(j==(n_y-1)/2 && i>=(n_x-1)/2){
// 				res[IDX(i,j)] = 0.0;
// 			}else{
				res[IDX(i,j)] =  // f-Au
               	  f[IDX(i,j)] -
               	    (1./(hx_local*hy_local))*
               	    (4.*u[IDX(i,j)]
               	     - u[IDX( i ,j-1)]
               	     - u[IDX( i ,j+1)]
               	     - u[IDX(i+1, j )]
               	     - u[IDX(i-1, j )]);
// 			}
		}
    }
#pragma omp parallel for
    for(int i = (n_x-1)/2; i<n_x-1;i++){
		res[IDX(i,(n_y-1)/2)] = 0.;
	}
}


double calcL2Norm(double *res, int n_x, int n_y){

    double norm = 0.;
//     #pragma omp parallel for reduction(+: norm)
    for(int j = 0; j<n_y; ++j){
        for(int i = 0 ; i<n_x; ++i){
            norm += res[IDX(i,j)]*res[IDX(i,j)];
        }
    }
    return sqrt(norm/(n_x*n_y));
}


void init_polar(double *u, const int n_x, const int n_y){

	double h_x = H;
	double h_y = H;
	//quadratisches Grid, also NX == NY
      //  #pragma omp parallel for
	//kein unterschied ob ohne pragma, mit 
	for(int i = 0; i<NX; ++i){  
		//0-te Spalte = linker Rand
		//compute with coordinates in x/y direction
		double x = -1.;
		double y = i * h_y - 1.;
		u[IDX(0,i)] = pow((x*x + y*y), (1./4.)) * sin(polar(x,y)/2.);
		
		//(NX-1)-te Spalte = rechter Rand
		x = (NX-1) * h_x - 1.;
		y = i * h_y -1.;
		u[IDX(NX-1,i)] = pow((x*x + y*y), (1./4.)) * sin(polar(x,y)/2.);

		//0-te Zeile = unterer Rand
		x = i * h_x - 1.;
		y = 0. - 1.;
		u[IDX(i,0)] = pow((x*x + y*y), (1./4.)) * sin(polar(x,y)/2.);
 
		//(NY-1)-te Zeile = oberer Rand
		x = i * h_x - 1.;
		y = (NY-1) * h_y - 1.;
		u[IDX(i, NY-1)] = pow((x*x + y*y), (1./4.)) * sin(polar(x,y)/2.);

		//mittlere Zeile bei x = 0;
		//NX und NY immer ungerade
		u[IDX(NX/2 + i%(NY/2 + 1), NY/2)] = 0.;
		//zweimal reinschreiben oder immer if-Abrage?
	}
}


double polar(const double x, const double y){
	
	double phi = atan2(y,x);
	if(y < 0)
		phi += (2*M_PI);
	return phi;
}

void solveMG(double *u, double *f, double *res){
// 	double l2norm = 1.;
// 	double tol = 9.18e-5;
// 	while(l2norm > tol){
	for(int i = 0; i<7; ++i){
		//multigrid steps
		mgm(u, f, 2, 1, NX, NY,2);
// 		residuum(res, f, u, NX, NY);
		// norm and convergence
// 		l2norm = calcL2Norm(res, NX, NY);
// 		cout<<"L2 Norm: "<<l2norm<<endl;
		//cout<<"Convergence rate: "<< l2norm / l2_old <<endl;
		//l2_old = l2norm;
	}
}
