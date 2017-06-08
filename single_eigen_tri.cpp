#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<cstring>
#include<eigen3/Eigen/Sparse>
#include<eigen3/Eigen/SparseLU>
using namespace std;
using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> Trip;

int N=10000;
int size;

struct paramStruct
{
	double dx;
	double Le,f,l,ua,beta,q,r;
};

paramStruct param;

void eulerStep(VectorXd x,VectorXd &u,VectorXd &y,double &t,double &dt);
void init(VectorXd x,VectorXd &u,VectorXd &y,double &t);
void output(VectorXd x,VectorXd u,VectorXd y,double t,const char fname[]);
void output(VectorXd x,VectorXd u,VectorXd y,double t,ofstream &out);
#include "triSolve.cpp"

#include "getSpeed_eigen.cpp"

double a2f(char str[])
{
        double num=0;
        int decLoc=-1;

        for (int i=0;i<strlen(str);i++)
                if (str[i]=='.')
                        decLoc=i;
                else
                        num=num*10+(str[i]-'0');

	if (decLoc==-1) // no decimal point
		return num;
	else
        	return num/pow(10.0,double(strlen(str)-decLoc-1));
}

int a2i(char str[])
{
        int num=0;

        for (int i=0;i<strlen(str);i++)
                        num=num*10+(str[i]-'0');

        return num;
}

int main(int argc,char ** argv)
{
	VectorXd u(N+1),y(N+1),x(N+1),yOld(N+1);

	double t0=0.0,t1=4.0e3,dt=0.1,t; // time
	double x0=0.0,dx=0.1,x1; // space variable

	double tOld;

	x1=N*dx;

        for (int i=0;i<=N;i++)
                x(i)=x0+i*dx;

	init(x,u,y,t);

        output(x,u,y,t,"init_eigen_tri");

        char fname[100];
/*
	string line;
	ifstream in;
	in.open("current.dat");
	if (in) // previous profile
	{
		getline(in,line); // get header
		cout << line << endl;
		for (int i=0;i<=N;i++)
			in >> x[i] >> u[i] >> y[i];
	}
	in.close();
*/
			

	param.dx=dx;
	t=t0;

	for (int i=0;i<=N;i++)
		yOld(i)=y(i);
	tOld=t;

	double facP=t+1;
	double facS=t+1;
	int index=1;

	sprintf(fname,"waveSpeed_eigen_tri.dat");
	ofstream outC(fname);

	double waveLoc;

	ofstream out("trace_eigen_tri.dat");
        output(x,u,y,t,out);
	
	while (t<t1)
	{
	dt=min(dt,t1-t);
	eulerStep(x,u,y,t,dt);
	if (t>=facP)
	{	
	cout << "start" << endl;
	cout << "dt = " << dt << "\tindex = " << index << "\t";
	output(x,u,y,t,"current_eigen_tri");
        output(x,u,y,t,out);
	cout << "end" << endl;
	if (facP<1000)
		facP=facP*2;
	else
		facP=facP+1000;
	index++;
	}
	if (t>=facS)
	{
		getSpeed(yOld,y,tOld,t,dx,outC,waveLoc);
		facS+=100;
	}
	}
	output(x,u,y,t,"final_eigen_tri");
        output(x,u,y,t,out);


	outC.close();
	out.close();

}	

void init(VectorXd x,VectorXd &u,VectorXd &y,double &t)
{
//      initialisation

	t=0.0;
        param.Le=2;
	param.l=5e-4;
	param.ua=1e-4;
	param.beta=3.25;
	param.f=3;
	param.q=5;
	param.r=25;

	double A=0.2,offset=100;
	for (int i=0;i<=N;i++)
	{
		u(i)=A*exp(-pow(x(i)-offset,2)*0.01)+param.ua;
		y(i)=1.0;
	}

}

void output(VectorXd x,VectorXd u,VectorXd y,double t,const char fname[])
{
	char fname2[100];
	sprintf(fname2,"%s.dat",fname);
	ofstream out(fname2);

	out << "#\t t = " << t << endl;
	for (int i=0;i<=N;i++)
	{
		out << x(i) << "\t" << u(i) << "\t" << y(i) << endl;
	}
	out << endl << endl;
	out.close();
	cout << "t = " << t << "\t u0 = " << u(0) << "\ty0 = " << y(0) << endl;
}

void output(VectorXd x,VectorXd u,VectorXd y,double t,ofstream &out)
{
        out << "#\t t = " << t << endl;
        for (int i=0;i<=N;i++)
        {
                out << x(i) << "\t" << u(i) << "\t" << y(i) << endl;
        }
        out << endl << endl;
}

void integ(VectorXd x,VectorXd u0,VectorXd &u1,VectorXd y0,VectorXd &y1,
			double dt)
{
	double Le=param.Le,l=param.l,beta=param.beta,ua=param.ua;
	double q=param.q,r=param.r,f=param.f;
	double Uxx,Yxx,rhs1,rhs2,rhs3,rhs4;
	double dx=x[1]-x[0];
	SpMat A(N+1,N+1);
	VectorXd rhs(N+1);
	vector<Trip> T;
	SparseLU <SpMat> solver;

	T.clear();A.setZero();
	for (int i=0;i<=N;i++)
	{
		double temp=max(u0[i],0.0001);
		if (i==0)
		{
			T.push_back(Trip(i,i,1.0));
			rhs(i)=ua;
		}
		else
		if (i==N)
		{
			T.push_back(Trip(i,i,1.0));
			rhs(i)=ua;
		}
		else
		{
			T.push_back(Trip(i,i-1,-0.5/pow(dx,2)));
			T.push_back(Trip(i,i,1.0/dt+1.0/pow(dx,2)));
			T.push_back(Trip(i,i+1,-0.5/pow(dx,2)));
			Uxx=(u0[i+1]-2.0*u0[i]+u0[i-1])/pow(dx,2);
			rhs(i)=0.5*Uxx+1.0/dt*u0[i]+y0[i]*exp(-1.0/temp)
				+q*r*y0[i]*exp(-f/u0[i])-l*(u0[i]-ua);
		}
	}
	A.setFromTriplets(T.begin(),T.end());
	solver.analyzePattern(A);
	solver.factorize(A);
	u1=solver.solve(rhs);

	T.clear();A.setZero();	
	for (int i=0;i<=N;i++)
	{
		double temp=max(u0[i],0.0001);
                if (i==0)
                {
			T.push_back(Trip(i,i,1.0/dt+1.0/pow(dx,2)/Le));
			T.push_back(Trip(i,i+1,-1.0/pow(dx,2)/Le));
			Yxx=2.0*(y0[1]-y0[0])/pow(dx,2);
			rhs(i)=0.5*Yxx/Le+1.0/dt*y0[i]
				-beta*y0[i]*(exp(-1.0/temp)+r*exp(-f/temp));
		}
                else
                if (i==N)
                {
			T.push_back(Trip(i,i,1.0));
			rhs(i)=1.0;
		}
		else
		{
			T.push_back(Trip(i,i-1,-0.5/pow(dx,2)/Le));
			T.push_back(Trip(i,i,1.0/dt+1.0/pow(dx,2)/Le));
			T.push_back(Trip(i,i+1,-0.5/pow(dx,2)/Le));
			Yxx=(y0[i+1]-2.0*y0[i]+y0[i-1])/pow(dx,2)/Le;
			rhs(i)=0.5*Yxx/Le+1.0/dt*y0[i]
				-beta*y0[i]*(exp(-1.0/temp)+r*exp(-f/temp));
		}
        }
	A.setFromTriplets(T.begin(),T.end());
	solver.analyzePattern(A);
	solver.factorize(A);
	y1=solver.solve(rhs);

}

void findError(VectorXd x0,VectorXd x1,VectorXd y0,VectorXd y1,double &error)
{
	error=0;
	for (int i=0;i<=N;i++)
		error+=abs(x0(i)-y0(i))+abs(x1(i)-y1(i));
	error/=N*2.0;
}
	
void eulerStep(VectorXd x,VectorXd &u,VectorXd &y,double &t,double &dt)
{
/*
	From initial point u(t) calculate u(t+dt)
	in one step Full and two steps Half1 and Half2
	and reduce time step until difference is below eps
*/

        VectorXd uHalf1(N+1),uHalf2(N+1),uFull(N+1);
        VectorXd yHalf1(N+1),yHalf2(N+1),yFull(N+1);
	double eps=1.0e-6;
	double error;

	do
	{
		integ(x,u,uHalf1,y,yHalf1,0.5*dt);
		integ(x,uHalf1,uHalf2,yHalf1,yHalf2,0.5*dt);
		integ(x,u,uFull,y,yFull,dt);
		findError(uFull,yFull,uHalf2,yHalf2,error);
		dt*=0.5;
	} while (error>eps);

	t+=dt*2;
//	dt=min(dt*4,1.0);
	dt=dt*4;

//	cout << "t = " << t << endl;

	for (int i=0;i<=N;i++)
	{
		u(i)=0.5*(uFull(i)+uHalf2(i));
                y(i)=0.5*(yFull(i)+yHalf2(i));
	}	

}

