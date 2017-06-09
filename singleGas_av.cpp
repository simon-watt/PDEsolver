#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<cstring>

using namespace std;

int N=10000;
int size;

struct paramStruct
{
	double dx;
	double Le,f,l,ua,beta,q,r;
};

paramStruct param;

void eulerStep(double x[],double u[],double y[],double &t,double &dt);
void init(double x[],double u[],double y[],double &t);
void output(double x[],double u[],double y[],double t,const char fname[]);
void output(double x[],double u[],double y[],double t,ofstream &out);
#include "triSolve.cpp"

#include "getSpeed.cpp"

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
	double *u=0,*y=0;

	double t0=0.0,t1=4.0e3,dt=0.1,t; // time
	double *x=0,x0=0.0,dx=0.1,x1; // space variable

	double *y_old=0,tOld;

	x1=N*dx;
	size=(N+1)*sizeof(double);

        u=(double *)malloc(size);
        y=(double *)malloc(size);
        x=(double *)malloc(size);
	y_old=(double *)malloc(size);
        if (u==0 || y==0 || x==0 || y_old==0)
        {
                cout << "Not enough memory" << endl;
                exit(0);
        }

        for (int i=0;i<=N;i++)
                x[i]=x0+i*dx;

	init(x,u,y,t);

        output(x,u,y,t,"init_av");

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
		y_old[i]=y[i];
	tOld=t;

	double facP=t+1;
	double facS=t+1;
	int index=1;

	sprintf(fname,"waveSpeed_av.dat");
	ofstream outC(fname);

	double waveLoc;

	ofstream out("trace_av.dat");
        output(x,u,y,t,out);
	
	while (t<t1)
	{
	dt=min(dt,t1-t);
	eulerStep(x,u,y,t,dt);
	if (t>=facP)
	{	
	cout << "start" << endl;
	cout << "dt = " << dt << "\tindex = " << index << "\t";
	output(x,u,y,t,"current_av");
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
		getSpeed(y_old,y,tOld,t,dx,outC,waveLoc);
		facS+=100;
	}
	}
	output(x,u,y,t,"final_av");
        output(x,u,y,t,out);


	outC.close();
	out.close();

	free(u);free(y);free(x);free(y_old);
}	

void init(double x[],double u[],double y[],double &t)
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
		u[i]=A*exp(-pow(x[i]-offset,2)*0.01)+param.ua;
		y[i]=1.0;
	}

}

void output(double x[],double u[],double y[],double t,const char fname[])
{
	char fname2[100];
	sprintf(fname2,"%s.dat",fname);
	ofstream out(fname2);

	out << "#\t t = " << t << endl;
	for (int i=0;i<=N;i++)
	{
		out << x[i] << "\t" << u[i] << "\t" << y[i] << endl;
	}
	out << endl << endl;
	out.close();
	cout << "t = " << t << "\t u0 = " << u[0] << "\ty = " << y[0] << endl;
}

void output(double x[],double u[],double y[],double t,ofstream &out)
{
        out << "#\t t = " << t << endl;
        for (int i=0;i<=N;i++)
        {
                out << x[i] << "\t" << u[i] << "\t" << y[i] << endl;
        }
        out << endl << endl;
}

double arr(double a,double b)
{
	if (b<=0)
		return 0.0;
	else
		return exp(-a/b);
}

double fcn1(double u,double y)
{
	double Le=param.Le,l=param.l,beta=param.beta,ua=param.ua;
        double q=param.q,r=param.r,f=param.f;

	return y*arr(1.0,u)+q*r*y*arr(f,u)-l*(u-ua);
}

double fcn2(double u,double y)
{
	double Le=param.Le,l=param.l,beta=param.beta,ua=param.ua;
        double q=param.q,r=param.r,f=param.f;

	return -beta*y*(arr(1.0,u)+r*arr(f,u));
}
 
void integ(double x[],double u0[],double u1[],double y0[],double y1[],
			double dt)
{
	double Le=param.Le,l=param.l,beta=param.beta,ua=param.ua;
	double q=param.q,r=param.r,f=param.f;
	double Uxx,Yxx;
	double dx=x[1]-x[0];

	double *a1=0,*b1=0,*c1=0,*rhs1=0;
	double *a2=0,*b2=0,*c2=0,*rhs2=0;

	a1=(double *)malloc(size);
	b1=(double *)malloc(size);
	c1=(double *)malloc(size);
	rhs1=(double *)malloc(size);
        a2=(double *)malloc(size);
        b2=(double *)malloc(size);
        c2=(double *)malloc(size);
        rhs2=(double *)malloc(size);

	for (int i=0;i<=N;i++)
	{
		u1[i]=u0[i];
		y1[i]=y0[i];
	}

	for (int i=0;i<=N;i++)
	{
		double temp=max(u0[i],0.0001);
		if (i==0)
		{
			a1[i]=0;
			b1[i]=1.0;
			c1[i]=0.0;
			rhs1[i]=ua;
		}
		else
		if (i==N)
		{
			a1[i]=0.0;
			b1[i]=1.0;
			c1[i]=0.0;
			rhs1[i]=ua;
		}
		else
		{
			a1[i]=-0.5/pow(dx,2);
			b1[i]=1.0/dt+1.0/pow(dx,2);
			c1[i]=-0.5/pow(dx,2);
			Uxx=(u0[i+1]-2.0*u0[i]+u0[i-1])/pow(dx,2);
			rhs1[i]=0.5*Uxx+1.0/dt*u0[i]
				+0.5*(fcn1(u1[i],y1[i])+fcn1(u0[i],y0[i]));
		}
	}
        triSolve(a1,b1,c1,rhs1,u1,N+1);

	for (int i=0;i<=N;i++)
	{
		double temp=max(u0[i],0.0001);
                if (i==0)
                {
			a2[i]=0.0;
                        b2[i]=1.0/dt+1.0/pow(dx,2)/Le;
                        c2[i]=-1.0/pow(dx,2)/Le;
			Yxx=2.0*(y0[1]-y0[0])/pow(dx,2);
			rhs2[i]=0.5*Yxx/Le+1.0/dt*y0[i]
				+0.5*(fcn2(u1[i],y1[i])+fcn2(u0[i],y0[i]));
		}
                else
                if (i==N)
                {
                        a2[i]=0.0;
                        b2[i]=1.0;
			c2[i]=0.0;
			rhs2[i]=1.0;
			Yxx=2.0*(y0[N-1]-y0[N])/pow(dx,2);
		}
		else
		{
			a2[i]=-0.5/pow(dx,2)/Le;
			b2[i]=1.0/dt+1.0/pow(dx,2)/Le;
			c2[i]=-0.5/pow(dx,2)/Le;
			Yxx=(y0[i+1]-2.0*y0[i]+y0[i-1])/pow(dx,2)/Le;
			rhs2[i]=0.5*Yxx/Le+1.0/dt*y0[i]
				+0.5*(fcn2(u1[i],y1[i])+fcn2(u0[i],y0[i]));
		}
        }
        triSolve(a2,b2,c2,rhs2,y1,N+1);

	free(a1);free(b1);free(c1);free(rhs1);
        free(a2);free(b2);free(c2);free(rhs2);

}

void findError(double x0[],double x1[],double y0[],double y1[],double &error)
{
	error=0;
	for (int i=0;i<=N;i++)
		error+=abs(x0[i]-y0[i])+abs(x1[i]-y1[i]);
	error/=N*2.0;
}
	
void eulerStep(double x[],double u[],double y[],double &t,double &dt)
{
/*
	From initial point u(t) calculate u(t+dt)
	in one step Full and two steps Half1 and Half2
	and reduce time step until difference is below eps
*/

	double *uHalf1=0,*uHalf2=0,*uFull=0;
	double *yHalf1=0,*yHalf2=0,*yFull=0;
	double eps=1.0e-6;
	double error;

	uHalf1=(double *)malloc(size);
	uHalf2=(double *)malloc(size);
	uFull=(double *)malloc(size);
        yHalf1=(double *)malloc(size);
        yHalf2=(double *)malloc(size);
        yFull=(double *)malloc(size);

	if (uHalf1==0 || uHalf2==0 || uFull==0 || 
		yHalf1==0 || yHalf2==0 || yFull==0)
	{
		cout << "Not enough memory" << endl;
		exit(0);
	}

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
		u[i]=0.5*(uFull[i]+uHalf2[i]);
                y[i]=0.5*(yFull[i]+yHalf2[i]);
	}	

	free(uHalf1);free(uHalf2);free(uFull);
        free(yHalf1);free(yHalf2);free(yFull);
}

