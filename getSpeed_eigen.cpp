void getSpeed(VectorXd &Cold,VectorXd C,VectorXd &xOld,VectorXd x,
		double &tOld,double t,double dx,ofstream &out)
{
	double *f=0;
	double sum1=0,sum2=0,sum;
	f=(double *)malloc((N+1)*sizeof(double));
	for (int i=0;i<=N;i++)
		f[i]=abs(Cold(i))/abs(t-tOld);

	sum1=f[0]+f[N];
	for (int i=1;i<=N/2-1;i++)
		sum1+=2*f[2*i];
	for (int i=1;i<=N/2;i++)
		sum1+=4*f[2*i-1];

	sum1=0.0;
	for (int i=1;i<=N;i++)
		sum1+=0.5*(x[i]-x[i-1])*(f[i]+f[i-1]);

        for (int i=0;i<=N;i++)
        	f[i]=abs(C(i))/abs(t-tOld);

        sum2=f[0]+f[N];
        for (int i=1;i<=N/2-1;i++)
                sum2+=2*f[2*i];
        for (int i=1;i<=N/2;i++)
                sum2+=4*f[2*i-1];

	sum2=0.0;
	for (int i=1;i<=N;i++)
		sum2+=0.5*(x[i]-x[i-1])*(f[i]+f[i-1]);
	
	sum=abs(sum1-sum2);

	out << t << " " << sum << endl;
	cout << "t = " << t << " speed = " << sum << endl;

	for (int i=0;i<=N;i++)
	{
		Cold(i)=C(i);
		xOld[i]=x(i);
	}
	tOld=t;
	free(f);
}

void getSpeed(VectorXd &Cold,VectorXd C,
		double &tOld,double t,double dx,ofstream &out, double &waveLoc)
{
	double *f=0;
	double sum1=0,sum2=0,sum;
	f=(double *)malloc((N+1)*sizeof(double));
	for (int i=0;i<=N;i++)
		f[i]=abs(Cold(i))/abs(t-tOld);

	sum1=0.0;
	for (int i=1;i<=N;i++)
		sum1+=0.5*dx*(f[i]+f[i-1]);

        for (int i=0;i<=N;i++)
        	f[i]=abs(C(i))/abs(t-tOld);

	sum2=0.0;
	for (int i=1;i<=N;i++)
		sum2+=0.5*dx*(f[i]+f[i-1]);
	
	sum=sum1-sum2;

	out << t << " " << sum << endl;
	cout << "t = " << t << " speed = " << sum << 
		" front loc = " << N*dx-sum1*(t-tOld) << endl;

	waveLoc=N*dx-sum1*(t-tOld);

	for (int i=0;i<=N;i++)
	{
		Cold(i)=C(i);
	}
	tOld=t;
	free(f);
}

