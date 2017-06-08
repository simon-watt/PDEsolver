void triSolve(double a[],double b[],double c[],double r[],double u[],int N)
{
	double beta,*gam;

	gam=(double *)malloc(N*sizeof(double));
	beta=b[0];
	u[0]=r[0]/beta;
	for (int j=1;j<N;j++)
        {
                gam[j]=c[j-1]/beta;
                beta=b[j]-a[j]*gam[j];
                u[j]=(r[j]-a[j]*u[j-1])/beta;
        }

        for (int j=N-2;j>=0;j--)
                u[j]-=gam[j+1]*u[j+1];

	free(gam);
}
