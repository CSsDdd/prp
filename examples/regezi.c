#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define M 128
#define N 128
#define Pr 0.71
#define g 0.981
#define H 128
#define rho0 1.0
#define Q 9
#define T0 300

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
int re[Q]={0,3,4,1,2,7,8,5,6};
int rr[M+1][N+1];
double rho[M+1][N+1],u[M+1][N+1],v[M+1][N+1],T[M+1][N+1],G[M+1][N+1];
double f[M+1][N+1][Q],f1[M+1][N+1][Q],q[M+1][N+1][Q],q1[M+1][N+1][Q];
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double T_tau,U_tau,niu,D,Ra,Th,Tc,d_T;

double U_feq(int k,double U,double V,double RHO);
double T_feq(int k,double U,double V,double T);
void evolution();
void output();

void main()
{
	printf("Ra=");
	scanf("%lf",&Ra);
	printf("U_tau=");
	scanf("%lf",&U_tau);
	U_tau=1.0/U_tau;
	niu=(U_tau-0.5)/3.0;
	D=niu/Pr;
	T_tau=D*3.0+0.5;
	printf("T_tau=%f\n",T_tau);
	d_T=Ra*T0*(niu*D)/(g*H*H*H);
	d_T=d_T/2.0;
	Th=T0+d_T;
	Tc=T0-d_T;
	printf("Th=%f,Tc=%f,niu=%f,D=%f.\n",Th,Tc,niu,D);

	
	//初始化
	int i,j,k;
	for(i=0;i<=M;i++)
		for(j=0;j<=N;j++)
		{
			u[i][j]=v[i][j]=0.0;
			rho[i][j]=rho0;
			T[i][j]=T0;
			G[i][j]=0.0;
		}
	for(i=0;i<=M;i++)
		for(j=0;j<=N;j++)
			for(k=0;k<Q;k++)
			{
				f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j]);
				q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j]);
			}

	int nmax;
	printf("nmax=");
	scanf("%d",&nmax);
	laber:;
	for(int n=1;n<=nmax;n++)
	{
		evolution();
		if(!(n%100))
		{
			printf("The conner U is %f, V is %f, T is %f, n=%d.\n",u[M/4][N/4],v[M/4][N/4],T[M/4][N/4],n);
		}
		output();
	}
	printf("Go on?\nnmax=");
	scanf("%d",&nmax);
	if(nmax)
		goto laber;
}

double U_feq(int k,double U,double V,double RHO)
{
	double eu=e[k][0]*U+e[k][1]*V;
	double uv=U*U+V*V;
	double x=w[k]*RHO*(1.0+3*eu+4.5*eu*eu-1.5*uv);
	return x;
}

double T_feq(int k,double U,double V,double T)
{
	double eu=e[k][0]*U+e[k][1]*V;
	double x=w[k]*T*(1.0+3*eu);
	return x;
}

void evolution()
{
	int i,j,k;
	double FM,FT;

	for(i=0;i<=M;i++)
		for(j=0;j<=N;j++)
			for(k=0;k<Q;k++)
			{
				FM=U_feq(k,u[i][j],v[i][j],rho[i][j]);
				f1[i][j][k]=f[i][j][k]-(f[i][j][k]-FM)/U_tau+3*w[k]*e[k][1]*G[i][j];
				FT=T_feq(k,u[i][j],v[i][j],T[i][j]);
				q1[i][j][k]=q[i][j][k]-(q[i][j][k]-FT)/T_tau;
			}
	for(i=1;i<M;i++)
		for(j=1;j<N;j++)
			for(k=0;k<Q;k++)
			{
				int ip=i-e[k][0];
				int jp=j-e[k][1];
				f[i][j][k]=f1[ip][jp][k];
				q[i][j][k]=q1[ip][jp][k];
			}

//边界条件
	for(j=0;j<=N;j++)
		for(k=0;k<Q;k++)
		{
			i=0;
			u[i][j]=0.0;
			v[i][j]=0.0;
			rho[i][j]=rho[i+1][j];
			f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j])+f[i+1][j][k]-U_feq(k,u[i+1][j],v[i+1][j],rho[i+1][j]);

			T[i][j]=Th;
			q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j])+q[i+1][j][k]-T_feq(k,u[i+1][j],v[i+1][j],T[i+1][j]);
			
			i=M;
			u[i][j]=0.0;
			v[i][j]=0.0;
			rho[i][j]=rho[i-1][j];
			f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j])+f[i-1][j][k]-U_feq(k,u[i-1][j],v[i-1][j],rho[i-1][j]);

			T[i][j]=Tc;
			q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j])+q[i-1][j][k]-T_feq(k,u[i-1][j],v[i-1][j],T[i-1][j]);
		}
	for(i=1;i<M;i++)
		for(k=0;k<Q;k++)
		{
			j=0;
			u[i][j]=0.0;
			v[i][j]=0.0;
			rho[i][j]=rho[i][j+1];
			f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j])+f[i][j+1][k]-U_feq(k,u[i][j+1],v[i][j+1],rho[i][j+1]);

			T[i][j]=T[i][j+1];
			q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j])+q[i][j+1][k]-T_feq(k,u[i][j+1],v[i][j+1],T[i][j+1]);

			j=N;
			u[i][j]=0.0;
			v[i][j]=0.0;
			rho[i][j]=rho[i][j-1];
			f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j])+f[i][j-1][k]-U_feq(k,u[i][j-1],v[i][j-1],rho[i][j-1]);

			T[i][j]=T[i][j-1];
			q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j])+q[i][j-1][k]-T_feq(k,u[i][j-1],v[i][j-1],T[i][j-1]);
		}

		//计算宏观量
	for(i=1;i<M;i++)
		for(j=1;j<N;j++)
		{
			u[i][j]=0.0;
			v[i][j]=0.0;
			T[i][j]=0.0;
			for(k=0;k<Q;k++)
			{
				u[i][j]+=f[i][j][k]*e[k][0];
				v[i][j]+=f[i][j][j]*e[k][1];
				T[i][j]+=q[i][j][k];
			}
			if(fabs(rho[i][j])>1.0e-9)
				{
					u[i][j]/=rho[i][j];
					v[i][j]/=rho[i][j];
				}
				else
					u[i][j]=v[i][j]=0.0;
			}
	for(i=0;i<=M;i++)
		for(j=0;j<=N;j++)
			G[i][j]=g*(T[i][j]-T0)/T0;
}

void output()
{
	FILE *fp;
	int i,j;
	if((fp=fopen("T.dat","w"))==NULL)
	{
		printf("File Open Error!\n");
		exit(1);
	}
	for(j=0;j<=N;j++)
	{
		for(i=0;i<=M;i++)
			fprintf(fp,"%e\t",T[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	if((fp=fopen("u.dat","w"))==NULL)
	{
		printf("File Open Error!\n");
		exit(1);
	}
	for(j=0;j<=N;j++)
	{
		for(i=0;i<=M;i++)
			fprintf(fp,"%e\t",u[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	if((fp=fopen("v.dat","w"))==NULL)
	{
		printf("File Open Error!\n");
		exit(1);
	}
	for(j=0;j<=N;j++)
	{
		for(i=0;i<=M;i++)
			fprintf(fp,"%e\t",v[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	if((fp=fopen("rho.dat","w"))==NULL)
	{
		printf("File Open Error!\n");
		exit(1);
	}
	for(j=0;j<=N;j++)
	{
		for(i=0;i<=M;i++)
			fprintf(fp,"%e\t",rho[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
}