#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define M 128//格子 x 轴 数目（格点数）
#define N 128//格子 y 轴 数目 （格点数）（是否需要和H对齐？）
#define Pr 0.71//普朗特数
#define g 0.981//重力加速度
#define H 128//上下边界距离
#define rho0 1.0//初始密度
#define Q 9//离散速度方向数
#define T0 300//初始温度

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};//离散速度方向
int re[Q]={0,3,4,1,2,7,8,5,6};//反弹方向
int rr[M+1][N+1];
double rho[M+1][N+1],u[M+1][N+1],v[M+1][N+1],T[M+1][N+1],G[M+1][N+1];//密度，速度（x,y），温度，局部浮力加速度
double f[M+1][N+1][Q],f1[M+1][N+1][Q],q[M+1][N+1][Q],q1[M+1][N+1][Q];//分布函数，碰撞后分布函数，温度分布函数，碰撞后温度分布函数
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};//权重系数
double T_tau,U_tau,niu,D,Ra,Th,Tc,d_T;//温度松弛时间，速度松弛时间，运动粘度，雷诺数，上边界温度，下边界温度，温差

double U_feq(int k,double U,double V,double RHO);//速度分布函数的平衡态
double T_feq(int k,double U,double V,double T);//温度分布函数的平衡态
void evolution();//演化函数
void output();//输出函数

void main()
{
	printf("Ra=");
	scanf("%lf",&Ra);//输入瑞利数，浮力驱动 vs 粘性/热扩散，描述自然对流强度（与实际情况对齐）
	printf("U_tau=");
	scanf("%lf",&U_tau);//输入碰撞的松弛频率
	//算式中由D2Q9 决定的值的说明：dt=dx=1（转换为格子了）,Cs=1/sqrt(3)（格子声速，分布作为约束会解出声速）
	U_tau=1.0/U_tau;//计算碰撞的松弛时间
	niu=(U_tau-0.5)/3.0;//计算运动粘度 利用 LBM 中的关系：ν = (τ - 0.5) * c_s^2 * dt 
	D=niu/Pr;//计算热扩散系数 α = ν / Pr 
	T_tau=D*3.0+0.5;// 计算温度松弛时间，利用 LBM 中的关系：α = (τ_T - 0.5) * c_s^2 * dt 
	printf("T_tau=%f\n",T_tau);//展示温度松弛时间
	d_T=Ra*T0*(niu*D)/(g*H*H*H);//计算冷热壁温差，利用 Ra = g * β * ΔT * H^3 / (ν * α) 以及 β = 1/T0 的关系
	d_T=d_T/2.0;// 准备计算热边界和冷边界温度
	Th=T0+d_T;//计算热边界温度
	Tc=T0-d_T;//计算冷边界温度
	printf("Th=%f,Tc=%f,niu=%f,D=%f.\n",Th,Tc,niu,D);//展示热边界温度，冷边界温度，运动粘度，热扩散系数

	
	//初始化
	int i,j,k;
	for(i=0;i<=M;i++)
		for(j=0;j<=N;j++)
		{
			u[i][j]=v[i][j]=0.0;//无初始速度
			rho[i][j]=rho0;//初始密度
			T[i][j]=T0;//初始温度
			G[i][j]=0.0;//初始浮力加速度
		}
	for(i=0;i<=M;i++)
		for(j=0;j<=N;j++)
			for(k=0;k<Q;k++)
			{
				f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j]);//计算速度分布函数的初始平衡态
				q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j]);//计算温度分布函数的初始平衡态
			}

	int nmax;
	printf("nmax=");
	scanf("%d",&nmax);//最大时间步数
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

double U_feq(int k,double U,double V,double RHO)//经典计算平衡态分布函数的公式
{
	double eu=e[k][0]*U+e[k][1]*V;
	double uv=U*U+V*V;
	double x=w[k]*RHO*(1.0+3*eu+4.5*eu*eu-1.5*uv);
	return x;
}

double T_feq(int k,double U,double V,double T)//温度分布函数的平衡态计算公式，但是返回的是温度分布函数的平衡态值（温度）
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
				f1[i][j][k]=f[i][j][k]-(f[i][j][k]-FM)/U_tau+3*w[k]*e[k][1]*G[i][j];//计算下一时刻的分布函数，增加了浮力项
				FT=T_feq(k,u[i][j],v[i][j],T[i][j]);
				q1[i][j][k]=q[i][j][k]-(q[i][j][k]-FT)/T_tau;//计算下一时刻的温度分布函数
			}
	for(i=1;i<M;i++)
		for(j=1;j<N;j++)
			for(k=0;k<Q;k++)
			{
				int ip=i-e[k][0];
				int jp=j-e[k][1];
				f[i][j][k]=f1[ip][jp][k];//如果是流体就正常扩散
				q[i][j][k]=q1[ip][jp][k];//温度分布函数正常扩散
			}

	for(j=0;j<=N;j++)//边界处理，这里处理的是i=0和i=M的边界（有温度壁的），j=0和j=N的边界在下面处理
		for(k=0;k<Q;k++)
		{
			i=0;
			u[i][j]=0.0;//流进来的流体没有x轴速度
			v[i][j]=0.0;//也没有y轴速度
			rho[i][j]=rho[i+1][j];//密度和内部流体一样
			f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j])+f[i+1][j][k]-U_feq(k,u[i+1][j],v[i+1][j],rho[i+1][j]);//利用对称性计算分布函数，保证边界处的分布函数满足平衡态的约束

			T[i][j]=Th;//热边界温度
			q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j])+q[i+1][j][k]-T_feq(k,u[i+1][j],v[i+1][j],T[i+1][j]);//利用对称性计算温度分布函数，保证边界处的温度分布函数满足平衡态的约束
			
			i=M;
			u[i][j]=0.0;//一样没有x轴速度
			v[i][j]=0.0;//和y轴速度
			rho[i][j]=rho[i-1][j];//密度和内部流体一样
			f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j])+f[i-1][j][k]-U_feq(k,u[i-1][j],v[i-1][j],rho[i-1][j]);//利用对称性计算分布函数，保证边界处的分布函数满足平衡态的约束

			T[i][j]=Tc;//冷边界温度
			q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j])+q[i-1][j][k]-T_feq(k,u[i-1][j],v[i-1][j],T[i-1][j]);//利用对称性计算温度分布函数，保证边界处的温度分布函数满足平衡态的约束
		}
	for(i=1;i<M;i++)//边界处理，处理j=0和j=N的边界（无温度壁的），这里处理的是y轴边界
		for(k=0;k<Q;k++)
		{
			j=0;
			u[i][j]=0.0;//流进来的流体没有x轴速度
			v[i][j]=0.0;//也没有y轴速度
			rho[i][j]=rho[i][j+1];//密度和内部流体一样
			f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j])+f[i][j+1][k]-U_feq(k,u[i][j+1],v[i][j+1],rho[i][j+1]);//利用对称性计算分布函数，保证边界处的分布函数满足平衡态的约束

			T[i][j]=T[i][j+1];//温度继承内部流体（不是直接取冷热壁了）
			q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j])+q[i][j+1][k]-T_feq(k,u[i][j+1],v[i][j+1],T[i][j+1]);//利用对称性计算温度分布函数，保证边界处的温度分布函数满足平衡态的约束

			j=N;
			u[i][j]=0.0;//流进来的流体没有x轴速度
			v[i][j]=0.0;//也没有y轴速度
			rho[i][j]=rho[i][j-1];//密度和内部流体一样
			f[i][j][k]=U_feq(k,u[i][j],v[i][j],rho[i][j])+f[i][j-1][k]-U_feq(k,u[i][j-1],v[i][j-1],rho[i][j-1]);//利用对称性计算分布函数，保证边界处的分布函数满足平衡态的约束

			T[i][j]=T[i][j-1];//温度继承内部流体（不是直接取冷热壁了）
			q[i][j][k]=T_feq(k,u[i][j],v[i][j],T[i][j])+q[i][j-1][k]-T_feq(k,u[i][j-1],v[i][j-1],T[i][j-1]);//利用对称性计算温度分布函数，保证边界处的温度分布函数满足平衡态的约束
		}

		//统计宏观量
	for(i=1;i<M;i++)
		for(j=1;j<N;j++)
		{
			u[i][j]=0.0;//计算宏观量之前先把速度和温度清零
			v[i][j]=0.0;
			T[i][j]=0.0;
			for(k=0;k<Q;k++)
			{
				u[i][j]+=f[i][j][k]*e[k][0];//x轴速度
				v[i][j]+=f[i][j][j]*e[k][1];//y轴速度
				T[i][j]+=q[i][j][k];//温度
			}
			if(fabs(rho[i][j])>1.0e-9)
				{
					u[i][j]/=rho[i][j];//算加权的
					v[i][j]/=rho[i][j];//算加权的
				}
				else
					u[i][j]=v[i][j]=0.0;//如果密度太小了就直接把速度设为0，避免除以一个很小的数导致数值不稳定（！）
			}
	for(i=0;i<=M;i++)
		for(j=0;j<=N;j++)
			G[i][j]=g*(T[i][j]-T0)/T0;//计算局部浮力加速度
}

void output()//输出函数，输出温度场，速度场，密度场
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