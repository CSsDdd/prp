#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string>

using namespace std;

#define M 420//格点数
#define N 210//格点数
#define Lx 420//物理长度
#define Ly 210//物理长度
#define Q 9//D2Q9模型
#define rho0 1.0//开始时格点总密度

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};//速度方向矢量
int re[Q]={0,3,4,1,2,7,8,5,6};//索引
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};//权重（似乎根据模型选取计算得来）
double f[M+1][N+1][Q],g[M+1][N+1][Q],u[M+1][N+1],v[M+1][N+1],rho[M+1][N+1];//f为分布函数，g为碰撞后分布函数，u,v为速度分量，rho为密度
int jud[M+1][N+1];//判断格点是否为障碍物，1为流体，0为障碍物
void evolution();//推进
void init();//初始化
double feq(int k,double u,double v,double rho);//平衡态分布函数，k为微观速度方向索引（0-8），u,v为宏观速度分量，rho为格点总密度
int output();//输出数据

double U,Re,Rx,Ry,D,tau;//宏观总速度，雷诺数，圆心坐标，圆半径，弛豫时间（碰撞用）
int n,nmax;//当前时间步数，总时间步数
int main()
{
	double dx=Lx/M;
	double dt=dx;
	D=Ly/7.0;//计算直径
	Rx=2.5*D;//计算圆心位置
	Ry=3.5*D;//计算圆心位置
	
	cout<<"Input U:";
	cin>>U;//初始x轴速度
	cout<<"Input Re:";
	cin>>Re;
	double niu=D*U/Re;
	tau=3*niu/dt+0.5;
	cout<<"tau="<<tau<<endl;//计算弛豫时间
	
	init();
	cout<<"nmax=";
	cin>>nmax;//总时长
	laber://重新开始标签
	for(n=1;n<=nmax;n++)
	{
		evolution();
		if(!(n%100))
			cout<<"The U of 10D is "<<u[10*int(D)][N/2]<<",n="<<n<<endl;
        if(!(n%1000))
            output();
	}
	output();
	cout<<"Go on?"<<endl<<"nmax=";
	cin>>nmax;
	if(nmax)
		goto laber;//返回上面的laber
	return 0;
}

void init()
{
	int i,j,k;
	for(i=0;i<=M;i++)
		for(j=0;j<=N;j++)
		{
			double RR=(i-Rx)*(i-Rx)+(j-Ry)*(j-Ry);//计算距离用的，R²
			if(RR<=D*D/4.0)//柱子
			{
				jud[i][j]=0;
				u[i][j]=0.0;
			}
			else//流体
			{
				jud[i][j]=1;
				u[i][j]=U;
			}
			v[i][j]=0.0;//初始没有y轴速度
			rho[i][j]=rho0;//初始每个格点的总密度都是1
			
			for(k=0;k<Q;k++)
				f[i][j][k]=feq(k,u[i][j],v[i][j],rho[i][j]);//计算每个格点每个方向的分布函数初值
		}
}

double feq(int k,double u,double v,double rho)// 数学计算
{
	double uv,eu,x;
  	eu=e[k][0]*u+e[k][1]*v;
  	uv=u*u+v*v;
  	x=w[k]*rho+w[k]*rho0*(3.0*eu+4.5*eu*eu-1.5*uv);
  	return x;
}

//修改1：增加了遍历
//修改2：增加了注释
//修改3：宏观速度改用了加权平均，之前的版本是直接加和的？？？
void evolution()
{
	int i,j,k;
	for(i=0;i<=M;i++)//碰撞
		for(j=0;j<=N;j++)
		{
			if(jud[i][j])//如果是流体在计算
			{
				for(k=0;k<Q;k++)
				{
				double FM=feq(k,u[i][j],v[i][j],rho[i][j]);//期望的平衡态分布函数
				g[i][j][k]=f[i][j][k]-(f[i][j][k]-FM)/tau;//先算碰撞
				}
			}
		}
	
	for(i=1;i<M;i++)//扩散
		for(j=1;j<N;j++)
		{
			if(jud[i][j])
			{			
				for(k=0;k<Q;k++)
				{
				int ip=i-e[k][0];
				int jp=j-e[k][1];
				if(jud[ip][jp])
					f[i][j][k]=g[ip][jp][k];//如果是流体就正常扩散
				else
					f[i][j][k]=g[i][j][re[k]];//如果是障碍物就反弹回去
				}
			}
		}
		
	for(i=1;i<M;i++)
		for(j=1;j<N;j++)
		{
			if(jud[i][j])
			{
				u[i][j]=v[i][j]=rho[i][j]=0.0;
				for(k=0;k<Q;k++)
				{
					u[i][j]+=f[i][j][k]*e[k][0];//计算宏观x轴速度，这里速度为什么是加和而不是加权平均啊
					v[i][j]+=f[i][j][k]*e[k][1];//计算宏观y轴速度
					rho[i][j]+=f[i][j][k];//计算宏观密度
				}
				u[i][j]/=rho[i][j];//计算平均速度
				v[i][j]/=rho[i][j];//
			}
		}
		
	for(i=1;i<M;i++)//处理y轴边界
	{
		u[i][0]=u[i][1];
		v[i][0]=0.0;
		rho[i][0]=rho[i][1];
		u[i][N]=u[i][N-1];
		v[i][N]=0.0;
		rho[i][N]=rho[i][N-1];
		for(k=0;k<Q;k++)
		{
			f[i][0][k]=feq(k,u[i][0],v[i][0],rho[i][0])+f[i][1][k]-feq(k,u[i][1],v[i][1],rho[i][1]);
			f[i][N][k]=feq(k,u[i][N],v[i][N],rho[i][N])+f[i][N-1][k]-feq(k,u[i][N-1],v[i][N-1],rho[i][N-1]);
		}
	}
	for(j=0;j<=N;j++)//处理x轴边界
	{
		u[0][j]=U;
		v[0][j]=0.0;
		rho[0][j]=rho[1][j];
		u[M][j]=u[M-1][j];
		v[M][j]=v[M-1][j];
		rho[M][j]=rho[M-1][j];
		for(k=0;k<Q;k++)
		{
			f[0][j][k]=feq(k,u[0][j],v[0][j],rho[0][j])+f[1][j][k]-feq(k,u[1][j],v[1][j],rho[1][j]);
			f[M][j][k]=f[M-1][j][k];
		}
	}
}

int output()
{
    int i,j;
	ostringstream name;
	const int a[] = { 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000 };
	int b[8];
	for (i = 0; i < 8; i++)
	{
		b[i] = n % a[i] / (a[i] / 10);
	}
	name << "tecplot"
		<< b[7] << b[6] << b[5] << b[4] << b[3] << b[2] << b[1] << b[0]
		<< ".dat";
	ofstream out(name.str().c_str());
	out << "TITLE = \"Tecplot Data\"" << "\n"
		<< "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\"\n"
		<< "ZONE T= \"Zone " << n / 1000 << "\", J=" << M+1 << ", I=" << N+1 << ", F=POINT\n";
	for (i = 0; i < M+1; i++)
		for (j = 0; j < N+1; j++)
		{
			out << i << "\t"
				<< j << "\t"
				<< rho[i][j] << "\t"
				<< u[i][j] << "\t"
				<< v[i][j] << "\t"
				<< endl;
		}
	return 0;
}
