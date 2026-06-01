#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string>

using namespace std;

#define M 40
#define N 20
#define Lx 40
#define Ly 20
#define Q 9
#define rho0 1.0

int e[Q][2] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};

int re[Q] = {0,3,4,1,2,7,8,5,6};

double w[Q] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};

// 流动分布函数
double f[M+1][N+1][Q], f_temp[M+1][N+1][Q];
// 温度分布函数
double g[M+1][N+1][Q], g_temp[M+1][N+1][Q];
// 宏观量
double u[M+1][N+1], v[M+1][N+1], rho[M+1][N+1], T[M+1][N+1];
// 标记数组：1=流体，0=固体（原jud）
int jud[M+1][N+1];
// 新增：热对流标记
int SorL[M+1][N+1];
// 新增：基础对流项
double gcol[M+1][N+1][Q];

// 物理参数
double U, Re, Rx, Ry, D, tau_f, tau_T;
int n, nmax;
// 新增：热参数
double Tw = 300.0;   // 壁面温度
double T_in = 280.0; // 入口温度

// 函数声明
void init();
void evolution();
double feq(int k, double u, double v, double rho);
double geq(int k, double u, double v, double T);
int output();

int main() {
    double dx = Lx / M;
    double dt = dx;
    
   
    D = Ly / 7.0;
    Rx = 2.5 * D;
    Ry = 3.5 * D;
    
    cout << "Input inlet velocity U: ";
    cin >> U;
    cout << "Input Reynolds number Re: ";
    cin >> Re;
    
   
    double nu = D * U / Re;
    tau_f = 3.0 * nu / dt + 0.5;
    
   
    double alpha = nu / 0.7;
    tau_T = 3.0 * alpha / dt + 0.5;
    
    cout << "Flow tau_f = " << tau_f << endl;
    cout << "Thermal tau_T = " << tau_T << endl;
    
    init();
    
    cout << "Input max steps nmax: ";
    cin >> nmax;
    
    for (n = 1; n <= nmax; n++) {
        evolution();
        
        if (!(n % 100)) {
            cout << "Step " << n << ": Centerline velocity = " 
                 << u[M/2][N/2] << ", Avg Temp = " << T[M/2][N/2] << endl;
        }
        
        if (!(n % 1000)) {
            output();
        }
    }
    
    output();
    return 0;
}

// 初始化
void init() {
    int i, j, k;
    
    // 初始化标记数组
    for (i = 0; i <= M; i++) {
        for (j = 0; j <= N; j++) {
            double dist = sqrt((i - Rx)*(i - Rx) + (j - Ry)*(j - Ry));
            
            // 圆柱内部为固体
            if (dist <= D/2.0) {
                jud[i][j] = 0;       
                SorL[i][j] = 0;  
                T[i][j] = Tw;       
            } else {
                jud[i][j] = 1;      
                SorL[i][j] = 1;    
                T[i][j] = T_in;     
            }
            
            u[i][j] = (jud[i][j] ? U : 0.0);
            v[i][j] = 0.0;
            rho[i][j] = rho0;
            
          
            for (k = 0; k < Q; k++) {
                f[i][j][k] = feq(k, u[i][j], v[i][j], rho[i][j]);
            }
            
           
            for (k = 0; k < Q; k++) {
                g[i][j][k] = geq(k, u[i][j], v[i][j], T[i][j]);
            }
            
           
            for (k = 0; k < Q; k++) {
                gcol[i][j][k] = 0.1 * w[k] * T[i][j];
            }
        }
    }
}

// 流动平衡分布函数
double feq(int k, double u, double v, double rho) {
    double eu = e[k][0]*u + e[k][1]*v;
    double uv = u*u + v*v;
    return w[k] * rho * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*uv);
}

// 温度平衡分布函数
double geq(int k, double u, double v, double T) {
    double eu = e[k][0]*u + e[k][1]*v;
    return w[k] * T * (1.0 + 3.0*eu);
}


void evolution() {
    int i, j, k;
    int ip, jp;
    
  
    for (i = 0; i <= M; i++) {
        for (j = 0; j <= N; j++) {
            if (jud[i][j]) {
                // 流动碰撞
                for (k = 0; k < Q; k++) {
                    double feqi = feq(k, u[i][j], v[i][j], rho[i][j]);
                    f_temp[i][j][k] = f[i][j][k] - (f[i][j][k] - feqi) / tau_f;
                }
                
                // 温度碰撞（热对流）
                for (k = 0; k < Q; k++) {
                    double geqi = geq(k, u[i][j], v[i][j], T[i][j]);
                    g_temp[i][j][k] = g[i][j][k] - (g[i][j][k] - geqi) / tau_T;
                }
            }
        }
    }
    
  
    for (i = 0; i <= M; i++) {
        for (j = 0; j <= N; j++) {
            if (!jud[i][j]) continue; // 跳过固体
            
            for (k = 0; k < Q; k++) {
                ip = i - e[k][0];
                jp = j - e[k][1];
                // 边界处理（嵌入你提供的热对流逻辑）
                if (jp < 0) { // 下边界
                    if (SorL[i][j]) {
                        // 热对流区：g = -gcol + 2*w*Tw
                        g[i][j][k] = -g_temp[i][j][re[k]] + 2.0 * w[k] * Tw;
                    } else {
                        // 导热区：g = gcol
                        g[i][j][k] = g_temp[i][j][re[k]];
                    }
                    f[i][j][k] = f_temp[i][j][re[k]]; // 流动边界
                    
                } else if (jp > N) { // 上边界
                    g[i][j][k] = g_temp[i][j][re[k]]; // 上边界恒温
                    f[i][j][k] = f_temp[i][j][re[k]];
                    
                } else if (ip < 0) { // 左边界（入口）
                    g[i][j][k] = geq(k, U, 0.0, T_in); // 入口温度
                    f[i][j][k] = feq(k, U, 0.0, rho[i][j]);
                    
                } else if (ip > M) { // 右边界（出口）
                    g[i][j][k] = g_temp[M-1][j][k]; // 零梯度
                    f[i][j][k] = f_temp[M-1][j][k];
                    
                }else if(!jud[ip][jp]){
                    if (SorL[i][j]) {
                        // 热对流区：g = -gcol + 2*w*Tw
                        g[i][j][k] = -g_temp[i][j][re[k]] + 2.0 * w[k] * Tw;
                    } else {
                        // 导热区：g = gcol
                        g[i][j][k] = g_temp[i][j][re[k]];
                    }
                    f[i][j][k] = f_temp[i][j][re[k]];
                } else { // 内部点
                    g[i][j][k] = g_temp[ip][jp][k];
                    f[i][j][k] = f_temp[ip][jp][k];
                }
            }
        }
    }
    
  
    for (i = 0; i <= M; i++) {
        for (j = 0; j <= N; j++) {
            if (!jud[i][j]) continue;
            
            // 流动宏观量
            rho[i][j] = 0.0;
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            for (k = 0; k < Q; k++) {
                rho[i][j] += f[i][j][k];
                u[i][j] += f[i][j][k] * e[k][0];
                v[i][j] += f[i][j][k] * e[k][1];
            }
            u[i][j] /= rho[i][j];
            v[i][j] /= rho[i][j];
            
            // 温度宏观量（热对流结果）
            T[i][j] = 0.0;
            for (k = 0; k < Q; k++) {
                T[i][j] += g[i][j][k];
            }
            if(T[i][j]>=400){
                //cout<<i<<' '<<j<<':'<<T[i][j]<<';';
            }
            //cout<<endl;
        }
    }
    /**/
}


int output() {
    ostringstream filename;
    filename << "flow_thermal_" << n << ".dat";
    ofstream out(filename.str().c_str());
    
    out << "TITLE = \"LBM Flow with Thermal Convection\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"Rho\", \"Temperature\"\n";
    out << "ZONE T=\"Step " << n << "\", I=" << M+1 << ", J=" << N+1 << ", F=POINT\n";
    
    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= N; j++) {
            out << i << " " << j << " "
                << u[i][j] << " " << v[i][j] << " "
                << rho[i][j] << " " << T[i][j] << "\n";
        }
    }
    
    out.close();
    cout << "Output saved: " << filename.str() << endl;
    return 0;
}