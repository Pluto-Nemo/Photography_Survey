#include<iostream>
#include <fstream>
#include <iomanip>
#include<string.h>
#include<math.h>
using namespace std;
double** caculateRotateMat(double** Mat, double a, double b, double c);
double** MatrixMutiply(double** Mat1, double** Mat2, int r, int n, int c);
double** MatrixTransport(double** Mat, int r, int c);
double** MatrixConverse(double** Mat, int a);
double** MatrixRemain(double** Mat, int a, int x, int y);
double MatrixDet(double** Mat, int a);

int main(){
    double M_1 = 40000;//1/m = 40000
    double f = 0.15324;//f = 153.24mm

    int pointNumber;
    double tmpxi, tmpyi, tmpXO, tmpYO, tmpZO;
    ifstream infile;

    infile.open("data1.txt", ios::in);

    if(!infile.is_open()){
        cout << "Open file failed!" << endl;
        system("pause");
    }
    infile >> pointNumber;
    double** imageControl = new double* [100]();//开辟图像控制点存储空间
    double** objectControl = new double* [100]();//开辟地物控制点存储空间
    for(int i = 0; i < pointNumber; i++){
        infile >> tmpxi >> tmpyi >> tmpXO >> tmpYO >> tmpZO;
        imageControl[i] = new double[2]();
        imageControl[i][0] = tmpxi/1000.0;
        imageControl[i][1] = tmpyi/1000.0;//换算单位mm->m
        objectControl[i] = new double[3]();
        objectControl[i][0] = tmpXO;
        objectControl[i][1] = tmpYO;
        objectControl[i][2] = tmpZO;
    }//为控制点数据赋值

    //到这里所有的测量数据和已知数据都已经存在内存中_input data part done!

    double Xs = 0.0, Ys = 0.0, Zs = 0.0, Phai = 0.0, Omega = 0.0, Ka = 0.0;
    //开辟外方位元素参数空间
    double** delta = new double* [6]();
    for(int i = 0; i < 6; i++){
        delta[i] = new double[1]();
    }
    //开辟外方位元素的改正数的空间
    double** Rotate;
    //开辟旋转矩阵的空间
    for(int i = 0; i < pointNumber; i++){
        Xs += objectControl[i][0];
        Ys += objectControl[i][1];
    }
    Xs /= pointNumber;
    Ys /= pointNumber;
    Zs = f * M_1;
    //求算三个外方位线元素初值
    double** x = new double* [pointNumber]();
    for(int i = 0; i < pointNumber; i++){
        x[i] = new double[1]();
    }
    double** y = new double* [pointNumber]();
    for(int i = 0; i < pointNumber; i++){
        y[i] = new double[1]();
    }
    double** A = new double* [2 * pointNumber]();
    for(int i = 0; i < 2 * pointNumber; i++){
        A[i] = new double[6]();
    }
    double** AT = new double* [6]();
    for(int i = 0; i < 6; i++){
        AT[i] = new double[2 * pointNumber]();
    }
    double** ATA = new double* [6]();
    for(int i = 0; i < 6; i++){
        ATA[i] = new double[6]();
    }
    double** ATAInverse = new double* [6]();
    for(int i = 0; i < 6; i++){
        ATAInverse[i] = new double[6]();
    }
    double** ATAInverseAT = new double* [6]();
    for(int i = 0; i < 6; i++){
        ATAInverseAT[i] = new double[2 * pointNumber]();
    }
    double** L = new double* [2 * pointNumber]();
    for(int i = 0; i < 2 * pointNumber; i++){
        L[i] = new double[1]();
    }
    double** ATAInverseATL = new double* [6]();
    for(int i = 0; i < 6; i++){
        ATAInverseATL[i] = new double[1]();
    }
    //开辟空间存储过程中的xy估计值，矩阵A，A转置AT，A乘A转置ATA，ATA的逆ATAInverse，L阵
    /*规模：
        x:n_1   y:n_1   A:2n_6  AT:6_2n     ATA:6_6     ATAInverse:6_6      ATAInverseAT:6_2n   L:2n_1      ATAInverseATL:6_1
    */
    int count = 0;//循环次数
    //循环次数限差
    while(1){
        printf("%.2lf %.2lf %.2lf %.6lf %.6lf %.6lf\n", Xs, Ys, Zs, Phai, Omega, Ka);

        Rotate = caculateRotateMat(Rotate, Phai, Omega, Ka);
for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
        cout << Rotate[i][j] << "  ";
    }
    cout << endl;
}
cout << "Rotate" << endl;        

        for(int i = 0; i < pointNumber; i++){
            x[i][0] = 0 - f * (Rotate[0][0] * (objectControl[i][0] - Xs) + Rotate[1][0] * (objectControl[i][1] - Ys) + Rotate[2][0] * (objectControl[i][2] - Zs))\
            /(Rotate[0][2] * (objectControl[i][0] - Xs) + Rotate[1][2] * (objectControl[i][1] - Ys) + Rotate[2][2] * (objectControl[i][2] - Zs));
            y[i][0] = 0 - f * (Rotate[0][1] * (objectControl[i][0] - Xs) + Rotate[1][1] * (objectControl[i][1] - Ys) + Rotate[2][1] * (objectControl[i][2] - Zs))\
            /(Rotate[0][2] * (objectControl[i][0] - Xs) + Rotate[1][2] * (objectControl[i][1] - Ys) + Rotate[2][2] * (objectControl[i][2] - Zs));
            //计算每个点的xy估计值
            L[2 * i][0] = imageControl[i][0] - x[i][0];
            L[2 * i + 1][0] = imageControl[i][1] - y[i][0];
            //计算每个点的lx，ly
            double H = Zs - objectControl[i][2];
            //计算航高H
            A[2 * i][0] = - f / H;
            A[2 * i][1] = 0;
            A[2 * i][2] = - x[i][0] / H;
            A[2 * i + 1][0] = 0;
            A[2 * i + 1][1] = - f / H;
            A[2 * i + 1][2] = - y[i][0] / H;
            A[2 * i][3] = - f * (1.0 + (x[i][0] * x[i][0]) / (f * f));
            A[2 * i][4] = - (x[i][0] * y[i][0]) / f;
            A[2 * i][5] = y[i][0];
            A[2 * i + 1][3] = - (x[i][0] * y[i][0]) / f;
            A[2 * i + 1][4] = - f * (1.0 + (y[i][0] * y[i][0]) / (f * f));
            A[2 * i + 1][5] = - x[i][0];
            //把该点这一轮的12个系数填进A矩阵
        }
        //求算每个点在最小二乘平差方程中的参数
for(int i = 0; i < 2 * pointNumber; i++){
    for(int j = 0; j < 6; j++){
        cout << A[i][j] << "  ";
    }
    cout << endl;
}
cout << "A" << endl;

        AT = MatrixTransport(A, 2 * pointNumber, 6);
for(int i = 0; i < 6; i++){
    for(int j = 0; j < 2 * pointNumber; j++){
        cout << AT[i][j] << "  ";
    }
    cout << endl;
}
cout << "AT" << endl;


        ATA = MatrixMutiply(AT, A, 6, 2 * pointNumber, 6);
for(int i = 0; i < 6; i++){
    for(int j = 0; j < 6; j++){
        cout << ATA[i][j] << "  ";
    }
    cout << endl;
}
cout << "ATA" << endl;

///////////////////////////////////////////////////////////////到这里没毛病

        ATAInverse = MatrixConverse(ATA, 6);                             //没写完这个玩意
        if(!ATAInverse){
            cout << "Inverse ERROR!";
            return 0;
        } 
for(int i = 0; i < 6; i++){
    for(int j = 0; j < 6; j++){
        cout << ATAInverse[i][j] << "  ";
    }
    cout << endl;
}
cout << "ATAINV" << endl;


        ATAInverseAT = MatrixMutiply(ATAInverse, AT, 6, 6, 2 * pointNumber);
for(int i = 0; i < 6; i++){
    for(int j = 0; j < 2 * pointNumber; j++){
        cout << ATAInverseAT[i][j] << "  ";
    }
    cout << endl;
}
cout << "ATAInverseAT" << endl;


        ATAInverseATL = MatrixMutiply(ATAInverseAT, L, 6, 2 * pointNumber, 1);
for(int i = 0; i < 6; i++){
    for(int j = 0; j < 1; j++){
        cout << ATAInverseAT[i][j] << "  ";
    }
    cout << endl;
}
cout << "ATAInverseATL" << endl;


        delta = ATAInverseATL;
        //解法方程，把解的6个参量存进delta矩阵

cout << delta[0][0] << delta[1][0] << delta[2][0] << delta[3][0] << delta[4][0] << delta[5][0] << endl;

        Xs += delta[0][0];
		Ys += delta[1][0];
		Zs += delta[2][0];
		Phai += delta[3][0];
		Omega += delta[4][0];
		Ka += delta[5][0];
        //改正原估计值
        count++;
        //计数器加一
        if((count > 99) || ((fabs(delta[3][0]) < 1e-6) && (fabs(delta[4][0]) < 1e-6) && (fabs(delta[5][0]) < 1e-6))){
            cout << "circulate done! count =" << count << endl;
            break;
            //如果循环超标或角元素已经在限差内，则跳出循环
        }
    } 

    if(count == 100){
        cout << "FlowOut!";
        return 0;
    }

    //如果到这里，说明已经求算完成了要求精度内的6个外方位元素， 所以接下来的任务是输出和释放内存
    
     printf("%.2lf %.2lf %.2lf %.6lf %.6lf %.6lf\n", Xs, Ys, Zs, Phai, Omega, Ka);

    delete[] imageControl, objectControl, delta, Rotate, x, y, A, AT, ATA, ATAInverse, ATAInverseAT, ATAInverseATL, L = NULL;

system("pause");

    return 0;
}

/////////////////函数定义////////////////

double** caculateRotateMat(double** Mat, double a, double b, double c){
    double** result = new double* [3]();
    for(int i = 0; i < 3; i++){
        result[i] = new double[3]();
    }
    result[0][0] = cos(a) * cos(c) - sin(a) * sin(b) * sin(c);
    result[0][1] = 0.0 - cos(a) * sin(c) - sin(a) * sin(b) * cos(c);
    result[0][2] = 0.0 - sin(a) * cos(b);
    result[1][0] = cos(b) * sin(c);
    result[1][1] = cos(b) * cos(c);
    result[1][2] = 0.0 - sin(b);
    result[2][0] = sin(a) * cos(c) + cos(a) * sin(b) * sin(c);
    result[2][1] = 0.0 - sin(a) * sin(c) + cos(a) * sin(b) * cos(c);
    result[2][2] = cos(a) * cos(b);
    //计算旋转矩阵
    return result;
}
double** MatrixMutiply(double** Mat1, double** Mat2, int r, int n, int c){
    double** result = new double* [r]();
    for(int i = 0; i < r; i++){
        result[i] = new double[c]();
    }
    //开辟内存存放结果矩阵，并规定好行列数
    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++){
            for(int k = 0; k < n; k++){
                result[i][j] += Mat1[i][k] * Mat2[k][j];
            }
        }
    }
    //矩阵乘法
    return result;
}
double** MatrixTransport(double** Mat, int r, int c){
    double** result = new double* [c]();
    for(int i = 0; i < r; i++){
        result[i] = new double[r]();
    }
    //开辟内存空间存放结果矩阵
    for(int i = 0; i < c; i++){
        for(int j = 0; j < r; j++){
            result[i][j] = Mat[j][i];
        }
    }
    //转置变换
    return result;
}


double** MatrixRemain(double** Mat, int a, int x, int y){
    double** result = new double* [a-1]();
    for(int i = 0; i < a-1; i++){
        result[i] = new double[a-1]();
    }
    //开辟空间为a-1*a-1的矩阵
    for(int i = 0; i < x-1; i++){
        for(int j = 0; j < y-1; j++){
            result[i][j] = Mat[i][j];
            //左上角的元素行列号不变
        }
    }
    for(int i = x-1; i < a-1; i++){
        for(int j = 0; j < y-1; j++){
            result[i][j] = Mat[i+1][j];
            //左下角跳过第x行元素
        }
    }
    for(int i = 0; i < x-1; i++){
        for(int j = y-1; j < a-1; j++){
            result[i][j] = Mat[i][j+1];
            //右上角跳过第y列元素
        }
    }
    for(int i = x-1; i < a-1; i++){
        for(int j = y-1; j < a-1; j++){
            result[i][j] = Mat[i+1][j+1];
            //右下角跳过第x行，第y列元素
        }
    }
    //完成余子矩阵
    return result;
}
double MatrixDet(double** Mat, int a){
    //递归法求矩阵行列式
    if(a == 1) return Mat[0][0];
    else{
        double result = 0.0;
        for(int i = 0; i < a; i++){
            result += Mat[0][i] * MatrixDet(MatrixRemain(Mat, a, 1, i+1), a-1) * pow(-1.0, (i + 2)*1.0);
        }
        return result;
    }
}
double** MatrixConverse(double** Mat, int a){
    double** result = new double* [a]();
    for(int i = 0; i < a; i++){
        result[i] = new double[a]();
    }
    double det = MatrixDet(Mat, a);
    //求算矩阵的行列式并判断
    if(det == 0){
        return NULL;
    }
    for(int i = 0; i < a; i++){
        for(int j = 0; j < a; j++){
            result[i][j] = MatrixDet(MatrixRemain(Mat, a, i+1, j+1), a-1) * pow(-1.0, 1.0 * (i + j + 2)) / det;
        }
    }
    return result;
}