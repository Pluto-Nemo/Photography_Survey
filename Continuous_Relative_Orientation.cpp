#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <math.h>
using namespace std;
double** caculateRotateMat(double a, double b, double c);
double** MatrixMutiply(double** Mat1, double** Mat2, int r, int n, int c);
double** MatrixTransport(double** Mat, int r, int c);
double** MatrixConverse(double** Mat, int a);
double** MatrixRemain(double** Mat, int a, int x, int y);
double MatrixDet(double** Mat, int a);
double** ALTOX(double** A, int a, int n, double** L, int b);

int main(){
    double f = 0.024;
    double Bx = 0.005185;
    //基础值
    int pointNumber;
    double tmpx1, tmpy1, tmpx2, tmpy2;
    ifstream infile;
    infile.open("data2.txt", ios::in);
    if(!infile.is_open()){
        cout << "Open file failed!" << endl;
        system("pause"); 
    }
    infile >> pointNumber;
    double** pointPairs = new double* [100]();
    for(int i = 0; i < pointNumber; i++){
        infile >> tmpx1 >> tmpy1 >> tmpx2 >> tmpy2;
        pointPairs[i] = new double[4]();
        pointPairs[i][0] = tmpx1 / 1000.0;
        pointPairs[i][1] = tmpy1 / 1000.0;
        pointPairs[i][2] = tmpx2 / 1000.0;
        pointPairs[i][3] = tmpy2 / 1000.0;
        printf("%.6lf %.6lf %.6lf %.6lf\n", pointPairs[i][0], pointPairs[i][1], pointPairs[i][2], pointPairs[i][3]);
    }
    //从文件赋值同名点对

    double  Miu = 0.0, Niu = 0.0, Phai = 0.0, Omega = 0.0, Ka = 0.0;
    //开辟相对定向元素参数空间
    double N1, N2;
    double* auxL = new double[3]();
    double* auxR = new double[3]();
    double* B = new double[3]{Bx, 0, 0};
    double** Rotate = new double* [3]();
    for(int i = 0; i < 3; i++){
        Rotate[i] = new double[3]();
    }

    double** delta = new double* [5]();
    for(int i = 0; i < 5; i++){
        delta[i] = new double[1]();
    }
    double** A = new double* [pointNumber]();
    for(int i = 0; i < pointNumber; i++){
        A[i] = new double[5]();
    }
    double** L = new double* [pointNumber]();
    for(int i = 0; i < pointNumber; i++){
        L[i] = new double[1]();
    }
    //开辟
int flag = 0;
    for(int count = 0; count < 99; count++){
        printf("%.6lf %.6lf %.6lf %.6lf %.6lf %d\n", Phai, Omega, Ka, Miu, Niu, count);
        system("pause");
        Rotate = caculateRotateMat(Phai, Omega, Ka);
        for(int i = 0; i < pointNumber; i++){
            auxL[0] = pointPairs[i][0]; auxL[1] = pointPairs[i][1]; auxL[2] = -f;
            auxR[0] = Rotate[0][0] * pointPairs[i][2] + Rotate[0][1] * pointPairs[i][3] - Rotate[0][2] * f;
            auxR[1] = Rotate[1][0] * pointPairs[i][2] + Rotate[1][1] * pointPairs[i][3] - Rotate[1][2] * f;
            auxR[2] = Rotate[2][0] * pointPairs[i][2] + Rotate[2][1] * pointPairs[i][3] - Rotate[2][2] * f;
// printf("%.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n", auxL[0], auxL[1], auxL[2], auxR[0], auxR[1], auxR[2]);
// system("pause");
            B[1] = Miu * B[0]; B[2] = Niu * B[0];
            N1 = (B[0] * auxR[2] - B[2] * auxR[0]) / (auxL[0] * auxR[2] - auxL[2] * auxR[0]);
            N2 = (B[0] * auxL[2] - B[2] * auxL[0]) / (auxL[0] * auxR[2] - auxL[2] * auxR[0]);
// printf("%.6lf %.6lf %.6lf %.6lf\n", B[1], B[2], N1, N2);
// system("pause");
            A[i][0] = -auxR[0] * auxR[1] / auxR[2] * N2;
            A[i][1] = -(auxR[2] + auxR[1] * auxR[1] / auxR[2]) * N2;
            A[i][2] = auxR[0] * N2;
            A[i][3] = B[0];
            A[i][4] = -auxR[1] / auxR[2] * B[0];
            L[i][0] = N1 * auxL[1] - N2 * auxR[1] - B[1];
// printf("%.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n", A[i][0], A[i][1], A[i][2], A[i][3], A[i][4], L[i][0]);
// system("pause");
        //求算A阵
        }
        delta = ALTOX(A, 5, pointNumber, L, 1);
        //法化
        Phai+=delta[0][0]; Omega+=delta[1][0]; Ka+=delta[2][0]; Miu+=delta[3][0]; Niu+=delta[4][0];
        if((fabs(delta[3][0]) < 1e-6) && 
        (fabs(delta[4][0]) < 1e-6) && 
        (fabs(delta[0][0]) < 1e-6) && 
        (fabs(delta[1][0]) < 1e-6) && 
        (fabs(delta[2][0]) < 1e-6)){
            cout << "circulate done! count =" << count << endl;
            flag = 1;
            break;
        }
    }
    //这里已经跳出循环
    if(flag == 0){
        cout << "FlowOut!";
        system("pause");
        return 0;
    }
    //考虑循环溢出的情况
    printf("Phai=%.6lf Omega=%.6lf Ka=%.6lf Miu=%.6lf Niu=%.6lf\n", Phai, Omega, Ka, Miu, Niu);
    delete[] pointPairs, Rotate, delta, A, L = NULL;
    system("pause");
    return 0;
}








//////////////////////函数定义///////////////////////

double** caculateRotateMat(double a, double b, double c){
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

double** ALTOX(double** A, int a, int n, double** L, int b = 1){
    double** AT = new double* [a]();
    for(int i = 0; i < a; i++){
        AT[i] = new double[n]();
    }
    double** ATA = new double* [a]();
    for(int i = 0; i < a; i++){
        ATA[i] = new double[a]();
    }
    double** ATAInverse = new double* [a]();
    for(int i = 0; i < a; i++){
        ATAInverse[i] = new double[a]();
    }
    double** ATAInverseAT = new double* [a]();
    for(int i = 0; i < a; i++){
        ATAInverseAT[i] = new double[n]();
    }
    double** ATAInverseATL = new double* [a]();
    for(int i = 0; i < a; i++){
        ATAInverseATL[i] = new double[1]();
    }
    AT = MatrixTransport(A, n, a);
    ATA = MatrixMutiply(AT, A, a, n, a);
    ATAInverse = MatrixConverse(ATA, a);
    ATAInverseAT = MatrixMutiply(ATAInverse, AT, a, a, n);
    ATAInverseATL = MatrixMutiply(ATAInverseAT, L, a, n, b);
    return ATAInverseATL;
}