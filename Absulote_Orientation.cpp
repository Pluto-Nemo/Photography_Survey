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
double** ALTOX(double** A, int a, int n, double** L, int b = 1);

int main(){
    int pointNumber;
    double X, Y, Z, Xtp, Ytp, Ztp;
    ifstream infile;
    infile.open("data3.txt", ios::in);
    if(!infile.is_open()){
        cout << "Open file failed!" << endl;
        system("pause"); 
    }
    infile >> pointNumber;
    double** Points = new double* [100]();
    for(int i = 0; i < pointNumber; i++){
        infile >> X >> Y >> Z >> Xtp >> Ytp >> Ztp;
        Points[i] = new double[6]();
        Points[i][0] = X;
        Points[i][1] = Y;
        Points[i][2] = Z;
        Points[i][3] = Xtp;
        Points[i][4] = Ytp;
        Points[i][5] = Ztp;
        printf("%.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n", Points[i][0], Points[i][1], Points[i][2], Points[i][3], Points[i][4], Points[i][5]);
    }

    double dX = 0, dY = 0, dZ = 0, Phai = 0, Omega = 0, Ka = 0, Lamda = 1;

    double** Rotate = new double* [3]();
    for(int i = 0; i < 3; i++){
        Rotate[i] = new double[3]();
    }
    double** delta = new double* [7]();
    for(int i = 0; i < 7; i++){
        delta[i] = new double[1]();
    }
    double** A = new double* [3 * pointNumber]();
    for(int i = 0; i < 3 * pointNumber; i++){
        A[i] = new double[7]();
    }
    double** L = new double* [3 * pointNumber]();
    for(int i = 0; i < 3 * pointNumber; i++){
        L[i] = new double[1]();
    }
int flag = 0;
    for(int count = 0; count < 99; count++){
        printf("%.3lf %.3lf %.3lf %.6lf %.6lf %.6lf %.6lf %d\n", dX, dY, dZ, Lamda, Phai, Omega, Ka, count);
        system("pause");
        Rotate = caculateRotateMat(Phai, Omega, Ka);
        for(int i = 0; i < pointNumber; i++){
            A[i * 3 + 0][0] = 1; A[i * 3 + 0][1] = 0; A[i * 3 + 0][2] = 0;
            A[i * 3 + 0][3] = Points[i][0];
            A[i * 3 + 0][4] = -Points[i][2];
            A[i * 3 + 0][5] = 0;
            A[i * 3 + 0][6] = -Points[i][1];
            A[i * 3 + 1][0] = 0; A[i * 3 + 1][1] = 1; A[i * 3 + 1][2] = 0;
            A[i * 3 + 1][3] = Points[i][1];
            A[i * 3 + 1][4] = 0;
            A[i * 3 + 1][5] = -Points[i][2];
            A[i * 3 + 1][6] = Points[i][0];
            A[i * 3 + 2][0] = 0; A[i * 3 + 2][1] = 0; A[i * 3 + 2][2] = 1;
            A[i * 3 + 2][3] = Points[i][2];
            A[i * 3 + 2][4] = Points[i][0];
            A[i * 3 + 2][5] = Points[i][1];
            A[i * 3 + 2][6] = 0;
            //这里赋值的时候因为角度比较小就直接赋值原来的值了
            L[i * 3 + 0][0] = Points[i][3] - Lamda * (Rotate[0][0] * Points[i][0] + Rotate[0][1] * Points[i][1] + Rotate[0][2] * Points[i][2]) - dX;
            L[i * 3 + 1][0] = Points[i][4] - Lamda * (Rotate[1][0] * Points[i][0] + Rotate[1][1] * Points[i][1] + Rotate[1][2] * Points[i][2]) - dY;
            L[i * 3 + 2][0] = Points[i][5] - Lamda * (Rotate[2][0] * Points[i][0] + Rotate[2][1] * Points[i][1] + Rotate[2][2] * Points[i][2]) - dZ;
        }
        delta = ALTOX(A, 7, 3 * pointNumber, L, 1);
        Phai+=delta[4][0]; Omega+=delta[5][0]; Ka+=delta[6][0]; dX+=delta[0][0]; dY+=delta[1][0]; dZ+=delta[2][0]; Lamda = Lamda * (delta[3][0] + 1);
        if((fabs(delta[3][0]) < 1e-6) && 
        (fabs(delta[4][0]) < 1e-6) && 
        (fabs(delta[0][0]) < 1e-6) && 
        (fabs(delta[1][0]) < 1e-6) && 
        (fabs(delta[2][0]) < 1e-6) &&
        (fabs(delta[5][0]) < 1e-6) && 
        (fabs(delta[6][0]) < 1e-6)){
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
    printf("dX=%.3lf dY=%.3lf dZ%.3lf Lamda=%.6lf Phai=%.6lf Omega=%.6lf Ka=%.6lf\n", dX, dY, dZ, Lamda, Phai, Omega, Ka);
    delete[] Points, Rotate, delta, A, L = NULL;
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

double** ALTOX(double** A, int a, int n, double** L, int b){
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