//Hi, as I was coding this, I realized I cannot handle all the exceptions nor have the time to.
//So be warned.
#include <bits/stdc++.h>
using namespace std;
class Matrix{
    vector<vector<double>> a;
    static Matrix dotProduct(const Matrix& l,const Matrix& r){
        if(l.dimX==r.dimX){
            return transpose(l)*r;
        } else if(l.dimY==r.dimY){
            return transpose(r)*l;
        } else {
            throw invalid_argument("Trying to gain dot product of matrices of unfit dimensions.");
        }
    }
    public:
    double dimX,dimY;
    Matrix(const vector<vector<double>> value): a(value), dimX(a.size()),dimY(a[0].size()){}
    Matrix(const double value):a(vector<vector<double>>(1,vector<double>(1,value))), dimX(a.size()),dimY(a[0].size()){}
    static Matrix transpose(const Matrix x){
        Matrix a(vector<vector<double>>(x.dimY,vector<double>(x.dimX)));
        for(int i=0; i<x.dimX;i++){
            for(int j=0;j<x.dimY;j++){
                a[j][i]=x[i][j];
            }
        }
        return a;
    }
    void transposeThis(){
        *this=transpose(*this);
        // vector<vector<double>> a(this->dimY,vector<double>(this->dimX));
        // for(int i=0;i++; i<this->dimX){
        //     for(int j=0;j++;j<this->dimY){
        //         a[j][i]=this->a[i][j];
        //     }
        // }
        // this->dimX=a.size();
        // this->dimY=a[0].size();
        // this->a=a;
    }
    bool containsNaN(){
        for(int i=0;i<this->dimX;i++){
            for(int j=0;j<this->dimY;j++){
                if(isinf(this->a[i][j])) return true;
            }
        }
        return false;
    }
    Matrix minor(int m,int n)const{
        if(m<0||n<0||m>=dimX||n>=dimY) return *this;
        int j=0;
        vector<vector<double>> a(dimX-1,vector<double>{});
        while(j<this->dimY){
            if(j==n){
                j++;
                continue;
            } 
            int i=0;
            while(i<this->dimX){
                if(m>i){
                    a[i].push_back(this->a[i][j]);
                }
                if(m<i){
                    a[i-1].push_back(this->a[i][j]);
                }
                i++;
            }
            j++;
        }
        return a;
    }
    static bool isSameDimension(const Matrix& x,const Matrix& y){
        return x.dimX==y.dimX&&x.dimY==y.dimY;
    }
    static double norm(const Matrix& x){
        return sqrt(normDoubled(x));
    }
    static double normDoubled(const Matrix& x){
        double sum=0;
        for(auto i:dotProduct(x,x)){
            for(auto j:i){
                sum+=j;
            }
        }
        return sum;
    }
    static Matrix inverse(Matrix x){
        if(!x.isSquare()) throw invalid_argument("Matrix: Trying to get inverse of nonsquare matrix");
        if(x.isSingular()) throw invalid_argument("Matrix: Trying to get inverse of singular matrix");
        vector<vector<double>> a(x.dimX,vector<double>{});
        double dt=Matrix::det(x);
        for(int i=0;i<x.dimX;i++){
            for(int j=0;j<x.dimX;j++){
                a[i].push_back((i+j%2==1?-1:1)*Matrix::det(x.minor(i,j))/dt);
            }
        }
        return transpose((Matrix)a);
    }
    static void printMatrix(const Matrix& x){
        for(auto i:x){
            for(auto j:i){
                cout<<j<<" ";
            }
            cout<<"\n";
        }
    }
    bool isSquare()const {
        return dimX==dimY;
    }
    static double det(const Matrix& x){
        if(x.dimX==0||x.dimY==0) return 0;
        if(x.dimX==1&&x.dimY==1){
            return x[0][0];
        }
        if(!x.isSquare()) throw invalid_argument("Matrix: Trying to get determinant of non square matrix");
        double sum=0;
        for(int i=0;i<x.dimY;i++){
            sum+=x[0][i]*(i%2==1?-1:1)*det(x.minor(0,i));
        }
        return sum;
    }
    bool isSingular() const{
        for(auto i:this->a){
            for(auto j:i){
                if(j!=0){
                    return false;
                }
            }
        }
        return true;
    }
    double operator=(const Matrix x){
        this->a=x.a;
        this->dimX=x.dimX;
        this->dimY=x.dimX;
        return 0;
    }
    vector<double>& operator[](const double n){
        if(n<0){
            throw out_of_range("Index is larger than matrix dimension.");
        }
        if(n>=a.size()){
            throw invalid_argument("Trying to access a negative index.");
        }
        return a[n];

    }
    const vector<double>& operator[](const double n)const{
        if(n<0){
            throw out_of_range("Trying to access a negative index.");
        }
        if(n>=a.size()){
            throw invalid_argument("Index is larger than matrix dimension.");
        }
        return a[n];
    }
    friend Matrix operator+(Matrix l,const Matrix& r){
        if(isSameDimension(l,r)){
            for(int x=0;x<l.dimX;x++){
                for(int y=0;y<l.dimY;y++){
                    l[x][y]+= r[x][y];
                }
            }
            return l;
        } else {
            throw invalid_argument("Trying to add together matrices of different dimensions.");
        }
    }
    friend Matrix operator-(Matrix l,const Matrix& r) {
        return l+(-r);
    }
    Matrix operator-() const{
        return (-1)*(*this);
    }
    friend Matrix operator*(const Matrix& l,const Matrix& r){
        //TODO: Replace this with more effective solutions
        if(l.dimY==r.dimX){
            vector<vector<double>> a(l.dimX,vector<double>(r.dimY,0));
            for(int x=0;x<l.dimX;x++){
                for(int y=0;y<l.dimY;y++){
                    for(int z=0;z<r.dimY;z++){
                        a[x][z]+=l[x][y]*r[y][z];
                    }
                }
            }
            return Matrix(a);
        } else {
            throw invalid_argument("Trying to add together matrices of different dimensions.");
        }
    }
    friend Matrix operator*(const double n,Matrix x){
        for(auto& i:x){
            for(auto& j:i){
                j*=n;
            }
        }
        return x;
    }
    friend Matrix operator*(const Matrix x,const double n){
        return n*x;
    }
    friend Matrix operator/(Matrix x, const double n){
        for(auto& i:x){
            for(auto& j:i){
                j/=n;
            }
        }
        return x;
    }
    vector<std::vector<double>>::iterator begin(){
        return this->a.begin();
    }
    vector<std::vector<double>>::const_iterator begin() const{
        return this->a.begin();
    }
    vector<std::vector<double>>::iterator end(){
        return this->a.end();
    }
    vector<std::vector<double>>::const_iterator end() const{
        return this->a.end();
    }
};
template <class T>
class Func{
    public:
    typedef T(*Function)(T);
    Function func;
    Function gradient;
    Function hessian;
    Func(Function func,Function gradient,Function hessian):func(func),gradient(gradient),hessian(hessian){}
};
// class Descent{
//     Func<Matrix> func;
//     Descent(Func<Matrix> func):func(func){}
//     bool isDescentDirection(Matrix x,Matrix d){
//         if(!Matrix::isSameDimension(x,d)) throw invalid_argument("Descent: Descent direction is not of the same dimension as x");
//         return Matrix::det(Matrix::transpose(x)*d)<0;
//     }
// };
class Descent{
    //WARNING: This might end up being not the correct minimizer t, but that's the limitation of Exact Line Search.
    double exactLineSearch_newtonMethod(Matrix x,double t, Matrix d, double precision){
        Matrix df=-Matrix::transpose(innerFunction.gradient(x-t*d))*d;
        if(Matrix::norm(df)<precision) {
            if(t<0) return 0;
            return t;
        }
        else {
            return exactLineSearch_newtonMethod(x,t-Matrix::det(innerFunction.func(x-t*d))/Matrix::det(df),d,precision);
        }
    }
    public:
    Func<Matrix> innerFunction;
    bool isDescentDirection(Matrix x,Matrix d){
        if(!Matrix::isSameDimension(x,d)) throw invalid_argument("Descent: Descent direction is not of the same dimension as x");
        return Matrix::det(Matrix::transpose(x)*d)<0;
    }
    Descent(Func<Matrix> func): innerFunction(func){}
    Matrix GradientDescent_exactLineSearch(Matrix x,double precision){
        // cout<<"OKAY\n";
        // Matrix::printMatrix(x);
        Matrix df=innerFunction.gradient(x);
        if(Matrix::norm(df)<precision) return x;
        return GradientDescent_exactLineSearch(x-exactLineSearch_newtonMethod(x,1,df,precision/1000)*df,precision);
    }
    Matrix GradientDescent_backtrack(Matrix x, double m,double alpha, double precision){
        Matrix df=innerFunction.gradient(x);
        double hDoubled=Matrix::normDoubled(innerFunction.hessian(x));
        double expected=m*hDoubled;
        double t=1;
        Loop:
        Matrix y=x-m*t*df;
        double diff=Matrix::normDoubled(innerFunction.func(y)-innerFunction.func(x));
        if(diff>=expected){
            t=alpha*t;
            goto Loop;
        }
        return diff>precision?GradientDescent_backtrack(y,m,alpha,precision):y;
    }
    Matrix newtonMethod_pure(Matrix x,double precision){
        Matrix df=innerFunction.gradient(x);
        if(Matrix::norm(df)<precision) return x;
        Matrix Hf=innerFunction.hessian(x);
        // //WARNING: Hf might be singular so maybe implement this?// Since this is pure newton, I won't add this.
        // if(Hf.isSingular()){
        //     Matrix x1=x-exactLineSearch_newtonMethod(x,1,df,precision/1000)*df;\
        //     if(Matrix::norm(innerFunction.gradient(x))<precision) return x;
        //     else return newtonMethodPure(x1,precision);
        // }
        Matrix x1=x-Matrix::inverse(Hf)*df;
        return newtonMethod_pure(x1,precision);
    }
    Matrix newtonMethod_lineSearch(Matrix x,double precision){
        cout<<"----\n";
        Matrix df=innerFunction.gradient(x);
        if(Matrix::norm(df)<precision) return x;
        Matrix Hf=innerFunction.hessian(x);
        if(Hf.isSingular()){
            Matrix x1=x-exactLineSearch_newtonMethod(x,1,df,precision/10)*df;\
            if(Matrix::norm(innerFunction.gradient(x))<precision) return x;
            else return newtonMethod_lineSearch(x1,precision);
        }
        Matrix direction=Matrix::inverse(Hf)*df;
        Matrix::printMatrix(Matrix::inverse(Hf));
        double t= exactLineSearch_newtonMethod(x,1,direction,precision/10);
        cout<<"OKE\n";
        Matrix x1=x-t*direction;
        return newtonMethod_lineSearch(x1,precision);
    }
    Matrix newtonMethod_backtrack(Matrix x,double alpha,double precision){
        Matrix df=innerFunction.gradient(x);
        if(Matrix::norm(df)<precision) return x;
        Matrix Hf=innerFunction.hessian(x);
        if(Hf.isSingular()){
            Matrix x1=x-exactLineSearch_newtonMethod(x,1,df,precision/10)*df;\
            if(Matrix::norm(innerFunction.gradient(x))<precision) return x;
            else return newtonMethod_backtrack(x1,alpha,precision);
        }
        Matrix direction=Matrix::inverse(Hf)*df;
        double t= 1;
        Loop:
        Matrix x1=x-t*direction;
        if(Matrix::det(innerFunction.func(x1))>=Matrix::det(innerFunction.func(x))){
            t=alpha*t;
            goto Loop;
        }
        return newtonMethod_backtrack(x1,alpha,precision);
    }
};
#define p 10
int main(){
    Descent g(Func<Matrix>(
        [](Matrix x)->Matrix{return -9*x[0][0]-10*x[1][0]+p*(-log(100-x[0][0]-x[1][0])-log(x[0][0])-log(x[1][0])-log(50-x[0][0]+x[1][0]));},
        [](Matrix x)->Matrix{return vector<vector<double>>{
            {-9+p*(1/(100-x[0][0]-x[1][0])-1/x[0][0]+1/(5-x[0][0]+x[1][0]))},
            {-10+p*(1/(100-x[0][0]-x[1][0]-1/x[1][0]-1/(5-x[0][0]+x[1][0])))}
        };},
        [](Matrix x)->Matrix{return p*Matrix(vector<vector<double>>{
            {1/pow(100-x[0][0]-x[1][0],2)+1/pow(x[1][0],2)+1/pow(50-x[0][0]+x[1][0],2),
            1/pow(100-x[0][0]-x[1][0],2)-1/pow(50-x[0][0]+x[1][0],2)},
            {1/pow(100-x[0][0]-x[1][0],2)-1/pow(50-x[0][0]+x[1][0],2),
            1/pow(100-x[0][0]-x[1][0],2)+1/pow(x[1][0],2)+1/pow(50-x[0][0]+x[1][0],2)}
        });
        }
    )
    );
    try{
    Matrix::printMatrix(g.newtonMethod_pure(Matrix(vector<vector<double>>{{8},{90}}),0.5));
    // Matrix::printMatrix(g.newtonMethod_backtrack(Matrix(vector<vector<double>>{{-1.333},{-11.2312}}),0.25,0.5));
    } catch(const exception& e){
        cout<<e.what();
    }
    // Descent g(Func<Matrix>(
    //     [](Matrix x)->Matrix{return pow((x[0][0]-4),2)+pow((x[1][0]-6),2);},
    //     [](Matrix x)->Matrix{return vector<vector<double>>{{2*(x[0][0]-4)},{2*(x[1][0]-6)}};},
    //     [](Matrix x)->Matrix{return vector<vector<double>>{{2,0},{0,2}};}
    // )
    // );
    // try{
    // // Matrix::printMatrix(g.GradientDescent_backtrack(Matrix(vector<vector<double>>{{-1.333},{-11.2312}}),0.5,0.75,0.1));
    // // Matrix::printMatrix(g.GradientDescent_exactLineSearch(Matrix(vector<vector<double>>{{-1.333},{-11.2312}}),0.1));
    // // Matrix::printMatrix(g.newtonMethod_pure(Matrix(vector<vector<double>>{{-1.333},{-11.2312}}),0.5));
    // Matrix::printMatrix(g.newtonMethod_lineSearch(Matrix(vector<vector<double>>{{-1.333},{-11.2312}}),0.5));
    // Matrix::printMatrix(g.newtonMethod_backtrack(Matrix(vector<vector<double>>{{-1.333},{-11.2312}}),0.25,0.5));
    // } catch(const exception& e){
    //     cout<<e.what();
    // }
}