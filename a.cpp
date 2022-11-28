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
    Matrix minor(int m,int n){
        if(m>=dimX||n>=dimY||m<0||n<0) return *this;
        int i=0,j=0;
        vector<vector<double>> a(m-1,vector<double>{});
        while(j<this->dimY){
            if(j==n) continue;
            while(i<this->dimX){
                if(m>i){
                    a[i].push_back(a[i][j]);
                }
                if(m<i){
                    a[i].push_back(a[i-1][j]);
                }
                i++;
            }
            j++;
        }
        return a;
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
    bool isSquare()const {
        return dimX==dimY;
    }
    static double det(const Matrix& x){
        if(x.dimX==0||x.dimY==0) return 0;
        if(x.dimX==1&&x.dimY==1){
            return x[0][0];
        }
        if(!x.isSquare()) throw invalid_argument("Matrix: Trying to get determinant of non square matrix");
        int sum=0;
        for(int i=0;i<x.dimY;i++){
            sum+=x[0][i]*(i%2==1?-1:1)*det(x.minor(0,i));
        }
        return sum;
    }
    double operator=(const Matrix x){
        this->a=x.a;
        this->dimX=x.dimX;
        this->dimY=x.dimX;
        return 0;
    }
    static void printMatrix(const Matrix& x){
        for(auto i:x){
            for(auto j:i){
                cout<<j<<" ";
            }
            cout<<"\n";
        }
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
    friend Matrix operator-(Matrix l,const Matrix& r){
        return l+(-1)*r;
        // if(l.dimX==r.dimX&&l.dimY==r.dimY){
        //     for(int x=0;x<l.dimX;x++){
        //         for(int y=0;y<l.dimY;y++){
        //             l[x][y]-= r[x][y];
        //         }
        //     }
        //     return l;
        // } else {
        //     throw invalid_argument("Trying to add together matrices of different dimensions.");
        // }
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
typedef Matrix(*Function)(Matrix);
class MatrixFunction{
    public:
    Function func;
    Function gradient;
    Function hessian;
    MatrixFunction(Function func,Function gradient,Function hessian):func(func),gradient(gradient),hessian(hessian){}
};
class Descent{
    MatrixFunction func;
    Descent(MatrixFunction func):func(func){}
    bool isDescentDirection(Matrix x,Matrix d){
        if(!Matrix::isSameDimension(x,d)) throw invalid_argument("Descent: Descent direction is not of the same dimension as x");
        return Matrix::det(Matrix::transpose(x)*d)<0;
    }
};
class GradientDescent{
    Function func;
    Function gradient;
    Function hessian;
    double dim=0;
    public:
    GradientDescent(Function func,Function gradient,Function hessian,double dim){
        this->func=func;
        this->gradient=gradient;
        this->hessian=hessian;
        this->dim=dim;
    }

    double exactLineSearch(Matrix x,double m, double alpha,double precision){
        Matrix df=gradient(x);
        double hDoubled=Matrix::normDoubled(hessian(x));
        double expected=m*hDoubled;
        throw logic_error("Not implemented yet");
        return 0;
    }
    Matrix backtrack(Matrix x, double m,double alpha, double precision){
        Matrix df=gradient(x);
        double hDoubled=Matrix::normDoubled(hessian(x));
        double expected=m*hDoubled;
        double t=1;
        Loop:
        Matrix y=x-m*t*df;
        double diff=Matrix::normDoubled(func(y)-func(x));
        if(diff>expected){
            t=alpha*t;
            goto Loop;
        }
        return diff>precision?backtrack(y,m,alpha,precision):y;
    }
};

int main(){
    vector<vector<double>> a{{1,2},{3,4}};
    Matrix x(a);
    Matrix::printMatrix(Matrix::det(x));
    // Matrix::printMatrix(2*x);
    // Matrix::printMatrix(x*2);
    // return 0;
    // GradientDescent g(
    //     [](Matrix x)->Matrix{return pow((x[0][0]-4),2)+pow((x[0][1]-6),2);},
    //     [](Matrix x)->Matrix{return vector<vector<double>>{{2*(x[0][0]-4),2*(x[0][1]-6)}};},
    //     [](Matrix x)->Matrix{return 4;},
    //     1
    // );
    // try{
    // Matrix::printMatrix(g.backtrack(Matrix(vector<vector<double>>(1,vector<double>{-1.333,-11.2312})),0.5,0.75,0.1));
    // } catch(exception e){
    //     cout<<e.what();
    // }
}