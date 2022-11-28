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
            throw out_of_range("Index is larger than matrix dimension.");
        }
        if(n>=a.size()){
            throw invalid_argument("Trying to access a negative index.");
        }
        return a[n];
    }
    friend Matrix operator+(Matrix l,const Matrix& r){
        if(l.dimX==r.dimX&&l.dimY==r.dimY){
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
// class MatrixFunction{
//     public:
//     Function func;
//     Function gradient;
//     Function hessian;
//     MatrixFunction(Function func,Function gradient,Function hessian):func(func),gradient(gradient),hessian(hessian){}
// };
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
    // vector<vector<double>> a{{1,2},{3,4}};
    // Matrix x(a);
    // Matrix::printMatrix((Matrix)1);
    // Matrix::printMatrix(2*x);
    // Matrix::printMatrix(x*2);
    // return 0;
    GradientDescent g(
        [](Matrix x)->Matrix{return pow((x[0][0]-4),2)+pow((x[0][1]-6),2);},
        [](Matrix x)->Matrix{return vector<vector<double>>{{2*(x[0][0]-4),2*(x[0][1]-6)}};},
        [](Matrix x)->Matrix{return 4;},
        1
    );
    try{
    Matrix::printMatrix(g.backtrack(Matrix(vector<vector<double>>(1,vector<double>{-1,-12})),0.5,0.75,0.1));
    } catch(exception e){
        cout<<e.what();
    }
}