#include<Fastor/Fastor.h>
#include<octave/oct.h>
#include<octave/octave.h>
#include<octave/parse.h>
#include<cmath>
using namespace Fastor;

template<int m,int n> Tensor<double,m,n>& OctMat2Ten(const Matrix& A){
    Tensor<double,m,n>* Ten = new Tensor<double,m,n>(0.0);
    for(int i=0;i<m;++i){
        for(int j=0;j<n;++j){
            (*Ten)(i,j) = A(i,j);
        }
    }
    return *Ten;
}

template<int m,int n> Matrix Ten2OctMat(const Tensor<double,m,n>& A){
    Matrix Mat(m,n);
    for(int i=0;i<m;++i){
        for(int j=0;j<n;++j){
            Mat(i,j) = A(i,j);
        }
    }
    return Mat;
}

template<int m,int n> Tensor<double,m,n>* gradient(const Tensor<double,m,n>& A, double dy,double dx){
    Tensor<double,m,n>* temp = nullptr;
    temp = new Tensor<double,m,n> [2];
    for(int i=0;i<m;++i){
        for(int j=0;j<n;++j){
            if(i==0){
                temp[0][i][j] = (A(i+1,j)-A(i,j))/dy;
            }else if(i==m-1){
                temp[0][i][j] = (A(i,j)-A(i-1,j))/dy;
            }else{
                temp[0][i][j] = (A(i+1,j)-A(i-1,j))/2/dy;
            }
            if(j==0){
                temp[1][i][j] = (A(i,j+1)-A(i,j))/dx;
            }else if(j==n-1){
                temp[1][i][j] = (A(i,j)-A(i,j-1))/dx;
            }else{
                temp[1][i][j] = (A(i,j+1)-A(i,j-1))/2/dx;
            }
        }
    }
    return temp;
}

template<int m,int n> Tensor<double,m,n>& ExtendBoundary(Tensor<double,m-2,n-2>& A){
    Tensor<double,m,n>* temp = new Tensor<double,m,n>;
    (*temp)(0,0) = A(0,0);
    (*temp)(0,n-1) = A(0,n-3);
    (*temp)(m-1,0) = A(m-3,0);
    (*temp)(m-1,n-1) = A(m-3,n-3);
    for(int i=1;i<m-1;++i){
        (*temp)(i,n-1) = A(i-1,n-3);
        (*temp)(i,0) = A(i-1,0);
    }
    for(int j=1;j<n-1;++j){
        (*temp)(0,j) = A(0,j-1);
        (*temp)(m-1,j) = A(m-3,j-1);
    }
    for(int i=1;i<m-1;++i){
        for(int j=1;j<n-1;++j){
            (*temp)(i,j) = A(i-1,j-1);          
        }
    }
    return *temp;
}

class OSCAR{
    private:
        Tensor<double,m,n>* U=nullptr;
        Tensor<double,m,n>* V=nullptr;
    public:
        OSCAR(){

        }
        template<int m,int n> Tensor<double,m,n>* get(int k){
            auto temp = new Tensor<double,m,n> [2];
            temp[0] = U[k];
            temp[1] = V[k];
            return temp;
        }
};

DEFMETHOD_DLD  (FasterSolver,interp,args,nargout,"This is a C++ implementation of the core PDE solver using the Fastor library."){
    //checking number of input parameters
    if(args.length()<12){
        octave_stdout << "Invalid number of input arguments";
        return NULL;
    }

    //parsing input
    const int K = args(8).int_value();
    const int r = args(9).int_value();
    double t0 = args(2).double_value();
    double finalt = args(3).double_value();
    double dx = args(4).double_value();
    double dy = args(5).double_value();
    double dt = args(6).double_value();
    Tensor<double,m,n> x = OctMat2Ten<m,n>(args(0).matrix_value());
    Tensor<double,m,n> y = OctMat2Ten<m,n>(args(1).matrix_value());
    Tensor<double,m,n> pdf0 = OctMat2Ten<m,n>(args(7).matrix_value());
    SparseMatrix M1 = args(10).sparse_matrix_value();
    SparseMatrix M2 = args(11).sparse_matrix_value();

    bool flag = false;  
    auto pdf = pdf0;
    int k=1;
    int c=1;
    int velo_track = 1;
    Tensor<double,m,n>* sol = nullptr;
    sol = new Tensor<double,m,n> [K] {0};
    sol[0] = pdf0;
    int jump = std::floor((finalt-t0)/dt)/(K-2);
    double t = t0;

    auto velo = OSCAR();

    while(t<finalt){
        auto u = velo.get(velo_track);
        auto u1 = velo.get(velo_track+1);
        auto vp = gradient(u[0],dy,dx);
        auto phi = evaluate(x-0.5*0.5*dt*(u[0]+u1[0])+0.5*(0.5*dt)*(0.5*dt)*(u1[0]*vp[1]+u1[1]*vp[0]));
        vp = gradient(u[1],dy,dx);
        auto psi = evaluate(y-0.5*0.5*dt*(u[1]+u1[1])+0.5*(0.5*dt)*(0.5*dt)*(u1[0]*vp[1]+u1[1]*vp[0]));
        //pdf = interpolation
        pdf = pdf/(sum(pdf)*dx*dy);

        auto pdf_t = Ten2OctMat(transpose(pdf));
        auto Vright = M2*Ten2OctMat(flatten(pdf_t));
        octave_value_list ml_arg;
        ml_arg(0) = octave_value(M1);
        ml_arg(1) = octave_value(Vright);
        auto retv = interp.feval("mldivide",ml_arg,1);
        auto pdf1 = transpose(reshape<n-2,m-2>(OctMat2Ten(retv(0).matrix_value())));
        pdf = ExtendBoundary<m,n>(pdf1);

        delete[] u;
        u=u1;
        u1 = velo.get(velo_track+2);

    }
}