#include<octave/oct.h>
#include<octave/octave.h>
#include<octave/parse.h>
#include<chrono>
#include"OSCAR.h"
#include <boost/interprocess/managed_shared_memory.hpp>

using namespace boost::interprocess;

NDArray LowRankApprox(octave::interpreter& interp,const Matrix& pdf,int r){
    octave_value_list arg_temp;
    arg_temp(0) = octave_value(pdf);
    arg_temp(1) = octave_value(r);
    return interp.feval("LowRankApprox",arg_temp,1)(0).array_value();
}

void SaveToSol(NDArray& sol,NDArray pdf,int page){
    Array<octave_idx_type> idx (dim_vector (3,1),0);
    idx(0) = page;
    int r = pdf.dims()(0);
    int l = pdf.dims()(1);
    pdf.resize(dim_vector (1,r,l));
    sol.insert(pdf,idx);
}

DEFMETHOD_DLD (SFinit,interp,args,nargout,"Initialize the shared memory"){
    managed_shared_memory segment(open_only, "OSCAR");
    OSCAR* wrapper = segment.construct<OSCAR>("OSCARwrapper")(interp);
    if(wrapper == 0){
        octave_stdout<<"Wrapper Creation Failed.\n";
    }
    return octave_value_list();
}

DEFMETHOD_DLD (SFclose,interp,args,nargout,"Inform the shared memory end of subprocess"){
    managed_shared_memory segment(open_only, "OSCAR");
    int* counter = segment.find<int>("PC").first;
    ++(*counter);
    if(*counter == 8){
        segment.destroy<OSCAR>("OSCARwrapper");
        segment.destroy<int>("PC");
    }
    return octave_value_list();
}

DEFMETHOD_DLD  (SFsolver,interp,args,nargout,"This is a slight improvement of the core solver by isolating the core part and move the iteration to C++."){
    octave_stdout << "SFsolver triggered\n";
    auto start_time = std::chrono::steady_clock::now();
    if(args.length()<12){
        octave_stdout << "Invalid number of input arguments";
        return octave_value_list();
    }

    //parsing input
    const int K = args(8).int_value();
    const int r = args(9).int_value();
    double t0 = args(2).double_value();
    double finalt = args(3).double_value();
    double dx = args(4).double_value();
    double dy = args(5).double_value();
    double dt = args(6).double_value();
    Matrix x = args(0).matrix_value();
    Matrix y = args(1).matrix_value();
    Matrix pdf0 = args(7).matrix_value();
    SparseMatrix M1 = args(10).sparse_matrix_value();
    SparseMatrix M2 = args(11).sparse_matrix_value();

    int m = x.dims()(0);
    int n = x.dims()(1);

    bool flag = false;  
    auto pdf = pdf0;
    int k=1;
    int c=1;
    int velo_track = 1;
    dim_vector dim(K,r,m+n+1);
    NDArray sol (dim,0);
    int jump = std::floor((finalt-t0)/dt)/(K-2);
    double t = t0;

    octave_value_list argsin;
    argsin(0) = octave_value(pdf);
    argsin(1) = octave_value(x);
    argsin(2) = octave_value(y);
    argsin(3) = octave_value(dx);
    argsin(4) = octave_value(dy);
    argsin(5) = octave_value(dt);
    argsin(6) = octave_value(M1);
    argsin(7) = octave_value(M2);

    managed_shared_memory segment(open_only, "OSCAR");
    //offset_ptr<MatVec> U = segment.find<MatVec>("OSCAR_U").first;
    //offset_ptr<MatVec> V = segment.find<MatVec>("OSCAR_V").first;

    //octave_stdout << (*U)[1];
    auto velo = segment.find<OSCAR>("OSCARwrapper").first;

    if(velo == 0||velo==nullptr){
        octave_stdout << "Fail to obtain the velo info\n";   
    }else{
        octave_stdout << velo->get(1);
        // auto temp = velo->get(1);
        // delete[] temp;
    }

    while(t<finalt){
        //const Matrix* temp[6];
        for(int k=0;k<3;++k){
            auto temp = velo->get(velo_track+k);
            argsin(8+2*k) = octave_value(temp[0]);
            argsin(8+2*k+1) = octave_value(temp[1]);
            delete[] temp;
        }
        auto res = interp.feval("SubSolver",argsin,1);
        
        pdf = res(0).matrix_value();
        argsin(0) = octave_value(pdf);

        t += dt;
        if(t>=finalt && !flag){
            t -= dt;
            dt = finalt - t;
            t = finalt;
            flag = true;
        }

        if((k-1)%jump == 0 && c<K-1){
            auto temp = LowRankApprox(interp,pdf,r);
            SaveToSol(sol,temp,c);
            c+=1;
        }
        k+=1;
    }
    auto temp = LowRankApprox(interp,pdf,r);
    SaveToSol(sol,temp,K-1);

    octave_value_list retr;
    retr(0) = octave_value(sol);

    auto end_time = std::chrono::steady_clock::now();
    octave_stdout<<"Solver Time takes (s): "<<std::chrono::duration_cast<std::chrono::seconds>(end_time-start_time).count()<<"\n";

    return retr;
}