#ifndef OSCAR_H_
#define OSCAR_H_

#include<octave/oct.h>
#include<octave/octave.h>
#include<octave/parse.h>

class OSCAR{
    private:
        Matrix* U=nullptr;
        Matrix* V=nullptr;
    public:
        OSCAR() = delete;
        OSCAR(octave::interpreter& interp){
            auto value_list = interp.feval("DataWrapper",octave_value_list(),2);
            NDArray uf = value_list(0).array_value();
            NDArray vf = value_list(1).array_value();
            int P = (uf.dims())(2);
            U = new Matrix [P];
            V = new Matrix [P];
            for(int i=0;i<P;++i){
                U[i] = uf.page(i).as_matrix();
                V[i] = vf.page(i).as_matrix();
            }
        }
        Matrix const ** get(int k){
            Matrix const ** temp = new Matrix const * [2];
            temp[0] = &(U[k]);
            temp[1] = &(V[k]);
            return temp;
        }
        ~OSCAR(){
            delete[] U;
            delete[] V;
        }
};
#endif