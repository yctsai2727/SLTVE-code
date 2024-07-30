#ifndef OSCAR_H_
#define OSCAR_H_

#include<octave/oct.h>
#include<octave/octave.h>
#include<octave/parse.h>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/offset_ptr.hpp>
#include <boost/interprocess/containers/vector.hpp>

using namespace boost::interprocess;

typedef allocator<Matrix, managed_shared_memory::segment_manager>  ShmemAllocator;
typedef vector<Matrix, ShmemAllocator> MatVec;

class OSCAR{
    private:
        //offset_ptr<MatVec> U=nullptr;
        //offset_ptr<MatVec> V=nullptr;
        int m,n,P;
    public:
        OSCAR() = delete;
        OSCAR(octave::interpreter& interp){
            octave_stdout << "Constructing OSCAR wrapper\n";
            auto value_list = interp.feval("DataWrapper",octave_value_list(),2);
            NDArray uf = value_list(0).array_value();
            NDArray vf = value_list(1).array_value();
            auto dim = uf.dims();
            m = dim(0);
            n = dim(1);
            P = dim(2);
            octave_stdout<< "Pages No.: " << P<<"\n";

            managed_shared_memory segment(open_only, "OSCAR");

            segment.construct<NDArray>("OSCAR_U")(uf);
            segment.construct<NDArray>("OSCAR_V")(vf);
        
            // const ShmemAllocator alloc_inst (segment.get_segment_manager());
            
            // offset_ptr<MatVec> U = segment.construct<MatVec>("OSCAR_U")(alloc_inst);
            // offset_ptr<MatVec> V = segment.construct<MatVec>("OSCAR_V")(alloc_inst);
            
            // for(int i=0;i<P;++i){
            //     U->push_back(uf.page(i).as_matrix());
            //     V->push_back(vf.page(i).as_matrix());
            // }
            
            octave_stdout << "Construction finished\n";
        }
        Matrix* get(int k){
            managed_shared_memory segment(open_only, "OSCAR");
            auto U = segment.find<NDArray>("OSCAR_U").first;
            auto V = segment.find<NDArray>("OSCAR_V").first;
            if(U==0 || V==0){
                octave_stdout << "Fail to obtain velo data.\n";
                return nullptr;
            }
            Matrix* temp = new Matrix[2];
            temp[0] = (U->page(k).as_matrix());
            temp[1] = (V->page(k).as_matrix());
            return temp;
        }
        const Matrix get_u(int k){
            managed_shared_memory segment(open_only, "OSCAR");
            auto U = segment.find<NDArray>("OSCAR_U").first;
            return U->page(k).as_matrix();
        }
        // const Matrix get_v(int k){
        //     managed_shared_memory segment(open_only, "OSCAR");
        //     offset_ptr<MatVec> V = segment.find<MatVec>("OSCAR_V").first;
        //     if(V==0){
        //         octave_stdout << "Fail to obtain velo data.\n";
        //         return Matrix();
        //     }
        //     return V->at(k);
        // }
        void check(){
            octave_stdout << "Checked.\n";
            return;
        }
        ~OSCAR(){
            managed_shared_memory segment(open_only, "OSCAR");
            segment.destroy<MatVec>("OSCAR_U");
            segment.destroy<MatVec>("OSCAR_V");
        }
};
#endif