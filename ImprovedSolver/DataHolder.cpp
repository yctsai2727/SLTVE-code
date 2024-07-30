#include <boost/interprocess/managed_shared_memory.hpp>
#include<stdlib.h>
#include<iostream>

using namespace boost::interprocess;

int main(int argc, char *argv[]){
    struct shm_remove{
         shm_remove() {  shared_memory_object::remove("OSCAR"); }
         ~shm_remove(){  shared_memory_object::remove("OSCAR"); }
    } remover;

    size_t kb = 1024;
    size_t gb_4 = kb*kb*kb*4;

    std::cout<<"Memory Allocated: "<<gb_4<<"bytes"<<std::endl;
    
    managed_shared_memory segment(create_only, "OSCAR", gb_4);

    system("cd .. & octave SF_init_route.m");
    
    system("Pause");
}