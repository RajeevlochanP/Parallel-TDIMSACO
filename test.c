#include <stdio.h>
#include <omp.h>
int main()
{

#pragma omp parallel num_threads(4)
    {
        int id = omp_get_thread_num();
        if (id == 0){
            return 0; 
        }
        #pragma omp barrier
        printf("After barrier %d\n", id);
    }
    return 0;
}