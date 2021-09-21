#include "fss.h"

void fss_interval(uint64_t a, uint64_t b){
    uint64_t lboundary = a;
    uint64_t rboundary = b;
    uint64_t l_lt_ans0, l_lt_ans1, l_lt_fin, r_lt_ans0, r_lt_ans1, r_lt_fin;
    Fss l_fClient, l_fServer, r_fClient, r_fServer;

    ServerKeyLt l_lt_k0;
    ServerKeyLt l_lt_k1;
    ServerKeyLt r_lt_k0;
    ServerKeyLt r_lt_k1;

    // Client does the key generation for the FSS

    // left interval (a, -1), x < a then - 1, x >= a then 0 
    initializeClient(&l_fClient, 64, 2);
    generateTreeLt(&l_fClient, &l_lt_k0, &l_lt_k1, a, -1);

    //right interval (b, 1), x < b then 1, x >= b then 0
    initializeClient(&r_fClient, 64, 2);
    generateTreeLt(&r_fClient, &r_lt_k0, &r_lt_k1, b, 1);

    //Server does the FSS evaluations
    initializeServer(&l_fServer, &l_fClient);
    initializeServer(&r_fServer, &r_fClient);

    // Checks to see if this works 

    l_lt_ans0 = evaluateLt(&l_fServer, &l_lt_k0, (a-1));
    l_lt_ans1 = evaluateLt(&l_fServer, &l_lt_k1, (a-1));
    l_lt_fin = l_lt_ans1 - l_lt_ans0;
    r_lt_ans0 = evaluateLt(&r_fServer, &r_lt_k0, (a-1));
    r_lt_ans1 = evaluateLt(&r_fServer, &r_lt_k1, (a-1));
    r_lt_fin = r_lt_ans0 - r_lt_ans1;
    
    cout << "FSS Lt Match (a - 1) " << l_lt_fin - r_lt_fin << endl;

    l_lt_ans0 = evaluateLt(&l_fServer, &l_lt_k0, a + 1);
    l_lt_ans1 = evaluateLt(&l_fServer, &l_lt_k1, a + 1);
    l_lt_fin = l_lt_ans0 - l_lt_ans1;
    r_lt_ans0 = evaluateLt(&r_fServer, &r_lt_k0, a + 1);
    r_lt_ans1 = evaluateLt(&r_fServer, &r_lt_k1, a + 1);
    r_lt_fin = r_lt_ans1 - r_lt_ans0;
    
    cout << "FSS Lt Match a + 1 " << l_lt_fin - r_lt_fin << endl;

    l_lt_ans0 = evaluateLt(&l_fServer, &l_lt_k0, b-1);
    l_lt_ans1 = evaluateLt(&l_fServer, &l_lt_k1, b-1);
    l_lt_fin = l_lt_ans0 - l_lt_ans1;
    r_lt_ans0 = evaluateLt(&r_fServer, &r_lt_k0, (b-1));
    r_lt_ans1 = evaluateLt(&r_fServer, &r_lt_k1, (b-1));
    r_lt_fin = r_lt_ans1 - r_lt_ans0;
    
    cout << "FSS Lt Match (b - 1) " << l_lt_fin - r_lt_fin << endl;

     l_lt_ans0 = evaluateLt(&l_fServer, &l_lt_k0, b+1);
    l_lt_ans1 = evaluateLt(&l_fServer, &l_lt_k1, b+1);
    l_lt_fin = l_lt_ans0 - l_lt_ans1;
    r_lt_ans0 = evaluateLt(&r_fServer, &r_lt_k0, (b+1));
    r_lt_ans1 = evaluateLt(&r_fServer, &r_lt_k1, (b+1));
    r_lt_fin = r_lt_ans1 - r_lt_ans0;
    
    cout << "FSS Lt Match (b + 1) " << l_lt_fin - r_lt_fin << endl;

}

void fss_parallel(int n, uint64_t l_boundary, uint64_t r_boundary){
    auto nt = std::thread::hardware_concurrency();
    std::cout << "concurrent threads allowed is " << nt << std::endl; 
    auto routine = [&](int threadIdx){
        auto begin = n * threadIdx / nt;
        auto end = n * (threadIdx +1) / nt;
        for(int i = begin; i < end; ++i)
            fss_interval(l_boundary, r_boundary);
    };
    std::vector<std::thread> thrds(nt);
    for(uint64_t i =0; i < nt; ++i)
       thrds[i] = std::thread(routine, i);
    for(uint64_t i =0; i < nt; ++i)
       thrds[i].join();
    }


void fss_interval2(uint64_t a, uint64_t b){
    uint64_t lboundary = a;
    uint64_t rboundary = b;
    uint64_t l_lt_ans0, l_lt_ans1, l_lt_fin, r_lt_ans0, r_lt_ans1, r_lt_fin;
    Fss l_fClient, l_fServer, r_fClient, r_fServer;

    ServerKeyLt l_lt_k0;
    ServerKeyLt l_lt_k1;
    ServerKeyLt r_lt_k0;
    ServerKeyLt r_lt_k1;

    
    // Client does the key generation for the FSS

    // left interval (a, -1), x < a then - 1, x >= a then 0 
    common_test();
    initializeClient(&l_fClient, 64, 2);
    generateTreeLt(&l_fClient, &r_lt_k0, &r_lt_k1, b, 1);
   
}