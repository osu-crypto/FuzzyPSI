#include "frontend/fss/fss-common.h"
#include "frontend/fss/fss-server.h"
#include "frontend/fss/fss-client.h"
#include<thread>
#include<vector>

void fss_interval(uint64_t a, uint64_t b);
void fss_parallel(int n, uint64_t l_boundary, uint64_t r_boundary);
void fss_interval2(uint64_t a, uint64_t b);