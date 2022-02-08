#include <iostream>
#include <unordered_map>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Network/Session.h>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Crypto/AES.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Crypto/RandomOracle.h>
#include <cryptoTools/Common/BitVector.h>
#include "cryptoTools/Common/BitIterator.h"
//adding for the transpose
#include <cryptoTools/Common/Matrix.h>
#include <cryptoTools/Common/MatrixView.h>
#include "libOTe/Tools/Tools.h"
//using encryption
#include <cryptoTools/Crypto/Rijndael256.h>


#include <algorithm>
#include "libOTe/Base/BaseOT.h"
#include "libOTe/Base/MasnyRindal.h"
/*
#include <fss-common.h>
#include <fss-client.h>
#include <fss-server.h>
*/
#include "fss.h"
#include "frontend/libPaXoS/ObliviousDictionary.h"


using namespace osuCrypto; 

void fss_psi(int nsquares, std::vector<uint64_t> inputs_x, std::vector<uint64_t> inputs_y, std::vector<uint64_t> inputs);
void Transpose_View_Test();
//void fssPSI(uint64_t sender_inputs); 
void fuzzyPSI(u64 keysize, u64 y_size, u64 x_volume); 
void fuzzyPSI_int(u64 keysize); 
void basic_transpose();