#include "okvsfss.h"

using namespace osuCrypto;

/*
    This implementation assumes 2-dimensions
*/
// Helper function for FSS Sharing
/* How we write a square and identify the grid cells
    1. Grid cell label is bottom left co-coordinate
    2. We can identify a square also by bottom left coord

    Note: each square is processed independent of the other
    this allows us to sample the input squares systematically
    For now, we will assume that all the squares are in a straight line

    How to generate input squares
    1. Sample the bottom left square/center of the square
    2. Sample the side length = 2 * delta
    3. First square Bottomleftlabel - (delta, delta)
    4. To get next square add (2*delta) along the the x-axis
*/
void return_grid(uint64_t x, uint64_t y, uint64_t &grd_x, uint64_t &grd_y, uint64_t sqrlen){
        grd_x = x / sqrlen;
        grd_y = y / sqrlen;
}

void PaxosEncode(std::vector<uint64_t> setKeys, const std::vector<block> setValues0, const std::vector<block> setValues1, std::vector<block>& okvs0, std::vector<block>& okvs1, uint64_t fieldSize)
{
    GF2E dhBitsVal;
    int hashSize=setKeys.size(), gamma = 60, v=20;
    double c1 = 1.5;
    //vector<uint64_t> keys;
    //keys.resize(hashSize);
    int fieldSizeBytes = fieldSize % 8 == 0 ? fieldSize/8 : fieldSize/8 + 1;
    int zeroBits = 8 - fieldSize % 8;

    //initialize 
    ObliviousDictionary * dict0 = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    ObliviousDictionary * dict1 = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dict0->init();
    dict1->init();

    // hashing keys
    /*for (int i=0; i < setKeys.size(); i++){
        keys[i] = XXH64(&setKeys[i], 64, 0);
    }*/

    vector<byte> values0, values1;
    values0.resize(setValues0.size() * sizeof(block));
    values1.resize(setValues1.size() * sizeof(block));
    memcpy(values0.data(), setValues0.data(), setValues0.size() * sizeof(block));
    memcpy(values1.data(), setValues1.data(), setValues1.size() * sizeof(block));
    
    dict0->setKeysAndVals(setKeys, values0);
    dict1->setKeysAndVals(setKeys, values1);
    dict0->encode();
    dict1->encode();

    vector<byte> x0 = dict0->getVariables();
    vector<byte> x1 = dict1->getVariables();
    
    // vector of **bytes to block** , vector<block> Okvs;
    okvs0.resize(x0.size() / fieldSizeBytes);
    okvs1.resize(x1.size() / fieldSizeBytes);
    memcpy(okvs0.data(), x0.data(), x0.size());
    memcpy(okvs1.data(), x1.data(), x1.size());

    //dict0->checkOutput();
    //dict1->checkOutput();
   
} 

vector<block> share_trivialFSS(uint64_t delta, uint64_t grid_x, uint64_t grid_y, uint64_t point_x, bool x, uint64_t point_y, bool y, block &rand0, block &rand1){

    // std::cout << " Let's write the trivial FSS in 2 Dimensions " << std::endl;
    BitVector FSS_keyzero, FSS_keyone;
    FSS_keyzero.assign(rand0);
    FSS_keyone.assign(rand1);

    // (grid_x, grid_y) = grid label
    grid_x = grid_x * (2 * delta);
    grid_y = grid_y * (2 * delta);
    uint64_t int_start = 64 - (2 * delta);
 
    // Let's process the X-coord which is the second half of the block, access block.mData[0]
    // in a block(y, x)

    if (x) { // fill right
    uint64_t start_point_x = int_start + (point_x - grid_x); 
    for (int i = start_point_x + 64; i < 128; i++)
        FSS_keyone[i] = FSS_keyzero[i];
    }
    else { // fill left 
     uint64_t end_point_x = int_start + (point_x - grid_x); 
        for (int i = int_start; i <= end_point_x; i++)
            FSS_keyone[64 + i] = FSS_keyzero[64 + i];    
    }

    if (y) { // fill top
    uint64_t start_point_y = int_start + (point_y - grid_y); 
    for (int i = start_point_y; i < 64; i++)
        FSS_keyone[i] = FSS_keyzero[i];
    }
    else { // fill bottom
     uint64_t end_point_y = int_start + (point_y - grid_y); 
        for (int i = int_start; i <= end_point_y; i++)
            FSS_keyone[i] = FSS_keyzero[i];
    }

    //  converting the BitVector -> Block
    u8 * u8_key0 = FSS_keyzero.data();
    u8 * u8_key1 = FSS_keyone.data();

    //  use span function to block type thing
    block k0(u8_key0[15],u8_key0[14],u8_key0[13],u8_key0[12],u8_key0[11],u8_key0[10],u8_key0[9],u8_key0[8],u8_key0[7],u8_key0[6],u8_key0[5],u8_key0[4],u8_key0[3],u8_key0[2],u8_key0[1],u8_key0[0]);
    block k1(u8_key1[15],u8_key1[14],u8_key1[13],u8_key1[12],u8_key1[11],u8_key1[10],u8_key1[9],u8_key1[8],u8_key1[7],u8_key1[6],u8_key1[5],u8_key1[4],u8_key1[3],u8_key1[2],u8_key1[1],u8_key1[0]);
    //block k00 = FSS_keyzero.getSpan<block>();
    //block k11 = FSS_keyone.getSpan<block>();

    //  comment below line later!!
    rand1 = k1;
    vector<block> fss_shares;
    fss_shares.push_back(k0);
    fss_shares.push_back(k1);
    return fss_shares;
}


vector<array<block, 440>> share_weakFSS(uint64_t delta, uint64_t grid_x, uint64_t grid_y, uint64_t point_x, bool x, uint64_t point_y, bool y){

    //  just sample random bitvector and then send it back na, what is this circus?
    PRNG sprng(toBlock(31));
    vector<BitVector> bRand0;
    vector<BitVector> bRand1;
    BitVector FSS_keyzero, FSS_keyone;
    FSS_keyzero.assign(toBlock(11));
    FSS_keyone.assign(toBlock(12));

    for (int i = 0; i < 440; i++){
        FSS_keyzero.randomize(sprng);
        bRand0.push_back(FSS_keyzero);
        FSS_keyone.randomize(sprng);
        bRand1.push_back(FSS_keyone);
    }

    //  (grid_x, grid_y) = grid label
    grid_x = grid_x * (2 * delta);
    grid_y = grid_y * (2 * delta);
    uint64_t int_start = 64 - (2 * delta);
 
    //  Let's process the X-coord which is the second half of the block, access block.mData[0]
    //  in a block (y, x)
    if (x) {    // fill right
    uint64_t start_point_x = int_start + (point_x - grid_x); 
    for (int i = start_point_x + 64; i < 128; i++)
        for (int j = 0; j < bRand0.size(); j++)
            bRand1[j][i] = bRand0[j][i];
    }

    else {     // fill left 
    uint64_t end_point_x = int_start + (point_x - grid_x); 
    for (int i = int_start; i <= end_point_x; i++)
        for (int j = 0; j < 440; j++)
            bRand1[j][64 + i] = bRand0[j][64 + i];    
    }

    if (y) {    // fill top
    uint64_t start_point_y = int_start + (point_y - grid_y); 
    for (int i = start_point_y; i < 64; i++)
        for (int j = 0; j < bRand0.size(); j++)
            bRand1[j][i] = bRand0[j][i];
    }

    else {    // fill bottom
    uint64_t end_point_y = int_start + (point_y - grid_y); 
    for (int i = int_start; i <= end_point_y; i++)
        for (int j = 0; j < bRand0.size(); j++)
            bRand1[j][i] = bRand0[j][i];
    }
    // better way to convert BitVector to Block??
    array<block, 440> zero_shares, one_shares;
    for (int k = 0; k < bRand0.size(); k++){
        u8 * u8_key0 = bRand0[k].data();
        block blk_k0(u8_key0[15],u8_key0[14],u8_key0[13],u8_key0[12],u8_key0[11],u8_key0[10],u8_key0[9],u8_key0[8],u8_key0[7],u8_key0[6],u8_key0[5],u8_key0[4],u8_key0[3],u8_key0[2],u8_key0[1],u8_key0[0]);
        zero_shares[k] = blk_k0;
    }
    for (int k = 0; k < bRand1.size(); k++){
        u8 * u8_key1 = bRand1[k].data();
        block blk_k1(u8_key1[15],u8_key1[14],u8_key1[13],u8_key1[12],u8_key1[11],u8_key1[10],u8_key1[9],u8_key1[8],u8_key1[7],u8_key1[6],u8_key1[5],u8_key1[4],u8_key1[3],u8_key1[2],u8_key1[1],u8_key1[0]);
        one_shares[k] = blk_k1;
    }
    vector<array<block, 440>> fss_shares;
    fss_shares.push_back(zero_shares);
    fss_shares.push_back(one_shares);
    return fss_shares;

}

void fulldomainEval(unordered_map<block, uint64_t> &recv_hash, uint64_t delta, array<block, 440> fss_block, uint64_t grd_x, uint64_t grd_y, 
uint64_t point_x, bool x, uint64_t point_y, bool y){
    
    array<block, 440> fss_transpose; 
    // TRANSPOSE
    MatrixView<u8> bitView((u8*)fss_block.data(), 440, 16);
	Matrix<u8> blockView(128, 55);
    transpose(bitView, blockView);

    uint64_t int_start = 64 - (2 * delta);
    uint64_t grd_key = 0;
    RandomOracle sha_fss(sizeof(block));
    //RandomOracle sha_fss2(sizeof(block));
    block hash_output, hash_output2; 
    uint64_t pos_x = point_x - (grd_x * 2 * delta);
    uint64_t pos_y = point_y - (grd_y * 2 * delta);
    uint64_t pt_y = point_y; 
    /*bool flag = false;
    if (point_x == 29 && point_y == 29){
        flag = true;
        grd_key = grd_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_x;
        std::cout << "fss sharing eval okvs key : " << grd_key << std::endl;
        grd_key = 0;
    }*/
    if(x){
        if(y){
            for (int i = 64 + int_start + pos_x; i < 128; i++){
                for (int j = int_start + pos_y; j < 64; j++){
                    grd_key = pt_y;
                    grd_key = grd_key << 32;
                    grd_key = grd_key + point_x;
                    pt_y = pt_y + 1; 
                    //sha_fss.Update(&blockView(i, 0), 55);
                    //sha_fss.Update(&blockView(j, 0), 55);
                    sha_fss.Update((u8*)blockView[i].data(), 55);
                    sha_fss.Update((u8*)blockView[j].data(), 55);
                    sha_fss.Final(hash_output);
                    //std::cout << "hash output quad 1 " << hash_output << std::endl;
                    recv_hash.insert({hash_output, grd_key});
                    /*if (i == 127 && j == 63) {
                        //std::cout << int(blockView(i, 50)) << " " << int(blockView(i, 40)) << " " << int(blockView(j, 20)) << std::endl;
                        std::cout << "receiver i " << i << "j " << j << " key " << grd_key << " hash " << hash_output << std::endl;
                        BitVector x;
                        x.append((u8*)blockView[i].data(), 440, 0);
                        std::cout << x << std::endl;
                        BitVector y;
                        y.append((u8*)blockView[j].data(), 440, 0);
                        std::cout << y << std::endl;
                    }*/
                    //std::cout << "key " << grd_key << " hash " << hash_output << std::endl;
                    grd_key = 0;
                    sha_fss.Reset();
                }
                point_x = point_x + 1; 
                pt_y = point_y; 
            }             
        }
        else {
            for (int i = 64 + int_start + pos_x; i < 128; i++){
                for (int j = int_start + pos_y; j >= int_start; j--){
                    grd_key = pt_y;
                    grd_key = grd_key << 32;
                    grd_key = grd_key + point_x;
                    pt_y = pt_y - 1; 
                    sha_fss.Update(&blockView(i, 0), 55);
                    sha_fss.Update(&blockView(j, 0), 55);
                    //sha_fss.Update(blockView[i].data(), 55);
                    //sha_fss.Update(blockView[j].data(), 55);
                    sha_fss.Final(hash_output);
                    recv_hash.insert({hash_output, grd_key});
                    //std::cout << "i " << i << "j " << j << " key " << grd_key << " hash " << hash_output << std::endl;
                    
                    /*if (j == 53){
                     //   std::cout << int(blockView(i, 50)) << " " << int(blockView(i, 40)) << " " << int(blockView(j, 20)) << std::endl;
                        std::cout << "i " << i << "j " << j << std::endl;
                        std::cout << "key " << grd_key << std::endl;
                        std::cout << "hash " << hash_output << std::endl;
                    }*/
                    //std::cout << "key " << grd_key << " hash " << hash_output << std::endl;
                    grd_key = 0;
                    sha_fss.Reset();
                }
                point_x = point_x + 1; 
                pt_y = point_y;
            }
        }

    }
    else{
        if (y){
           for (int i = 64 + int_start + pos_x; i >= 64 + int_start; i--){
                for (int j = int_start + pos_y; j < 64; j++){
                    grd_key = pt_y;
                    grd_key = grd_key << 32;
                    grd_key = grd_key + point_x;
                    pt_y = pt_y + 1;
                    sha_fss.Update(&blockView(i, 0), 55);
                    sha_fss.Update(&blockView(j, 0), 55);
                    //sha_fss.Update(blockView[i].data(), 55);
                    //sha_fss.Update(blockView[j].data(), 55);
                    sha_fss.Final(hash_output);
                    recv_hash.insert({hash_output, grd_key});
                    /*if (i == 117){
                     //   std::cout << int(blockView(i, 50)) << " " << int(blockView(i, 40)) << " " << int(blockView(j, 20)) << std::endl;
                        std::cout << "i " << i << "j " << j << std::endl;
                        std::cout << "key " << grd_key << std::endl;
                        std::cout << "hash " << hash_output << std::endl;
                    }*/
                    grd_key = 0;
                    sha_fss.Reset();
                }
                point_x = point_x - 1; 
                pt_y = point_y;
            } 
        }
        else {
            for (int i = 64 + int_start + pos_x; i >= 64 + int_start; i--){
                for (int j = int_start + pos_y; j >= int_start; j--){
                    //std::cout << "x " << point_x << "y " << pt_y << std::endl;
                    grd_key = pt_y;
                    grd_key = grd_key << 32;
                    grd_key = grd_key + point_x; 
                    pt_y = pt_y - 1;
                    sha_fss.Update(&blockView(i, 0), 55);
                    sha_fss.Update(&blockView(j, 0), 55);
                    //sha_fss.Update(blockView[i].data(), 55);
                    //sha_fss.Update(blockView[j].data(), 55);
                    sha_fss.Final(hash_output);
                    recv_hash.insert({hash_output, grd_key});
                    
                    /*if (j == 53 && i == 117 && flag == true){
                       // std::cout << int(blockView(i, 50)) << " " << int(blockView(i, 40)) << " " << int(blockView(j, 20)) << std::endl;
                        std::cout << "x " << point_x << "y " << pt_y << std::endl;
                        std::cout << "i " << i << "j " << j << std::endl;
                        //std::cout << "FSS shareEv : key " << grd_key << std::endl;
                        std::cout << "hash " << hash_output << std::endl;
                        BitVector rx;
                        rx.append((u8*)blockView[i].data(), 440, 0);
                        //std::cout << "rx " << rx << std::endl;
                        BitVector ry;
                        ry.append((u8*)blockView[j].data(), 440, 0);
                        //std::cout << "ry " << ry << std::endl;
                    }*/
                    //std::cout << "key " << grd_key << " hash " << hash_output << std::endl;
                    grd_key = 0;
                    sha_fss.Reset();
                }
                point_x = point_x - 1;
                pt_y = point_y; 
            } 
        }

    }
}

// FSS_Share for PSI Receiver 
// output : gives you 440 instances of okvs0, okvs1 which are inputs to OT messages
uint64_t psi_FssShareEval(std::unordered_map<block, uint64_t> &recv_hash, uint64_t delta, int nSquares, array<vector<block>, 440> &okvs0, array<vector<block>, 440> &okvs1){
    //Full domain evaluation of the FSS

    //initialize variables
    uint64_t len_sqr = 2 * delta; 
    array<vector<block>, 440> okvsVal0, okvsVal1;  /// NOT SURE ABOUT THIS!!!!
    vector<uint64_t> okvsKeys; // vals0 -> fsskeys0, vals1 -> fsskeys1
    //below are variables to process a single square 
    uint64_t bl_x, bl_y, br_x, br_y, tl_x, tl_y, tr_x, tr_y, grd_bl_x, grd_bl_y, grd_br_x, grd_br_y, grd_tl_x, grd_tl_y, grd_tr_x, grd_tr_y, grd_key;
    vector<block> Okvs0, Okvs1;
    PRNG mprng(toBlock(31));

    // let's process the squares: identify the 2^dim = 4 grid cells -> compute fsskey0, fsskey1 
    // input: square1 (delta, delta), square2 (delta + 4*delta, delta) ..... square_n(delta + 4(n -1)delta, delta) ...
    bl_x = delta;
    bl_y = delta; 

    for (int i = 0; i < nSquares; i++){
        //BL 
        bl_x = delta + (i*2*len_sqr);
        bl_y = delta; // make bl_y = delta + (i*2*len_sqr)
        return_grid(bl_x, bl_y, grd_bl_x, grd_bl_y, len_sqr);
        vector<array<block, 440>> weakfss_bl = share_weakFSS(delta, grd_bl_x, grd_bl_y, bl_x, true, bl_y, true);
        fulldomainEval(recv_hash, delta, weakfss_bl[0], grd_bl_x, grd_bl_y, bl_x, true, bl_y, true);
        for (int j = 0; j < weakfss_bl[0].size(); j++){
            okvsVal0[j].push_back(weakfss_bl[0][j]); // this is where the transpose happens
            okvsVal1[j].push_back(weakfss_bl[1][j]); // do this step 440 times!
        }
        //fulldomainEval(recv_hash, delta, weakfss_bl[1], grd_bl_x, grd_bl_y, bl_x, true, bl_y, true);
        grd_key = grd_bl_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_bl_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;
        
        //TL
        tl_x = bl_x;
        tl_y = bl_y + len_sqr - 1;
        return_grid(tl_x, tl_y, grd_tl_x, grd_tl_y, len_sqr);
        vector<array<block, 440>> weakfss_tl = share_weakFSS(delta, grd_tl_x, grd_tl_y, tl_x, true, tl_y, false);
        fulldomainEval(recv_hash, delta, weakfss_tl[0], grd_tl_x, grd_tl_y, tl_x, true, tl_y, false);        
        for (int j = 0; j < weakfss_tl[0].size(); j++){
            okvsVal0[j].push_back(weakfss_tl[0][j]);
            okvsVal1[j].push_back(weakfss_tl[1][j]);
        }
        //fulldomainEval(recv_hash, delta, weakfss_tl[1], grd_tl_x, grd_tl_y, tl_x, true, tl_y, false);
        grd_key = grd_tl_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_tl_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;
        
        //BR
        br_x = bl_x + len_sqr - 1;
        br_y = bl_y;
        return_grid(br_x, br_y, grd_br_x, grd_br_y, len_sqr);
        vector<array<block, 440>> weakfss_br = share_weakFSS(delta, grd_br_x, grd_br_y, br_x, false, br_y, true);
        fulldomainEval(recv_hash, delta, weakfss_br[0], grd_br_x, grd_br_y, br_x, false, br_y, true);
        for (int j = 0; j < weakfss_br[0].size(); j++){
            okvsVal0[j].push_back(weakfss_br[0][j]);
            okvsVal1[j].push_back(weakfss_br[1][j]);
        }
        //fulldomainEval(recv_hash, delta, weakfss_br[1], grd_br_x, grd_br_y, br_x, false, br_y, true);
        grd_key = grd_br_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_br_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;

        //TR
        tr_x = bl_x + len_sqr - 1;
        tr_y = bl_y + len_sqr - 1;
        return_grid(tr_x, tr_y, grd_tr_x, grd_tr_y, len_sqr);
        vector<array<block, 440>> weakfss_tr = share_weakFSS(delta, grd_tr_x, grd_tr_y, tr_x, false, tr_y, false);
        fulldomainEval(recv_hash, delta, weakfss_tr[0], grd_tr_x, grd_tr_y, tr_x, false, tr_y, false);
        for (int j = 0; j < weakfss_tr[0].size(); j++){
            okvsVal0[j].push_back(weakfss_tr[0][j]);
            okvsVal1[j].push_back(weakfss_tr[1][j]);
        }
        //fulldomainEval(recv_hash, delta, weakfss_tr[1], grd_tr_x, grd_tr_y, tr_x, false, tr_y, false);
        grd_key = grd_tr_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_tr_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;
        
    }

    for (int i = 0; i < okvsVal0.size(); i++){
        PaxosEncode(okvsKeys, okvsVal0[i], okvsVal1[i], okvs0[i], okvs1[i], 128);
    }
    // BELOW CHECKS DECODE WORKS, for 100 squares and c = 2.4
    /*std::cout << "key " << okvsKeys[3] << std::endl;
    std::cout << "<key, value pair> " << okvsKeys[3] << " " << okvsVal0[10][3] << " " << okvsVal1[10][3] << std::endl;
    BitVector a;
    a.assign(okvsVal1[10][3]);
    std::cout << a << std::endl;
    uint64_t dec_key = okvsKeys[3];
    int gamma = 60, v=20, fieldSizeBytes = 16, fieldSize = 128; // for okvs
    double c1 = 1.5; // for okvs
    vector<BitVector> okvsbits0, okvsbits1;
    ObliviousDictionary * dict = new OBD3Tables(okvsKeys.size(), c1, fieldSize, gamma, v);
    dict->init();
    auto indices = dict->dec(dec_key);
    std::cout << "indices " << indices[0] << std::endl;
    BitVector bits1, bits3, bits4;
    //for (int k = 0; k < 10; k++){
    bits1.assign(okvs1[10][indices[0]]);
    BitVector bits = bits1;
    bits3.assign(okvs1[10][indices[1]]);
    bits ^= bits3;
    bits4.assign(okvs1[10][indices[2]]);
    bits ^= bits4;
    //std::cout << "bits1 " << bits1 << std::endl;
    //std::cout << "bits2 " << bits3 << std::endl;
    //std::cout << "bits3 " << bits4 << std::endl;
    std::cout << bits << std::endl;
    if (bits == a)
        std::cout << "decode works " << std::endl;
    */
    return okvsKeys.size();
}

//FSS_Eval for the PSI Sender, #baseOT OKVS instances, returns the SHA evaluation, without using GF2E
/*void psiSender_FssEval(vector<uint64_t> inputs_x, vector<uint64_t> inputs_y, vector<uint64_t> inputs, Matrix<block> fss_keys, uint64_t delta, uint64_t hashSize){ 
    
    uint64_t len_sqr = 2 * delta;
    // setting OKVS parameters
    int gamma = 60, v=20, fieldSizeBytes = 16, fieldSize = 128; // for okvs
    double c1 = 1.3; // for okvs
    ObliviousDictionary * dict = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dict->init();
    std::cout << "OKVS size " << okvs.size() << std::endl;
    
    uint64_t okvs_key, yx_share, grd_x, grd_y, x, y; 
    x = x_coord;
    y = y_coord;
    return_grid(x, y, grd_x, grd_y, len_sqr);
    okvs_key = grd_y;
    okvs_key = okvs_key << 32;
    okvs_key = okvs_key + grd_x;
    std::cout << "OKVS key " << okvs_key << std::endl;
    
    // start by learning which indices we need for a key
    // decoding from checkoutput() OKVS
    //initialize just to be able to dec() -- below 3 calls only set up the hash functions
    //auto t1 = high_resolution_clock::now();

    
    vector<BitVector> okvsbits;
        //BitVector eval;
    auto indices = dict->dec(okvs_key);
    std::cout << "INDICES " << indices[0] << " " << indices[1] << " " << indices[2] << std::endl;
    //std::cout << "each vector size " << okvs[indices[0]].size() << " " << okvs[indices[1]].size() << " " << okvs[indices[2]].size() << std::endl;
    for (int i = 0; i < okvs.size(); i++)
        okvsbits.push_back(okvs[i][indices[0]]);
    //okvsbits = okvs[indices[0]];
    
    for (int j = 1; j < 3; j++){
        for (int i = 0; i < okvs.size(); i++){
            okvsbits[i] = okvsbits[i] ^ okvs[i][indices[j]];
        }
    }
    std::cout << "# of OKVS " << okvsbits.size() << std::endl;

    RandomOracle sha_fss(sizeof(block));
    block fss_output;

    // computing OKVS key
    //computing the starting indices within block to access the shares
    int int_start_y = 63 - (2 * delta);
    int int_start_x = 64 + int_start_y;
    grd_y = grd_y * 2 * delta;
    grd_x = grd_x * 2 * delta; 
    std::cout << "y " << int_start_y + (y - grd_y) << std::endl;
    std::cout << "x " << int_start_x + (x - grd_x) << std::endl;
    for (int i = 0; i < okvsbits.size(); i++){ // hardcoding in the size of okvsbits = # of baseOt

        u8 y_share = okvsbits[i][int_start_y + (y - grd_y)];
        u8 x_share = okvsbits[i][int_start_x + (x - grd_x)];
        sha_fss.Update(x_share);
        sha_fss.Update(y_share);
    }
    sha_fss.Final(fss_output);
       //sha_fss.Reset();
    
    std::cout << "Eval for points " << fss_output << std::endl;
}*/

//FSS_Eval for PSI Sender, #baseOT OKVS instances, returns the SHA evaluation    
void psi_FssEval(uint64_t x_coord, uint64_t y_coord, array<vector<osuCrypto::block>, 440> okvs, uint64_t delta, uint64_t hashSize){ 
    
    uint64_t nbaseOT = 2; // this will be 440, can be sent a parameter
    uint64_t len_sqr = 2 * delta;
    // setting OKVS parameters
    int gamma = 60, v=20, fieldSizeBytes = 16, fieldSize = 128; // for okvs
    double c1 = 1.3; // for okvs

    uint64_t okvs_key, yx_share, grd_x, grd_y, x, y; 
    x = x_coord;
    y = y_coord;
    for (int l = 0; l < 1; l++){ 
        return_grid(x, y, grd_x, grd_y, len_sqr);
        okvs_key = grd_y;
        okvs_key = okvs_key << 32;
        okvs_key = okvs_key + grd_x;
    }
    // start by learning which indices we need for a key
    // decoding from checkoutput() OKVS
    //initialize just to be able to dec() -- below 3 calls only set up the hash functions
    ObliviousDictionary * dict = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dict->init();
    

    // handling okvs matrix
    // converting block to bytes
    vector<GF2EVector> okvs_gf2e;
    for (int k = 0; k < okvs.size(); k++) {
        vector<byte> okvs_idx;
        GF2EVector tmp_vec;
        GF2X temp;
        okvs_idx.resize(okvs[k].size() * sizeof(block)); // size can also be 440 * sizeof(block)
        memcpy(okvs_idx.data(), okvs[k].data(), okvs[k].size() * sizeof(block));
        for (int i = 0; i < okvs[k].size(); i++){
            GF2XFromBytes(temp, okvs_idx.data() + i*fieldSizeBytes, fieldSizeBytes);
            tmp_vec.push_back(to_GF2E(temp)); 
        }
        okvs_gf2e.push_back(tmp_vec);
        //cout << " okvs size " << k << " " << okvs_gf2e[k].size() << std::endl;
    }

    // not sure if we can add two variables of type of GF2EVector, might have to iterate and add
    if (okvs_gf2e.size() == 0)
        std::cout << "ERROR: okvs_gf2e is empty!" << std::endl;
    GF2EVector dhBitsVals;
    GF2E eval;
    eval = 0; 
    for (int l = 0; l < 1; l++){
        dhBitsVals.clear(); 
        auto indices = dict->dec(okvs_key);
        for (int i = 0; i < okvs_gf2e[indices[0]].size(); i++)
            dhBitsVals.push_back(okvs_gf2e[indices[0]][i]);
        for (int i = 1; i < 3; i++){
            for (int j = 0; j < dhBitsVals.size(); j++){
                dhBitsVals[j] += okvs_gf2e[indices[i]][j];
            }
        }
        // this is just for simulation purposes, remove later for a single evaluation
    }
    RandomOracle sha_fss(sizeof(block));
    block fss_output;
    // computing OKVS key

    
    //computing the starting indices within block to access the shares
    int int_start_y = 63 - (2 * delta);
    int int_start_x = 64 + int_start_y;
    grd_y = grd_y * 2 * delta;
    grd_x = grd_x * 2 * delta; 
    //BitVector val_in_bits
    //std::cout << "dhbitsvals size " << dhBitsVals.size() << std::endl;
    for (int l = 0; l < 1; l++){
    
        for (int i = 0; i < 440; i++){ // hardcoding in the size of dhBitsVals
            vector<byte> valBytes(fieldSizeBytes);
            BytesFromGF2X(valBytes.data(), rep(dhBitsVals[i]), fieldSizeBytes);
            BitIterator iter_y(valBytes.data(), int_start_y + (y - grd_y));
            BitIterator iter_x(valBytes.data(), int_start_x + (x - grd_x));
            u8 y_share, x_share;
            y_share = *iter_y;
            x_share = *iter_x;
            sha_fss.Update(y_share);
            sha_fss.Update(x_share);
        }
        sha_fss.Final(fss_output);
        sha_fss.Reset();
    }
    std::cout << "Eval for points " << fss_output << std::endl;
}

vector<vector<BitVector>> blocks_to_bits(array<vector<block>, 440> okvs){
    vector<vector<BitVector>> bitVector_okvs;
    for (int k = 0; k < okvs.size(); k++) {
        ///vector<byte> okvs_idx;
        vector<BitVector> tmp_vec;
        BitVector temp;
        for (int i = 0; i < okvs[k].size(); i++){
            temp.assign(okvs[k][i]);
            tmp_vec.push_back(temp); 
        }
        bitVector_okvs.push_back(tmp_vec);
    }
    return bitVector_okvs;
}


// FSS_Eval optimized for Full Domain Evalutation by the PSI Receiver 
/*
void psiRecver_FssEval(uint64_t nSquares, vector<vector<osuCrypto::block>> okvs, uint64_t delta, uint64_t nkeys){
    
    uint64_t len_sqr = 2 * delta;
    std::unordered_map<block, uint64_t> recv_hash; 

    // setting OKVS parameters
    int gamma = 60, v=20, fieldSizeBytes = 16, fieldSize = 128; // for okvs
    double c1 = 1.3; // for okvs

    // handling okvs matrix
    // converting block to bitvector
    vector<vector<BitVector>> bitVector_okvs;
    for (int k = 0; k < okvs.size(); k++) {
        vector<BitVector> tmp_vec;
        BitVector temp;
        for (int i = 0; i < okvs[k].size(); i++){
            temp.assign(okvs[k][i]);
            tmp_vec.push_back(temp); 
        }
        bitVector_okvs.push_back(tmp_vec);
    }
    std::cout << "bitvector okvs param " << bitVector_okvs.size() << "  " << bitVector_okvs[0].size() << std::endl;
   
    //initialize just to be able to dec() -- below 3 calls only set up the hash functions
    ObliviousDictionary * dict = new OBD3Tables(nkeys, c1, fieldSize, gamma, v);
    dict->init();
    
    // below we are exploiting the known structure of the input set; we know the bottom left label for each square is at the centre of a grid cell
    // TODO: can be generalized to any input set
    uint64_t bl_x, bl_y, br_x, br_y, tl_x, tl_y, tr_x, tr_y, grd_bl_x, grd_bl_y, grd_br_x, grd_br_y, grd_tl_x, grd_tl_y, grd_tr_x, grd_tr_y, grd_key;
    int start_y = 63 - (2 * delta);
    int start_x = 64 + start_y;
    int end_y = 64;
    int end_x = 128;
    for (int i = 0; i < nSquares; i++){
        //BL 
        bl_x = delta + (i*2*len_sqr);
        bl_y = delta; // make bl_y = delta + (i*2*len_sqr)
        return_grid(bl_x, bl_y, grd_bl_x, grd_bl_y, len_sqr);
        grd_key = grd_bl_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_bl_x;
        // single okvs key -> set of fixed indices 
        // 1. here we have okvs decryption call using <grd_key>
        auto indices = dict->dec(okvs_key);
        vector<BitVector> okvsbits;
        okvsbits = okvs[indices[0]];
        for (int i = 1; i < 3; i++)
            for (int j = 0; j < okvs[indices[i]].size(); j++)
                okvsbits[j] = okvsbits[j] ^ okvs[indices[i]][j];

        RandomOracle sha_fss(sizeof(block));
        block fss_output;
    
        //computing the starting indices within block to access the shares
        grd_bl_y = grd_bl_y * 2 * delta;
        grd_bl_x = grd_bl_x * 2 * delta; 
        for (){
            for (){
                for (int i = 0; i < 440; i++){ // hardcoding in the size of okvsbits = # of baseOt
                    u8 y_share = okvsbits[i][start_y + (bl_y - grd_bl_y)];
                    u8 x_share = okvsbits[i][start_x + (bl_x - grd_bl_x)];
                    sha_fss.Update(y_share);
                    sha_fss.Update(x_share);
                }
                sha_fss.Final(fss_output);
            }
        }
        grd_key = 0;

        //TL
        tl_x = bl_x;
        tl_y = bl_y + len_sqr;
        return_grid(tl_x, tl_y, grd_tl_x, grd_tl_y, len_sqr);
        grd_key = grd_tl_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_tl_x;
        // single okvs key -> set of fixed indices 
        // 1. here we have okvs decryption call using <grd_key>
        grd_key = 0;

        //BR
        br_x = bl_x + len_sqr;
        br_y = bl_y;
        return_grid(br_x, br_y, grd_br_x, grd_br_y, len_sqr);
        grd_key = grd_br_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_br_x;
        // single okvs key -> set of fixed indices 
        // 1. here we have okvs decryption call using <grd_key>
        grd_key = 0;

        
        //TR
        tr_x = bl_x + len_sqr;
        tr_y = bl_y + len_sqr;
        return_grid(tr_x, tr_y, grd_tr_x, grd_tr_y, len_sqr);
        grd_key = grd_tr_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_tr_x;
        // single okvs key -> set of fixed indices 
        // 1. here we have okvs decryption call using <grd_key>
        grd_key = 0;

    } 
     
    
}
*/

// FSS for a single OKVS and reading the trivial FSS share
void far_apart_FssEval(uint64_t x_coord, uint64_t y_coord, vector<block> okvs, uint64_t delta, uint64_t hashSize){
    //std::cout << "we are now evaluating the FSS, far apart " << std::endl;  
    uint64_t len_sqr = 2* delta;
    uint64_t grd_x, grd_y, x, y; 
    int gamma = 60, v=20, fieldSizeBytes = 16, fieldSize = 128; // for okvs
    double c1 = 1.3; // for okvs

    uint64_t okvs_key, yx_share;
    
     //handling the key
    x = x_coord;
    y = y_coord;
    return_grid(x, y, grd_x, grd_y, len_sqr);
  
    okvs_key = grd_y;
    okvs_key = okvs_key << 32;
    okvs_key = okvs_key + grd_x;
    
    //std::cout << "eval side key " << okvs_key << std::endl;
    //okvsKeys.push_back(grd_key);
    //    grd_key = 0;
    //block okvskey(toBlock(grd_y, grd_x));
    //std::cout << "okvskey" << okvs_key << std::endl;
   // uint64_t okvs_key = XXH64(&okvskey, 64, 0);


    //handling okvs object 
    // converting block to bytes
    vector<byte> okvs_bytes;
    okvs_bytes.resize(okvs.size() * sizeof(block));
    memcpy(okvs_bytes.data(), okvs.data(), okvs.size() * sizeof(block));
    // from bytes to GF2E!! yahooo 
    GF2EVector okvs_gf2e(okvs.size());
    GF2X temp;
    for (int i=0; i < okvs_gf2e.size(); i++){
        GF2XFromBytes(temp, okvs_bytes.data() + i*fieldSizeBytes ,fieldSizeBytes);
        okvs_gf2e[i] = to_GF2E(temp); 
        
    }
    std::cout << "trying out some tricks " << std::endl;
    BitVector tmp0, tmp1, tmp; 
    GF2E gf_tmp;
    tmp0.assign(okvs[0]);
    tmp1.assign(okvs[5]);
    tmp = tmp0 ^ tmp1;
    gf_tmp += okvs_gf2e[0];
    gf_tmp += okvs_gf2e[5];
    //std::cout << "block " << okvs[0] << "  " << okvs[5] << std::endl;
    std::cout << "tmp1 " << tmp << std::endl;
    std::cout << "dhbits " << gf_tmp << std::endl;
    //std::cout << "dhbits " << okvs_gf2e[5] << std::endl;

    // decoding from checkoutput()

    //initialize just to be able to dec() -- below 3 calls only set up the hash functions
    ObliviousDictionary * dict = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dict->init();
    BitVector val_in_bits;
    vector<byte> valBytes(fieldSizeBytes);
    block returndecode;
    int int_start_y;
    int int_start_x;
    u8 y_share, x_share;
    
    auto indices = dict->dec(okvs_key);
    //cout << "ENCODE " << okvs_key << " " << okvs_key << " " << indices[0] << " " << indices[1] << std::endl;
    GF2E dhBitsVals;
    for (int j=0; j<indices.size(); j++){
        //std::cout << indices[j] << std::endl;
        dhBitsVals += okvs_gf2e[indices[j]];
    }
    //std::cout << "dhBitsVals  " << dhBitsVals << std::endl;
    // can generalize this later
    
    BytesFromGF2X(valBytes.data(), rep(dhBitsVals), fieldSizeBytes);
    
    //std::cout << "bit15 " << int(bit15) << std::endl;
    //block returndecode(valBytes[15],valBytes[14],valBytes[13],valBytes[12],valBytes[11],valBytes[10],valBytes[9],valBytes[8],valBytes[7],valBytes[6],valBytes[5],valBytes[4],valBytes[3],valBytes[2],valBytes[1],valBytes[0]);
    //std::cout <<  "return block from okvs " << returndecode << std::endl;
    //val_in_bits.assign(returndecode); 
    //std::cout << "blk " << val_in_bits << std::endl;
    
    // UGH, below stuff somehow works out, phew
    int_start_y = 63 - (2 * delta);
    int_start_x = 64 + int_start_y;
    grd_y = grd_y * 2 * delta;
    grd_x = grd_x * 2 * delta; 
    
    
    //std::cout << "x, y " << " " << int_start_x << " " << int_start_y << std::endl;
    //std::cout << "x, y " << " " << (x - grd_x) << " " << (y - grd_y) << std::endl;

    //y_share = u8(val_in_bits[int_start_y + (y - grd_y)]);
    //x_share = u8(val_in_bits[int_start_x + (x - grd_x)]);

    //std::cout << "y, x" << y_share << x_share << std::endl;
  //  std::cout << "y share " << val_in_bits[int_start_y + (y - grd_y)] << std::endl;
  //  std::cout << "x share " << val_in_bits[int_start_x + (x - grd_x)] << std::endl;

    BitIterator iter_y(valBytes.data(), int_start_y + (y - grd_y));
    BitIterator iter_x(valBytes.data(), int_start_x + (x - grd_x));
    y_share = *iter_y;
    x_share = *iter_x;
    
    //std::cout << "(y, x) " << int(y_share) << "  " << int(x_share) << std::endl;
    
}

// UNDO THIS : need to modify fullDomain eval to work wtih a single block 
/*
void far_apart_FssShare(uint64_t delta, int nSquares, vector<block> &okvs0, vector<block> &okvs1){

    //std::cout << "balls are pairwise 3*delta apart " << std::endl;
    auto t1 = high_resolution_clock::now();

    //Full domain evaluation of the FSS
    std::unordered_map<block, uint64_t> recv_hash;

    //initialize variables
    uint64_t len_sqr = 2 * delta;
    vector<block> okvsVal0, okvsVal1, shares; // TODO:  change to vec<vec> block
    vector<uint64_t> okvsKeys; // vals0 -> fsskeys0, vals1 -> fsskeys1
    //below are variables to process a single square 
    uint64_t bl_x, bl_y, br_x, br_y, tl_x, tl_y, tr_x, tr_y, grd_bl_x, grd_bl_y, grd_br_x, grd_br_y, grd_tl_x, grd_tl_y, grd_tr_x, grd_tr_y, grd_key;
    block rand0, rand1; // TODO: vec<block> of size 440
    PRNG mprng(toBlock(31));

    // let's process the squares: identify the 2^dim = 4 grid cells -> compute fsskey0, fsskey1 
    // input: square1 (delta, delta), square2 (delta + 4*delta, delta) ..... square_n(delta + 4(n -1)delta, delta) ...
    bl_x = delta;
    bl_y = delta; 

    for (int i = 0; i < nSquares; i++){
        //BL 
        bl_x = delta + (i*2*len_sqr);
        bl_y = delta; // make bl_y = delta + (i*2*len_sqr)
        return_grid(bl_x, bl_y, grd_bl_x, grd_bl_y, len_sqr);
        rand0 = mprng.get(); //sample 440 such 
        rand1 = mprng.get(); //sample 440 such - easy
        uint64_t pos_bl_x = bl_x - (grd_bl_x * len_sqr);
        uint64_t pos_bl_y = bl_y - (grd_bl_y * len_sqr);
        fulldomainEval(recv_hash, delta, rand0, pos_bl_x, true, pos_bl_y, true); // modify the full domain eval to work for vector<block> rand0
        shares = share_trivialFSS(delta, grd_bl_x, grd_bl_y, bl_x, true, bl_y, true, rand0, rand1); // do this 440 times, once for each pair of <rand0, rand1>
        okvsVal0.push_back(shares[0]);
        okvsVal1.push_back(shares[1]); // do this step 440 times!
        grd_key = grd_bl_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_bl_x;
        okvsKeys.push_back(grd_key);
        //std::cout << "keys!! " << grd_key << " " << okvsVal0.back()  << "  " << okvsVal1.back() << std::endl;
        grd_key = 0;
        
        
        //TL
        tl_x = bl_x;
        tl_y = bl_y + len_sqr;
        return_grid(tl_x, tl_y, grd_tl_x, grd_tl_y, len_sqr);
        rand0 = mprng.get();
        rand1 = mprng.get();
        //fulldomainEval(recv_hash, rand0, true, false);
        shares = share_trivialFSS(delta, grd_tl_x, grd_tl_y, tl_x, true, tl_y, false, rand0, rand1);
        okvsVal0.push_back(shares[0]);
        okvsVal1.push_back(shares[1]);
        grd_key = grd_tl_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_tl_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;
        //okvsKeys.push_back(toBlock(grd_tl_y, grd_tl_x));
        //std::cout << "keys " << okvsKeys.back() << " " << okvsVal0.back()  << "  " << okvsVal1.back() << std::endl;
        
        //BR
        br_x = bl_x + len_sqr;
        br_y = bl_y;
        return_grid(br_x, br_y, grd_br_x, grd_br_y, len_sqr);
        rand0 = mprng.get();
        rand1 = mprng.get();
        //fulldomainEval(recv_hash, rand0, false, true);
        shares = share_trivialFSS(delta, grd_br_x, grd_br_y, br_x, false, br_y, true, rand0, rand1);
        okvsVal0.push_back(shares[0]);
        okvsVal1.push_back(shares[1]);
        grd_key = grd_br_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_br_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;
        //okvsKeys.push_back(toBlock(grd_br_y, grd_br_x));
        //std::cout << "keys " << okvsKeys.back() << " " << okvsVal0.back()  << "  " << okvsVal1.back() << std::endl;
        //std::cout << "grd id " << grd_br_x << "  " << grd_br_y << std::endl;

        //TR
        tr_x = bl_x + len_sqr;
        tr_y = bl_y + len_sqr;
        return_grid(tr_x, tr_y, grd_tr_x, grd_tr_y, len_sqr);
        rand0 = mprng.get();
        rand1 = mprng.get();
        //fulldomainEval(recv_hash, rand0, false, false);
        shares = share_trivialFSS(delta, grd_tr_x, grd_tr_y, tr_x, false, tr_y, false, rand0, rand1);
        okvsVal0.push_back(shares[0]);
        okvsVal1.push_back(shares[1]);
        grd_key = grd_tr_y;
        grd_key = grd_key << 32;
        grd_key = grd_key + grd_tr_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;
        //okvsKeys.push_back(toBlock(grd_tr_y, grd_tr_x));
        //std::cout << "keys " << okvsKeys.back() << " " << okvsVal0.back()  << "  " << okvsVal1.back() << std::endl;
        //std::cout << "grd id " << grd_tr_x << "  " << grd_tr_y << std::endl;
        

    }
    std::cout << "recv hash " << recv_hash.size() << std::endl;
    PaxosEncode(okvsKeys, okvsVal0, okvsVal1, okvs0, okvs1, 128); //
}
*/

