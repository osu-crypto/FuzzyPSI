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
    double c1 = 2.4;
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
    /*std::cout << "Okvs.size() " << Okvs.size() << std::endl;
    for (int i = 0; i < (Okvs.size() - 50); i++)
        std::cout << "okvs  " << i << "  " << Okvs[i] << std::endl;*/
    
    //below is function call to check output!!

    //dict0->checkOutput();
    //dict1->checkOutput();
    
    //auto indices = dict0->dec(setKeys[0]);
    //cout << "ENCODE " << setKeys[0] << " " << indices[0] << " " << indices[1] << std::endl;
   
} 
vector<block> share_trivialFSS(uint64_t delta, uint64_t grid_x, uint64_t grid_y, uint64_t point_x, bool x, uint64_t point_y, bool y, block rand0, block rand1){
    //std::cout << "Let's write the trivial FSS in 2 Dimensions " << std::endl;
    BitVector FSS_keyzero, FSS_keyone;
    FSS_keyzero.assign(rand0);
    FSS_keyone.assign(rand1);

    //(grid_x, grid_y) = grid label
    grid_x = grid_x * (2 * delta);
    grid_y = grid_y * (2 * delta);
    uint64_t int_start = 64 - (2 * delta);
 
    //Let's process the X-coord which is the second half of the block, access block.mData[0]
    // in a block(y, x)

    if (x) { // fill right
    uint64_t start_point_x = int_start + (point_x - grid_x); 
    for (int i = start_point_x; i < 64; i++)
        FSS_keyone[64 + i] = FSS_keyzero[64 + i];
    }
    else { // fill left 
     uint64_t end_point_x = int_start + (point_x - grid_x); 
        for (int i = int_start; i < end_point_x; i++)
            FSS_keyone[64 + i] = FSS_keyzero[64 + i];    
    }

    if (y) { // fill top
    uint64_t start_point_y = int_start + (point_y - grid_y); 
    for (int i = start_point_y; i < 64; i++)
        FSS_keyone[i] = FSS_keyzero[i];
    }
    else { // fill bottom
     uint64_t end_point_y = int_start + (point_y - grid_y); 
        for (int i = int_start; i < end_point_y; i++)
            FSS_keyone[i] = FSS_keyzero[i];
    }

    // converting the BitVector -> Block
    u8 * u8_key0 = FSS_keyzero.data();
    u8 * u8_key1 = FSS_keyone.data();

    // use span function to block type thing
    block k0(u8_key0[15],u8_key0[14],u8_key0[13],u8_key0[12],u8_key0[11],u8_key0[10],u8_key0[9],u8_key0[8],u8_key0[7],u8_key0[6],u8_key0[5],u8_key0[4],u8_key0[3],u8_key0[2],u8_key0[1],u8_key0[0]);
    block k1(u8_key1[15],u8_key1[14],u8_key1[13],u8_key1[12],u8_key1[11],u8_key1[10],u8_key1[9],u8_key1[8],u8_key1[7],u8_key1[6],u8_key1[5],u8_key1[4],u8_key1[3],u8_key1[2],u8_key1[1],u8_key1[0]);
    vector<block> fss_shares;
    fss_shares.push_back(k0);
    fss_shares.push_back(k1);
    return fss_shares;
}

void fulldomainEval(unordered_map<block, uint64_t> &recv_hash, uint64_t delta, block fss_block, uint64_t point_x, bool x, uint64_t point_y, bool y){
    BitVector fss_bits;
    fss_bits.assign(fss_block);
    uint64_t int_start = 64 - (2 * delta);
    uint64_t grd_key = 0;
    RandomOracle sha_fss(sizeof(block));
    block hash_output; 
    if(x){
        if(y){
            for (int i = 64 + int_start + point_x; i < 128; i++){
                for (int j = int_start + point_y; j < 64; j++){
                    grd_key = j;
                    grd_key = j << 32;
                    grd_key = grd_key + i;
                    u8 x_share = fss_bits[i];
                    u8 y_share = fss_bits[j];
                    //std::cout << "x share , y share " << int(x_share) << int(y_share) << std::endl;
                    sha_fss.Update(x_share);
                    sha_fss.Update(y_share);
                    sha_fss.Final(hash_output);
                    //std::cout << hash_output << std::endl;
                    recv_hash.insert({hash_output, grd_key});
                    grd_key = 0;
                    sha_fss.Reset();
                }
            }                
        }
        else {
            for (int i = 64 + int_start + point_x; i < 128; i++){
                for (int j = int_start; j < int_start + point_y; j++){
                    grd_key = j;
                    grd_key = j << 32;
                    grd_key = grd_key + i;
                    u8 x_share = fss_bits[i];
                    u8 y_share = fss_bits[j];
                    sha_fss.Update(x_share);
                    sha_fss.Update(y_share);
                    sha_fss.Final(hash_output);
                    recv_hash.insert({hash_output, grd_key});
                    grd_key = 0;
                    sha_fss.Reset();
                }
            }
        }

    }
    else{
        if (y){
           for (int i = 64 + int_start; i < 64 + int_start + point_x; i++){
                for (int j = int_start + point_y; j < 64; j++){
                    grd_key = j;
                    grd_key = j << 32;
                    grd_key = grd_key + i;
                    u8 x_share = fss_bits[i];
                    u8 y_share = fss_bits[j];
                    sha_fss.Update(x_share);
                    sha_fss.Update(y_share);
                    sha_fss.Final(hash_output);
                    recv_hash.insert({hash_output, grd_key});
                    grd_key = 0;
                    sha_fss.Reset();
                }
            } 
        }
        else {
            for (int i = 64 + int_start; i < 64 + int_start + point_x; i++){
                for (int j = int_start; j < int_start + point_y; j++){
                    grd_key = j;
                    grd_key = j << 32;
                    grd_key = grd_key + i;
                    u8 x_share = fss_bits[i];
                    u8 y_share = fss_bits[j];
                    sha_fss.Update(x_share);
                    sha_fss.Update(y_share);
                    sha_fss.Final(hash_output);
                    recv_hash.insert({hash_output, grd_key});
                    grd_key = 0;
                    sha_fss.Reset();
                }
            } 
        }

    }
}

// FSS_Share for PSI Receiver 
void far_apart_FssShare(uint64_t delta, int nSquares, vector<block> &okvs0, vector<block> &okvs1){

    //std::cout << "balls are pairwise 3*delta apart " << std::endl;
    auto t1 = high_resolution_clock::now();

    //Full domain evaluation of the FSS
    std::unordered_map<block, uint64_t> recv_hash;

    //initialize variables
    uint64_t len_sqr = 2 * delta;
    vector<block> okvsVal0, okvsVal1, shares;
    vector<uint64_t> okvsKeys; // vals0 -> fsskeys0, vals1 -> fsskeys1
    //below are variables to process a single square 
    uint64_t bl_x, bl_y, br_x, br_y, tl_x, tl_y, tr_x, tr_y, grd_bl_x, grd_bl_y, grd_br_x, grd_br_y, grd_tl_x, grd_tl_y, grd_tr_x, grd_tr_y, grd_key;
    block rand0, rand1;
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
        rand0 = mprng.get();
        rand1 = mprng.get();
        //grd_bl_x = grd_bl_x * len_sqr; 
        //grd_bl_y = grd_bl_y * len_sqr; 
        uint64_t pos_bl_x = bl_x - (grd_bl_x * len_sqr);
        uint64_t pos_bl_y = bl_y - (grd_bl_y * len_sqr);
        fulldomainEval(recv_hash, delta, rand0, pos_bl_x, true, pos_bl_y, true);
        shares = share_trivialFSS(delta, grd_bl_x, grd_bl_y, bl_x, true, bl_y, true, rand0, rand1);
        okvsVal0.push_back(shares[0]);
        okvsVal1.push_back(shares[1]);
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
    //std::cout << "sizes " << okvsKeys[0] << "  " << okvsVal0[0]  << "  " << okvsVal1[0] << std::endl;

    // OKVS (keys, values = fsskeys0) || (keys, values = fsskeys1) UNCOMMENT BELOW
    PaxosEncode(okvsKeys, okvsVal0, okvsVal1, okvs0, okvs1, 128);
    
    //auto t2 = high_resolution_clock::now();
    //auto duration = duration_cast<milliseconds>(t2-t1).count();
    //cout << "FSS_Share took in milliseconds: " << duration << endl;

    //std::cout << "key okvsVals0[0] " << okvsKeys[3] << "  " << okvsVal1[3] << std::endl;
    
}

//FSS_Eval for PSI Sender, #baseOT OKVS instances, returns the SHA evaluation    
void psi_FssEval(uint64_t x_coord, uint64_t y_coord, vector<vector<osuCrypto::block>> okvs, uint64_t delta, uint64_t hashSize){ 
    
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

vector<vector<BitVector>> blocks_to_bits(vector<vector<block>> okvs){
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

//FSS_Eval for the PSI Sender, #baseOT OKVS instances, returns the SHA evaluation, without using GF2E
void psiSender_FssEval(uint64_t x_coord, uint64_t y_coord, vector<vector<osuCrypto::BitVector>> okvs, uint64_t delta, uint64_t hashSize){ 
    
    uint64_t nbaseOT = 2; // this will be 440, can be sent a parameter
    uint64_t len_sqr = 2 * delta;
    // setting OKVS parameters
    int gamma = 60, v=20, fieldSizeBytes = 16, fieldSize = 128; // for okvs
    double c1 = 1.3; // for okvs

    uint64_t okvs_key, yx_share, grd_x, grd_y, x, y; 
    x = x_coord;
    y = y_coord;
    return_grid(x, y, grd_x, grd_y, len_sqr);
    okvs_key = grd_y;
    okvs_key = okvs_key << 32;
    okvs_key = okvs_key + grd_x;
    
    //cout << " okvs size " << " " << bitVector_okvs.size() << "  " << bitVector_okvs[0].size()<< std::endl;

    // start by learning which indices we need for a key
    // decoding from checkoutput() OKVS
    //initialize just to be able to dec() -- below 3 calls only set up the hash functions
    //auto t1 = high_resolution_clock::now();

    ObliviousDictionary * dict = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dict->init();
    //auto t2 = high_resolution_clock::now();
    //auto duration = duration_cast<milliseconds>(t2-t1).count();
    //cout << "okvs initialization time: " << duration << endl;

    vector<BitVector> okvsbits;
        //BitVector eval;
    auto indices = dict->dec(okvs_key);
    okvsbits = okvs[indices[0]];
    //for (int i = 0; i < okvs[indices[0]].size(); i++){
        //eval.assign(okvs[indices[0]][i]);
      //  okvsbits.pu
    //}
            
    for (int i = 1; i < 3; i++){
        for (int j = 0; j < okvs[indices[i]].size(); j++){
            //eval.assign(okvs[indices[i]][j]);
            okvsbits[j] = okvsbits[j] ^ okvs[indices[i]][j];
        }
    }
        

    RandomOracle sha_fss(sizeof(block));
    block fss_output;
    // computing OKVS key

    //computing the starting indices within block to access the shares
    int int_start_y = 63 - (2 * delta);
    int int_start_x = 64 + int_start_y;
    grd_y = grd_y * 2 * delta;
    grd_x = grd_x * 2 * delta; 
    for (int i = 0; i < 440; i++){ // hardcoding in the size of okvsbits = # of baseOt
        //BitIterator iter_y(valBytes.data(), int_start_y + (y - grd_y));
        //BitIterator iter_x(valBytes.data(), int_start_x + (x - grd_x));
        u8 y_share = okvsbits[i][int_start_y + (y - grd_y)];
        u8 x_share = okvsbits[i][int_start_x + (x - grd_x)];
        sha_fss.Update(y_share);
        sha_fss.Update(x_share);
    }
    sha_fss.Final(fss_output);
       //sha_fss.Reset();
    
    std::cout << "Eval for points " << fss_output << std::endl;
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
    double c1 = 2.4; // for okvs

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

/*
void batchFssEval(vector<uint64_t> x_coord, vector<uint64_t> y_coord, vector<osuCrypto::block> okvs, uint64_t delta, uint64_t hashSize){
    uint64_t len_sqr = 2* delta;
    uint64_t grd_x, grd_y, x, y; 
    int gamma = 60, v=20, fieldSizeBytes = 16, fieldSize = 128; // for okvs
    double c1 = 2.4; // for okvs
    //vector<uint64_t> okvs_keys;
    uint64_t okvs_key;

    if (x_coord.size() != y_coord.size()){
        std::cout << "inconsistent point set " << std::endl;
    }
    // decoding from checkoutput()
    //initialize just to be able to dec() -- below 3 calls only set up the hash functions
    ObliviousDictionary * dict = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dict->init();
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

    vector<byte> valBytes(fieldSizeBytes);
    int int_start_y = 63 - (2 * delta);
    int int_start_x = 64 + int_start_y;
    u8 y_share, x_share;
    GF2E dhBitsVals;

     //handling the key
    for (int i = 0; i < x_coord.size(); i++){
        x = x_coord[i];
        y = y_coord[i];
        return_grid(x, y, grd_x, grd_y, len_sqr);
        okvs_key = grd_y;
        okvs_key = okvs_key << 32;
        okvs_key = okvs_key + grd_x;
        //okvs_keys.push_back(okvs_key);
        auto indices = dict->dec(okvs_key);
        for (int j=0; j<indices.size(); j++){
            dhBitsVals += okvs_gf2e[indices[j]];
        }
        BytesFromGF2X(valBytes.data(), rep(dhBitsVals), fieldSizeBytes);
        grd_y = grd_y * 2 * delta;
        grd_x = grd_x * 2 * delta; 
        BitIterator iter_y(valBytes.data(), int_start_y + (y - grd_y));
        BitIterator iter_x(valBytes.data(), int_start_x + (x - grd_x));
        y_share = *iter_y;
        x_share = *iter_x;
    }
}

void FssShare_Enumerate(uint64_t delta, int nSquares, vector<vector<block>> &okvs0, vector<vector<block>> &okvs1){

    //std::cout << "balls are pairwise 3*delta apart " << std::endl;
    auto t1 = high_resolution_clock::now();

    //initialize variables
    uint64_t len_sqr = 2 * delta;
    vector<vector<block>> okvsVal0, okvsVal1, shares;
    vector<uint64_t> okvsKeys; // vals0 - fsskeys0, vals1 - fsskeys1
    uint64_t bl_x, bl_y, br_x, br_y, tl_x, tl_y, tr_x, tr_y, grd_bl_x, grd_bl_y, grd_br_x, grd_br_y, grd_tl_x, grd_tl_y, grd_tr_x, grd_tr_y, grd_key;
    block rand0, rand1;
    PRNG mprng(toBlock(31));

    // let's process the squares: identify the 2^dim = 4 grid cells -> compute fsskey0, fsskey1 
    // Input: square1 (delta, delta), square2 (delta + 4*delta, delta) ..... square_n(delta + 4(n -1)delta, delta) ...
    bl_x = delta;
    bl_y = delta; 

    for (int i = 0; i < nSquares; i++){
        //BL 
        bl_x = delta + (i*2*len_sqr);
        bl_y = delta; // make bl_y = delta + (i*2*len_sqr)
        return_grid(bl_x, bl_y, grd_bl_x, grd_bl_y, len_sqr);
        rand0 = mprng.get();
        rand1 = mprng.get();
        shares = share_trivialFSS2(delta, grd_bl_x, grd_bl_y, bl_x, true, bl_y, true, rand0, rand1);
        for (int k = 0; k < 128; k++){
                okvsVal0[k].push_back(shares[0][k]);
                okvsVal1[k].push_back(shares[1][k]);
        }
        
        grd_key = grd_bl_y;
        grd_key = grd_bl_y << 32;
        grd_key = grd_key + grd_bl_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;
        
        
        //TL
        tl_x = bl_x;
        tl_y = bl_y + len_sqr;
        return_grid(tl_x, tl_y, grd_tl_x, grd_tl_y, len_sqr);
        rand0 = mprng.get();
        rand1 = mprng.get();
        shares = share_trivialFSS2(delta, grd_tl_x, grd_tl_y, tl_x, true, tl_y, false, rand0, rand1);
        for (int k = 0; k < 128; k++){
                okvsVal0[k].push_back(shares[0][k]);
                okvsVal1[k].push_back(shares[1][k]);
        }
        grd_key = grd_tl_y;
        grd_key = grd_tl_y << 32;
        grd_key = grd_key + grd_tl_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;
        
        
        //BR
        br_x = bl_x + len_sqr;
        br_y = bl_y;
        return_grid(br_x, br_y, grd_br_x, grd_br_y, len_sqr);
        rand0 = mprng.get();
        rand1 = mprng.get();
        shares = share_trivialFSS2(delta, grd_br_x, grd_br_y, br_x, false, br_y, true, rand0, rand1);
        for (int k = 0; k < 128; k++){
                okvsVal0[k].push_back(shares[0][k]);
                okvsVal1[k].push_back(shares[1][k]);
        }
        grd_key = grd_br_y;
        grd_key = grd_br_y << 32;
        grd_key = grd_key + grd_br_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;

        //TR
        tr_x = bl_x + len_sqr;
        tr_y = bl_y + len_sqr;
        return_grid(tr_x, tr_y, grd_tr_x, grd_tr_y, len_sqr);
        rand0 = mprng.get();
        rand1 = mprng.get();
        shares = share_trivialFSS2(delta, grd_tr_x, grd_tr_y, tr_x, false, tr_y, false, rand0, rand1);
        for (int k = 0; k < 128; k++){
                okvsVal0[k].push_back(shares[0][k]);
                okvsVal1[k].push_back(shares[1][k]);
        }
        grd_key = grd_tr_y;
        grd_key = grd_tr_y << 32;
        grd_key = grd_key + grd_tr_x;
        okvsKeys.push_back(grd_key);
        grd_key = 0;

    }

    // OKVS (keys, values = fsskeys0) || (keys, values = fsskeys1)
    vector<block> Okvs0, Okvs1;
    for(int i = 0; i < 128; i++){
        PaxosEncode(okvsKeys, okvsVal0[i], okvsVal1[i], Okvs0, Okvs1, 128);
        okvs0.push_back(Okvs0);
        okvs1.push_back(Okvs1);
        Okvs0.clear();
        Okvs1.clear();
    }
    
    //auto t2 = high_resolution_clock::now();
    //auto duration = duration_cast<milliseconds>(t2-t1).count();
    //cout << "FSS_Share took in milliseconds: " << duration << endl;

    //testing SHA
    
} 

void OkvsEncode(std::vector<uint64_t> setKeys, const std::vector<array<osuCrypto::block, 2>> setValues0, std::vector<array<osuCrypto::block, 12>>& okvs0, uint64_t fieldSize)
{
    GF2E dhBitsVal;
    int hashSize=setKeys.size(), gamma = 60, v=20;
    double c1 = 2.4;
    //vector<uint64_t> keys;
    //keys.resize(hashSize);
    int fieldSizeBytes = fieldSize % 8 == 0 ? fieldSize/8 : fieldSize/8 + 1;
    int zeroBits = 8 - fieldSize % 8;

    //initialize 
    ObliviousDictionary * dict0 = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    //ObliviousDictionary * dict1 = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dict0->init();
    //dict1->init();

    
    vector<byte> values0;
    values0.resize(setValues0.size() * sizeof(block)*12);
    memcpy(values0.data(), setValues0.data(), setValues0.size() * sizeof(block)*12);
    
    dict0->setKeysAndVals(setKeys, values0);
    dict0->encode();

    vector<byte> x0 = dict0->getVariables();
    
    // vector of **bytes to block** , vector<block> Okvs;
    okvs0.resize(x0.size() / fieldSizeBytes);
    //okvs1.resize(x1.size() / fieldSizeBytes);
    memcpy(okvs0.data(), x0.data(), x0.size());
    //memcpy(okvs1.data(), x1.data(), x1.size());
    std::cout << "Okvs.size() " << okvs0.size() << std::endl;
    for (int i = 0; i < (okvs0.size() - 60); i++)
        std::cout << "okvs  " << i << "  " << okvs0[i][0] << "  " << okvs0[i][3] << std::endl;
    
    //below is function call to check output!!

    //dict0->checkOutput();
    //dict1->checkOutput();
    
    //auto indices = dict0->dec(setKeys[0]);
    //cout << "ENCODE " << setKeys[0] << " " << indices[0] << " " << indices[1] << std::endl;
   
} 

void OKVSDecode(vector<block> okvs, block key){
    cout << "Let's write decode " << std::endl;
    // 1. convert block key -> to unint64_t keys using xxH64
    // nothing changes in setKeysandValues()
  // converting block to bytes
    vector<byte> okvs_bytes;
    okvs_bytes.resize(Okvs.size() * sizeof(block));
    memcpy(okvs_bytes.data(), Okvs.data(), Okvs.size() * sizeof(block));

    // from bytes to GF2E!! yahooo 
    GF2EVector okvs_gf2e(Okvs.size());
    GF2X temp;
    for (int i=0; i < okvs_gf2e.size(); i++){
        GF2XFromBytes(temp, okvs_bytes.data() + i*fieldSizeBytes ,fieldSizeBytes);
        okvs_gf2e[i] = to_GF2E(temp); 
        //std::cout << okvs_gf2e[i] << std::endl;
    }
    
    // decoding a single key
    std::cout << "setKey[4] " << setKeys[4] << std::endl;

    // decoding from checkoutput()
    GF2E dhBitsVals;
    auto indices = dic->dec(keys[5]);
    for (int j=0; j<indices.size(); j++){
        dhBitsVals += okvs_gf2e[indices[j]];
    }
    //std::cout << "dhBitsVals  " << dhBitsVals << std::endl;
    vector<byte> valBytes(fieldSizeBytes);
    BytesFromGF2X(valBytes.data(), rep(dhBitsVals), fieldSizeBytes);
    block returndecode(valBytes[15],valBytes[14],valBytes[13],valBytes[12],valBytes[11],valBytes[10],valBytes[9],valBytes[8],valBytes[7],valBytes[6],valBytes[5],valBytes[4],valBytes[3],valBytes[2],valBytes[1],valBytes[0]);
    std::cout << "setVals[4] " << setValues[5] << "return " << returndecode << std::endl;

    // 2. 
}*/
   
   
/*
    for (int i = 0; i < okvs_gf2e[0].size(); i++){
        for (int j = 0; j < indices.size(); j++){
            eval += okvs_gf2e[indices[j]][i];
            dhBitsVals.push_back(eval);
        }
    }
    

    for (int j = 0; j < indices.size(); j++){
        if (j == 0){
            dhBitsVals = okvs_gf2e[indices[0]];
            std::cout << "dhbitsvals[0]" << dhBitsVals[0] << std::endl;
        }
        else {
            for (int i = 0; i < okvs_gf2e[j].size(); i++){
                dhBitsVals[i] = dhBitsVals[i] + okvs_gf2e[indices[j]][i];
                //if (j < 5)
                //    std::cout << "dhbitsvals[0]" << dhBitsVals[0] << std::endl;
            }
        }
    }
*/

