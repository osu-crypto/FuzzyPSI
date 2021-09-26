#include "okvsfss.h"

using namespace osuCrypto;

void PaxosEncode(const std::vector<block> setKeys, const std::vector<block> setValues, std::vector<block>& okvs, uint64_t fieldSize)
{
    int hashSize=setKeys.size(), gamma = 60, v=20;
    double c1 = 2.4;
    vector<uint64_t> keys;
//    vector<unsigned char> values;
    keys.resize(hashSize);
    //okvs.resize(c1*hashSize);
    int fieldSizeBytes = fieldSize % 8 == 0 ? fieldSize/8 : fieldSize/8 + 1;
    int zeroBits = 8 - fieldSize % 8;
//    values.resize(hashSize*fieldSizeBytes);
    ObliviousDictionary * dic = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dic->init();
    for (int i=0; i < setKeys.size(); i++){
        keys[i] = XXH64(&setKeys[i], 64, 0);
    }
    cout << "Done with keys" << endl;
    vector<byte> values;
    values.resize(setValues.size() * sizeof(block));
    memcpy(values.data(), setValues.data(), setValues.size() * sizeof(block));
    cout << "setValue size " << setValues.size() << std::endl;
    cout << "block size " << sizeof(block) << "  fieldSize " << fieldSizeBytes << std::endl;
    cout << "printting the keys, values size " << keys.size() << "    " << values.size() << std::endl;
    dic->setKeysAndVals(keys, values);
    dic->encode();

    vector<byte> x = dic->getVariables();
    std::cout << "getTableSize() " << dic->getTableSize() << std::endl;
    //std::cout << "x.size() " << x.size() << std::endl;
    //vector<block> okvs_ret(c1*setKeys.size());
    //memcpy(&okvs_ret, x.data(), x.size());
    //cout << "printting the x size " << okvs_ret.size() << std::endl;
    dic->checkOutput();
    cout << "here" << endl;

}
vector<uint64_t> return_grid(vector<uint64_t> coords, uint64_t delta){
    vector<uint64_t> grid_id;
    uint64_t temp;
    for (int i = 0; i < coords.size(); i++){
        temp = coords[i] / delta;
        grid_id.push_back(temp);
    }
    
    return grid_id;
}

block d_trivialFSS(uint64_t delta, uint64_t grid_x, uint64_t grid_y, uint64_t point_x, bool x, uint64_t point_y, bool y, block rand){
    std::cout << "Let's write the trivial FSS in 2 Dimensions " << std::endl;
    BitVector FSS_keyzero, FSS_keyone;
    FSS_keyzero.assign(rand);
    FSS_keyone.assign(rand);

    //(grid_x, grid_y) = grid label
    grid_x = grid_x * (2 * delta);
    grid_y = grid_y * (2 * delta);
    uint64_t int_start = 64 - (2 * delta);
 

    //Let's process the X-coord which is the second half of the block, access block.mData[0]
   
    if (x) { // fill right
    uint64_t start_point_x = int_start + (point_x - grid_x); 
    for (int i = start_point_x; i < 64; i++)
        if(FSS_keyone[64 + i])
            FSS_keyone[64 + i] = 0;
        else 
            FSS_keyone[64 + i] = 1;
    }
    else { // fill left 
     uint64_t end_point_x = int_start + (point_x - grid_x); 
        for (int i = int_start; i < end_point_x; i++)
        if(FSS_keyone[64 + i])
            FSS_keyone[64 + i] = 0;
        else 
            FSS_keyone[64 + i] = 1;
        
    }

    if (y) { // fill right
    uint64_t start_point_y = int_start + (point_y - grid_y); 
    for (int i = start_point_y; i < 64; i++)
        if(FSS_keyone[i])
            FSS_keyone[i] = 0;
        else 
            FSS_keyone[i] = 1;
    }
    else { // fill left 
     uint64_t end_point_y = int_start + (point_y - grid_y); 
        for (int i = int_start; i < end_point_y; i++)
        if(FSS_keyone[i])
            FSS_keyone[i] = 0;
        else 
            FSS_keyone[i] = 1;
        
    }
    std::cout << "FSS interval  " << (FSS_keyzero ^ FSS_keyone) << std::endl;
    //std::cout << "key 1  " << FSS_keyone << std::endl;
    // Let's process the y-coord, first/left half of block, access by block.mData[1]
    return rand;
}

/* How we write a square and identify the grid cells
    1. Grid cell label is bottom left co-coordinate
    2. We can identify a square also by bottom left coo

    Note: each square is processed independent of the other
    this allows us to sample the input squares systematically
    For now, we will assume that all the squares are in a straight line

    How to generate input squares
    1. Sample the bottom left square/center of the square
    2. Sample the side length = 2 * delta
    3. First square Bottomleftlabel - (delta, delta)
*/
void ShareFss_farapart(uint64_t delta, int nSquares){
    std::cout << "balls are pairwise 3*delta apart " << std::endl;
    PRNG mprng(toBlock(31));
    vector<uint64_t> coord;
    uint64_t a = 12345;
    uint64_t b = 5678;
    coord.push_back(a);
    coord.push_back(b);
    vector<uint64_t> grid = return_grid(coord, delta);
    
    // trivial FSS
    block test = toBlock(a, b);
     std::cout << "a " << a << " b " << b << std::endl;
    std::cout << "test " << test << std::endl;
    std::cout << "a " << test.mData[1] << "  " << test.mData[0] << std::endl;

    // testing bitvector of size 10 where first half is 1 and second half is 0
    BitVector test1;

    //test1.randomize(mprng);
    //std::cout << "random bitvector " << test1  << std::endl;

    //reserve(10);
    /*
    for (int i = 0; i < 5; i++)
        test1.pushBack(1);
    for (int i = 0; i < 5; i++)
        test1.pushBack(0);*/
    test1.assign(test);
    std:: cout << "printing bitvector size  " << test1 << std::endl;
    test1[0] = 1;
    test1[1] = 1;
    std:: cout << "printing bitvector size  " << test1 << std::endl;
    block test2 = d_trivialFSS(10, 0, 0, 15, true, 5, true, test);

   
}

/* useful notes

*/