#include "fuzzyPSI.h"
//#include "paxos.h"

using namespace osuCrypto;
//using namespace cryptoTools;

/* experiment, issue is reusing te pointers, fix later
array<vector<block>, 2> s_psi_interval(Fss* f_sender, int n, uint64_t a, uint64_t b){

    array<vector<block>, 2> return_blks;
    vector<block> int_keys0, int_keys1;
    block * key_in_blks; 
    //pair of key
    ServerKeyLt l_lt_k0, l_lt_k1, r_lt_k0, r_lt_k1;
    for (int i = 0; i < n; i++){
        generateTreeLt(f_sender, &r_lt_k0, &r_lt_k1, b, 1);  
        key_in_blks = reinterpret_cast<block *>(&r_lt_k1); 
        for (int i = 0; i < 11; i++){
            int_keys1.push_back(key_in_blks[i]);
            std::cout << "one keys  " << key_in_blks[i] << std::endl;
        }
        key_in_blks = reinterpret_cast<block *>(&r_lt_k0); 
        for (int i = 0; i < 11; i++){
            int_keys0.push_back(key_in_blks[i]);
            std::cout << "zero keys  " << key_in_blks[i] << std::endl;
        }
    }
    return_blks[0] = int_keys0;
    return_blks[1] = int_keys1;
    return return_blks; 

    
}*/

// helper function call for testing transpose, printing matrix
void printMtx(std::array<block, 1>& data)
	{
		for (auto& d : data)
		{
			std::cout << d << std::endl;
		}
}

void basic_transpose(){

    PRNG r_prng(toBlock(11));
    //-------------------basic test---------------------------------

    std::array<block, 128> data, data2, data3;
    r_prng.get(data.data(), data.size());

    MatrixView<block> dataView(data.begin(), data.end(), 1);
	MatrixView<block> data2View(data2.begin(), data2.end(), 1);
    MatrixView<block> data3View(data3.begin(), data3.end(), 1);
    
    transpose(dataView, data2View);
    transpose(data2View, data3View);
    std::cout << "data before matrix  " << data[5] << std::endl;
    std::cout << "data after tranpose " << data3[5] << std::endl;

    //------------------next step----------------------------------

    std::array<std::array<block, 4>, 128> t, t3;
    r_prng.get((u8*)t.data(), sizeof(block) * 4 * 128);
    std::array<block, 128 * 4> t2;

	MatrixView<block> tView((block*)t.data(), 128, 4);
	MatrixView<block> t2View((block*)t2.data(), 128 * 4, 1);
    MatrixView<block> t3View((block*)t3.data(), 128, 4);
	transpose(tView, t2View);
    transpose(t2View, t3View);

    // why no index out of bounds question here?
    std::cout << "transposed matrix " << std::endl;
    std::cout << "data before " << t[0][0] << std::endl;
    std::cout << "data after " << t3[0][0] << std::endl;

    //--------transpose needed for PSI-----------------------

    // does not work for 440 and breaking it into blocks, 
    std::array<block, 440> f, f2, f3;
    r_prng.get(f.data(), f.size());

    MatrixView<u8> fView((u8*)f.data(), 440, 16);
	MatrixView<u8> f2View((u8*)f2.data(), 128, 55);
    MatrixView<u8> f3View((u8*)f3.data(), 440, 16);

    transpose(fView, f2View);

    transpose(f2View, f3View);
   
    std::cout << "f2 size " << std::endl;
    std::cout << "data before " << f[400] << std::endl;
    std::cout << "data after " << f3[400] << std::endl;

    //--------------let's try with Matrix---------------------------
    Matrix<u8> in(440, 16);
	r_prng.get((u8*)in.data(), sizeof(u8) *in.bounds()[0] * in.stride());

	Matrix<u8> out2(128, 55);
	transpose(in, out2); 
}

void Transpose_View_Test() 
{
        PRNG prng1(toBlock(12));
        std::array<block, 2> input; 
        prng1.get(input.data(), input.size());
        std::cout << "input " << input[0] << "   " << input[1] << std::endl;
        

/*

		{

			PRNG prng(ZeroBlock);

			std::array<block, 128> data, data2;
			prng.get(data.data(), data.size());
			//std::array<block, 128> data2;

			MatrixView<block> dataView(data.begin(), data.end(), 1);
			MatrixView<block> data2View(data2.begin(), data2.end(), 1);

			transpose(dataView, data2View);
            //printMtx();
			transpose128(data);




			for (u64 i = 0; i < 128; ++i)
			{
				if (neq(data[i], data2[i]))
				{
					std::cout << i << "\n";
					printMtx(data);
					std::cout << "\n";
					printMtx(data2);

					//throw UnitTestFail();
				}
			}
		}


		{
			PRNG prng(ZeroBlock);

			std::array<std::array<block, 8>, 128> data;

			prng.get((u8*)data.data(), sizeof(block) * 8 * 128);


			std::array<std::array<block, 8>, 128> data2;

			MatrixView<block> dataView((block*)data.data(), 128, 8);
			MatrixView<block> data2View((block*)data2.data(), 128 * 8, 1);
			transpose(dataView, data2View);


			for (u64 i = 0; i < 8; ++i)
			{
				std::array<block, 128> data128;

				for (u64 j = 0; j < 128; ++j)
				{
					data128[j] = data[j][i];
				}

				transpose128(data128);


				for (u64 j = 0; j < 128; ++j)
				{
					if (neq(data128[j], data2View[i * 128 + j][0]))
						//throw UnitTestFail();
                        std::cout << "issue" << std::endl;
				}
			}

		}
*/

		{
			PRNG prng(ZeroBlock);

			std::array<std::array<std::array<block, 8>, 128>, 2> data;

			Matrix<block> dataView(208, 8);
			prng.get((u8*)dataView.data(), sizeof(block) *dataView.bounds()[0] * dataView.stride());

			Matrix<block> data2View(1024, 2);
			memset(data2View.data(), 0, data2View.bounds()[0] * data2View.stride() * sizeof(block));
			transpose(dataView, data2View);

			for (u64 b = 0; b < 2; ++b)
			{

				for (u64 i = 0; i < 8; ++i)
				{
					std::array<block, 128> data128;

					for (u64 j = 0; j < 128; ++j)
					{
						if (dataView.bounds()[0] > 128 * b + j)
							data128[j] = dataView[128 * b + j][i];
						else
							data128[j] = ZeroBlock;
					}

					transpose128(data128);

					for (u64 j = 0; j < 128; ++j)
					{
						if (neq(data128[j], data2View[i * 128 + j][b]))
						{
							std::cout << "failed " << i << "  " << j << "  " << b << std::endl;
							std::cout << "exp: " << data128[j] << "\nact: " << data2View[i * 128 + j][b] << std::endl;
							//throw UnitTestFail();
						}
					}
				}
			}
		}
/*
		{
			PRNG prng(ZeroBlock);

			Matrix<u8> in(16, 8);
			prng.get((u8*)in.data(), sizeof(u8) *in.bounds()[0] * in.stride());

			Matrix<u8> out(63, 2);
			transpose(in, out);


			Matrix<u8> out2(64, 2);
			transpose(in, out2);

			for (u64 i = 0; i < out.bounds()[0]; ++i)
			{
				if (memcmp(out[i].data(), out2[i].data(), out[i].size()))
				{
					std::cout << "bad " << i << std::endl;
					//throw UnitTestFail();
				}
			}
		}

		{
			PRNG prng(ZeroBlock);

			//std::array<std::array<std::array<block, 8>, 128>, 2> data;

			Matrix<u8> in(25, 9);
			Matrix<u8> in2(32, 9);

			prng.get((u8*)in.data(), sizeof(u8) *in.bounds()[0] * in.stride());
			memset(in2.data(), 0, in2.bounds()[0] * in2.stride());

			for (u64 i = 0; i < in.bounds()[0]; ++i)
			{
				for (u64 j = 0; j < in.stride(); ++j)
				{
					in2[i][j] = in[i][j];
				}
			}

			Matrix<u8> out(72, 4);
			Matrix<u8> out2(72, 4);

			transpose(in, out);
			transpose(in2, out2);

			for (u64 i = 0; i < out.bounds()[0]; ++i)
			{
				for (u64 j = 0; j < out.stride(); ++j)
				{
					if (out[i][j] != out2[i][j])
					{
						std::cout << (u32)out[i][j] << " != " << (u32)out2[i][j] << std::endl;
						//throw UnitTestFail();
                        std::cout << "issue" << std::endl;
					}
				}
			}
		}

        */
}

void fss_interval3(uint64_t a, uint64_t b){
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


    block* struct_in_blk = reinterpret_cast<block *>(&r_lt_k1); 
    vector<block> struct_in_blocks;
    for (int i = 0; i < 11; i++){
        std::cout << "blocks in struct " << struct_in_blk[i] << std::endl; 
        struct_in_blocks.push_back(struct_in_blk[i]);
    }
    generateTreeLt(&l_fClient, &r_lt_k0, &r_lt_k1, b, 1);
    block* struct_in_blk2 = reinterpret_cast<block *>(&r_lt_k1); 
    vector<block> struct_in_blocks2;
    for (int i = 0; i < 11; i++){
        std::cout << "blocks in struct2 " << struct_in_blk2[i] << std::endl; 
        struct_in_blocks2.push_back(struct_in_blk2[i]);
    }

    ServerKeyLt * original_struct = reinterpret_cast<ServerKeyLt*>(struct_in_blocks.data());
    l_lt_ans0 = evaluateLt(&l_fClient, &r_lt_k1, (a+1));
    r_lt_ans0 = evaluateLt(&l_fClient, &r_lt_k0, (a+1));
    l_lt_ans1 = evaluateLt(&l_fClient, original_struct, (a+1));

    std::cout << l_lt_ans0 << "   " << r_lt_ans0 << "   " << l_lt_ans1 << std::endl; 
}

//Below is full simulated PSI protocol!
// 1. using SHA instead of AES for hashing
// 2. which base OT do we use? we need semi-honest -- cost to performance is so low that using malicious baseOT is fine
// 3. use unordered_map / hash_map to lookup items in the end
        // map<uint64_t, block> (key is represented as unint64_t, final evaluation is a hash so using sha it will be a block)
// 4. Assuming 512 base OTs - 4 blocks (might be easy for the transpose notation)

void fuzzyPSI(u64 keysize, u64 y_input_size, u64 x_volume)
{

    // Setup networking. Setting up channels for PSI, sender and recver according to OT not PSI 
    IOService ios;
    Channel senderChl = Session(ios, "localhost:1212", SessionMode::Server).addChannel();
    Channel recverChl = Session(ios, "localhost:1212", SessionMode::Client).addChannel();

    // key length = # of blocks needed to represent the FSS keys
    // determine this according to the FSS function we calling PSI on
    // for a single interval the **keysize is 22**
    u64 fsskeySize =  110; //  12000 works
    u64 x_input_volume = x_volume;

    // The number of OTs - we need \kappa OTs
    int baseCount = 440; // for  BFSS - baseCount = hamming distance
    
    // **PSI receiver** has 128 choices and OT received blocks 
    // **PSI receiver** has an unstructured set r_inputs as inputs
    // we just call the Hash * r_inputs 
    std::vector<osuCrypto::block> baseRecv(baseCount);
    BitVector choices(baseCount);
    for (int i = 0; i < baseCount; i++)
        choices[i] = 0;
    u64 r_input_size = y_input_size;

    // **PSI sender** operates out of the receiver thread
    auto recverThread = std::thread([&]() {

        //initializing
        block ciphertxt, hash, theirHashseed, r_psiHash;
        
        // sampling the l - bit vector or choices as the OT receiver 
        PRNG r_prng(toBlock(14));
        //NOTE: we set choices to all 0 to make FSSrecv and FSSkeys0 match and check rest of
        //choices.randomize(r_prng); 
        MasnyRindal recver;
        
        // Receive the messages
        recver.receiveChosen(choices, baseRecv, r_prng, recverChl);
        
        // basically receiver gets ciphertexts of form [# base OT] * [enc(k1, m1), enc(k2, m2)] here
        std::vector<block> recv_ciphertxt0, recv_ciphertxt1, fsskeysRecv;
        block recv_ciphertxt;
        recverChl.recv(recv_ciphertxt0);
        recverChl.recv(recv_ciphertxt1);

        //now we use the baseRecv as decryption key and choices 
        details::AESTypes::Portable;
        AESDec aeskey;
            for (int i = 0; i < baseCount; i++){
                aeskey.setKey(baseRecv[i]);
                if (choices[i] == 0){
                    for (int j = 0; j < fsskeySize; j++){
                        recv_ciphertxt = aeskey.ecbDecBlock(recv_ciphertxt0[(i*fsskeySize) + j]);
                        fsskeysRecv.push_back(recv_ciphertxt);
                    }  
                }
                else{ 
                    for (int j = 0; j < fsskeySize; j++){
                        recv_ciphertxt = aeskey.ecbDecBlock(recv_ciphertxt1[(i*fsskeySize) + j]);
                        fsskeysRecv.push_back(recv_ciphertxt);
                    }
                }

            }
       

        
        // initializing variables, getting common hash seed
        std::vector<block> hashed_inputs, hashed_item(fsskeySize*baseCount+1);
        block r_HashSeed = r_prng.get<block>();
        details::AESTypes::Portable;
        recverChl.asyncSend(r_HashSeed);
        recverChl.recv(theirHashseed);
        r_psiHash = r_HashSeed ^ theirHashseed;
        //AES hasher(r_psiHash);
        RandomOracle r_sha(sizeof(block));
       
        // hash each of the received keys with the appropriate inputs
        // ^^insert *FSS eval for each of the r_inputs happens here*
        // we Hash size (fsskeySize * baseCount )  ---> size (baseCount) for a single inputs!

        // suppose we have r_inputs_size, append input and hash (for now)
        
        for (u64 r = 0; r < r_input_size; r++){
            /*
            fsskeysRecv.push_back(toBlock(r));
            hasher.ecbEncBlocks(fsskeysRecv.data(), fsskeySize*baseCount+1, hashed_item.data());
            fsskeysRecv.pop_back();
            hash = hashed_item[0];
            for (int k = 1; k < hashed_item.size(); k++){ // here we xor all blocks associated with "one fsskey", there are baseCount such keys
                hash = hash ^ hashed_item[k];
            }
            */
            r_sha.Update(r_psiHash);
            r_sha.Update(toBlock(r));
            for (int k = 0; k < fsskeysRecv.size(); k++)
                r_sha.Update(fsskeysRecv[k]);
            //r_sha.Update(fsskeysRecv.data(), fsskeysRecv.size() * sizeof(block));
            r_sha.Final(hash);
            r_sha.Reset();
            hashed_inputs.push_back(hash);
        }
        //sender sending hash of FSS eval of his inputs
        recverChl.asyncSend(hashed_inputs); 
        //for (int i = 0; i < hashed_inputs.size(); i++)
        //    std::cout << "recver hash  " << hashed_inputs[i] << std::endl;  

       
        });

    // *PSI receiver*
    block theirHash, s_HashSeed, s_psiHash, s_ciphertxt;
    std::vector<block> recvd_outputs;
    
    //initializing *OT sender*
    PRNG s_prng(toBlock(12));
    MasnyRindal sender;
    std::vector<std::array<osuCrypto::block, 2>> baseSend(baseCount);
    s_prng.get(baseSend.data(), baseSend.size());
    

    // Send the messages, that is, short keys (k0, k1) that will encrypt the actual FSS keys (FSSkey0, FSSkey1)
    sender.sendChosen(baseSend, s_prng, senderChl);
   
    //sampling uniformly random FSS keys ---> later this should come from some function call
    //Fss f_sender;
    vector<block> fsskeys0(baseCount*fsskeySize), fsskeys1(baseCount*fsskeySize);
    //= s_psi_interval(&f_sender, 1, 5, 10);
    std::vector<block> keys_ciphertxt0, keys_ciphertxt1;
    s_prng.get(fsskeys0.data(), fsskeys0.size());
    s_prng.get(fsskeys1.data(), fsskeys1.size());
   
    AES aeskey0, aeskey1;
    block ciphertxt0, ciphertxt1;  
    //encrypt the fsskeys using baseSend(encryption keys)
    for (int i = 0; i < baseCount; i++){
        aeskey0.setKey(baseSend[i][0]);
        aeskey1.setKey(baseSend[i][1]);
        for (int j = 0; j < fsskeySize; j++){
            aeskey0.ecbEncBlock(fsskeys0[(i*fsskeySize) + j], ciphertxt0);
            aeskey1.ecbEncBlock(fsskeys1[(i*fsskeySize) + j], ciphertxt1); 
            keys_ciphertxt0.push_back(ciphertxt0);
            keys_ciphertxt1.push_back(ciphertxt1);
        }
    }
    senderChl.asyncSend(keys_ciphertxt0);
    senderChl.asyncSend(keys_ciphertxt1);
    

    // establishing common hash seed for the last step 
    s_HashSeed = s_prng.get<block>();
    senderChl.recv(theirHash);
    senderChl.asyncSend(s_HashSeed);
    s_psiHash = s_HashSeed ^ theirHash; 
    //AES s_hasher(s_psiHash);
    RandomOracle s_sha(sizeof(block));
    

    // now we recv the hash of the each of the PSI sender's i nputs 
    vector<block> recvHashinputs(r_input_size), recver_hashes, s_hashed_item(fsskeySize*baseCount+1);
    block s_hash; 
    senderChl.recv(recvHashinputs);
    recverThread.join();
    // **PSI receiver learns the actual output here**
    // # idea 
    // use fsskeys0 and compute the hashes, |hashes| = x_volume (#structured input) and sort the vector hashinput
    std::unordered_map<block, uint64_t> recv_hash; 

    for (uint64_t r = 0; r < x_input_volume; r++){
            /*
            fsskeys0.push_back(toBlock(r));
            s_hasher.ecbEncBlocks(fsskeys0.data(), fsskeySize*baseCount+1, s_hashed_item.data());
            fsskeys0.pop_back();
            s_hash = s_hashed_item[0];
            for (int k = 1; k < s_hashed_item.size(); k++) // here we xor all blocks associated with "one input"
                s_hash = s_hash ^ s_hashed_item[k];
            */
            s_sha.Update(s_psiHash);
            s_sha.Update(toBlock(r));
            //s_sha.Update(fsskeys0.data(), fsskeys0.size() * sizeof(block));
            for (int k = 0; k < fsskeys0.size(); k++) 
                s_sha.Update(fsskeys0[k]);
            s_sha.Final(s_hash);
            s_sha.Reset();
            //recver_hashes.push_back(s_hash);
            recv_hash.insert({s_hash, r});
        }

    vector<uint64_t> psi_outputs;
    for (uint64_t i = 0; i < recvHashinputs.size(); i++){
        if(recv_hash.find(recvHashinputs[i]) != recv_hash.end()){
            psi_outputs.push_back(recv_hash[recvHashinputs[i]]);
        }
        else {
            cout << "problematic point is " << recv_hash[recvHashinputs[i]] << std::endl; 
        }

    } 
    
    //for (int i = 0; i < recver_hashes.size(); i++)
    //    std::cout << "sender keys  " << recver_hashes[i] << std::endl;
    

   
    if(psi_outputs.size() == recvHashinputs.size())
        cout << " PSI WORKS " << std::endl;
    else
        cout << " PSI FAILS " << std::endl;
    cout << " output size " << psi_outputs.size() << std::endl;
    
   
}

/*
    RandomOracle sha(sizeof(block)), sha1(sizeof(block));
    u8 a = 1;
    u8 b = 0;
    sha.Update(a);
    sha.Update(b);
    block r0, r1, r2;
    sha.Final(r0);
    std::cout << "sha " << r0 << std::endl;
    sha.Reset();
    sha.Update(b);
    sha.Final(r1);
    std::cout << "sha " << r1 << std::endl;
    sha1.Update(b);
    sha1.Final(r2);
    std::cout << "sha " << r2 << std::endl;


    s_sha.Update(toBlock(r));
            //for (int k = 1; k < fsskeys0.size(); k++)
                s_sha.Update((u8*)(fsskeys0.data()), fsskeys0.size() * sizeof(block));
            s_sha.Final(s_buff);
            memcpy(&s_hash, s_buff, std::min<u64>(RandomOracle::HashSize, sizeof(block)));
            //s_hash = *(block*)s_buff;
*/
 /*
    std::sort(recver_hashes.begin(), recver_hashes.end());
    
    // search rechHashinputs against recver_hashes
    
    int output_size = 0;
    for (u64 i = 0; i < recvHashinputs.size(); i++){
        if(binary_search(recver_hashes.begin(), recver_hashes.end(), recvHashinputs[i])){
            psi_outputs.push_back(recvHashinputs[i]);
            output_size = output_size + 1;
        }
        else 
            cout << "problematic point is " << i << std::endl;
    }
    
*/

//incomplete things!!
void fuzzyPSI_int(u64 keysize)
{
    //space to test out dependencies
    //fss_interval3(1, 10);

    // Setup networking. Setting up channels for PSI, sender and recver according to OT not PSI 
    IOService ios;
    Channel senderChl = Session(ios, "localhost:1212", SessionMode::Server).addChannel();
    Channel recverChl = Session(ios, "localhost:1212", SessionMode::Client).addChannel();

    // key length = # of blocks needed to represent the FSS keys
    // determine this according to the FSS function we calling PSI on
    // for a single interval the **keysize is 22**
    u64 fsskeySize =  11; 

    // The number of OTs - we need \kappa OTs
    int baseCount = 2; // for  BFSS - baseCount = hamming distance
    
    // **PSI receiver** has 128 choices and OT received blocks 
    std::vector<osuCrypto::block> baseRecv(baseCount);
    BitVector choices(baseCount);

    // **PSI sender** operates out of the receiver thread
    auto recverThread = std::thread([&]() {

        //initializing
        block ciphertxt, hash, theirHashseed, r_psiHash;
        std::vector<block> hashed_inputs, hashed_item(fsskeySize);

        // sampling the l - bit vector or choices as the OT receiver 
        PRNG r_prng(toBlock(14));
        choices.randomize(r_prng);
        MasnyRindal recver;
        
        // Receive the messages, establish **common hash key**
        recver.receiveChosen(choices, baseRecv, r_prng, recverChl);
        block r_HashSeed = r_prng.get<block>();
        details::AESTypes::Portable;
        recverChl.asyncSend(r_HashSeed);
        recverChl.recv(theirHashseed);
        r_psiHash = r_HashSeed ^ theirHashseed;
        AES hasher(r_psiHash);

        // basically receiver ciphertexts of form # base OT * [enc(k1, m1), enc(k2, m2)] here
        std::vector<block> recv_ciphertxt0, recv_ciphertxt1, fsskeysRecv;
        block recv_ciphertxt;
        recverChl.recv(recv_ciphertxt0);
        recverChl.recv(recv_ciphertxt1);

        //now we use the baseRecv as decryption key and choices 
        details::AESTypes::Portable;
        AESDec aeskey;
            for (int i = 0; i < baseCount; i++){
                aeskey.setKey(baseRecv[i]);
                for (int j = 0; j < fsskeySize; j++){
                        if (choices[i] == 0)
                            recv_ciphertxt = aeskey.ecbDecBlock(recv_ciphertxt0[(i*fsskeySize) + j]);  
                        else 
                            recv_ciphertxt = aeskey.ecbDecBlock(recv_ciphertxt1[(i*fsskeySize) + j]);
                        fsskeysRecv.push_back(recv_ciphertxt);
                }
            }

        // hash each of the received keys with the appropriate inputs
       
        });

    // *PSI receiver*
    block theirHash, s_HashSeed, s_psiHash, s_ciphertxt;
    std::vector<block> recvd_outputs;
    
    //initializing *OT sender*
    PRNG s_prng(toBlock(12));
    MasnyRindal sender;
    std::vector<std::array<osuCrypto::block, 2>> baseSend(baseCount);
    s_prng.get(baseSend.data(), baseSend.size());

    // Send the messages, that is, short keys (k0, k1) that will encrypt the actual FSS keys (FSSkey0, FSSkey1)
    sender.sendChosen(baseSend, s_prng, senderChl);
    s_HashSeed = s_prng.get<block>();
    senderChl.recv(theirHash);
    senderChl.asyncSend(s_HashSeed);
    s_psiHash = s_HashSeed ^ theirHash; 
    AES s_hasher(s_psiHash);
    
    //sampling uniformly random FSS keys ---> later this should come from some function call
    //Fss f_sender;
    vector<block> fsskeys0(baseCount*fsskeySize), fsskeys1(baseCount*fsskeySize);
    //= s_psi_interval(&f_sender, 1, 5, 10);
    std::vector<block> keys_ciphertxt0, keys_ciphertxt1;
    s_prng.get(fsskeys0.data(), fsskeys0.size());
    s_prng.get(fsskeys1.data(), fsskeys1.size());
   
    AES aeskey0, aeskey1;
    block ciphertxt0, ciphertxt1;  
    //encrypt the fsskeys using baseSend(encryption keys)
    for (int i = 0; i < baseCount; i++){
        aeskey0.setKey(baseSend[i][0]);
        aeskey1.setKey(baseSend[i][1]);
        for (int j = 0; j < fsskeySize; j++){
            aeskey0.ecbEncBlock(fsskeys0[(i*fsskeySize) + j], ciphertxt0);
            aeskey1.ecbEncBlock(fsskeys1[(i*fsskeySize) + j], ciphertxt1); 
            keys_ciphertxt0.push_back(ciphertxt0);
            keys_ciphertxt1.push_back(ciphertxt1);
        }
    }
    senderChl.asyncSend(keys_ciphertxt0);
    senderChl.asyncSend(keys_ciphertxt1);
    
    
    recverThread.join();


}


