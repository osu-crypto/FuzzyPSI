#include "fuzzyPSI.h"

using namespace osuCrypto;
//using namespace cryptoTools;
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
   
}

void fuzzyPSI(u64 keysize)
{
  
    fss_interval(1, 10);
    // Setup networking. See cryptoTools\frontend_cryptoTools\Tutorials\Network.cpp
    IOService ios;
    Channel senderChl = Session(ios, "localhost:1212", SessionMode::Server).addChannel();
    Channel recverChl = Session(ios, "localhost:1212", SessionMode::Client).addChannel();

    // key length = # of blocks needed to represent the FSS keys
    u64 fsskeySize =  keysize; 

    // The number of OTs - we need \kappa OTs
    int baseCount = 3; // for  BFSS - baseCount = hamming distance 
    std::vector<osuCrypto::block> baseRecv(baseCount*fsskeySize);
    BitVector choices(baseCount*fsskeySize);

    // The code to be run by the OT receiver.
    auto recverThread = std::thread([&]() {

        PRNG r_prng(toBlock(14));
        MasnyRindal recver;

        // Choose which messages should be received.
        choices.randomize(r_prng);
        //making fsskeys blocks of choices same value
        for (u64 i = 0; i < baseCount; i++){
            for (u64 j = 1; j < fsskeySize; j++)
                choices[(i*fsskeySize) + j] = choices[(i*fsskeySize) + 0];
        }

        // Receive the messages
        recver.receiveChosen(choices, baseRecv, r_prng, recverChl);

        //AES.cpp has hash block, twoblocks, many blocks, woohoo 
        std::cout << "plain block " << baseRecv[0] << std::endl; 
        block ciphertxt, hash, theirHashseed, r_psiHash;
        std::vector<block> hashed_inputs, hashed_item(fsskeySize);
        block r_HashSeed = r_prng.get<block>();
        details::AESTypes::NI;
       
        std::cout << "I am in the receiver side sending my hash seed " << r_HashSeed << std::endl; 
        recverChl.asyncSend(r_HashSeed);
        recverChl.recv(theirHashseed);
        std::cout << "I am in the receiver side receiving my hash seed " << theirHashseed << std::endl; 
        r_psiHash = r_HashSeed ^ theirHashseed;
        AES hasher(r_psiHash);
        hasher.ecbEncBlock(toBlock(1), ciphertxt);
        std::cout << "receiver checking hashing " << ciphertxt << std::endl;
        //hasher.AES(rand);
        for (int i = 0; i < baseCount; i++){
            for (int j = 0; j < fsskeySize; j++){
                hasher.ecbEncBlock(baseRecv[(i*fsskeySize) + j], ciphertxt);
                hashed_item.push_back(ciphertxt);
                //hasher.ecbEncBlocks(baseRecv.data(), baseRecv.size(), ciphertxts.data());
            }
            for (int k = 0; k < hashed_item.size(); k++)
                hash = hash ^ hashed_item[k];
            hashed_inputs.push_back(hash);
            hashed_item.clear();
        }
        recverChl.asyncSend(hashed_inputs);      
        for (int i = 0; i < hashed_inputs.size(); i++)
          std::cout << " hash values " << hashed_inputs[i] << std::endl; 
        });

    PRNG s_prng(toBlock(12));
    MasnyRindal sender;

    // Choose which messages should be sent.
    std::vector<std::array<osuCrypto::block, 2>> baseSend(baseCount*fsskeySize);
    s_prng.get(baseSend.data(), baseSend.size());

    // Send the messages.
    sender.sendChosen(baseSend, s_prng, senderChl);
    block theirHash, s_HashSeed, s_psiHash, s_ciphertxt;
    std::vector<block> recvd_outputs;
    s_HashSeed = s_prng.get<block>();
    senderChl.recv(theirHash);
    std::cout << "I am in the sender side and receiving a hash seed " << theirHash << std::endl; 
    senderChl.asyncSend(s_HashSeed);
    std::cout << "I am in the sender side and sending a hash seed " << s_HashSeed << std::endl; 
    s_psiHash = s_HashSeed ^ theirHash; 
    AES s_hasher(s_psiHash);
    //s_hasher.ecbEncBlock(toBlock(1), s_ciphertxt);
    //std::cout << "sender checking hashing " << s_ciphertxt << std::endl;
    senderChl.recv(recvd_outputs);
    for (int i = 0; i < recvd_outputs.size(); i++)
          std::cout << " hash values " << recvd_outputs[i] << std::endl; 
    recverThread.join();


}


/*    std::cout << "testing chosen base OT " << std::endl; 
    for (u64 i = 0; i < (baseCount*fsskeySize); ++i)
        {
            
            if (neq(baseRecv[i], baseSend[i][choices[i]]))
            {
                std::cout << "failed " << i << " exp = m[" << int(choices[i]) << "], act = " << baseRecv[i] << " true = " << baseSend[i][0] << ", " << baseSend[i][1] << std::endl;
            }
            
        }*/