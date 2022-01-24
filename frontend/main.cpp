#include <iostream>

//using namespace std;
#include "tests_cryptoTools/UnitTests.h"
#include "libOTe_Tests/UnitTests.h"
#include <cryptoTools/Common/Defines.h>

using namespace osuCrypto;

#include <string.h>
#include <stdio.h>

#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Network/Session.h>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Common/BitVector.h>
#include <numeric>
#include <cryptoTools/Common/Timer.h>

#include <iomanip>
#include "util.h"

#include "ExampleBase.h"
#include "ExampleTwoChooseOne.h"
#include "ExampleNChooseOne.h"
#include "ExampleSilent.h"
#include "ExampleVole.h"
#include "libOTe/Tools/LDPC/LdpcImpulseDist.h"

#include "fuzzyPSI.h"
#include "okvsfss.h"
#include "frontend/libPaXoS/ObliviousDictionary.h"

static const std::vector<std::string>
f{"f", "fuzzyPSI"},
unitTestTag{ "u", "unitTest" },
kos{ "k", "kos" },
dkos{ "d", "dkos" },
kkrt{ "kk", "kkrt" },
iknp{ "i", "iknp" },
diknp{ "diknp" },
oos{ "o", "oos" },
moellerpopf{ "p", "moellerpopf" },
ristrettopopf{ "r", "ristrettopopf" },
mr{ "mr" },
mrb{ "mrb" },
Silent{ "s", "Silent" },
vole{ "vole" },
akn{ "a", "akn" },
np{ "np" },
simple{ "simplest" },
simpleasm{ "simplest-asm" };

#ifdef ENABLE_IKNP
void minimal()
{
    // Setup networking. See cryptoTools\frontend_cryptoTools\Tutorials\Network.cpp
    IOService ios;
    Channel senderChl = Session(ios, "localhost:1212", SessionMode::Server).addChannel();
    Channel recverChl = Session(ios, "localhost:1212", SessionMode::Client).addChannel();

    // The number of OTs.
    int n = 100;

    // The code to be run by the OT receiver.
    auto recverThread = std::thread([&]() {
        PRNG prng(sysRandomSeed());
        IknpOtExtReceiver recver;

        // Choose which messages should be received.
        BitVector choices(n);
        choices[0] = 1;
        //...

        // Receive the messages
        std::vector<block> messages(n);
        recver.receiveChosen(choices, messages, prng, recverChl);

        // messages[i] = sendMessages[i][choices[i]];
        });

    PRNG prng(sysRandomSeed());
    IknpOtExtSender sender;

    // Choose which messages should be sent.
    std::vector<std::array<block, 2>> sendMessages(n);
    sendMessages[0] = { toBlock(54), toBlock(33) };
    //...

    // Send the messages.
    sender.sendChosen(sendMessages, prng, senderChl);
    recverThread.join();
}
#endif



#include "cryptoTools/Crypto/RandomOracle.h"
int main(int argc, char** argv)
{

    CLP cmd;
    cmd.parse(argc, argv);
    bool flagSet = false;

    //if (cmd.isSet("triang"))
    //{
    //    ldpc(cmd);
    //    return 0;
    //}


    //if (cmd.isSet("encode"))
    //{
    //    encodeBench(cmd);
    //    return 0;
    //}

    if (cmd.isSet("ldpc"))
    {
        LdpcDecode_impulse(cmd);
        return 0;
    }

    // add calls here to test the fuzzyPSI code 
    if (cmd.isSet("f"))
    {
        std::cout << "Let's do fuzzy PSI " << std::endl; 
        //fuzzyPSI(11, 2000, 1000000);
        //fuzzyPSI(11, 3, 10);

        //here we are testing a basic okvs call 
        /*cout << "Let's do OKVS " << std::endl;
        PRNG pprng(toBlock(123));
        //vector<block> okvs; 
        vector<uint64_t> okvskeys(25);
        vector<array<block, 12>> okvsvals(25);
        vector<array<block, 12>> okvs;
        pprng.get(okvskeys.data(), okvskeys.size());
        pprng.get(okvsvals.data(), okvsvals.size());
        std::cout << okvsvals[0][0] << " " << okvsvals[0][1] << " " << okvsvals[0][11] << std::endl;
        uint64_t fieldsize = 128*12;
        auto t1 = high_resolution_clock::now();
        OkvsEncode(okvskeys, okvsvals, okvs, fieldsize);  
        auto t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "OKVS encode + decode  " << duration << endl;
        */

        //here we are testing a bas ic share FSS for far apart 
        
        cout << "OKVS FSS + share + eval by the PSI receiver " << std::endl;
        fss_psi();
        /*uint64_t delta = 10;
        uint64_t nsquares = 1;
        uint64_t nkeys = nsquares * 4;
        array<vector<block>, 440> okvs_fsskey0, okvs_fsskey1;
        auto t1 = high_resolution_clock::now();
        psi_FssShareEval(delta, nsquares, okvs_fsskey0, okvs_fsskey1);
        auto t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "FSS_Eval simulation in milliseconds: " << duration << endl;
        */
        
         cout << "TESTING PSI STUFF " << std::endl;
         //basic_transpose();

         // 1. figure out the parameters of OKVS, OT messages size, check run time of 440 instances of OT 
         // 2. try a basic transpose for a matrix of arbitrary size - DONE
         // 3. Share+Eval - go through every point in the square - DONE
         // 4. Share+Eval - go through the indices within block - DONE check if the hash values match for both shares *inside* square - ALMOST DONE, check again, TR 117
         // 5. Share+Eval - write code outside full domain eval and compute the hash outside the square for some arbit points
         // 6. PSI_sender - (write the code + no need to test it now) -> Transpose, adding over paxos indices and SHA calls  
         // 7. Share + Eval - after the PaxosEncode() -> write the transpose after the encode and compute the right format OT messages
         // 8. PSI_sender - testing!!
         // 9. Does the Paxos encode and decode correctly with xor of just 3 hash functions??

         

        
        return 0;
    }


    if (cmd.isSet(unitTestTag))
    {
        flagSet = true;
        auto tests = tests_cryptoTools::Tests;
        tests += tests_libOTe::Tests;

        auto r = tests.runIf(cmd);
        return r == TestCollection::Result::passed ? 0 : -1;
    }

    if (cmd.isSet("latency"))
    {
        getLatency(cmd);
        flagSet = true;
    }

#ifdef ENABLE_SIMPLESTOT
    flagSet |= runIf(baseOT_example<SimplestOT>, cmd, simple);
#endif

#ifdef ENABLE_SIMPLESTOT_ASM
    flagSet |= runIf(baseOT_example<AsmSimplestOT>, cmd, simpleasm);
#endif

#ifdef ENABLE_MRR_TWIST
#ifdef ENABLE_SSE
    flagSet |= runIf([&](Role role, int totalOTs, int numThreads, std::string ip, std::string tag, CLP& clp) {
        DomainSepEKEPopf factory;
        const char* domain = "EKE POPF OT example";
        factory.Update(domain, std::strlen(domain));
        baseOT_example_from_ot(role, totalOTs, numThreads, ip, tag, clp, McRosRoyTwist(factory));
    }, cmd, moellerpopf, {"eke"});
#endif

    flagSet |= runIf([&](Role role, int totalOTs, int numThreads, std::string ip, std::string tag, CLP& clp) {
        DomainSepMRPopf factory;
        const char* domain = "MR POPF OT example";
        factory.Update(domain, std::strlen(domain));
        baseOT_example_from_ot(role, totalOTs, numThreads, ip, tag, clp, McRosRoyTwistMR(factory));
    }, cmd, moellerpopf, {"mrPopf"});

    flagSet |= runIf([&](Role role, int totalOTs, int numThreads, std::string ip, std::string tag, CLP& clp) {
        DomainSepFeistelPopf factory;
        const char* domain = "Feistel POPF OT example";
        factory.Update(domain, std::strlen(domain));
        baseOT_example_from_ot(role, totalOTs, numThreads, ip, tag, clp, McRosRoyTwistFeistel(factory));
    }, cmd, moellerpopf, {"feistel"});

    flagSet |= runIf([&](Role role, int totalOTs, int numThreads, std::string ip, std::string tag, CLP& clp) {
        DomainSepFeistelMulPopf factory;
        const char* domain = "Feistel With Multiplication POPF OT example";
        factory.Update(domain, std::strlen(domain));
        baseOT_example_from_ot(role, totalOTs, numThreads, ip, tag, clp, McRosRoyTwistMul(factory));
    }, cmd, moellerpopf, {"feistelMul"});
#endif

#ifdef ENABLE_MRR
    flagSet |= runIf([&](Role role, int totalOTs, int numThreads, std::string ip, std::string tag, CLP& clp) {
        DomainSepFeistelRistPopf factory;
        const char* domain = "Feistel POPF OT example (Risretto)";
        factory.Update(domain, std::strlen(domain));
        baseOT_example_from_ot(role, totalOTs, numThreads, ip, tag, clp, McRosRoy(factory));
    }, cmd, ristrettopopf, {"feistel"});

    flagSet |= runIf([&](Role role, int totalOTs, int numThreads, std::string ip, std::string tag, CLP& clp) {
        DomainSepFeistelMulRistPopf factory;
        const char* domain = "Feistel With Multiplication POPF OT example (Risretto)";
        factory.Update(domain, std::strlen(domain));
        baseOT_example_from_ot(role, totalOTs, numThreads, ip, tag, clp, McRosRoyMul(factory));
    }, cmd, ristrettopopf, {"feistelMul"});
#endif

#ifdef ENABLE_FUZZY
    std::cout << "enable fuzzy " << std::endl; 
    flagSet |= runIf(baseOT_example<MasnyRindal>, cmd, f);
#endif

#ifdef ENABLE_MR
    flagSet |= runIf(baseOT_example<MasnyRindal>, cmd, mr);
#endif

#ifdef ENABLE_NP
    flagSet |= runIf(baseOT_example<NaorPinkas>, cmd, np);
#endif

#ifdef ENABLE_IKNP
    flagSet |= runIf(TwoChooseOne_example<IknpOtExtSender, IknpOtExtReceiver>, cmd, iknp);
#endif

#ifdef ENABLE_KOS
    flagSet |= runIf(TwoChooseOne_example<KosOtExtSender, KosOtExtReceiver>, cmd, kos);
#endif

#ifdef ENABLE_DELTA_KOS
    flagSet |= runIf(TwoChooseOne_example<KosDotExtSender, KosDotExtReceiver>, cmd, dkos);
#endif

#ifdef ENABLE_KKRT
    flagSet |= runIf(NChooseOne_example<KkrtNcoOtSender, KkrtNcoOtReceiver>, cmd, kkrt);
#endif

#ifdef ENABLE_OOS
    flagSet |= runIf(NChooseOne_example<OosNcoOtSender, OosNcoOtReceiver>, cmd, oos);
#endif

    flagSet |= runIf(Silent_example, cmd, Silent);
    flagSet |= runIf(Vole_example, cmd, vole);



    if (flagSet == false)
    {

        std::cout
            << "#######################################################\n"
            << "#                      - libOTe -                     #\n"
            << "#               A library for performing              #\n"
            << "#                  oblivious transfer.                #\n"
            << "#                     Peter Rindal                    #\n"
            << "#######################################################\n" << std::endl;


        std::cout
            << "Protocols:\n"
            << Color::Green << "  -simplest-asm " << Color::Default << "  : to run the ASM-SimplestOT  active secure  1-out-of-2  base OT      " << Color::Red << (spaEnabled ? "" : "(disabled)")             << "\n"   << Color::Default
            << Color::Green << "  -simplest     " << Color::Default << "  : to run the SimplestOT      active secure  1-out-of-2  base OT      " << Color::Red << (spEnabled ? "" : "(disabled)")              << "\n"   << Color::Default
            << Color::Green << "  -moellerpopf  " << Color::Default << "  : to run the McRosRoyTwist   active secure  1-out-of-2  base OT      " << Color::Red << (popfotMoellerEnabled ? "" : "(disabled)")   << "\n"   << Color::Default
            << Color::Green << "  -ristrettopopf" << Color::Default << "  : to run the McRosRoy active secure  1-out-of-2  base OT      " << Color::Red << (popfotRistrettoEnabled ? "" : "(disabled)") << "\n"   << Color::Default
            << Color::Green << "  -mr           " << Color::Default << "  : to run the MasnyRindal     active secure  1-out-of-2  base OT      " << Color::Red << (mrEnabled ? "" : "(disabled)")              << "\n"   << Color::Default
            << Color::Green << "  -np           " << Color::Default << "  : to run the NaorPinkas      active secure  1-out-of-2  base OT      " << Color::Red << (npEnabled ? "" : "(disabled)")              << "\n"   << Color::Default
            << Color::Green << "  -iknp         " << Color::Default << "  : to run the IKNP            passive secure 1-out-of-2       OT      " << Color::Red << (iknpEnabled ? "" : "(disabled)")            << "\n"   << Color::Default
            << Color::Green << "  -diknp        " << Color::Default << "  : to run the IKNP            passive secure 1-out-of-2 Delta-OT      " << Color::Red << (diknpEnabled ? "" : "(disabled)")           << "\n"   << Color::Default
            << Color::Green << "  -Silent       " << Color::Default << "  : to run the Silent          passive secure 1-out-of-2       OT      " << Color::Red << (silentEnabled ? "" : "(disabled)")          << "\n"   << Color::Default
            << Color::Green << "  -kos          " << Color::Default << "  : to run the KOS             active secure  1-out-of-2       OT      " << Color::Red << (kosEnabled ? "" : "(disabled)")             << "\n"   << Color::Default
            << Color::Green << "  -dkos         " << Color::Default << "  : to run the KOS             active secure  1-out-of-2 Delta-OT      " << Color::Red << (dkosEnabled ? "" : "(disabled)")            << "\n"   << Color::Default
            << Color::Green << "  -oos          " << Color::Default << "  : to run the OOS             active secure  1-out-of-N OT for N=2^76 " << Color::Red << (oosEnabled ? "" : "(disabled)")             << "\n"   << Color::Default
            << Color::Green << "  -kkrt         " << Color::Default << "  : to run the KKRT            passive secure 1-out-of-N OT for N=2^128" << Color::Red << (kkrtEnabled ? "" : "(disabled)")            << "\n\n" << Color::Default

            << "POPF Options:\n"
            << Color::Green << "  -eke          " << Color::Default << "  : to run the EKE POPF (Moeller only)                                  " << "\n"<< Color::Default
            << Color::Green << "  -mrPopf       " << Color::Default << "  : to run the MasnyRindal POPF (Moeller only)                          " << "\n"<< Color::Default
            << Color::Green << "  -feistel      " << Color::Default << "  : to run the Feistel POPF                                             " << "\n"<< Color::Default
            << Color::Green << "  -feistelMul   " << Color::Default << "  : to run the Feistel With Multiplication POPF                         " << "\n\n"<< Color::Default

            << "Other Options:\n"
            << Color::Green << "  -n            " << Color::Default << ": the number of OTs to perform\n"
            << Color::Green << "  -r 0/1        " << Color::Default << ": Do not play both OT roles. r 1 -> OT sender and network server. r 0 -> OT receiver and network client.\n"
            << Color::Green << "  -ip           " << Color::Default << ": the IP and port of the netowrk server, default = localhost:1212\n"
            << Color::Green << "  -t            " << Color::Default << ": the number of threads that should be used\n"
            << Color::Green << "  -u            " << Color::Default << ": to run the unit tests\n"
            << Color::Green << "  -u -list      " << Color::Default << ": to list the unit tests\n"

            << Color::Green << "  -u 1 2 15     " << Color::Default << ": to run the unit tests indexed by {1, 2, 15}.\n"
            << std::endl;
    }

    return 0;
}
/*
        auto t1 = high_resolution_clock::now();
        far_apart_FssEval(86, 25, okvs_fsskey1, delta, nkeys); 
        auto t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "FSS_Eval * 1M times took in milliseconds: " << duration << endl;
        
        far_apart_FssEval(96, 15, okvs_fsskey1, delta, nkeys); 
        far_apart_FssEval(96, 15, okvs_fsskey0, delta, nkeys); 

        
        far_apart_FssEval(27, 24, okvs_fsskey0, delta, nkeys);
        far_apart_FssEval(27, 24, okvs_fsskey1, delta, nkeys);


        far_apart_FssEval(66, 25, okvs_fsskey0, delta, nkeys); 
        far_apart_FssEval(66, 25, okvs_fsskey1, delta, nkeys); 

        far_apart_FssEval(52, 24, okvs_fsskey0, delta, nkeys);
        far_apart_FssEval(52, 24, okvs_fsskey1, delta, nkeys);

        far_apart_FssEval(30, 22, okvs_fsskey1, delta, nkeys);
        far_apart_FssEval(30, 22, okvs_fsskey0, delta, nkeys);
        
        auto t1 = high_resolution_clock::now();
        far_apart_FssEval(15, 15, okvs_fsskey0, delta, nkeys); 
        auto t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "FSS_Eval times took in milliseconds: " << duration << endl;
        //far_apart_FssEval(95, 29, okvs_fsskey1, delta, nkeys); 
        */