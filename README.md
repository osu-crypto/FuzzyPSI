

A  C++17 implementation of Fuzzy PSI protocol written using the Oblivious Transfer extension [(OTe) library](https://github.com/osu-crypto/libOTe). 
 
 
## Introduction
 
This is an implementation of our new, fast [FuzzyPSI](https://eprint.iacr.org/2022/1011.pdf) protocol accepted at [CRYPTO'22](https://crypto.iacr.org/2022/). 
As a representative example we consider the case where Alice has a structured set A with 10
million total points, and Bob has an unstructured set B of 1.2 million points. We hold the total
cardinality of Alice’s set constant and consider two different ways that her points could be arranged
into balls:
• 6250 balls of radius δ = 20, in 2 dimensions
• 2778 balls of radius δ = 30, in 2 dimensions
We assign dimensions d = 2, our choice of bFSS is a (p = 0.75, k = 2)-bFSS instantiated within the structure-aware PSI framework. 

Our protocol involves cryptographic components (1) base oblivious transfers (2) hamming correlation hash function (3) encryption/decryption functionalities and a general communication framework. For the hamming-correlation robust hash function we use
SHA256. For base OTs and general framework we use the [libOTe library](https://github.com/osu-crypto/libOTe). We implement the bFSS recipe spatial hash ◦ concat ◦ tt, for balls of pairwise distance > 4δ, and use the OKVS implementation from [NTY21] for spatial hashing. 


## Build
 
We implement our protocol within the [LibOTe](https://github.com/osu-crypto/libOTe) library. Please follow the libOTe library to first build the repository and install the necessary dependencies. 

The library can be built as
```
git clone --recursive https://github.com/osu-crypto/libOTe.git
cd libOTe
python build.py --setup --boost --relic
python build.py -- -D ENABLE_RELIC=ON -D ENABLE_ALL_OT=ON
```
The main executable with examples is `frontend` and is located in the build directory, eg `out/build/linux/frontend/frontend.exe, out/build/x64-Release/frontend/Release/frontend.exe` depending on the OS. 

### Dependancies

Dependancies can be managed via the `build.py` script or or installed via an external tool. If an external tool is used install to system location or set  `-D CMAKE_PREFIX_PATH=path/to/install`.


## Install

libOTe can be installed and linked the same way as other cmake projects. By default the dependancies are not installed. To install all dependancies, run the following
```
python build.py --setup --boost --relic --sodium --install
```
You can also selectively install the dependancies. The install location can be specifying as `--install=path/to/installation`. Otherwise the system default is used.

The main library is similarly installed as
```
python build.py --install 
```

By default, sudo is not used. If installation requires sudo access, then add `--sudo` to the `build.py` script arguments. See `python build.py --help` for full details.



## Help
 
Contact Gayathri Garimella garimelg@oregonstate.edu for any assistance on building 
or running the library.

 
 ## License
 
This project has been placed in the public domain and/or MIT license. As such, you are unrestricted in how 
you use it, commercial or otherwise. However, no warranty of fitness is provided. If you 
found this project helpful, feel free to spread the word and cite us.

