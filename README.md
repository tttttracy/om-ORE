# om-ORE
## Prerequisites ##
Required environment
- [OpenSSL-1.1.1](https://www.openssl.org/source/)
- [GMP-6.2.0](https://gmplib.org/)
- [PBC-0.5.14](https://crypto.stanford.edu/pbc/download.html)
## Installation ##
``` shell
git clone git@github.com:tttttracy/om-ORE.git
cd om-ORE
make
```
## Run the test ##
Run the correctness check by 
``` shell
# Requires type-d parameter of PBC library as input to generate asymmetric pairing
./tests/test_ore < location_of_your_pbc_library/pbc-0.5.14/param/d159.param
``` 
Run the benchmark by
``` shell
./tests/time_ore < location_of_your_pbc_library/pbc-0.5.14/param/d159.param
``` 

