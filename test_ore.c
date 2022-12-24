#include "ore.h"
#include "errors.h"

#include <stdint.h>
#include <stdio.h>

#define ERR_CHECK(x) if((err = x) != ERROR_NONE) { return err; }

//./tests/test_ore < ../library/pbc-0.5.14/param/d159.param

/**
 * Generates a random 64-bit integers and encrypts it, then generates a
 * 64-bit integers and use it to generate a token
 * The encrypted integers are chosen randomly.
 *
 * @return 0 on success, -1 on failure, and an error if it occurred during the
 * encryption or comparison phase
 */
static int check_ore(pairing_t pairing, element_t g1, element_t g2, int err) {
    // length of plaintext
    uint32_t nbits = 32;

    // Randomly generate plaintext
    uint64_t n1 = rand() % ((uint64_t)1 << nbits);
    uint64_t n2 = rand() % ((uint64_t)1 << nbits);

    // 0 is equal, -1 is n1 < n2, 1 is n1 > n2
    int cmp = (n1 < n2) ? -1 : 1;
    if (n1 == n2) {
        cmp = 0;
    }

    ore_params params;
    ERR_CHECK(init_ore_params(params, nbits));

    element_t k1, k21, k22;
    element_init_Zr(k1, pairing);
    element_init_Zr(k21, pairing);
    element_init_Zr(k22, pairing);
    element_random(k1);
    element_random(k21);
    element_random(k22);

    ore_ciphertext ctxt;
    element_t xi, xj;//secret value
    // element_t kFi, kFj;//key for prf F
    mpz_t xjm_inv, xjm;//xj^(-1)
    mpz_t kFi_m;//key for prf F
    ERR_CHECK(init_ore_ciphertext(ctxt, params, pairing, g1, g2, xi, xj, k21, xjm_inv, xjm, kFi_m));
    // gmp_printf("xjm_inv:%Zd\n", xjm_inv);

    ore_token token, token_test;
    mpz_t mim[nbits], mim_inv[nbits];
    ore_exp_tmp exp_tmp;
    ERR_CHECK(init_ore_token(token, params, pairing, g2, k22, xj, mim, mim_inv, exp_tmp));
    // for(int i = 0; i < nbits; i++){
    //     gmp_printf("mim_inv2:%Zd\n", mim_inv[i]);
    // }

    mpz_t q;
    mpz_init_set_str(q, "625852803282871856053922297323874661378036491717", 10);
    mpz_t r;
    mpz_init_set_str(r, "208617601094290618684641029477488665211553761021", 10);

    ERR_CHECK(ore_encryption(ctxt, n1, pairing, k1, kFi_m, q, r));

    ERR_CHECK(ore_token_gen(token, n2, pairing, k1, xj, xi, g1, g2, q, r, xjm_inv, mim, mim_inv, exp_tmp));

    int ret = 0;
    int res;
    ERR_CHECK(ore_compare(&res, ctxt, token, pairing));
    if (res == cmp) {
        ret = 0;  // success
    }
    else {
        ret = -1; // fail
    }

    ERR_CHECK(clear_ore_ciphertext(ctxt));
    ERR_CHECK(clear_ore_token(token));

    element_clear(k1);
    element_clear(k21);
    element_clear(k22);
    return ret;
}

int main(int argc, char **argv) {
    srand((unsigned)time(NULL));

    printf("Testing ORE...\n");

    fflush(stdout);

    int err = 0;
    pairing_t pairing;
    element_t g1, g2;
    ERR_CHECK(init_pairing(pairing, g1, g2));

    int test_round = 10;//round
    for (int i = 0; i < test_round; i++) {
        printf("round %d\n", i + 1);

        if (check_ore(pairing, g1, g2, err) != ERROR_NONE) {
            printf("FAIL\n");
            return -1;
        }
    }

    printf("PASS\n");
    ERR_CHECK(clear_pairing(pairing, g1, g2));
    return 0;
}
