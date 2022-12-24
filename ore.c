#include "crypto.h"
#include "errors.h"
#include "ore.h"

#include <stdio.h>
#include <stdlib.h>
#include <pbc/pbc.h>

// For permutation
typedef struct {
  uint32_t index;
  uint32_t rando;
} rand_permute[1];

// For qsort funtion
int comp(const void *p1, const void *p2) {
  rand_permute *c = (rand_permute *)p1;
  rand_permute *d = (rand_permute *)p2;
  return ((*(rand_permute *)c)->rando - (*(rand_permute *)d)->rando);
}

int init_ore_params(ore_params params, uint32_t nbits) {
  params->initialized = true;
  params->nbits = nbits;

  return ERROR_NONE;
}

int init_pairing(pairing_t pairing, element_t g1, element_t g2) {
  char param[1024];
  size_t count = fread(param, 1, 1024, stdin);
  if (!count) {
    pbc_die("input error");
    return ERROR_PAIRING_NOT_INITIALIZED;
    }
  pairing_init_set_buf(pairing, param, count);
  if (pairing_is_symmetric(pairing)) {return ERROR_PAIRING_IS_SYMMETRIC;}

  element_init_G1(g1, pairing);
  element_random(g1);
  element_init_G2(g2, pairing);
  element_random(g2);

  return ERROR_NONE;
}

int init_ore_ciphertext(ore_ciphertext ctxt, ore_params params,\
                        pairing_t pairing, element_t g1, element_t g2,\
                        element_t xi, element_t xj, element_t k21,\
                        mpz_t xjm_inv, mpz_t xjm, mpz_t kFi_m) {
  if (ctxt == NULL || params == NULL) {
    return ERROR_NULL_POINTER;
  }

  memcpy(ctxt->params, params, sizeof(ore_params));
  uint32_t nbits = ctxt->params->nbits;
  for (uint32_t i = 0; i < nbits; i++) {
    element_init_G1(ctxt->bit_ctxt[i], pairing);
  }
  
  element_t r1;
  element_init_Zr(r1, pairing);
  element_random(r1);
  element_init_G1(ctxt->g1r1k21, pairing);
  element_pow_zn(ctxt->g1r1k21, g1, r1);
  element_pow_zn(ctxt->g1r1k21, ctxt->g1r1k21, k21);//c_0=g1_r1k21
  element_clear(r1);

  element_init_Zr(xi, pairing);
  element_random(xi);

  element_init_Zr(xj, pairing);

  //calculate xj^(-1)
  mpz_t q, faiq, gcd;
  mpz_init(xjm_inv);
  mpz_init_set_str(faiq, "625852803282871856053922297323874661378036491716", 10);
  mpz_init_set_str(q, "625852803282871856053922297323874661378036491717", 10);
  mpz_init(xjm);
  mpz_init(gcd);

  do{
        element_random(xj);
        element_to_mpz(xjm, xj);//xj
        mpz_gcd(gcd, xjm, faiq);
    }while(mpz_cmp_ui(gcd, 1) != 0);

  mpz_invert(xjm_inv, xjm, faiq);//xjm_inv = xj^(-1)

  mpz_t g1m;
  mpz_init(kFi_m);
  mpz_init(g1m);
  element_to_mpz(g1m, (pairing->G1->get_x(g1)));
  mpz_t xim;
  mpz_init(xim);
  element_to_mpz(xim, xi);
  mpz_powm(kFi_m, g1m, xim, q);//kFi_m = g1m^xim = (g1.x)^xi mod q

  ctxt->initialized = true;
  return ERROR_NONE;
}

int init_ore_token(ore_token token, ore_params params, pairing_t pairing,\
                   element_t g2, element_t k22, element_t xj, mpz_t* mim, mpz_t* mim_inv, ore_exp_tmp exp_tmp) {
  if (token == NULL || params == NULL) {
    return ERROR_NULL_POINTER;
  }

  memcpy(token->params, params, sizeof(ore_params));
  uint32_t nbits = token->params->nbits;
  for (uint32_t i = 0; i < nbits; i++) {
    element_init_G2(token->token_bit[i]->add_one, pairing);
    element_init_G2(token->token_bit[i]->minus_one, pairing);
    mpz_init(exp_tmp->exp_bit[i]->add_one);
    mpz_init(exp_tmp->exp_bit[i]->minus_one);
  }
  
  element_t r2;
  element_init_Zr(r2, pairing);
  element_init_G2(token->g2r2k22, pairing);
  element_random(r2);
  element_pow_zn(token->g2r2k22, g2, r2);
  element_pow_zn(token->g2r2k22, token->g2r2k22, k22);//t_0=g2_r2k22
  

  //j calculate n random number mi^(-1)
  mpz_t q, faiq;
  mpz_init_set_str(faiq, "625852803282871856053922297323874661378036491716", 10);
  mpz_init_set_str(q, "625852803282871856053922297323874661378036491717", 10);

  element_t mi[nbits];  
  
  mpz_t gcd;
  mpz_init(gcd);
  for(int i = 0; i < nbits; i++) {
    element_init_Zr(mi[i], pairing);
    mpz_init(mim[i]);
    mpz_init(mim_inv[i]);
    do{
        element_random(mi[i]);
        element_to_mpz(mim[i], mi[i]);//mi
        mpz_gcd(gcd, mim[i], faiq);
    }while(mpz_cmp_ui(gcd, 1) != 0);
    mpz_invert(mim_inv[i], mim[i], faiq);//mim_inv = mi^(-1)
  }

  //cal k22r2
  mpz_init(exp_tmp->k22r2);
  mpz_t k22m, r2m;
  mpz_init(k22m);
  mpz_init(r2m);
  element_to_mpz(r2m, r2);
  element_to_mpz(k22m, k22);
  mpz_mul(exp_tmp->k22r2, k22m, r2m);
  mpz_mod(exp_tmp->k22r2, exp_tmp->k22r2, q);

  element_clear(r2);
  token->initialized = true;
  exp_tmp->initialized = true;
  return ERROR_NONE;
}

/**
 * Function which performs the encryption of an input, storing the result
 * in a ciphertext.
 *
 * @param ctxt    The ciphertext to store the encryption
 * @param buf     The input in a byte array, encoded in little-endian
 * @param buflen  The length of the byte array input
 * @param pairing pairing
 * @param k1      The prf key
 *
 * @return ERROR_NONE on success
 */
static int _ore_encryption(ore_ciphertext ctxt, byte *buf, uint32_t buflen,\
                           pairing_t pairing, element_t k1, mpz_t kFi_m, mpz_t q, mpz_t r) {
  if (!ctxt->initialized) {
    return ERROR_CTXT_NOT_INITIALIZED;
  }

  uint32_t nbits = ctxt->params->nbits;

  // byte length of ctxt
  uint32_t nbytes = (nbits + 7) / 8;
  
  // initial the input buffer for prf 
  byte prf_input_buf[sizeof(uint32_t) + nbytes];//4 + byte length of ctxt
  memset(prf_input_buf, 0, sizeof(prf_input_buf));

  // message buffer
  byte msgbuf[nbytes];



  // initial the output buffer for prf H(k1,x)
  byte prf_output_buf[SHA256_OUTPUT_BYTES];
  byte prf_input_buf_2[SHA256_OUTPUT_BYTES];
  byte key[32];
  memset(key, 0, sizeof(key));
  mpz_export(key, NULL, 1, 1, -1, 0, kFi_m);
  mpz_mod(kFi_m, kFi_m, q);

  // drop any extra bytes that have been provided
  if (buflen >= nbytes) {
    memcpy(msgbuf, buf, nbytes);
  }
  else {
    memcpy(msgbuf, buf, buflen);
  }

  mpz_t SHA_to_mpz;
  mpz_init(SHA_to_mpz);
  mpz_t SHA_256_MAX;
  mpz_init(SHA_256_MAX);
  mpz_init_set_str(SHA_256_MAX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\
                                 FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",16);


  element_t prf_tmp, prf_result;
  element_init_Zr(prf_result, pairing);
  element_init_Zr(prf_tmp, pairing);


  // index is the same as prf_input_buf[0]
  uint32_t *index = (uint32_t *)prf_input_buf;


  byte *value = &prf_input_buf[sizeof(uint32_t)];

  uint32_t offset = (8 - (nbits % 8)) % 8;

  // Generate a random array
  rand_permute *permute = (rand_permute *)malloc(sizeof(rand_permute) * nbits);
  for (uint32_t i = 0; i < nbits; i++) {
    permute[i]->index = i;
    permute[i]->rando = rand();
  }
  qsort(permute, nbits, sizeof(permute), comp);//升序排序

  mpz_t k1m;
  mpz_init(k1m);
  element_to_mpz(k1m, k1);
  for (uint32_t i = 0; i < nbits; i++) {//0->n-1
    // get the current bit of the message
    // little-endian
    uint32_t byteind = nbytes - 1 - (i + offset) / 8;//高位开始->低位

    // mask indicates whether this bit is 1 or 0
    byte mask = msgbuf[byteind] & (1 << ((7 - (i + offset)) % 8));

    HMAC_SHA256(prf_output_buf, sizeof(prf_output_buf), key,
            prf_input_buf, sizeof(prf_input_buf));//prf_output_buf = F(k1,i||b1b2...) = F(k1,x)    初始i = 0

    // first output as the second inout
    memcpy(prf_input_buf_2, prf_output_buf, 32);

    // If this bit is 1, it will be converted to mpz_t and then added 1 and modulo 2^256-1
    // and then converted to byte type as the second input
    if (mask > 0) {//bi=1
      mpz_import(SHA_to_mpz, 32, 1, 1, -1, 0, prf_output_buf);//converted to mpz_t: order=1 for most significant word first;endian=-1 for least significant first;The most significant nails=0 bits of each word are skipped,this can be 0 to use the full words.
      mpz_add_ui(SHA_to_mpz, SHA_to_mpz, 1);//F(k1,x)+1
      mpz_and(SHA_to_mpz, SHA_to_mpz, SHA_256_MAX);//ui=(F(k1,x)+1) & FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF = F(k1,x)+1 mod 2^256   lamda=256
      mpz_export(prf_input_buf_2, NULL, 1, 1, -1, 0, SHA_to_mpz);//converted to char : prf_input_buf_2=ui
    }

    // Compute the second hash and use the result to generate elements in Zr 
    // to be used in the subsequent modular exponential operation, covering the result of the first prf operation
    sha_256(prf_output_buf, sizeof(prf_output_buf), prf_input_buf_2, sizeof(prf_input_buf_2));//prf_output_buf=h(ui)

    mpz_import(SHA_to_mpz, 32, 1, 1, -1, 0, prf_output_buf);//converted to mpz_t :h(ui)
    mpz_mod(SHA_to_mpz, SHA_to_mpz, q);
    mpz_powm(SHA_to_mpz, SHA_to_mpz, k1m, q);  //h(ui)^k1 mod r
    mpz_mod(SHA_to_mpz, SHA_to_mpz, r);

    element_set_mpz(prf_result, SHA_to_mpz);
    element_pow_zn(ctxt->bit_ctxt[permute[i]->index], ctxt->g1r1k21, prf_result);//h1 = g1_r1k21H(k1,ui)

    // add the current bit of the message to the running prefix
    value[byteind] |= mask;//b1b1b2...||0^(n-i+1)

    // increment the index for the next iteration of the loop
    (*index)++;//i=i+1
  }

  free(permute);
  mpz_clear(SHA_to_mpz);
  mpz_clear(SHA_256_MAX);
  element_clear(prf_result);
  element_clear(prf_tmp);

  return ERROR_NONE;
}

int ore_encryption(ore_ciphertext ctxt, uint64_t msg, pairing_t pairing,\
                   element_t k1, mpz_t kFi_m, mpz_t q, mpz_t r) {
  return _ore_encryption(ctxt, (byte *)&msg, sizeof(msg), pairing, k1, kFi_m, q, r);
}

/**
 * Real main function which performs the token gen of an input, storing the result
 * in token.
 *
 * This function implements the encrypt algorithm for order revealing
 * encryption, using the secret key and input passed in and storing the
 * resulting ciphertext in ctxt.
 *
 * @param token    The ciphertext to store the encryption
 * @param buf      The input in a byte array, encoded in little-endian
 * @param buflen   The length of the byte array input
 * @param k1       The prf key
 *
 * @return ERROR_NONE on success
 */
static int _ore_token_gen(ore_token token, byte *buf, uint32_t buflen,\
                          pairing_t pairing, element_t k1, element_t xj,\
                          element_t xi, element_t g1, element_t g2, mpz_t q, mpz_t r, mpz_t xjm_inv,\
                          mpz_t* mim, mpz_t* mim_inv, ore_exp_tmp exp_tmp) {
  if (!token->initialized || !exp_tmp->initialized) {
    return ERROR_TOKEN_NOT_INITIALIZED;
  }

  uint32_t nbits = token->params->nbits;

  uint32_t nbytes = (nbits + 7) / 8;

  byte prf_input_buf[sizeof(uint32_t) + nbytes];
  memset(prf_input_buf, 0, sizeof(prf_input_buf));

  byte msgbuf[nbytes];
  byte prf_output_buf[SHA256_OUTPUT_BYTES];
  byte prf_input_buf_2[SHA256_OUTPUT_BYTES];

  // drop any extra bytes that have been provided
  if (buflen >= nbytes) {
    memcpy(msgbuf, buf, nbytes);
  }
  else {
    memcpy(msgbuf, buf, buflen);
  }

/******************step1 obtaining kF(optional)**************************/
  mpz_t kFi_m;
  mpz_init(kFi_m);
  mpz_t a, b, g1m, xim, xjm;
  mpz_init(a);
  mpz_init(b);
  mpz_init(g1m);
  mpz_init(xim);
  mpz_init(xjm);
  element_to_mpz(g1m,(pairing->G1->get_x(g1)));
  element_to_mpz(xim, xi);
  element_to_mpz(xjm, xj);


  mpz_powm(a, g1m, xjm, q);      //a = (g1.x)^xj
  mpz_powm(b, a, xim, q);        //b = a^xi
  mpz_powm(kFi_m,b, xjm_inv, q); //kFi = b^(1/xj)

/**********************step1 end(optional)********************************/

/**********************step2 generating token*****************************/
  mpz_t tmp_mpz;
  mpz_init(tmp_mpz);
  mpz_t SHA_to_mpz;
  mpz_init(SHA_to_mpz);
  mpz_t SHA_256_MAX;
  mpz_init(SHA_256_MAX);
  mpz_init_set_str(SHA_256_MAX, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\
                                 FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16);

  element_t prf_tmp, prf_result;
  element_init_Zr(prf_result, pairing);
  element_init_Zr(prf_tmp, pairing);

  uint32_t *index = (uint32_t *)prf_input_buf;

  byte *value = &prf_input_buf[sizeof(uint32_t)];

  uint32_t offset = (8 - (nbits % 8)) % 8;

  rand_permute *permute = (rand_permute *)malloc(sizeof(rand_permute) * nbits);
  for (uint32_t i = 0; i < nbits; i++) {
    permute[i]->index = i;
    permute[i]->rando = rand();
  }
  qsort(permute, nbits, sizeof(permute), comp);


  //client j
  byte key[32];
  memset(key, 0, sizeof(key));
  mpz_export(key, NULL, 1, 1, -1, 0, kFi_m);
  mpz_mod(kFi_m, kFi_m, q);

  mpz_t k1m, kFi_e_mi;
  mpz_init(k1m);
  mpz_init(kFi_e_mi);
  element_to_mpz(k1m, k1);

  for (uint32_t i = 0; i < nbits; i++) {
    // get the current bit of the message
    uint32_t byteind = nbytes - 1 - (i + offset) / 8;

    byte mask = msgbuf[byteind] & (1 << ((7 - (i + offset)) % 8));

    HMAC_SHA256(prf_output_buf, sizeof(prf_output_buf), key,
            prf_input_buf, sizeof(prf_input_buf));//F(kFi,x)

    mpz_powm(kFi_e_mi, kFi_m, mim[i], q);// kFi^mi mod q  --const

    element_t dr;
    element_init_Zr(dr,pairing);
    element_random(dr);
    mpz_t d, drm, g1_k1m, g1_k1_inv, kfi_r, kfi_r_inv, g1xm;
    mpz_init(d);
    mpz_init(drm);
    mpz_init(g1_k1m);//g1^k1
    mpz_init(g1_k1_inv);
    mpz_init(kfi_r);//kFi^r
    mpz_init(kfi_r_inv);
    mpz_init(g1xm);
    element_to_mpz(g1xm, (pairing->G1->get_x(g1)));//g1.x
    element_to_mpz(drm, dr);//r

    mpz_powm(g1_k1m, g1xm, k1m, q);//(g1.x)^k1
    mpz_powm(kfi_r, kFi_m, drm, q);//kFi^r
    mpz_powm(d, g1xm, drm, q);//d = (g1.x)^r
    // gmp_printf("d 1: %Zd\n", d);
    
    mpz_invert(g1_k1_inv, g1_k1m, q);//((g1.x)^k1)^(-1)
    mpz_invert(kfi_r_inv, kfi_r, q);//(kFi^r)^(-1)

    mpz_t con_e;
    mpz_init(con_e);
    mpz_mul(d, d, g1_k1_inv);//  d/((g1.x)^k1) mod q
    mpz_mod(d, d, q);
    mpz_powm(con_e, d, xim, q);// e = (d/((g1.x)^k1))^xi mod q  --const

    if (mask > 0) {//bi=1
      memcpy(prf_input_buf_2, prf_output_buf, 32);

      sha_256(prf_output_buf, sizeof(prf_output_buf), prf_input_buf_2, sizeof(prf_input_buf_2));//prf_output_buf=h(ui-1)

      mpz_import(tmp_mpz, 32, 1, 1, -1, 0, prf_output_buf);
      mpz_mod(tmp_mpz, tmp_mpz, q);

      mpz_powm(tmp_mpz, tmp_mpz, mim[i], q);// h(ui-1)^mi mod q

      /****************online****************/
      mpz_mul(exp_tmp->exp_bit[i]->minus_one, tmp_mpz, kFi_e_mi);//a = (kFi^mi)*(h(ui-1)^mi) mod q
      mpz_mod(exp_tmp->exp_bit[i]->minus_one, exp_tmp->exp_bit[i]->minus_one, q);

      mpz_powm(tmp_mpz, exp_tmp->exp_bit[i]->minus_one, k1m, q); //b=a^k1 mod q

      mpz_powm(tmp_mpz, tmp_mpz, mim_inv[i], q);//c = b^(1/mi) mod q
     
      mpz_t e;
      mpz_init(e);
      mpz_mul(e, con_e, tmp_mpz);//f = e*c mod q 
      mpz_mod(e, e, q);
      // gmp_printf("e 1: %Zd\n", e);//f = (g^rxi)*(h^k1) mod q

      mpz_mul(e, e, kfi_r_inv);//h = f*(kFi^r)^(-1) mod q
      mpz_mod(e, e, q);//h mod q
      mpz_mod(e, e, r);//h mod r

      /****************online****************/
      element_set_mpz(prf_result, e);
      element_pow_zn(token->token_bit[permute[i]->index]->minus_one, token->g2r2k22, prf_result);//g2_r2k22H(k1,ui-1)

      //ui+1
      mpz_import(SHA_to_mpz, 32, 1, 1, -1, 0, prf_input_buf_2);
      mpz_add_ui(SHA_to_mpz, SHA_to_mpz, 2);//F(k1,x)+2 = F(k1,x)+1+1 = ui+1
      mpz_and(SHA_to_mpz, SHA_to_mpz, SHA_256_MAX);//ui+1=(F(k1,x)+1)+1 mod 2^256
      mpz_export(prf_input_buf_2, NULL, 1, 1, -1, 0, SHA_to_mpz);
     
      sha_256(prf_output_buf, sizeof(prf_output_buf), prf_input_buf_2, sizeof(prf_input_buf_2));//prf_output_buf=h(ui+1)

      mpz_import(tmp_mpz, 32, 1, 1, -1, 0, prf_output_buf);
      mpz_mod(tmp_mpz, tmp_mpz, q);
      mpz_powm(tmp_mpz, tmp_mpz, mim[i], q);//h(ui+1)^mi mod q

      /****************online****************/
      mpz_mul(exp_tmp->exp_bit[i]->add_one, tmp_mpz, kFi_e_mi);//a = (kFi^mi)*(h(ui+1)^mi) mod q
      mpz_mod(exp_tmp->exp_bit[i]->add_one, exp_tmp->exp_bit[i]->add_one, q);

      mpz_powm(tmp_mpz, exp_tmp->exp_bit[i]->add_one, k1m, q); //b=a^k1 mod q

      mpz_powm(tmp_mpz, tmp_mpz, mim_inv[i], q);//c = b^(1/mi) mod q

      mpz_mul(e, con_e, tmp_mpz);//f = e*c mod q 
      mpz_mod(e, e, q);

      mpz_mul(e, e, kfi_r_inv);//h = f*(kFi^r)^(-1) mod q
      mpz_mod(e, e, q);//h mod q
      mpz_mod(e, e, r);//h mod r
      /****************online****************/

      //mpz_powm(tmp_mpz, tmp_mpz, mim[i], q);//H(k1,ui+1) = h(ui+1)^k1 mod r
      //mpz_powm(tmp_mpz, tmp_mpz, mim_inv[i], q);//test
      //mpz_powm(tmp_mpz, tmp_mpz, k1m, q);//test
      //mpz_mod(tmp_mpz, tmp_mpz, r);//test
      //element_set_mpz(prf_result, tmp_mpz);
      element_set_mpz(prf_result, e);
      element_pow_zn(token->token_bit[permute[i]->index]->add_one, token->g2r2k22, prf_result);//g2_r2k22H(k1,ui+1)
    } else {//bi=0
      mpz_import(SHA_to_mpz, 32, 1, 1, -1, 0, prf_output_buf);
      mpz_add_ui(SHA_to_mpz, SHA_to_mpz, 1);//SHA_to_mpz=F(k1,x)+1
      mpz_and(SHA_to_mpz, SHA_to_mpz, SHA_256_MAX);//ui+1=(F(k1,x)+0)+1 mod 2^256
      mpz_export(prf_input_buf_2, NULL, 1, 1, -1, 0, SHA_to_mpz);
     
      sha_256(prf_output_buf, sizeof(prf_output_buf), prf_input_buf_2, sizeof(prf_input_buf_2));//prf_output_buf=h(ui)

      mpz_import(tmp_mpz, 32, 1, 1, -1, 0, prf_output_buf);
      mpz_mod(tmp_mpz, tmp_mpz, q);
      mpz_powm(tmp_mpz, tmp_mpz, mim[i], q);//H(k1,ui+1) = h(ui+1)^mi mod q

      /****************online****************/
      mpz_mul(exp_tmp->exp_bit[i]->add_one, tmp_mpz, kFi_e_mi);//a = (kFi^mi)*(h(ui+1)^mi) mod q
      mpz_mod(exp_tmp->exp_bit[i]->add_one, exp_tmp->exp_bit[i]->add_one, q);

      mpz_powm(tmp_mpz, exp_tmp->exp_bit[i]->add_one, k1m, q); //b=a^k1 mod q

      mpz_powm(tmp_mpz, tmp_mpz, mim_inv[i], q);//c = b^(1/mi) mod q

      mpz_t e;
      mpz_init(e);
      mpz_mul(e, con_e, tmp_mpz);//f = e*c mod q 
      mpz_mod(e, e, q);

      mpz_mul(e, e, kfi_r_inv);//h = f*(kFi^r)^(-1) mod q
      mpz_mod(e, e, q);//h mod q
      mpz_mod(e, e, r);//h mod r
      /****************online****************/

      //mpz_powm(tmp_mpz, tmp_mpz, mim_inv[i], q);//test
      //mpz_powm(tmp_mpz, tmp_mpz, k1m, q);//H(k1,ui+1) = h(ui+1)^k1 mod q
      //mpz_mod(tmp_mpz, tmp_mpz, r);//H(k1,ui+1) = h(ui+1)^k1 mod r
      //element_set_mpz(prf_result, tmp_mpz);
      element_set_mpz(prf_result, e);
      element_pow_zn(token->token_bit[permute[i]->index]->add_one , token->g2r2k22, prf_result);//g2_r2k22H(k1,ui+1)

      mpz_sub_ui(SHA_to_mpz, SHA_to_mpz, 2);//F(k1,x)+1-2 = F(k1,x)-1 
      mpz_and(SHA_to_mpz, SHA_to_mpz, SHA_256_MAX);//ui-1=F(k1,x)-1 mod 2^256
      mpz_export(prf_input_buf_2, NULL, 1, 1, -1, 0, SHA_to_mpz);
      
      sha_256(prf_output_buf, sizeof(prf_output_buf), prf_input_buf_2, sizeof(prf_input_buf_2));//prf_output_buf=h(ui)

      mpz_import(tmp_mpz, 32, 1, 1, -1, 0, prf_output_buf);
      mpz_mod(tmp_mpz, tmp_mpz, q);
      mpz_powm(tmp_mpz, tmp_mpz, mim[i], q);

      /****************online****************/
      mpz_mul(exp_tmp->exp_bit[i]->minus_one, tmp_mpz, kFi_e_mi);//a = (kFi^mi)*(h(ui-1)^mi) mod q
      mpz_mod(exp_tmp->exp_bit[i]->minus_one, exp_tmp->exp_bit[i]->minus_one, q);

      mpz_powm(tmp_mpz, exp_tmp->exp_bit[i]->minus_one, k1m, q); //b=a^k1 mod q

      mpz_powm(tmp_mpz, tmp_mpz, mim_inv[i], q);//c = b^(1/mi) mod q
     
      mpz_mul(e, con_e, tmp_mpz);//f = e*c mod q 
      mpz_mod(e, e, q);
      // gmp_printf("e 1: %Zd\n", e);//f = (g^rxi)*(h^k1) mod q

      mpz_mul(e, e, kfi_r_inv);//h = f*(kFi^r)^(-1) mod q
      mpz_mod(e, e, q);//h mod q
      mpz_mod(e, e, r);//h mod r
      /****************online****************/

      // mpz_powm(tmp_mpz, tmp_mpz, mim_inv[i], q);//test
      //mpz_powm(tmp_mpz, tmp_mpz, k1m, q);//test
      //mpz_mod(tmp_mpz, tmp_mpz, r);//test//H(k1,ui+1) = h(ui+1)^k1 mod r
      //element_set_mpz(prf_result, tmp_mpz);
      element_set_mpz(prf_result, e);
      element_pow_zn(token->token_bit[permute[i]->index]->minus_one , token->g2r2k22, prf_result);//g2_r2k22H(k1,ui-1)
    }
    
    // add the current bit of the message to the running prefix
    value[byteind] |= mask;

    // increment the index for the next iteration of the loop
    (*index)++;
  }

  
  free(permute);
  mpz_clear(SHA_to_mpz);
  mpz_clear(SHA_256_MAX);
  element_clear(prf_result);
  element_clear(prf_tmp);

  return ERROR_NONE;
}

int ore_token_gen(ore_token token, uint64_t msg, pairing_t pairing, element_t k1, element_t xj,\
                  element_t xi, element_t g1, element_t g2, mpz_t q, mpz_t r, mpz_t xjm_inv, \
                  mpz_t* mim, mpz_t* mim_inv, ore_exp_tmp exp_tmp) {
  return _ore_token_gen(token, (byte *)&msg, sizeof(msg), pairing, k1, xj, xi, g1, g2, q, r, xjm_inv, mim, mim_inv, exp_tmp);
}

int ore_compare(int *result_p, ore_ciphertext ctxt, ore_token token, pairing_t pairing) {

  if ((ctxt->params->initialized != token->params->initialized) ||
      (ctxt->params->nbits != token->params->nbits)) {
    return ERROR_PARAMS_MISMATCH;
  }
  if(!ctxt->initialized) {
    return ERROR_CTXT_NOT_INITIALIZED;
  }
  if (!token->initialized) {
    return ERROR_TOKEN_NOT_INITIALIZED;
  }

  int res = 0;

  element_t temp1, temp2, temp3;
  element_init_GT(temp1, pairing);
  element_init_GT(temp2, pairing);
  element_init_GT(temp3, pairing);

  uint32_t nbits = token->params->nbits;

  // Because the ciphertext and token have been randomly permuted,
  // this test function requires O(n^2) compare complexity, but
  // it only need 3n pairing because we can compute and store
  // pairing results for token before perform the real compare.

  // preprocessing
  pairing_pp_t pp;
  pairing_pp_init(pp, ctxt->g1r1k21, pairing);
  // temporarily store pairing results for token and ctxt.
  ore_token_bit token_result[nbits];
  for (uint32_t i = 0; i < nbits; i++) {
    element_init_GT(token_result[i]->add_one, pairing);
    element_init_GT(token_result[i]->minus_one, pairing);
    pairing_pp_apply(token_result[i]->add_one, token->token_bit[i]->add_one, pp);//e(g1_r1k21,h2)
    pairing_pp_apply(token_result[i]->minus_one, token->token_bit[i]->minus_one, pp);//e(g1_r1k21,h3)
  }
  pairing_pp_clear(pp);

  bool break_flag = false;
  for (uint32_t i = 0; i < nbits; i++) {
    pairing_apply(temp1, ctxt->bit_ctxt[i], token->g2r2k22, pairing);//temp1 = e(h1,g2_r2k22)
    for (uint32_t j = 0; j < nbits; j++) {
      if (!element_cmp(temp1, token_result[j]->add_one)) {//if e(h1,g2_r2k22) = e(g1_r1k21,h2)
        res = 1;//n1 > n2
        break_flag = true;
        break;
      }
      if (!element_cmp(temp1, token_result[j]->minus_one)) {//if e(h1,g2_r2k22) = e(g1_r1k21,h3)
        res = -1;//n1 < n2
        break_flag = true;
        break;
      }
    }
    if (break_flag)
      break;
  }

  // clear pairing results for token.
  for (uint32_t i = 0; i < nbits; i++) {
    element_clear(token_result[i]->add_one);
    element_clear(token_result[i]->minus_one);
  }

  *result_p = res;
  element_clear(temp1);
  element_clear(temp2);
  element_clear(temp3);
  return ERROR_NONE;
}

int clear_ore_ciphertext(ore_ciphertext ctxt) {
  if (ctxt == NULL) {
    return ERROR_NONE;
  }
  uint32_t nbits = ctxt->params->nbits;
  for (uint32_t i = 0; i < nbits; i++) {
    element_clear(ctxt->bit_ctxt[i]);
  }
  element_clear(ctxt->g1r1k21);

  return ERROR_NONE;
}

int clear_ore_token(ore_token token) {
  if (token == NULL) {
    return ERROR_NONE;
  }
  uint32_t nbits = token->params->nbits;
  for (uint32_t i = 0; i < nbits; i++) {
    element_clear(token->token_bit[i]->add_one);
    element_clear(token->token_bit[i]->minus_one);
  }
  element_clear(token->g2r2k22);

  return ERROR_NONE;
}

int clear_pairing(pairing_t pairing, element_t g1, element_t g2) {
  // Pairing must be the last one to be cleared.
  element_clear(g1);
  element_clear(g2);
  pairing_clear(pairing);

  return ERROR_NONE;
}