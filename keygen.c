// gcc -march=native -mtune=native -funroll-loops -O3 -std=c99 -o keygen keygen.c ntt.c packing.c poly.c shake.c
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include "dilithium.h"
#include "shake.h"


void randombytes(uint8_t *out, size_t outlen) {
  static int fd = -1;
  ssize_t ret;

  while(fd == -1) {
    fd = open("/dev/urandom", O_RDONLY);
    if(fd == -1 && errno == EINTR)
      continue;
    else if(fd == -1)
      abort();
  }

  while(outlen > 0) {
    ret = read(fd, out, outlen);
    if(ret == -1 && errno == EINTR)
      continue;
    else if(ret == -1)
      abort();

    out += ret;
    outlen -= ret;
  }
}

static unsigned int rej_eta(int32_t *a, unsigned int len, const uint8_t *buf, unsigned int buflen)
{
  unsigned int ctr, pos;
  uint32_t t0, t1;

  ctr = pos = 0;
  while(ctr < len && pos < buflen) {
    t0 = buf[pos] & 0x0F;
    t1 = buf[pos++] >> 4;

#if ETA == 2
    if(t0 < 15) {
      t0 = t0 - (205*t0 >> 10)*5;
      a[ctr++] = 2 - t0;
    }
    if(t1 < 15 && ctr < len) {
      t1 = t1 - (205*t1 >> 10)*5;
      a[ctr++] = 2 - t1;
    }
#elif ETA == 4
    if(t0 < 9)
      a[ctr++] = 4 - t0;
    if(t1 < 9 && ctr < len)
      a[ctr++] = 4 - t1;
#endif
  }

  return ctr;
}

/* Sample polynomial with uniformly random coefficients in [-ETA,ETA] by performing 
rejection sampling on the output stream from SHAKE256(seed|nonce) */
#if ETA == 2
#define POLY_UNIFORM_ETA_NBLOCKS ((136 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#elif ETA == 4
#define POLY_UNIFORM_ETA_NBLOCKS ((227 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#endif
void poly_uniform_eta(poly *a, const uint8_t seed[CRHBYTES], uint16_t nonce)
{
  unsigned int ctr;
  unsigned int buflen = POLY_UNIFORM_ETA_NBLOCKS*STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_ETA_NBLOCKS*STREAM128_BLOCKBYTES];
  stream128_state state;

  shake128_stream_init(&state, seed, nonce);
  shake128_squeezeblocks(buf, POLY_UNIFORM_ETA_NBLOCKS, &state);

  ctr = rej_eta(a->coeffs, N, buf, buflen);

  while(ctr < N) {
    shake128_squeezeblocks(buf, 1, &state);
    ctr += rej_eta(a->coeffs + ctr, N - ctr, buf, STREAM128_BLOCKBYTES);
  }
}

void polyvecl_uniform_eta(polyvecl *v, const uint8_t seed[HASHBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < L; i++)
    poly_uniform_eta(&v->vec[i], seed, nonce++);
}

void polyveck_uniform_eta(polyveck *v, const uint8_t seed[HASHBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < K; i++)
    poly_uniform_eta(&v->vec[i], seed, nonce++);
}

int main(int argc,char** argv)
{
  if(argc < 3)
  {
    printf("Please enter the Private key and Public key file name....\n");
    exit(1);
  }
    
  uint8_t hashbuf[96]={0};  // 3*HASHBYTES
  uint8_t *rho, *rhoprime, *key;
  uint8_t tr[CRHBYTES]={0};
  FILE *publickeyfile, *privatekeyfile;
  polyvecl A[K]={0};
  uint8_t packed_pk[CRYPTO_PUBLICKEYBYTES]={0};
  uint8_t packed_sk[CRYPTO_SECRETKEYBYTES]={0};
  polyvecl s1, s1hat;
  polyveck s2, t, t0, t1;   

  randombytes(hashbuf, HASHBYTES);  // Generate randomness for the first 32 bytes   
  
  /* Test seed  */  
  /*
  char *test = "7c9935a0b07694aa0c6d10e4db6b1add2fd81a25ccb148032dcd739936737f2d";    
  char *temp=test;
  int ll=strlen(temp); 
  for (int i=0;i<ll/2;i++)
  {
    sscanf(temp, "%2x", (unsigned int *)&hashbuf[i]);   //convert hex byte to integer
    temp += 2;   
  }
  */
  
  shake256(hashbuf, 3*HASHBYTES, hashbuf, HASHBYTES);
  rho = hashbuf;
  rhoprime = hashbuf + HASHBYTES;
  key = hashbuf + 2*HASHBYTES;  

  clock_t start =clock();
  /* Expand matrix */ 
   expandA(A,rho);   
  //  for(int i=0;i<K;i++)           
  //      polyvecl_ntt(&A[i]); 
  
  /* Sample short vectors s1 and s2 */
  polyvecl_uniform_eta(&s1, rhoprime, 0);
  polyveck_uniform_eta(&s2, rhoprime, L);  
    
  s1hat = s1;
  polyvecl_ntt(&s1hat); 
  
  /* Matrix-vector multiplication */
  polyvec_matrix_mul(&t, A, &s1hat);  // t = NTT(A) * NTT(S1)

  polyveck_intt(&t); 
 
  polyveck_add(&t, &t, &s2);   // t = A*s1 + s2

  poly_modq(&t);

  polyveck_power2round(&t1, &t0, &t);

  /* Packing public key*/
  pack_pk(packed_pk,rho,&t1); 

  shake256(tr,CRHBYTES,packed_pk,CRYPTO_PUBLICKEYBYTES);
  
  /* packing private key */
  pack_sk(packed_sk,rho,tr,key,&t0,&s1,&s2);  

  printf("Key generated successfully!!\n");
  printf("Time taken for keygen:%f\n",(((double)(clock()-start)/CLOCKS_PER_SEC)*1000.0)); 

  publickeyfile = fopen(argv[2],"w+");
  
  for(int i=0;i<CRYPTO_PUBLICKEYBYTES;i++){
    fprintf(publickeyfile,"%02X",packed_pk[i]);
  }

  fclose(publickeyfile);

  privatekeyfile = fopen(argv[1],"w+");

  for(int i=0;i<CRYPTO_SECRETKEYBYTES;i++){
    fprintf(privatekeyfile,"%02X",packed_sk[i]);
  }
  
  fclose(privatekeyfile);

  return 0;
}
