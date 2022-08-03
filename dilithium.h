#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define mod_q(x,y) ((x % y) + y) % y

static inline int modulo(int64_t a,int b)
{
     int r = a/b;
     int64_t moding  = (a - r*b);  
     if(moding<0){ moding = b + moding; } 
     return (int)moding; 
}

#define HASHBYTES 32
#define CRHBYTES 48
#define N 256
#define Q 8380417
#define D 13

#define DILITHIUM_MODE 2

#if DILITHIUM_MODE == 2
#define K 4
#define L 4
#define ETA 2
#define TAU 39
#define BETA 78
#define GAMMA1 (1 << 17)
#define GAMMA2 ((Q-1)/88)
#define ALPHA (2*GAMMA2)
#define OMEGA 80

#elif DILITHIUM_MODE == 3
#define K 6
#define L 5
#define ETA 4
#define TAU 49
#define BETA 196
#define GAMMA1 (1 << 19)
#define GAMMA2 ((Q-1)/32)
#define ALPHA (2*GAMMA2)
#define OMEGA 55

#elif DILITHIUM_MODE == 5
#define K 8
#define L 7
#define ETA 2
#define TAU 60
#define BETA 120
#define GAMMA1 (1 << 19)
#define GAMMA2 ((Q-1)/32)
#define ALPHA (2*GAMMA2)
#define OMEGA 75

#endif

#define POLYT1_PACKEDBYTES  320
#define POLYT0_PACKEDBYTES  416
#define POLYVECH_PACKEDBYTES (OMEGA + K)

#if GAMMA1 == (1 << 17)
#define POLYZ_PACKEDBYTES   576
#elif GAMMA1 == (1 << 19)
#define POLYZ_PACKEDBYTES   640
#endif

#if GAMMA2 == (Q-1)/88
#define POLYW1_PACKEDBYTES  192
#elif GAMMA2 == (Q-1)/32
#define POLYW1_PACKEDBYTES  128
#endif

#if ETA == 2
#define POLYETA_PACKEDBYTES  96
#elif ETA == 4
#define POLYETA_PACKEDBYTES 128
#endif
#define CRYPTO_PUBLICKEYBYTES (HASHBYTES + K*POLYT1_PACKEDBYTES) 

#define CRYPTO_SECRETKEYBYTES (2*HASHBYTES + CRHBYTES  \
                               + L*POLYETA_PACKEDBYTES \
                               + K*POLYETA_PACKEDBYTES \
                               + K*POLYT0_PACKEDBYTES)

#define CRYPTO_BYTES (HASHBYTES + L*POLYZ_PACKEDBYTES + POLYVECH_PACKEDBYTES)

#define STREAM128_BLOCKBYTES 168
#define STREAM256_BLOCKBYTES 136

#ifndef SHA3_H
#define SHA3_H
#endif

typedef struct {
  int32_t coeffs[N];
} poly;

/* Vectors of polynomials of length K */
typedef struct {
  poly vec[K];
} polyveck;

/* Vectors of polynomials of length L */
typedef struct {
  poly vec[L];
} polyvecl;

void poly_uniform(poly *a,const uint8_t seed[HASHBYTES],uint16_t nonce);
void polyvecl_uniform_eta(polyvecl *v, const uint8_t seed[HASHBYTES], uint16_t nonce);
void polyveck_uniform_eta(polyveck *v, const uint8_t seed[HASHBYTES], uint16_t nonce);

int32_t power2round(int32_t *a0, int32_t a) ;
void poly_power2round(poly *a1, poly *a0, poly *a);
void polyveck_power2round(polyveck *v1, polyveck *v0, polyveck *v);
void poly_uniform_gamma1(poly *a, const uint8_t seed[CRHBYTES], uint16_t nonce);
void decompose(poly* r1, poly* r0, poly* r);
// int modulus(int32_t a,int32_t b);
void poly_modq(polyveck* w);
void poly_challenge(poly *c, const uint8_t seed[HASHBYTES]);
int polyvecl_checknorm(polyvecl *v, int32_t bound);
int polyveck_checknorm(polyveck *v, int32_t bound);

void lowbits(poly* r0, poly* r);
void highbits(poly* r1, poly* r);
void poly_neg(poly* neg, poly* a);
unsigned int poly_compare(poly* h,poly* a, poly* b);
unsigned int make_hint(polyveck* h,polyveck* z, polyveck* r);
void poly_shiftl(poly *tval);
void polyveck_shiftl(polyveck *a);


void ntt(int32_t a[N]);
void poly_ntt(poly *a);
void polyvecl_ntt(polyvecl *v);
void polyveck_ntt(polyveck *v);

void intt(int32_t a[N]);
void poly_intt(poly *a);
void polyvecl_intt(polyvecl *v);
void polyveck_intt(polyveck *v);

void poly_sub(poly* result, poly* a, poly* b);
void poly_add(poly* result, poly* a, poly* b);
void polyvecl_add(polyvecl *w, polyvecl *u, polyvecl *v);
void polyveck_add(polyveck *w, polyveck *u, polyveck *v);
void polyveck_sub(polyveck *w, polyveck *u, polyveck *v);
void poly_mul(poly *c, poly *a, poly *b);
void polyvecl_mul(polyvecl *r, poly *a, polyvecl *v);
void polyveck_mul(polyveck *r, poly *a, polyveck *v);
void polyvec_matrix_mul(polyveck *t, polyvecl mat[K], polyvecl *v);

void polyveck_decompose(polyveck *v1, polyveck *v0, polyveck *v);

void expandA(polyvecl *matrix, const uint8_t *rho);
void usehint(polyveck* w1,polyveck* h, polyveck* r);

void polyz_unpack(poly *r, const uint8_t *a);
void polyt1_pack(uint8_t *r, poly *a);
void polyt1_unpack(poly *r, const uint8_t *a);
void polyt0_pack(uint8_t *r, poly *a);
void polyt0_unpack(poly *r, const uint8_t *a);
void polyeta_pack(uint8_t *r, poly *a);
void polyeta_unpack(poly *r, const uint8_t *a);
void pack_pk(uint8_t pk[CRYPTO_PUBLICKEYBYTES],const uint8_t rho[32],polyveck *t1);
void unpack_pk(uint8_t rho[HASHBYTES],polyveck *t1,const uint8_t pk[CRYPTO_PUBLICKEYBYTES]);
void pack_sk(uint8_t sk[CRYPTO_SECRETKEYBYTES], const uint8_t rho[HASHBYTES],const uint8_t tr[CRHBYTES],const uint8_t key[HASHBYTES],polyveck *t0,polyvecl *s1,polyveck *s2);
void unpack_sk(uint8_t rho[HASHBYTES],uint8_t tr[CRHBYTES],uint8_t key[HASHBYTES],polyveck *t0,polyvecl *s1,polyveck *s2,const uint8_t sk[CRYPTO_SECRETKEYBYTES]);
void polyw1_pack(uint8_t *r, poly *a);
void polyveck_pack_w1(uint8_t r[K*POLYW1_PACKEDBYTES], polyveck *w1);
void pack_sig(uint8_t sig[CRYPTO_BYTES],const uint8_t c[HASHBYTES],polyvecl *z,polyveck *h);
int unpack_sig(uint8_t c[HASHBYTES],polyvecl *z,polyveck *h,const uint8_t sig[CRYPTO_BYTES]);

