
#include "dilithium.h"
#include "shake.h"

// void expandA(polyvecl *matrix, const uint8_t *rho){
//     keccak_state expand_state;
//     uint8_t ob[K*L*N];
//     uint32_t t;
//     shake256_init(&expand_state);
//     shake256_absorb(&expand_state, rho, HASHBYTES);
//     shake256_finalize(&expand_state);    
//     shake256_squeeze(ob,(K*L*N), &expand_state);
//     for(int i=0;i<K;i++){
//       for(int j=0;j<L;j++){
//         for(int m = 0;m<N;m++){         
//           t  = ob[(i*L)+j] << 16;
//           t |= (uint32_t)ob[(i*L)+j+1] << 8;
//           t |= (uint32_t)ob[(i*L)+j+2];         
//           matrix[i].vec[j].coeffs[m] = mod_q(t,Q);
//           t = 0;
//         }
//       }
//     }

// }

void expandA(polyvecl *matrix, const uint8_t *rho){
    keccak_state expand_state;
    uint32_t t;
    unsigned int pos =0;
    unsigned ctr = K*L*N;
    uint8_t buf[ctr];
    shake128_init(&expand_state);
    shake128_absorb(&expand_state, rho, HASHBYTES);
    shake128_finalize(&expand_state);    
    shake128_squeeze(buf, ctr, &expand_state);
 
    for(int i=0;i<K;i++){
      for(int j=0;j<L;j++){
        for(int m = 0;m<N;){    
          t  = buf[pos++];
          t |= (uint32_t)buf[pos++] << 8;
          t |= (uint32_t)buf[pos++] << 16;          
          t &= 0x7FFFFF;
          if((t < Q) && (pos < ctr))
            matrix[i].vec[j].coeffs[m++] = t; 
          else if((pos >= ctr)&& m<N)  
          {          
            shake128_squeeze(buf,N, &expand_state);
            pos=0;
          }
            
        }
      }
    }

}

/* Compute mod without using % operator */    
// int modulus(int32_t a,int32_t b)
// {
//   int32_t r = a/b;
//   int32_t moding  = (a - r*b);
//   if(moding<0)
//     moding = b + moding;
//   return moding;
// }

void decompose(poly* r1, poly* r0, poly* r)
{
  poly r_temp;
  for(int i=0;i<N;++i)
  {
      r_temp.coeffs[i] = mod_q(r->coeffs[i],Q);
      r0->coeffs[i]  = mod_q((r_temp.coeffs[i] + (ALPHA >> 1)),ALPHA) - (ALPHA >> 1);
      if((r_temp.coeffs[i] - r0->coeffs[i]) == (Q-1))
      {
        r1->coeffs[i] = 0;
        r0->coeffs[i] = r0->coeffs[i]-1;
      }
      else
      {
        r1->coeffs[i] = ((r_temp.coeffs[i]-r0->coeffs[i])/ALPHA);
      }
  }  

}

void highbits(poly* r1, poly* r)
{
  poly r0;
  decompose(r1,&r0,r);
}

void lowbits(poly* r0, poly* r)
{
  poly r1;
  decompose(&r1,r0,r);
}

void poly_shiftl(poly *tval) { 
  for(int i = 0; i < N; i++)
    tval->coeffs[i] <<= D; 
}

void polyveck_shiftl(polyveck *a)
{
   for (int i = 0; i < K; i++)
     poly_shiftl(&a->vec[i]);
}

int32_t power2round(int32_t *a0, int32_t a)  {
  int32_t a1;
  a1 = (a + (1 << (D-1)) - 1) >> D;
  *a0 = a - (a1 << D);  
  return a1;
}

void poly_power2round(poly *a1, poly *a0, poly *a) {  

  for(int i = 0; i < N; i++)
    a1->coeffs[i] = power2round(&a0->coeffs[i], a->coeffs[i]); 
}

void polyveck_power2round(polyveck *v1, polyveck *v0, polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; i++)
    poly_power2round(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}


#if GAMMA1 == (1 << 17)
#define POLY_UNIFORM_GAMMA1_NBLOCKS ((576 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#elif GAMMA1 == (1 << 19)
#define POLY_UNIFORM_GAMMA1_NBLOCKS ((640 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#endif

/* Sample polynomial with uniformly random coefficients in [-(GAMMA1 - 1), GAMMA1] by unpacking output stream of SHAKE256(seed|nonce) */
void poly_uniform_gamma1(poly *a,
                         const uint8_t seed[CRHBYTES],
                         uint16_t nonce)
{
  uint8_t buf[POLY_UNIFORM_GAMMA1_NBLOCKS*STREAM256_BLOCKBYTES];
  stream256_state state;

  shake256_stream_init(&state, seed, nonce);
  shake256_squeezeblocks(buf, POLY_UNIFORM_GAMMA1_NBLOCKS, &state);  
  polyz_unpack(a, buf);  
}

/* Samples polynomial with TAU nonzero coefficients in {-1,1} using the output stream of SHAKE256 */
void poly_challenge(poly *c, const uint8_t seed[HASHBYTES]) {
  unsigned int i, b, pos;
  uint64_t signs;
  uint8_t buf[STREAM256_BLOCKBYTES];
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, seed, HASHBYTES);
  shake256_finalize(&state);
  shake256_squeezeblocks(buf, 1, &state);

  signs = 0;
  for(i = 0; i < 8; i++)
    signs |= (uint64_t)buf[i] << 8*i;
  pos = 8;

  for(i = 0; i < N; i++)
    c->coeffs[i] = 0;
  for(i = N-TAU; i < N; i++) {
    do {
      if(pos >= STREAM256_BLOCKBYTES) {
        shake256_squeezeblocks(buf, 1, &state);
        pos = 0;
      }

      b = buf[pos++];
    } while(b > i);

    c->coeffs[i] = c->coeffs[b];
    c->coeffs[b] = 1 - 2*(signs & 1);
    signs >>= 1;
  }
}

int poly_checknorm(poly *a, int32_t B)
{   
  int32_t max_element;
  max_element = abs(a->coeffs[0]);

  for (int i = 1; i < N; i++) {    
    if (max_element < abs(a->coeffs[i])) {
      max_element = abs(a->coeffs[i]);
    }
  }   
  
  // for (int i = 0; i < N; i+=2) {    
  //   if(max_element < abs(a->coeffs[i])) {
  //     max_element = abs(a->coeffs[i]);
  //   }  
  //   if(max_element < abs(a->coeffs[i+1])) {
  //     max_element = abs(a->coeffs[i+1]);
  //   }    
  // }
    if(max_element >= B) { 
      return 1;
    }

  return 0;
}

int polyvecl_checknorm(polyvecl *v, int32_t bound)  {
  unsigned int i;

  for(i = 0; i < L; i++)
    if(poly_checknorm(&v->vec[i], bound))
      return 1;

  return 0;
}

int polyveck_checknorm(polyveck *v, int32_t bound)  {
  unsigned int i;

  for(i = 0; i < K; i++)
    if(poly_checknorm(&v->vec[i], bound))
      return 1;

  return 0;
}

/* Compute neg */    
void poly_neg(poly* neg, poly* a)
{
  for(int i=0;i<N;i++)
  {
    neg->coeffs[i] = (-1)*a->coeffs[i];
  }

}

/*Addition of two polynomials  */    
void poly_sub(poly* result, poly* a, poly* b)
{
  for(int i=0;i<N;i++)
  {
    result->coeffs[i] = a->coeffs[i]-b->coeffs[i];
  }

}

/*Addition of two polynomials  */    
void poly_add(poly* result, poly* a, poly* b)
{
  for(int i=0;i<N;i++)  
    result->coeffs[i] = a->coeffs[i]+b->coeffs[i];     

}
/* Compare two polynomials */    
unsigned int poly_compare(poly* h,poly* a, poly* b)
{
  unsigned int s=0;
  for(int i=0;i<N;i++)
  {
    if(a->coeffs[i]==b->coeffs[i]){
      h->coeffs[i] = 0;
    }  
    else{       
       h->coeffs[i] = 1;
       ++s;
    }
  }   
  return s;
}

/* Compute hint bit indicating whether the low bits of the input element overflow into the high bits */
// unsigned int make_hint(polyveck* h,polyveck* z, polyveck* r)
unsigned int make_hint(polyveck* h,polyveck* r, polyveck* z)
{
  unsigned int hint=0;
  poly r1,v1; 

  for(int i=0;i<K;i++)
  {
      highbits(&r1,&r->vec[i]);
      highbits(&v1,&z->vec[i]);      
      hint+=poly_compare(&h->vec[i],&r1,&v1);  
      // printf("\nh:%d\t",hint);   
  }
  
  return hint;
    
}

void usehint(polyveck* w1,polyveck* h, polyveck* r)
{
    int m = (Q-1)/(2*GAMMA2);
    polyveck r0;    
    for(int i = 0; i < K; i++)
    { 
      decompose(&w1->vec[i],&r0.vec[i],&r->vec[i]);
      for(int j = 0; j < N; j++)
      {   
          if((h->vec[i].coeffs[j]==1) && (r0.vec[i].coeffs[j]>0))
              w1->vec[i].coeffs[j] = mod_q((w1->vec[i].coeffs[j]+1),m);
              
          if((h->vec[i].coeffs[j]==1) && (r0.vec[i].coeffs[j]<=0))          
              w1->vec[i].coeffs[j] = mod_q((w1->vec[i].coeffs[j]-1),m);
              
      }
    }
    
    // for(int i = 0; i < K; i++)
    // {
    //   for(int j = 0; j < N; j++)
    //   {   
    //       if((h->vec[i].coeffs[j]==1) && (r0.vec[i].coeffs[j]>0))
    //           w1->vec[i].coeffs[j] = mod_q((w1->vec[i].coeffs[j]+1),m);
              
    //       if((h->vec[i].coeffs[j]==1) && (r0.vec[i].coeffs[j]<=0))          
    //           w1->vec[i].coeffs[j] = mod_q((w1->vec[i].coeffs[j]-1),m);
              
    //   }
    // }
}


void poly_ntt(poly *a) { 
  ntt(a->coeffs);
}

void polyvecl_ntt(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; i++)
    poly_ntt(&v->vec[i]);
}

void polyveck_ntt(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; i++)
    poly_ntt(&v->vec[i]);
}


void poly_intt(poly *a) {
  intt(a->coeffs); 
}

void polyvecl_intt(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; i++)
    poly_intt(&v->vec[i]);
}

void polyveck_intt(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; i++)
    poly_intt(&v->vec[i]);
}


void polyveck_decompose(polyveck *v1, polyveck *v0, polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; i++)
    decompose(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}

void polyvecl_add(polyvecl *w, polyvecl *u, polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; i++)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyveck_add(polyveck *w, polyveck *u, polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; i++)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}


void polyveck_sub(polyveck *w, polyveck *u, polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; i++)
    poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyveck_mul(polyveck *r, poly *a, polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; i++)
    poly_mul(&r->vec[i], a, &v->vec[i]);
}


void poly_mul(poly *c, poly *a, poly *b) {
  unsigned int i;

  for(i = 0; i < N; i++)
     c->coeffs[i] = mod_q(((int64_t)a->coeffs[i] * b->coeffs[i]) + (Q >> 1), Q) - (Q >> 1);
}

void polyvecl_pointwise_acc(poly *w, polyvecl *u,polyvecl *v)
{
  unsigned int i;
  poly t;

  poly_mul(w, &u->vec[0], &v->vec[0]);
  for(i = 1; i < L; i++) {
    poly_mul(&t, &u->vec[i], &v->vec[i]);
    poly_add(w, w, &t);
  }
}

void polyvec_matrix_mul(polyveck *t, polyvecl mat[K], polyvecl *v) {
  unsigned int i;

  for(i = 0; i < K; i++)
    polyvecl_pointwise_acc(&t->vec[i], &mat[i], v);
}

void polyvecl_mul(polyvecl *r, poly *a, polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; i++)
    poly_mul(&r->vec[i], a, &v->vec[i]);
}

void poly_modq(polyveck* w)
{
         for(int i=0;i<K;i++)         
            for(int j=0;j<N;j++)
                w->vec[i].coeffs[j] = mod_q(w->vec[i].coeffs[j],Q);       
          
}


void shake128_stream_init(keccak_state *state,
                                    const uint8_t seed[HASHBYTES],
                                    uint16_t nonce)
{
  uint8_t t[2]={0};
  t[0] = nonce;
  t[1] = nonce >> 8;  

  shake128_init(state);
  shake128_absorb(state, seed, HASHBYTES);
  shake128_absorb(state, t, 2);
  shake128_finalize(state);
}

void shake256_stream_init(keccak_state *state,
                                    const uint8_t seed[CRHBYTES],
                                    uint16_t nonce)
{
  uint8_t t[2]={0};
  t[0] = nonce;
  t[1] = nonce >> 8;

  shake256_init(state);
  shake256_absorb(state, seed, CRHBYTES);
  shake256_absorb(state, t, 2);
  shake256_finalize(state);
}


