// gcc -march=native -mtune=native -funroll-loops -O3 -std=c99 -o sign sign.c ntt.c  packing.c poly.c shake.c
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include "dilithium.h"
#include "shake.h"

void polyvecl_uniform_gamma1(polyvecl *v, const uint8_t seed[CRHBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < L; i++)
    poly_uniform_gamma1(&v->vec[i], seed, L*nonce + i);
}

void polyvecl_modq2(polyvecl* w)
{
         for(int i=0;i<L;i++)
            for(int j=0;j<N;j++)
                w->vec[i].coeffs[j] = mod_q((w->vec[i].coeffs[j]+ (Q >> 1)),Q) - (Q >> 1);
}

void polyveck_modq2(polyveck* w)
{
          for(int i=0;i<K;i++)
            for(int j=0;j<N;j++)
                w->vec[i].coeffs[j] = mod_q((w->vec[i].coeffs[j]+ (Q >> 1)),Q) - (Q >> 1);
}

int main(int argc, char **argv)
{    
    if(argc < 4)
  {
    printf("Please enter the Private key, document and signature file name....\n");
    exit(1);
  }

    uint8_t hashbuf[2*HASHBYTES + 3*CRHBYTES]={0}, w1_pack[K*POLYW1_PACKEDBYTES]={0},c_tild[HASHBYTES]={0};  
    uint8_t *rho, *tr, *key, *rhoprime,*mu;
    // uint8_t rho[HASHBYTES], key[HASHBYTES], tr[CRHBYTES], rhoprime[CRHBYTES], mu[CRHBYTES];
    polyvecl A[K]={0}, s1, y,  z;  
    polyveck s2,  t0,  w, w1, w0, h;
    uint16_t nonce = 0;
    poly c;
    unsigned int hint=0;

    char *privatekey_buffer = malloc((CRYPTO_SECRETKEYBYTES*2)*sizeof(*privatekey_buffer));   
    uint8_t privatekey[CRYPTO_SECRETKEYBYTES]={0};
   
    long long int bytes;
    char *document_buffer;
    keccak_state shake_state, rhoprime_state,w1_state;
    int sign_flag=0;
    uint8_t sig[CRYPTO_BYTES]={0};
    
    FILE *skey, *document, *signature;
    
    skey = fopen(argv[1],"r");
    int ret_val = fread(privatekey_buffer,sizeof(char),CRYPTO_SECRETKEYBYTES*2,skey);    
    fclose(skey);
    if(ret_val == 0){exit(1); }

    for (int i=0;i<CRYPTO_SECRETKEYBYTES;i++)
    {
      sscanf(privatekey_buffer, "%2x", (unsigned int *)&privatekey[i]);   //convert hex byte to integer
      privatekey_buffer += 2;
    }

    document = fopen(argv[2], "rb");  // Open the file in binary mode
    fseek(document, 0, SEEK_END);          // Jump to the end of the file
    bytes = ftell(document);             // Get the current byte offset in the file
    rewind(document);                      // Jump back to the beginning of the file
    document_buffer = (char *)malloc(bytes * sizeof(char)); // Enough memory for the file
    ret_val = fread(document_buffer, bytes, 1, document); // Read in the entire file
    fclose(document);
    if(ret_val == 0){exit(1); }

    rho = hashbuf;
    tr = rho + HASHBYTES;
    key = tr + CRHBYTES;
    mu = key + HASHBYTES;
    rhoprime = mu + CRHBYTES;
    
    unpack_sk(rho,tr,key,&t0,&s1,&s2,privatekey);
  
    clock_t start =clock();
    
    expandA(A,rho);
    // for(int i=0;i<K;i++)           
    //   polyvecl_ntt(&A[i]); 

    /* Computing mu */ 
    shake256_init(&shake_state);
    shake256_absorb(&shake_state, tr, CRHBYTES);
    shake256_absorb(&shake_state, document_buffer, bytes);
    shake256_finalize(&shake_state);
    shake256_squeeze(mu, CRHBYTES, &shake_state);
    
    /* Computing rho prime */   
    // shake256_init(&rhoprime_state);
    // shake256_absorb(&rhoprime_state, key, HASHBYTES);
    // shake256_absorb(&rhoprime_state, mu, CRHBYTES);
    // shake256_finalize(&rhoprime_state);
    // shake256_squeeze(rhoprime, CRHBYTES, &rhoprime_state); 
    shake256(rhoprime, CRHBYTES, key, HASHBYTES + CRHBYTES);

    /* Computing NTT for s1, s2 and t0 */ 
    polyvecl_ntt(&s1);  
    polyveck_ntt(&s2);
    polyveck_ntt(&t0);  

    while(sign_flag==0){
      
      polyvecl_uniform_gamma1(&y, rhoprime, nonce++);       

      z = y;
      polyvecl_ntt(&z);      

     /* Matrix-vector multiplication */
      polyvec_matrix_mul(&w, A, &z);   
      polyveck_intt(&w);        
      poly_modq(&w);

      /* Decompose w for w0 and w1*/  
      polyveck_decompose(&w1, &w0, &w);
      
      /* Packing w1 */ 
      polyveck_pack_w1(w1_pack,&w1);

      /* Computing c_tild */ 
      shake256_init(&w1_state);
      shake256_absorb(&w1_state, mu, CRHBYTES);
      shake256_absorb(&w1_state, w1_pack, K*POLYW1_PACKEDBYTES);
      shake256_finalize(&w1_state);      
      shake256_squeeze(c_tild, HASHBYTES, &w1_state);          
      
      poly_challenge(&c, c_tild);        
      poly_ntt(&c);    

      polyvecl_mul(&z, &c, &s1);    // z = c*s1
      polyvecl_intt(&z);
      polyvecl_add(&z, &z, &y);     // z = y + c*s1
      polyvecl_modq2(&z);

      polyveck_mul(&h, &c, &s2);    // h = c*s2
      polyveck_intt(&h);

      polyveck_sub(&w, &w, &h);   // w = w - h = w - c*s2      
      polyveck_decompose(&w1, &w0, &w);   
     
      if(polyvecl_checknorm(&z, GAMMA1 - BETA) || polyveck_checknorm(&w0, GAMMA2 - BETA))
      {
        // printf("1. Norm not satisfied\n");       
        sign_flag = 0;
      }
      else{
        
        polyveck_mul(&h, &c, &t0);   // h= c*t0
        polyveck_intt(&h);        
        
        polyveck_add(&w0, &w, &h);    // w = w + c*t0           
          
        /* Compute hint bit indicating whether the low bits of the input element overflow into the high bits */
        hint = make_hint(&h, &w0, &w);          

        // printf("\nhint:%d",hint);

        // printf("\nh:\n");
        //   for(int i=0;i<K;i++)
        //   {
        //     for(int j=0;j<N;j++)
        //     {
        //         printf("%d, ",h.vec[i].coeffs[j]);
        //     }
        //     printf("\n");
        //   }

        if(polyveck_checknorm(&h, GAMMA2)|| (hint > OMEGA))
        {
          // printf("2. Norm not satisfied\n");
          sign_flag = 0;            
        }       
        else
        {             
          sign_flag = 1;  
          pack_sig(sig, c_tild, &z, &h);
          printf("Signature generated successfully!!\n");
          printf("Time taken for signature:%f\n",(double)(clock() - start)/CLOCKS_PER_SEC*1000.0); 
          
          // printf("w1:\n");
          // for(int i = 0; i <K;i++)
          // {
          //   for(int n = 0; n <N;n++)
          //   {
          //         printf("%d, ",w1.vec[i].coeffs[n]);
          //   }
          //   printf("\n");
          // }

          signature = fopen(argv[3],"w+");
          
          for(int i=0;i<CRYPTO_BYTES;i++)
            fprintf(signature,"%02X",sig[i]);
                    
          fclose(signature);       
        }
        
      }

    }
  return 0; 
}