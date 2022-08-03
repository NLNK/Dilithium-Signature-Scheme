// gcc -march=native -mtune=native -funroll-loops -O3 -o verify verify.c ntt.c packing.c poly.c shake.c
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

int main(int argc, char **argv)
{    
    if(argc < 4)
  {
    printf("Please enter the Public key, document and signature file name....\n");
    exit(1);
  }  
    
  uint8_t rho[HASHBYTES],mu[CRHBYTES];    
  uint8_t c_tild[HASHBYTES], c2[HASHBYTES];
  poly c;
  polyvecl A[K], z;   
  polyveck t1, w1, h;
  char *publickey_buffer = malloc((CRYPTO_PUBLICKEYBYTES*2)*sizeof(*publickey_buffer));
  
  uint8_t unpack_publickey[CRYPTO_PUBLICKEYBYTES]={0};
  uint8_t sig[CRYPTO_BYTES];
  keccak_state state; //, state2,state3;
  char *publickey_value; long long int bytes;
  char *document_buffer;
  FILE *publickeyfile, *document, *signature;
  char *signature_buffer = malloc((CRYPTO_BYTES*2)*sizeof(*signature_buffer));
  uint8_t buf[K*POLYW1_PACKEDBYTES];    

  publickeyfile = fopen(argv[1],"r");

  int ret_val = fread(publickey_buffer,sizeof(char),CRYPTO_PUBLICKEYBYTES*2,publickeyfile);
  fclose(publickeyfile);
  if(ret_val == 0){exit(1); }

  for (int i=0;i<CRYPTO_PUBLICKEYBYTES;i++)
  {
      sscanf(publickey_buffer, "%2x", (unsigned int *)&unpack_publickey[i]);   //convert hex byte to integer
      publickey_buffer += 2;
  }
  /*  unpacking rho , t1 from public key */
  unpack_pk(rho,&t1,unpack_publickey);

  signature = fopen(argv[3],"r");
  ret_val = fread(signature_buffer,sizeof(char),CRYPTO_BYTES*2,signature);
  fclose(signature);
  if(ret_val == 0){exit(1); }

  for (int i=0;i<CRYPTO_BYTES;i++)
  {
    sscanf(signature_buffer, "%2x", (unsigned int *)&sig[i]);   //convert hex byte to integer
    signature_buffer += 2;
  }    

  document = fopen(argv[2], "rb");  // Open the file in binary mode
  fseek(document, 0, SEEK_END);          // Jump to the end of the file
  bytes = ftell(document);             // Get the current byte offset in the file
  rewind(document);                      // Jump back to the beginning of the file
  document_buffer = (char *)malloc(bytes * sizeof(char)); // Enough memory for the file
  ret_val = fread(document_buffer, bytes, 1, document); // Read in the entire file
  fclose(document);
  if(ret_val == 0){exit(1); }
  
  /* unpacking signature */
  if(unpack_sig(c_tild, &z, &h, sig))
    return -1;        

  if(polyvecl_checknorm(&z, GAMMA1 - BETA))
  return -1;   
   
  clock_t start =clock();  

  /* Expanding matrix A */  
  expandA(A,rho);
  
  // for(int i=0;i<K;i++)           
  //     polyvecl_ntt(&A[i]);   
  
  /* computing mu from publikey */
  shake256(mu, CRHBYTES, unpack_publickey, CRYPTO_PUBLICKEYBYTES);    
  
  shake256_init(&state);
  shake256_absorb(&state,mu,CRHBYTES);
  shake256_absorb(&state,document_buffer,bytes);
  shake256_finalize(&state);
  shake256_squeeze(mu,CRHBYTES,&state);    

  poly_challenge(&c,c_tild); 
 
  polyvecl_ntt(&z);  

  /* Reconstruct w1 */
  polyvec_matrix_mul(&w1, A, &z); 
  poly_ntt(&c);

  polyveck_shiftl(&t1);

  polyveck_ntt(&t1);
  polyveck_mul(&t1, &c, &t1);   // t1 = c*t1

  polyveck_sub(&w1, &w1, &t1);  // w1 = w1 - c*t1

  polyveck_intt(&w1);

  poly_modq(&w1);
  
  usehint(&w1,&h,&w1);

  polyveck_pack_w1(buf, &w1);
  
  /* Call random oracle and verify challenge */
  shake256_init(&state);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_absorb(&state, buf, K*POLYW1_PACKEDBYTES);
  shake256_finalize(&state);
  shake256_squeeze(c2, HASHBYTES, &state);
  
  // printf("\nc_tild value: \n");
  //   for (int i = 0; i < 32 ; i++)
  //   printf("%d, ",c_tild[i]);
  //   printf("\n");

  // printf("\nc2:\n");
  // for(int n = 0; n <HASHBYTES;n++)
  // {
  //     printf("%d, ",c2[n]);
  // }
  
  for(int i = 0; i < HASHBYTES; i++)
    if(c_tild[i] != c2[i])
    {
      printf("Verification failed\n");
      return -1;
    }
    else
    {
      printf("Verified Signature!!\n");
      printf("Time taken for verification:%f\n",(double)(clock() - start)/CLOCKS_PER_SEC*1000.0); 
      return 0;
    }
}

