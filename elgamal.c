#include <unistd.h>
#include "elgamal.h"
#include <time.h>


/*
  Sets r to a random GMP integer with the specified number
  of bits.
*/
void get_random_n_bits(mpz_t r, size_t bits)
{
	size_t size = (size_t) ceilf(bits/8);
	char *buffer = (char*) malloc(sizeof(char)*size);
	int prg = open("/dev/random", O_RDONLY);
	read(prg, buffer, size);
	close(prg);
	mpz_import (r, size, 1, sizeof(char), 0, 0, buffer);
	free(buffer);
}


/*
  Sets r to a random GMP *prime* integer, smaller than max.
*/
void get_random_n_prime(mpz_t r, mpz_t max) 
{
	do {
		get_random_n_bits(r, mpz_sizeinbase(max, 2));
		mpz_nextprime(r, r);
	} while (mpz_cmp(r, max) >= 0);
}


/*
  Sets r to a random GMP integer smaller than max.
*/
void get_random_n(mpz_t r, mpz_t max) 
{
	do {
		get_random_n_bits(r, mpz_sizeinbase(max, 2));
	} while (mpz_cmp(r, max) >= 0);
}


/*
 Init structure. Set domain parameters p, q and g
 */
void init_elgam(elgam_ctx **ectx, size_t bits) 
{
	*ectx = (elgam_ctx*) malloc(sizeof(elgam_ctx));
	// 1. find large prime p for domain parameter
	mpz_t p, g, x, h;
	mpz_init((*ectx)->dom_par_p);
	mpz_init((*ectx)->dom_par_g);
	mpz_init((*ectx)->priv_x);
	mpz_init((*ectx)->pub_h);
	mpz_init((*ectx)->eph_k);
	get_random_n_bits((*ectx)->dom_par_p, bits);
	mpz_nextprime((*ectx)->dom_par_p, (*ectx)->dom_par_p);
	gmp_printf("\n\np = %Zd\n", (*ectx)->dom_par_p);

	get_random_n_prime((*ectx)->dom_par_g, (*ectx)->dom_par_p);
	gmp_printf("g = %Zd\n", (*ectx)->dom_par_g);

	get_random_n((*ectx)->priv_x, (*ectx)->dom_par_p);
	gmp_printf("x = %Zd\n", (*ectx)->priv_x);
	/* h = g^x (mod n) */
	mpz_powm_sec((*ectx)->pub_h, (*ectx)->dom_par_g, (*ectx)->priv_x, (*ectx)->dom_par_p);
	gmp_printf("h = %Zd\n\n", (*ectx)->pub_h);
}



void destroy_elgam(elgam_ctx *ectx) 
{
	if (ectx) {
		mpz_clears(ectx->dom_par_p, ectx->dom_par_g, ectx->dom_par_q, NULL);
		mpz_clears(ectx->priv_x, ectx->pub_h, ectx->eph_k, NULL);
		free(ectx);
		ectx = NULL;
	}
}


void destroy_ciphertxt(ciphertext *ct) 
{
	if (ct) {
		mpz_clears(ct->c1, ct->c2, NULL);
		free(ct);
		ct = NULL;
	}
}


ciphertext* encrypt(mpz_t m, elgam_ctx *ectx) 
{
	ectx->eph_k;
	get_random_n(ectx->eph_k, ectx->dom_par_p);
	ciphertext *ct = malloc(sizeof(ciphertext));
	mpz_init(ct->c1);
	mpz_init(ct->c2);
	mpz_powm_sec(ct->c1, ectx->dom_par_g, ectx->eph_k, ectx->dom_par_p);
	mpz_powm_sec(ct->c2, ectx->pub_h, ectx->eph_k, ectx->dom_par_p);
	mpz_mul(ct->c2, m, ct->c2);
	mpz_mod(ct->c2, ct->c2, ectx->dom_par_p);
	gmp_printf("c1 = %Zd\n", ct->c1);
	gmp_printf("c2 = %Zd\n\n", ct->c2);
	return ct;
}


void decrypt(mpz_t msg, ciphertext *ct, elgam_ctx *ectx) 
{
	mpz_powm_sec(ct->c1, ct->c1, ectx->priv_x, ectx->dom_par_p);
	mpz_invert(ct->c1, ct->c1, ectx->dom_par_p);
	mpz_mul(msg, ct->c2, ct->c1);
	mpz_mod(msg, msg, ectx->dom_par_p);
}


/* setup elliptic curve, public and private key
 Using the brainpoolP160r1 - EC domain parameters
 http://www.ecc-brainpool.org/download/Domain-parameters.pdf
 */
void init_elgam_ec(elgam_ec_ctx **eec_ctx, char* curve_name)
{
    *eec_ctx = (elgam_ec_ctx*) malloc(sizeof(elgam_ec_ctx));
    elliptic_curve *ecc = malloc(sizeof(elliptic_curve));
    (*eec_ctx)->ec = ecc;

	init_point(&(ecc->base));

	// select curve
	if(strcmp(curve_name, "secp224k1")==0){
		mpz_set_str(ecc->a, "00000000000000000000000000000000000000000000000000000000", 16); 
		mpz_set_str(ecc->b, "00000000000000000000000000000000000000000000000000000005", 16); 
		mpz_set_str(ecc->p, "fffffffffffffffffffffffffffffffffffffffffffffffeffffe56d", 16); 
		mpz_set_str(ecc->base->x, "a1455b334df099df30fc28a169a467e9e47075a90f7e650eb6b7a45c", 16); 
		mpz_set_str(ecc->base->y, "7e089fed7fba344282cafbd6f7e319f7c0b0bd59e2ca4bdb556d61a5", 16); 
	}else if(strcmp(curve_name, "brainpoolP160r1")==0){
		mpz_set_str(ecc->a, "340E7BE2A280EB74E2BE61BADA745D97E8F7C300", 16); 
		mpz_set_str(ecc->b, "1E589A8595423412134FAA2DBDEC95C8D8675E58", 16); 
		mpz_set_str(ecc->p, "E95E4A5F737059DC60DFC7AD95B3D8139515620F", 16); 
		mpz_set_str(ecc->base->x, "BED5AF16EA3F6A4F62938C4631EB5AF7BDBCDBC3", 16); 
		mpz_set_str(ecc->base->y, "1667CB477A1A8EC338F94741669C976316DA6321", 16); 		
	}
	printf("use the curve: %s\n", curve_name);

}
void key_gen(elgam_ec_ctx **eec_ctx){

	mpz_init((*eec_ctx)->priv_key);
	init_point(&((*eec_ctx)->pub_key));
	get_random_n((*eec_ctx)->priv_key, (*eec_ctx)->ec->p);
	mpz_t tmp;
	mpz_init_set(tmp, (*eec_ctx)->priv_key);
	(*eec_ctx)->pub_key = ecc_scalar_mul((*eec_ctx)->ec, tmp, (*eec_ctx)->ec->base);


	// gmp_printf("\np = %Zd\n", (*eec_ctx)->ec->p);
	// gmp_printf("x = %Zd\n", (*eec_ctx)->priv_key);

	mpz_clears(tmp, NULL);
	// gmp_printf("Base point P = (%Zd,%Zd)\n", (*eec_ctx)->ec->base->x, (*eec_ctx)->ec->base->y);
	// gmp_printf("Public key xP =  (%Zd,%Zd)\n\n", ((*eec_ctx)->pub_key)->x, ((*eec_ctx)->pub_key)->y);
}

void test_init_elgam_ec(elgam_ec_ctx **eec_ctx)
{
    *eec_ctx = (elgam_ec_ctx*) malloc(sizeof(elgam_ec_ctx));
    elliptic_curve *ecc = malloc(sizeof(elliptic_curve));
    (*eec_ctx)->ec = ecc;

	mpz_init_set_ui(ecc->a, 1);
	mpz_init_set_ui(ecc->b, 3);
	mpz_init_set_ui(ecc->p, 23);

	mpz_init((*eec_ctx)->priv_key);
	init_point(&(ecc->base));
	init_point(&((*eec_ctx)->pub_key));

	//mpz_init_set_ui(ecc->base->x, 21);
	//mpz_init_set_ui(ecc->base->y, 1);
}


void destroy_elgam_ec(elgam_ec_ctx *eec_ctx) 
{
	if (eec_ctx) {
 		mpz_clears(eec_ctx->priv_key, eec_ctx->eph_k, NULL);
		mpz_clears(eec_ctx->ec->a, eec_ctx->ec->b, eec_ctx->ec->p, NULL);
		destroy_point(eec_ctx->ec->base);
		destroy_point(eec_ctx->pub_key);
		if (eec_ctx->ec) {
			free(eec_ctx->ec);
			eec_ctx->ec = NULL;
		}
		free(eec_ctx);
		eec_ctx = NULL;
	}
}


cipherec* encrypt_ec(elgam_ec_ctx *eec, point *pm)
{
	// gmp_printf("Encrypted: (%Zd,%Zd)\n", pm->x, pm->y);  

	mpz_init(eec->eph_k);
	get_random_n(eec->eph_k, eec->ec->p);
	// gmp_printf("\nEphemeral key = %Zd\n", eec->eph_k);

	cipherec *cipher = malloc(sizeof(cipherec));
	init_point(&cipher->c1);
	init_point(&cipher->c2);
	mpz_t tmp;
	mpz_init_set(tmp, eec->eph_k);
	cipher->c1 = ecc_scalar_mul(eec->ec, tmp, eec->ec->base);
	mpz_clears(tmp, NULL);

	mpz_init_set(tmp, eec->eph_k);
	cipher->c2 = ecc_scalar_mul(eec->ec, tmp, eec->pub_key);
	mpz_clears(tmp, NULL);
	// gmp_printf("Cipher c1: (%Zd,%Zd)\n", cipher->c1->x, cipher->c1->y);
	// gmp_printf("Cipher c2 without msg: (%Zd,%Zd)\n", cipher->c2->x, cipher->c2->y);
	cipher->c2 = ecc_addition(eec->ec, cipher->c2, pm);
	// gmp_printf("Cipher c2 with msg: (%Zd,%Zd)\n", cipher->c2->x, cipher->c2->y);
	mpz_clears(eec->eph_k, NULL);
	return cipher;
}


point* decrypt_ec(elgam_ec_ctx *eec, cipherec *c)
{
  	point *d1, *d2;
  	init_point(&d1);
  	init_point(&d2);
	mpz_t tmp;
  	mpz_init_set(tmp, eec->priv_key);
  	d1 = ecc_scalar_mul(eec->ec, tmp, c->c1);

  	mpz_clears(tmp, NULL);
  	// gmp_printf("D1=(%Zd,%Zd)\n", d1->x, d1->y);
	// gmp_printf("Before neg: (%Zd,%Zd)\n", d1->x, d1->y);
	mpz_neg(d1->y, d1->y);
  	// gmp_printf("After neg: (%Zd,%Zd)\n", d1->x, d1->y);
  	d2 = ecc_addition(eec->ec, c->c2, d1);
  	// gmp_printf("Decrypted: (%Zd,%Zd)\n", d2->x, d2->y);
	destroy_point(d1);
	return d2;
}


void destroy_cipherec(cipherec *c)
{
	if (c) {
		destroy_point(c->c1);
		destroy_point(c->c2);
		free(c);
		c = NULL;
	}
}


void test() {

	elgam_ec_ctx *eec;
	test_init_elgam_ec(&eec);
	// P + Q = R = (5,15)  
	point *p, *q;
	init_point(&p);
	init_point(&q);
	mpz_set_ui(p->x, 10);
	mpz_set_ui(p->y, 1);
	mpz_set_ui(q->x, 21);
	mpz_set_ui(q->y, 4);

	point *r;
	init_point(&r);
	point *r2;
	init_point(&r);	
	ecc_addition2(eec->ec, p, q, r);
	gmp_printf("Addition R=(%Zd, %Zd)\n", r->x, r->y);

	ecc_doubling2(eec->ec, p, r);
	// 2P = R = (4,5)
	gmp_printf("R=(%Zd, %Zd)\n", r->x, r->y);
	// 4P = R = (19,2)
	// 3P = R = (12,8)
	mpz_t m;
	mpz_init(m);
	mpz_set_ui(m, 4);

	r2 = ecc_scalar_mul2(eec->ec, m, p); 
	gmp_printf("R=(%Zd, %Zd)\n", r2->x, r2->y);

	// free up stuff
	mpz_clears(m, NULL);

	destroy_point(p);
	destroy_point(q);
	destroy_point(r);
	destroy_point(r2);
	destroy_elgam_ec(eec);

}


int main() 
{
	//test();
	// ElGamal-EC
	int run_times = 20, total_time = 0;
    clock_t  start, stop;
    double duration;
	elgam_ec_ctx *eec;
	// char* curve = "brainpoolP160r1";
	char* curve = "secp224k1";
	init_elgam_ec(&eec, curve);

	// generate key pair
	for(int i = 0; i < run_times; i++){
		start = clock();
		key_gen(&eec);
		stop = clock(); 
		duration = ((double)(stop - start)*1000.0)/CLOCKS_PER_SEC;
		total_time += duration;
		printf( "%d KeyGen %f ms\n", i+1, duration );  
	}
	printf( "Average Run Time: KeyGen %f ms\n", total_time/run_times ); 

	point *p;
	cipherec *c;
	init_point(&p);
	mpz_init_set_ui(p->x, 666);
	mpz_init_set_ui(p->y, 123);
	gmp_printf("Encrypted: (%Zd,%Zd)\n", p->x, p->y);  
	// c1 = kP (rand k * base point)
	// c2 = xkP + Pm (public key xP * rand k) + point on curve (secret message)

	// encrypt
	total_time = 0;
	for(int i = 0; i < run_times; i++){
		start = clock();
		c = encrypt_ec(eec, p);
		stop = clock(); 
		duration = ((double)(stop - start)*1000.0)/CLOCKS_PER_SEC;
		total_time += duration;
		printf( "%d Enc %f ms\n", i+1, duration );  
	}
	printf( "Average Run Time: Enc %f ms\n", total_time/run_times ); 

	destroy_point(p);
	init_point(&p);
	// c1 * x = c1' = xkP
	// Pm = c1' - c2 = xkP - xkP + Pm

	// decrypt
	total_time = 0;
	for(int i = 0; i < run_times; i++){
		start = clock();
		p = decrypt_ec(eec, c);
		stop = clock(); 
		duration = ((double)(stop - start)*1000.0)/CLOCKS_PER_SEC;
		total_time += duration;
		printf( "%d Dec %f ms\n", i+1, duration );  
	}
	printf( "Average Run Time: Dec %f ms\n", total_time/run_times ); 

	curve = "brainpoolP160r1";
	init_elgam_ec(&eec, curve);

	// generate key pair
	for(int i = 0; i < run_times; i++){
		start = clock();
		key_gen(&eec);
		stop = clock(); 
		duration = ((double)(stop - start)*1000.0)/CLOCKS_PER_SEC;
		total_time += duration;
		printf( "%d KeyGen %f ms\n", i+1, duration );  
	}
	printf( "Average Run Time: KeyGen %f ms\n", total_time/run_times ); 

	gmp_printf("Encrypted: (%Zd,%Zd)\n", p->x, p->y);  
	// c1 = kP (rand k * base point)
	// c2 = xkP + Pm (public key xP * rand k) + point on curve (secret message)

	// encrypt
	total_time = 0;
	for(int i = 0; i < run_times; i++){
		start = clock();
		c = encrypt_ec(eec, p);
		stop = clock(); 
		duration = ((double)(stop - start)*1000.0)/CLOCKS_PER_SEC;
		total_time += duration;
		printf( "%d Enc %f ms\n", i+1, duration );  
	}
	printf( "Average Run Time: Enc %f ms\n", total_time/run_times ); 

	destroy_point(p);
	init_point(&p);
	// c1 * x = c1' = xkP
	// Pm = c1' - c2 = xkP - xkP + Pm

	// decrypt
	total_time = 0;
	for(int i = 0; i < run_times; i++){
		start = clock();
		p = decrypt_ec(eec, c);
		stop = clock(); 
		duration = ((double)(stop - start)*1000.0)/CLOCKS_PER_SEC;
		total_time += duration;
		printf( "%d Dec %f ms\n", i+1, duration );  
	}
	printf( "Average Run Time: Dec %f ms\n", total_time/run_times ); 

	destroy_point(p);
	destroy_cipherec(c);
	destroy_elgam_ec(eec);

}



