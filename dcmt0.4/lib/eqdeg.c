/* eqdeg.c */

/* Copyright (C) 2001-2007 Makoto Matsumoto and Takuji Nishimura.  */
/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */ 
/* 02111-1307  USA                                                 */

#include <stdio.h>
#include <stdlib.h>
#include "dci.h"

/**** mt19937 *****/
void _sgenrand_dc(uint32_t);
uint32_t _genrand_dc(void);
/******************/

/**************************************/
#define SSS 7
#define TTT 15
/* #define S00 11 */
#define S00 12 
#define S01 18
/**************************************/

/** for get_tempering_parameter_hard **/
#define LIMIT_V_BEST_OPT 15 
/**************************************/

#define WORD_LEN 32
#define MIN_INFINITE (-2147483647-1)

#define TRNSTMP  tmp ^= ( (tmp>>shift_0) & greal_mask )
#define MASKTMP  tmp ^= ( (tmp<<shift_s) & mask_b);tmp ^= ( (tmp<<shift_t) & mask_c)
#define LSB(x) ( ((x) >> ggap) & 0x1 )

typedef struct {
    uint32_t z;    /* integer part */
    uint32_t *cf;  /* fraction part */
    int start;     /* beginning of fraction part */
    int degree;    /* maximum degree */
    uint32_t bp;   /* bitpatterm (shifted&bitmasked) at the maximum degree */
}Vector;

typedef struct mask_node{
    uint32_t b,c;
    int v,leng;
    struct mask_node *next;
}MaskNode;
	
void _get_tempering_parameter_dc(mt_struct *mts);
void _get_tempering_parameter_hard_dc(mt_struct *mts);

static void show_distrib(mt_struct *mts);
static int push_stack(uint32_t b, uint32_t c, int v, uint32_t *bbb, uint32_t *ccc);
static int push_mask(int l, int v, uint32_t b, uint32_t c, uint32_t *bbb, uint32_t *ccc);
static int lenstra(int v);
static int degree_of_vector( Vector *v);
static void init_tempering(mt_struct *mts);
static void free_Vector( Vector *v );
static void free_lattice( Vector **lattice, int v);
static void add( Vector *u, Vector *v);
static void add_i_to_j( Vector **lattice, int i, int j, int k);
static void pull_max( int i, int v, Vector *vec);
static void pull_min_row( int i, int v, Vector **lattice);
static void hakidasi( int k, Vector **lattice);
static void optimize_v(uint32_t b, uint32_t c, int v);
static MaskNode *optimize_v_hard(int v, MaskNode *prev);
static Vector *new_Vector(void);
static Vector **make_lattice(int v);
static Vector *mult( Vector *v, int k);

static void delete_MaskNodes(MaskNode *head);
static MaskNode *delete_lower_MaskNodes(MaskNode *head, int l);
static MaskNode *cons_MaskNode(MaskNode *head, uint32_t b, uint32_t c, int leng);
/* static void count_MaskNodes(MaskNode *head); */

/***********************************/
/*******  global variables  ********/
/***********************************/
static uint32_t bitmask[WORD_LEN];
static uint32_t cur_bitmask[WORD_LEN];
static uint32_t mask_b, mask_c;
static uint32_t upper_v_bits;
static int shift_0, shift_1, shift_s, shift_t;
static int mmm, nnn, rrr, www;
static uint32_t aaa[2];
/** most significant  (WWW - RRR) bits **/
static uint32_t gupper_mask;
/** least significant RRR bits **/
static uint32_t glower_mask;
/** upper WWW bitmask **/
static uint32_t greal_mask;
/** difference between machine wordsize and dest wordsize **/
static int ggap;
/** for optimize_v_hard **/
static int gcur_maxlengs[WORD_LEN];
static uint32_t gmax_b, gmax_c;
/*************************************/

#if 0
int main(int argc, char **argv)
{
    mt_struct mt = {AAA,MMM,NNN,RRR,WWW,0,0,0,S00,SSS,TTT,0,0,0,NULL};
    
    get_tempering_parameter(&mt);

    return 0;
}
#endif

void _get_tempering_parameter_dc(mt_struct *mts)
{
    init_tempering(mts);
    optimize_v(0, 0, 0);
    mts->shift0 = shift_0;
    mts->shift1 = shift_1;
    mts->shiftB = shift_s;
    mts->shiftC = shift_t;
    mts->maskB = mask_b >> ggap;
    mts->maskC = mask_c >> ggap;
}

void _get_tempering_parameter_hard_dc(mt_struct *mts)
{
    int i;
    MaskNode mn0, *cur, *next;

    init_tempering(mts);

    for (i=0; i<www; i++) 
	gcur_maxlengs[i] = -1;

    mn0.b = mn0.c = mn0.leng = 0;
    mn0.next = NULL;
	
    cur = &mn0;
    for (i=0; i<LIMIT_V_BEST_OPT; i++) {
	next = optimize_v_hard(i, cur);
	if (i > 0) 
	    delete_MaskNodes(cur);
	cur = next;
    }
    delete_MaskNodes(cur);

    optimize_v(gmax_b, gmax_c,i);
    mts->shift0 = shift_0;
    mts->shift1 = shift_1;
    mts->shiftB = shift_s;
    mts->shiftC = shift_t;
    mts->maskB = mask_b >> ggap;
    mts->maskC = mask_c >> ggap;

    /* show_distrib(mts); */
}

static void init_tempering(mt_struct *mts)
{
    int i;

    mmm = mts->mm;
    nnn = mts->nn;
    rrr = mts->rr; 
    www = mts->ww;
    shift_0 = S00;
    shift_1 = S01;
    shift_s = SSS;
    shift_t = TTT; 
    ggap = WORD_LEN - www;
    /* bits are filled in mts->aaa from MSB */
    aaa[0] = 0; aaa[1] = (mts->aaa) << ggap;


    for( i=0; i<WORD_LEN; i++)
        bitmask[i] = cur_bitmask[i] = UINT32_C(0x80000000) >> i;

	for( i=0, glower_mask=0; i<rrr; i++)
		glower_mask = (glower_mask<<1)| 0x1;
	gupper_mask = ~glower_mask;
	gupper_mask <<= ggap;
	glower_mask <<= ggap;

	greal_mask = (gupper_mask | glower_mask);

#if 0
	printf ("n=%d m=%d r=%d w=%d\n", nnn, mmm, rrr, www);
	printf ("nw-r=%d\n", nnn*www-rrr);
	printf ("a=%x(%x << %d)\n", aaa[1],mts->aaa,ggap); 
	printf ("upper (w-r) bit mask = %x\n", gupper_mask);
	printf ("lower r bit mask     = %x\n", glower_mask);
	printf ("w bit mask           = %x\n", greal_mask);
	fflush(stdout);
#endif
}

/* (v-1) bitmasks of b,c */
static MaskNode *optimize_v_hard(int v, MaskNode *prev_masks)
{
    int i, ll, t;
    uint32_t bbb[8], ccc[8];
    MaskNode *cur_masks;

    cur_masks = NULL;

    while (prev_masks != NULL) {

	ll = push_stack(prev_masks->b,prev_masks->c,v,bbb,ccc);

	for (i=0; i<ll; ++i) {
	    mask_b = bbb[i];
	    mask_c = ccc[i];
	    t = lenstra(v+1);
	    if (t >= gcur_maxlengs[v]) {
		gcur_maxlengs[v] = t;
		gmax_b = mask_b;
		gmax_c = mask_c;
		cur_masks = cons_MaskNode(cur_masks, mask_b, mask_c, t);
	    }
	}
	prev_masks = prev_masks->next;
    }

    cur_masks = delete_lower_MaskNodes(cur_masks, gcur_maxlengs[v]);

    return cur_masks;
}


/* (v-1) bitmasks of b,c */
static void optimize_v(uint32_t b, uint32_t c, int v)
{
    int i, max_len, max_i, ll, t;
    uint32_t bbb[8], ccc[8];

    ll = push_stack(b,c,v,bbb,ccc);

    max_len = max_i = 0;
    if (ll > 1) { 
	for (i=0; i<ll; ++i) {
	    mask_b = bbb[i];
	    mask_c = ccc[i];
	    t = lenstra(v+1);
	    if (t > max_len) {
		max_len = t;
		max_i = i;
	    }
	}
    }

    if ( v >= www-1 ) {
	mask_b = bbb[max_i];
	mask_c = ccc[max_i];
	return;
    }

    optimize_v(bbb[max_i], ccc[max_i], v+1);
}

static int push_stack(uint32_t b, uint32_t c, int v, uint32_t *bbb, uint32_t *ccc)
{
    int i, ll, ncv;
    uint32_t cv_buf[2];

    ll = 0;

    if( (v+shift_t) < www ){
        ncv = 2; cv_buf[0] = c | bitmask[v]; cv_buf[1] = c;
    }
    else {
        ncv = 1; cv_buf[0] = c;
    }

    for( i=0; i<ncv; ++i)
        ll += push_mask( ll, v, b, cv_buf[i], bbb, ccc);

    return ll;
}

static int push_mask(int l, int v, uint32_t b, uint32_t c, uint32_t *bbb, uint32_t *ccc)
{
    int i, j, k, nbv, nbvt;
    uint32_t bmask, bv_buf[2], bvt_buf[2];

    k = l;
    if( (shift_s+v) >= www ){
        nbv = 1; bv_buf[0] = 0;
    }
    else if( (v>=shift_t) && (c&bitmask[v-shift_t] ) ){
        nbv = 1; bv_buf[0] = b&bitmask[v];
    }
    else {
        nbv = 2; bv_buf[0] = bitmask[v]; bv_buf[1] = 0;
    }

    if( ((v+shift_t+shift_s) < www) && (c&bitmask[v]) ){
        nbvt = 2; bvt_buf[0] = bitmask[v+shift_t]; bvt_buf[1] = 0;
    }
    else {
        nbvt = 1; bvt_buf[0] = 0;
    }

    bmask = bitmask[v];
    if( (v+shift_t) < www )
        bmask |= bitmask[v+shift_t];
    bmask = ~bmask;
    for( i=0; i<nbvt; ++i){
        for( j=0; j<nbv; ++j){
            bbb[k] = (b&bmask) | bv_buf[j] | bvt_buf[i];
            ccc[k] = c;
            ++k;
        }
    }

    return k-l;
}


/**********************************/
/****  subroutines for lattice ****/
/**********************************/
static int lenstra(int v)
{
    Vector **lattice, *ltmp;
    int i, j, deg, max_row;

    upper_v_bits = 0;
    for( i=0; i<v; i++) {
        cur_bitmask[i] = bitmask[i];
        upper_v_bits |= bitmask[i];
    }

    lattice = make_lattice( v );

    i = -1; max_row=v;
    while( i<max_row ){ /* normalized until i-th row */

        pull_min_row( i+1, max_row, lattice );
        hakidasi( i+1, lattice);

        if( lattice[i+1]->bp & upper_v_bits ) {
            pull_max( i+1, v, lattice[i+1] );
            ++i;
        }
        else {
            deg = degree_of_vector( lattice[i+1]); 

            if(deg==MIN_INFINITE){
            /* if deg==MIN_INFINITE, */ 
            /* exchange (i+1)-th row and v-th row */
                ltmp = lattice[i+1]; lattice[i+1] = lattice[v];
                lattice[v] = ltmp;
                --max_row; 
            }
            else {
                for( j=i; j>=0; j--){
                    if( lattice[j]->degree <= deg ) break;
                }
                i = j;
            }
        }

    }

    i = lattice[max_row]->degree;
    free_lattice( lattice, v );
    return -i;
}



/********************************/
/** allocate momory for Vector **/
/********************************/
static Vector *new_Vector(void)
{
    Vector *v;

    v = (Vector *)malloc( sizeof( Vector ) );
    if( v == NULL ){
        printf("malloc error in \"new_Vector()\"\n");
        exit(1);
    }

    v->cf = (uint32_t *)calloc( nnn, sizeof( uint32_t ) );
    if( v->cf == NULL ){
        printf("calloc error in \"new_Vector()\"\n");
        exit(1);
    }

    v->start = 0;

    return v;
}


/************************************************/
/* frees *v which was allocated by new_Vector() */
/************************************************/
static void free_Vector( Vector *v )
{
    if( NULL != v->cf ) free( v->cf );
    if( NULL != v ) free( v );
}

static void free_lattice( Vector **lattice, int v)
{
    int i;

    for( i=0; i<=v; i++)
        free_Vector( lattice[i] );
    free( lattice );
}

/* Vector v is multiplied by t^k */
static Vector *mult( Vector *v, int k)  
{
    int i, j, jmmm; /* jmmm = j + mmm */
    uint32_t tmp;
    Vector *u;

    u = new_Vector();
    for( i=0; i<nnn; i++) u->cf[i] = v->cf[i]; /* copy v to u */
    u->start = v->start;
    u->degree = v->degree;
    u->bp = v->bp;
    u->z = v->z;

    if( k == 0 ) return u; /* if k==0, only v is copied to u  */

    j=u->start; 
    jmmm = j+mmm; /* note : jmmm < nnn+mmm < 2nnn */
    for(i=0; i<k; i++){

	if (j>=nnn) j -= nnn; /* same as j%=nnn */
	if (jmmm>=nnn) jmmm -= nnn; /* same as jmmm %= nnn */

	u->z = u->cf[j];
	tmp =  (u->cf[j]&gupper_mask) | (u->cf[(j+1)%nnn]&glower_mask) ;
	tmp = u->cf[jmmm] ^ ( (tmp>>1) ^ aaa[LSB(tmp)] );
	tmp &= greal_mask;
	u->cf[j] =  tmp;
	
	++j;
	++jmmm;
    }

    u->start += k; u->start %= nnn;

    /* integer part is shifted and bitmasked */
    tmp = u->z;
    TRNSTMP;
    MASKTMP;
    u->z = tmp;
    u->degree += k;

    return u;
}

/* adds v to u (then u will change) */
static void add( Vector *u, Vector *v) 
{
    int i, stu, stv;

    stu = u->start; stv = v->start;
    for( i=0; i<nnn; i++){

	/**  0 <= stu,stv < 2*nnn always holds          **/
	/** so, modulo nnn can be calculated as follows **/
	if (stu>=nnn) stu -= nnn; /* same as stu %= nnn  */
	if (stv>=nnn) stv -= nnn; /* same as stv %= nnn  */

	u->cf[stu] ^= v->cf[stv];
	stu++; stv++;
    }

    u->z ^=  v->z;
}

/* returns the max degree of v */
static int degree_of_vector( Vector *v)
{
    int i,j,k;
    int immm; /* immm = i + mmm */
    uint32_t tmp;
    Vector *u;

    if( v->z & upper_v_bits ){ /* if integer part is not 0 */
        v->bp = v->z;
        v->degree = 0;
        return 0;
    }

    for(k=v->start, j=0; j<nnn; j++,k++){

	if (k>=nnn) k -= nnn; /* same as k %= nnn (note 0<= k < 2*nnn) */

        tmp = v->cf[k];
	if (tmp) {
	    TRNSTMP;
	    MASKTMP;
	    if( tmp & upper_v_bits ) {
		v->bp = tmp;
		v->degree = -(j+1);
		return -(j+1);
	    }
	}
    }

    u = new_Vector(); /* copy v to u */
    for( j=0; j<nnn; j++) u->cf[j] = v->cf[j];
    u->z = v->z;
    u->start = v->start;


    k = nnn * (www-1) - rrr;
    i=u->start; 
    immm = i + mmm; /* note : immm < nnn+mmm < 2nnn */
    for(j=0; j<k; j++){ 

	  /* i = (u->start + j) % nnn */
	  if (i>=nnn) i -= nnn; /* same as i%=nnn: note always 0<=i<2*nnn */
	  if (immm>=nnn) immm -= nnn; /* same as immm %= nnn */

	  tmp = (u->cf[i]&gupper_mask) | (u->cf[(i+1)%nnn]&glower_mask);
	  tmp = u->cf[immm] ^ ( (tmp>>1) ^ aaa[LSB(tmp)] );
	  tmp  &=  greal_mask;
	  u->cf[i] = tmp;

	  if (tmp) {
	      TRNSTMP;
	      MASKTMP;
	      if( tmp & upper_v_bits ) {
		  v->bp = tmp;
		  v->degree = -(j+nnn+1);
		  free_Vector(u);
		  return -(j+nnn+1);
	      }
	  }
	  
	  ++i;
	  ++immm;
    }

    free_Vector(u);
    v->bp = 0;
    v->degree = MIN_INFINITE;
    return MIN_INFINITE; /* if 0 inspite of  (nw-r) times of generation */
}


/* add t^k*(i-th row) to j-th row */
static void add_i_to_j( Vector **lattice, int i, int j, int k)
{
    Vector *ith;

    ith = mult( lattice[i], k);
    add( lattice[j], ith);
    free_Vector( ith );
}

/* exchange columns so that i-th element of variable vec */
/* gives the norm of vec */
static void pull_max( int i, int v, Vector *vec)
{
    int j;
    uint32_t tmp;

    if( vec->bp & cur_bitmask[i] ) return;

    for( j=i+1; j<v; j++){
        if( vec->bp & cur_bitmask[j] ){
            tmp = cur_bitmask[i];
            cur_bitmask[i] = cur_bitmask[j];
            cur_bitmask[j] = tmp;
            break;
        }
    }

}

/* puts i-th row be the minimum one in i-th row ... v-th row */
static void pull_min_row( int i, int v, Vector **lattice)
{
    int j, min_deg, min_j;
    Vector *vtmp;

    min_deg = lattice[i]->degree;
    min_j = i;
    for( j=i+1; j<=v; j++){
        if( min_deg > lattice[j]->degree ){
            min_deg = lattice[j]->degree;
            min_j = j;
        }
    }

    vtmp = lattice[min_j]; lattice[min_j] = lattice[i];
    lattice[i] = vtmp;
}

/* sweeps out k-th row with 0th ... (k-1)th rows */
static void hakidasi( int k, Vector **lattice)
{
  int i;
  
  for( i=0; i<k; i++){
	if( lattice[k]->bp & cur_bitmask[i] ){
	  add_i_to_j( lattice, i, k, lattice[k]->degree - lattice[i]->degree);
	  lattice[k]->bp ^= lattice[i]->bp;
        }
    }
}

/* makes a initial lattice */
static Vector **make_lattice(int v)
{
    int i;
    uint32_t tmp;
    Vector **lattice, *top;

    lattice = (Vector **)malloc( (v+1) * sizeof( Vector *) );
    if( NULL == lattice ){
        printf("malloc error in \"make_lattice\"\n");
        exit(1);
    }

    for( i=1; i<=v; i++){ /* from 1st row to v-th row */
        lattice[i] = new_Vector();
        lattice[i]->z = bitmask[i-1];
        lattice[i]->start = 0;
        lattice[i]->bp = bitmask[i-1];
        lattice[i]->degree = 0;
    }


    top = new_Vector(); /* 0th row */
    for(i=0; i<nnn; i++) {
	  top->cf[i] = _genrand_dc();
	  top->cf[i] &= greal_mask;
	}
    tmp = ( top->cf[0] & gupper_mask ) | ( top->cf[1] & glower_mask );
    top->cf[0] = top->cf[mmm] ^ ( (tmp>>1) ^ aaa[LSB(tmp)] );
	top->cf[0] &= greal_mask;
    top->z = 0; top->start = 1; 
    degree_of_vector( top );
    lattice[0] = top;

    return lattice;
}

/***********/
static MaskNode *cons_MaskNode(MaskNode *head, uint32_t b, uint32_t c, int leng)
{
    MaskNode *t;

    t = (MaskNode*)malloc(sizeof(MaskNode));
    if (t == NULL) {
	printf("malloc error in \"cons_MaskNode\"\n");
        exit(1);
    }

    t->b = b;
    t->c = c;
    t->leng = leng;
    t->next = head;

    return t;
}

static void delete_MaskNodes(MaskNode *head)
{
    MaskNode *t;

    while(head != NULL) {
	t = head->next;
	free(head);
	head = t;
    }
}

static MaskNode *delete_lower_MaskNodes(MaskNode *head, int l)
{
    MaskNode *s, *t, *tail;

    s = head;
    while(1) { /* heading */
	if (s == NULL)
	    return NULL;
	if (s->leng >= l)
	    break;
	t = s->next;
	free(s);
	s = t;
    }

    head = tail = s;
    
    while (head != NULL) {
	t = head->next;
	if (head->leng < l) {
	    free(head);
	}
	else {
	    tail->next = head;
	    tail = head;
	}
	head = t;
    }
	
    tail->next = NULL;
    return s;
}

#if 0
static void count_MaskNodes(MaskNode *head)
{
    int c; 
    
    c = 0;
    while(head != NULL) {
	head = head->next;
	c++;
    }
    printf ("---> number of nodes = %d\n",c);
}
#endif

static void show_distrib(mt_struct *mts)
{
    int i, lim, diff, t;
    double per;

    init_tempering(mts);
	
    mask_b = (mts->maskB) << ggap;
    mask_c = (mts->maskC) << ggap;
    for (i=0; i<www; i++) {
	t = lenstra(i+1);
	lim = (nnn*www-rrr)/(i+1);
	diff = lim  - t;
	per = (double)t / (double)lim;
	printf ("%d %d %d %d %4.2f\n", i+1, t,  lim, diff, per);
    }
}
