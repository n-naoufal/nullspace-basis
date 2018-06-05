
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */

#include "slu_ddefs.h"

#define HANDLE_SIZE  8
/* kind of integer to hold a pointer.  Use int.
   This might need to be changed on 64-bit systems. */
typedef long long int fptr;  /* 32-bit by default */

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;

void
c_fortran_dgssv_(double *values, int *rowind, int *colptr,int *iopt, 
                  int *m, int *n, int *nnz, int *nrhs, int *PermC,
                  int *PermL, double *b, int *ldb,
		 fptr *f_factors, /* a handle containing the address
				     pointing to the factored matrices */
		 int *info)

{
/* 
 * This routine can be called from Fortran.
 *int *PermL,
                 int *PermC,
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs LU decomposition for the first time
 *      = 2, performs triangular solve
 *      = 3, free all the storage in the end
 *
 * f_factors (input/output) fptr* 
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 */
 
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      panel_size, permc_spec, relax;
    int *col_to_sup, *sup_to_col, *Row_Ind, *rowind_colptr;
    int nbnz, indcol, ind_diag;
    register int i, j, k, c, d, nsup;
    double       *dp;
    trans_t  trans;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    factors_t *LUfactors;
    double *at;
    int *colind;
    int *rowptr;

  /* Gestion de la transposition */

    trans = NOTRANS;

    if (*iopt < 0) {
      *iopt=-*iopt;
      trans=TRANS;
    }

  if ( *iopt == 1 ) { /* LU decomposition */

 /* Set the default input options. */
  set_default_options(&options);

	/* Initialize the statistics variables. */
	StatInit(&stat);


	/* Adjust to 0-based indexing */
	for (i = 0; i < *nnz; ++i) --rowind[i];
	for (i = 0; i <= *n; ++i) --colptr[i];



	dCreate_CompCol_Matrix(&A, *m, *n, *nnz, values, rowind, colptr,
			       SLU_NC, SLU_D, SLU_GE);



	L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	if ( !(perm_r = intMalloc(*m)) ) ABORT("Malloc fails for perm_r[].");
	 if ( *nrhs < 4 ) {
  if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");
	}
  if ( !(etree = intMalloc(*m)) ) ABORT("Malloc fails for etree[].");

//  dPrint_CompCol_Matrix("A", &A);
  
	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree on structure of A'*A
	 *   permc_spec = 2: minimum degree on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
    if ( *nrhs < 4 ) {
    /* permc_spec = options.ColPerm; */
      permc_spec = *nrhs;                      
    get_perm_c(permc_spec, &A, perm_c);
    sp_preorder(&options, &A, perm_c, etree, &AC);
      }


	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
   
	dgstrf(&options, &AC, relax, panel_size, etree,
                NULL, 0, perm_c, perm_r, L, U, &stat, info);

    printf(" \n");
  //  printf("La factorisation LU retourne INFO= %d\n", *info);
	if ( *info <= n ) { 
             printf(" \n");
             printf(" *** Factorization LU *** \n");
            
	     dQuerySpace(L, U, &mem_usage);
	     printf(" - taille mémoire L\\U MB %.3f \n",mem_usage.for_lu/1e6); 
      printf(" - total MB needed %.3f \n", mem_usage.total_needed/1e6);
 
	}

        /*
        printf(" Recuperation des infos - Matrice L \n");
        printf("Stype %d, Dtype %d, Mtype %d\n", L->Stype,L->Dtype,L->Mtype);
        printf("  \n");
        printf(" Taille de la matrice : %d x %d\n",L->nrow,L->ncol);
        printf(" nb. valeurs non nulles :%d\n",Lstore->nnz);
        */

	Lstore = (SCformat *) L->Store;
	Ustore = (NCformat *) U->Store;

  dp = (double *) Lstore->nzval;
  col_to_sup = Lstore->col_to_sup;
  sup_to_col = Lstore->sup_to_col;
  rowind_colptr = Lstore->rowind_colptr;
  Row_Ind = Lstore->rowind;
  nbnz=0;

  /* Comptage pour utilisation exterieure */

  for (k = 0; k <= Lstore->nsuper; ++k) {
    c = sup_to_col[k];
    nsup = sup_to_col[k+1] - c;
    for (j = c; j < c + nsup; ++j) {
d = Lstore->nzval_colptr[j];
      nbnz=nbnz+rowind_colptr[c+1]-rowind_colptr[c];
    }
  }

  *nrhs= (long long) nbnz;
  /*
  printf(" Nombre de valeurs non nulles trouve : %d \n",nbnz);
  dPrint_SuperNode_Matrix("L", L);    
  */    
        
/*	dPrint_SuperNode_Matrix("L", L); */
  
	/* Restore to 1-based indexing */
	for (i = 0; i < *nnz; ++i) ++rowind[i];
	for (i = 0; i <= *n; ++i) ++colptr[i];

	/* Save the LU factors in the factors handle */
	LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
	LUfactors->L = L;
	LUfactors->U = U;
	LUfactors->perm_c = perm_c;
	LUfactors->perm_r = perm_r;
	*f_factors = (fptr) LUfactors;

	/* Free un-wanted storage */
	SUPERLU_FREE(etree);
	Destroy_SuperMatrix_Store(&A);
	Destroy_CompCol_Permuted(&AC);
	StatFree(&stat);

} else if ( *iopt == 2 ) { /* Triangular solve */
  /* Initialize the statistics variables. */
  StatInit(&stat);

  /* Extract the LU factors in the factors handle */
  LUfactors = (factors_t*) *f_factors;
  L = LUfactors->L;
  U = LUfactors->U;
  perm_c = LUfactors->perm_c;
  perm_r = LUfactors->perm_r;
  


  double *Bmat = (double*) ((DNformat*) B.Store)->nzval; 

  dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_D, SLU_GE);

 //  for (i = 0; i < *ldb; ++i) printf("%f\n",b[i] ); 

  /* Solve the system A*X=B, overwriting B with X. */
  dgstrs (trans, L, U, perm_c, perm_r, &B, &stat, info);

  // for (i = 0; i < *ldb; ++i) printf("%f\n",b[i] ); 

  Destroy_SuperMatrix_Store(&B);
  StatFree(&stat);

	} else if ( *iopt == 4 ) { /* Copie les valeurs vers le fortran */
	/* Initialize the statistics variables. */
	
   StatInit(&stat);

	/* Extract the LU factors in the factors handle */
	LUfactors = (factors_t*) *f_factors;
	L = LUfactors->L;
        Lstore = (SCformat *) L->Store;
	perm_c = LUfactors->perm_c;
	perm_r = LUfactors->perm_r;

  dp = (double *) Lstore->nzval;
  col_to_sup = Lstore->col_to_sup;
  sup_to_col = Lstore->sup_to_col;
  rowind_colptr = Lstore->rowind_colptr;
  Row_Ind = Lstore->rowind;
  indcol=0;
  ind_diag=0;
  colptr[indcol]=1;
  /* Recopie de la matrice */
  for (k = 0; k <= Lstore->nsuper; ++k) {
    c = sup_to_col[k];
    nsup = sup_to_col[k+1] - c;
    for (j = c; j < c + nsup; ++j) {
	      d = Lstore->nzval_colptr[j];
        indcol++;
        colptr[indcol]=colptr[indcol-1]+rowind_colptr[c+1]-rowind_colptr[c];
        for (i = rowind_colptr[c]; i < rowind_colptr[c+1]; ++i) {
          values[d]=dp[d];
          rowind[d]=Row_Ind[i]+1;
          //printf("%d\t%d\t%e\n", Row_Ind[i], j, dp[d++]); 
          if ( Row_Ind[i] == j) {
            //printf("%d\n",123456789);
            //printf("%d\t%d\t%e\n", Row_Ind[i], j, dp[d]); 
            b[ind_diag]=dp[d];
            //printf("%d\t%d\t%e\n", ind_diag, d, b[ind_diag]); 
            ind_diag++;
          }
          d++;     
        }
    }
  }
    
    for ( i = 0; i < L->nrow ; i++) {
      PermL[i]=++perm_r[i];
    }
    
    for ( i = 0; i < L->ncol ; i++) {
      PermC[i]=++perm_c[i];
    } 
    

         printf(" \n");
         printf(" *** Extraction des facteurs LU  *** ");
         printf(" \n");
StatFree(&stat);
		

 } else if ( *iopt == 5) { /* Transpose une matrice  */
  /* Initialize the statistics variables. */
  StatInit(&stat);
        
        /* Adjust to 0-based indexing */
  for (i = 0; i < *nnz; ++i) --rowind[i];

  for (i = 0; i <= *n; ++i) --colptr[i]; 

        at=doubleMalloc(*nnz);
        colind=intMalloc(*nnz);
        rowptr=intMalloc(*m+1);
        
        dCompRow_to_CompCol(*n, *m, *nnz, 
      values,rowind, colptr,
      &at, &colind, &rowptr);
        for (i = 0; i < *nnz; ++i) {printf("%d\n",colind[i]);} 
        for (i = 0; i < *nnz ; i++) {
          values[i]=at[i];
          rowind[i]=colind[i];
          
        }
        for (i = 0; i < *m+1 ; i++) {
          colptr[i]=rowptr[i];
        }
        
        /* Back to 1-based indexing */
  for (i = 0; i < *nnz; ++i) {++rowind[i];printf("%d\n",rowind[i]);}
  for (i = 0; i <= *m; ++i) ++colptr[i];
            

        SUPERLU_FREE(at);
        SUPERLU_FREE(colind);
        SUPERLU_FREE(rowptr);
        
  StatFree(&stat);



    } else if ( *iopt == 3 ) { /* Free storage */
	/* Free the LU factors in the factors handle */

	LUfactors = (factors_t*) *f_factors;
	SUPERLU_FREE (LUfactors->perm_r);
	SUPERLU_FREE (LUfactors->perm_c);
	Destroy_SuperNode_Matrix(LUfactors->L);
	Destroy_CompCol_Matrix(LUfactors->U);
  SUPERLU_FREE (LUfactors->L);
  SUPERLU_FREE (LUfactors->U);
	SUPERLU_FREE (LUfactors);
    printf(" *** Nettoyage des facteurs LU *** ");
  printf(" \n");

    } else {
	fprintf(stderr,"Invalid iopt=%d passed to c_fortran_dgssv()\n",*iopt);
	exit(-1);
    }
}

