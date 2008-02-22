#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include "polynomial.h"
#include "polyDiags.h"
//#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

/* Variables for tracking internal polynomial memory usage. */
unsigned long maxHashListLength = HASH_TABLE_INCREASE;
unsigned long constantPs = 0, constantHashPs = 0, variablePs = 0, variableHashPs = 0, sumPs = 0, 
  sumHashPs = 0, productPs = 0, productHashPs = 0, functionPs = 0, functionHashPs = 0;
unsigned long peakConstantPs = 0, peakVariablePs = 0, peakSumPs = 0, peakProductPs = 0, peakFunctionPs = 0;
unsigned long constantPLExpansions = 0, variablePLExpansions = 0, sumPCollectExpansions = 0,
  sumPTermMergeExpansions = 0, sumPListExpansions = 0, productPCollectExpansions = 0, 
  productPTermMergeExpansions = 0, productPListExpansions = 0;
unsigned long constantPsSize = 0, variablePsSize = 0, variablePsExpSize = 0, sumPsSize = 0, sumPColExpSize = 0, 
  sumPTrmMrgExpSize = 0, productPsSize = 0, productPColExpSize = 0, productPTrmMrgExpSize = 0;

/* Display the internal polynomial memory usage statistics. */
void dumpPStats(char *stateDescription) {
  fprintf(stderr, "State: %s...\n", stateDescription);
  fprintf(stderr, "Max hash list len is initially %d, now %lu (list growth allocs not reported)\n",
	  HASH_TABLE_INCREASE, maxHashListLength);
  fprintf(stderr, "Constant poly count:%lu, hits:%lu, peak:%lu, list expansions (@%d):%lu\n",
	  constantPs, constantHashPs, peakConstantPs, CONSTANT_LIST_INCREASE, constantPLExpansions);
  fprintf(stderr, "Variable poly count:%lu, hits:%lu, peak:%lu, list expansions (@%d):%lu\n",
	  variablePs, variableHashPs, peakVariablePs, VARIABLE_LIST_INCREASE, variablePLExpansions);
  fprintf(stderr, "Sum poly count:%lu, hits:%lu, ratio:%d-to-1, peak:%lu, size:%lu\n", 
	  sumPs, sumHashPs, (int)(sumHashPs/sumPs), peakSumPs, sumPsSize);
  fprintf(stderr,
	  "...term coll realloc sz:%lu, term merge expansions:%lu, realloc sz:%lu, list expansions (@%d):%lu\n",
	 sumPColExpSize, sumPTermMergeExpansions, sumPTrmMrgExpSize, SUM_LIST_INCREASE, sumPListExpansions);
  fprintf(stderr, "Product poly count:%lu, hits:%lu, ratio:%d-to-1, peak:%lu, sz:%lu\n",
	 productPs, productHashPs, (int)(productHashPs/productPs), peakProductPs, productPsSize);
  fprintf(stderr,
	  "...term coll realloc sz:%lu, term merge expansions:%lu, realloc sz:%lu, list expansions (@%d):%lu\n",
	  productPColExpSize, productPTermMergeExpansions, productPTrmMrgExpSize, PRODUCT_LIST_INCREASE, 
	  productPListExpansions);
}

#define MAXPOLYTIERS 13
int polyTiers[MAXPOLYTIERS][5];
int peakPolyTiers;
char *polyTypes[] = {"constant", "variable", "sum", "product", "function"};
/* It might not be a sumPoly, but we can treat it as one. */
void traversePoly(struct polynomial *p, int currentTier)
{
  struct sumPoly *pS;
  struct functionPoly *pF;
  struct polynomial *q;
  int i;
  if (currentTier > peakPolyTiers) peakPolyTiers = currentTier;
  if (currentTier >= MAXPOLYTIERS) return;
  polyTiers[currentTier][p->eType]++;
  switch(p->eType){
  case T_CONSTANT:
  case T_VARIABLE:
    break;
  case T_SUM:
  case T_PRODUCT:
    pS = p->e.s;
    for (i=0; i<pS->num; i++) {
      q = pS->sum[i];
      traversePoly(q, currentTier+1);
    }
    break;
  case T_FUNCTIONCALL:
    pF = p->e.f;
    for (i=0; i<pF->paraNum; i++) {
      q = pF->para[i];
      traversePoly(q, currentTier+1);
    }
    break;
  default:
    fprintf(stderr,"Error, Unknown expression type: [%d], exiting!\n", p->eType);
    exit(1);
  }
}
void printSummaryPoly(struct polynomial *p)
{
  int i, j;
  printf("Summary of Polynomial:\n");
  memset(polyTiers, 0, sizeof(polyTiers));
  peakPolyTiers = 0;
  traversePoly(p, 0);
  for (i=0; i<=peakPolyTiers; i++) {
    printf("Tier %d:",i);
    int firstPrint = 1;
    for (j=0; j<5; j++) {
      if (polyTiers[i][j] != 0) {
	if (!firstPrint) printf(",");
	firstPrint = 0;
	printf(" %d %s", polyTiers[i][j], polyTypes[j]);
	if (polyTiers[i][j] > 1) printf("s");
      }
    }
    printf("\n");
  }
  printf("---\n");
}

void expPrint(struct polynomial *p) {
  int i;

  switch(p->eType){
  case T_CONSTANT:
    fprintf(stderr, "%f", p->value);
    break;
  case T_VARIABLE:
    fprintf(stderr,"%s",p->e.v->vName);
    break;
  case T_SUM:
    if(p->e.s->num>1) fprintf(stderr,"(");
    if (p->e.s->factor[0] != 1)
      fprintf(stderr, "%f*", p->e.s->factor[0]);
    expPrinting(p->e.s->sum[0]);
    for(i=0;i<p->e.s->num;i++) {
      fprintf(stderr, "+");
      if (p->e.s->factor[0] != 1)
	fprintf(stderr, "%f*", p->e.s->factor[i]);
      expPrinting(p->e.s->sum[i]);
    }
    if(p->e.s->num>1) fprintf(stderr,")");
    break;
  case T_PRODUCT:
    expPrinting(p->e.p->product[0]);
    if(p->e.p->exponent[0] != 1)
      fprintf(stderr,"^%d",p->e.p->exponent[0]);
    for(i=1;i<p->e.s->num;i++) {
      fprintf(stderr, "*");
      expPrinting(p->e.p->product[i]);
      if(p->e.p->exponent[i] != 1)
	fprintf(stderr,"^%d",p->e.p->exponent[i]);
    }
    break;
  case T_FUNCTIONCALL:
    fprintf(stderr,"%s(",p->e.f->name);
    for(i=0;i<p->e.f->paraNum-1;i++) {
      expPrinting(p->e.f->para[i]);
      fprintf(stderr,", ");
    }
    expPrinting(p->e.f->para[p->e.f->paraNum-1]);
    fprintf(stderr,")");
    break;
  default:
    fprintf(stderr,"Unknown expression type %d exit!!!!!, exit(5)\n",p->eType);
    exit(1);
  }
}

void printAllPolynomials()
{
   int i,j;
   if(constantCount>0)
   {
           fprintf(stderr,"All %d constants:\n",constantCount );
           for(i=0;i<CONSTANT_HASH_SIZE;i++)
           {
              if(constantHash[i].num<=0)
                 continue;
              for(j=0;j<constantHash[i].num;j++)
              {
                fprintf(stderr,"(%d  %d) index=%d key=%d valid=%d constant: ",i,j,constantHash[i].index[j],constantHash[i].key[j], constantList[constantHash[i].index[j]]->valid);
                expPrinting(constantList[constantHash[i].index[j]]);
                fprintf(stderr,"\n");
              }
              fprintf(stderr,"\n");
           }
           fprintf(stderr,"\n");
   }
   if(variableCount>0)
   {
           fprintf(stderr,"All %d variables:\n",variableCount);
           for(i=0;i<VARIABLE_HASH_SIZE;i++)
           {
              if(variableHash[i].num<=0)
                 continue;
              for(j=0;j<variableHash[i].num;j++)
              {
                fprintf(stderr,"(%d  %d) index=%d key=%d variable: ",i,j,variableHash[i].index[j],variableHash[i].key[j]);
                expPrinting(variableList[variableHash[i].index[j]]);
                fprintf(stderr,"\n");
              }
           }
           fprintf(stderr,"\n");
   }
   if(sumCount>0)
   {
           fprintf(stderr,"All %d sums:\n",sumCount);
           for(i=0;i<SUM_HASH_SIZE;i++)
           {
              if(sumHash[i].num<=0)
                 continue;
              for(j=0;j<sumHash[i].num;j++)
              {
                fprintf(stderr,"(%d  %d) index=%d Key=%d sum: ",i,j,sumHash[i].index[j],sumHash[i].key[j]);
                expPrinting(sumList[sumHash[i].index[j]]);
                fprintf(stderr,"\n");
              }
              fprintf(stderr,"\n");
           }
           fprintf(stderr,"\n");
   }
   if(productCount>0)
   {
           fprintf(stderr,"All %d products:\n",productCount);
           for(i=0;i<PRODUCT_HASH_SIZE;i++)
           {
              if(productHash[i].num<=0)
                 continue;
              for(j=0;j<productHash[i].num;j++)
              {
                fprintf(stderr,"(%d  %d) index=%d Key=%d product: ",i,j,productHash[i].index[j],productHash[i].key[j]);
                expPrinting(productList[productHash[i].index[j]]);
                fprintf(stderr,"\n");
              }
              fprintf(stderr,"\n");
           }
           fprintf(stderr,"\n");
   }

   if(functionCallCount>0)
   {
        fprintf(stderr,"All %d function calls:\n",functionCallCount);
           for(i=0;i<FUNCTIONCALL_HASH_SIZE;i++)
           {
              if(functionCallHash[i].num<=0)
                 continue;
              for(j=0;j<functionCallHash[i].num;j++)
              {
                fprintf(stderr,"(%d  %d) index=%d Key=%d functionCall: ",i,j,functionCallHash[i].index[j],functionCallHash[i].key[j]);
                expPrinting(functionCallList[functionCallHash[i].index[j]]);
                fprintf(stderr,"\n");
              }
              fprintf(stderr,"\n");
           }
           fprintf(stderr,"\n");
   }
}
