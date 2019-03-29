/**********************************************************************
 * Copyright 2007, Columbus Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include "polynomial.h"
//#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

//Recursively evaluate a polynomial
//This version of polynomail evaluation doesn't use
//polynomial sorting list.  It just compute the values
// of each sub polynomial recursively.
double evaluateValue(struct polynomial *p)
{
   int i;
   double result;
   struct sumPoly *sP;
   struct productPoly *pP;
   struct functionPoly *fp;
   double value0,value1;

   switch(p->eType){
     //If a sub polynomial is a contant, return the value
     //of this constant
     case T_CONSTANT:
         return p->value;
     //If a sub polynomial is a variable, return the value
     //of this variable         
     case T_VARIABLE:
         if(p->e.v->vType=='D')
         {
            p->value= *(p->e.v->vAddrD);
            return *(p->e.v->vAddrD);
         }
         else if(p->e.v->vType=='I')
         {
            p->value= *(p->e.v->vAddrI);
            return *(p->e.v->vAddrI);
         }
         else
         {
            fprintf(stderr,"Wrong variable type, exit!\n");
            exit(1);
         }
     
     //If a sub polynomial is a sum, evaluate the values of all the polynomials
     //in the items of this sum and then the sum of the items
     case T_SUM:
         result=0;
         sP=p->e.s;
         for(i=0;i<sP->num;i++)
           result += evaluateValue(sP->sum[i])*sP->factor[i];
         p->value=result;
         return result;

     //If a sub polynomial is a product, evaluate the values of all the polynomials
     //in the items of this product and then the product of items
     case T_PRODUCT:
         result=1;
         pP=p->e.p;
         for(i=0;i<pP->num;i++)
         {
           result *= pow(evaluateValue(pP->product[i]),pP->exponent[i]);

         }
         p->value=result;
         return result;

     //If a sub polynomial is a function call, evaluate the values of all the parameters
     //and then the function
     case T_FUNCTIONCALL:
           fp=p->e.f;
           if(strcmp(fp->name,"log10")==0)
              result=log10(evaluateValue(fp->para[0]));
           else if(strcmp(fp->name,"gsl_ran_tdist_pdf")==0)
           {
              value0=evaluateValue(fp->para[0]);
              value1=evaluateValue(fp->para[1]);
              result=gsl_ran_tdist_pdf(value0, value1);
           }
           else if(strcmp(fp->name,"gsl_cdf_tdist_Q")==0)
           {
              value0=evaluateValue(fp->para[0]);
              value1=evaluateValue(fp->para[1]);
              result=gsl_cdf_tdist_Q(value0, value1);
           }
           else if(strcmp(fp->name,"gsl_cdf_tdist_P")==0)
           {
              value0=evaluateValue(fp->para[0]);
              value1=evaluateValue(fp->para[1]);
              result=gsl_cdf_tdist_P(value0, value1);
           }
           else if(strcmp(fp->name,"gsl_ran_ugaussian_pdf")==0)
           {
              value0=evaluateValue(fp->para[0]);
              result=gsl_ran_ugaussian_pdf(value0);
           }
           else if(strcmp(fp->name,"gsl_cdf_ugaussian_Q")==0)
           {
              value0=evaluateValue(fp->para[0]);
              result=gsl_cdf_ugaussian_Q(value0);
           }
           else if(strcmp(fp->name,"gsl_cdf_ugaussian_P")==0)
           {
              value0=evaluateValue(fp->para[0]);
              result=gsl_cdf_ugaussian_P(value0);
           }
           else if(strcmp(fp->name,"gsl_cdf_chisq_P")==0)
           {
              value0=evaluateValue(fp->para[0]);
              value1=evaluateValue(fp->para[1]);
              result=gsl_cdf_chisq_P(value0, value1);
           }
           else if(strcmp(fp->name,"gsl_cdf_chisq_Q")==0)
           {
              value0=evaluateValue(fp->para[0]);
              value1=evaluateValue(fp->para[1]);
              result=gsl_cdf_chisq_Q(value0, value1);
           }
           else if(strcmp(fp->name,"gsl_ran_chisq_pdf")==0)
           {
              value0=evaluateValue(fp->para[0]);
              value1=evaluateValue(fp->para[1]);
              result=gsl_ran_chisq_pdf(value0, value1);
           }
           else if(strcmp(fp->name,"pow")==0)
           {
              value0=evaluateValue(fp->para[0]);
              value1=evaluateValue(fp->para[1]);
              result=pow(value0,value1);
                        
           }
           else if(strcmp(fp->name,"exp")==0)
           {
//              fprintf(stderr,"exp parameter:\n");
//              expPrinting(fp->para[0]);
//              fprintf(stderr,"  exp(%f)=%f\n",evaluateValue(fp->para[0]),exp(evaluateValue(fp->para[0])));
              result=exp(evaluateValue(fp->para[0]));
            
           }
           else if(strcmp(fp->name,"sqrt")==0)
           {
              result=sqrt(fp->para[0]->value);
           }
           else
           {
              fprintf(stderr,"unknown function name %s in polynomials\n",fp->name);
              exit(1);
           }
           p->value=result;
           return result;

     default:
         fprintf(stderr,"Error, Unknown expression type!!!!, exit(1)");
         exit(1);
    }
};

//Count the numbers of polynomials belonging to different categories
//in a polynomial sorting list.  This function gives us a quantitative
//measure of computational complexity of a polynomial
void countPoly(struct polyList *l, int *cCounter,
               int *vCounter,int *sCounter,int *pCounter, int *fCounter)
{
   int j;
   struct polynomial *p;
   if(l->listNext==0)
   {
       return;
   }
   for(j=l->listNext-1;j>=0;j--)
   {
      p=l->pList[j];
      switch(p->eType){
        case T_CONSTANT:
             (*cCounter)++;
        case T_VARIABLE:
             (*vCounter)++;
             break;
        case T_SUM:
           (*sCounter)++;
           break;
        case T_PRODUCT:
           (*pCounter)++;
           break;
        case T_FUNCTIONCALL:
           (*fCounter)++;
           break;
        default:
           fprintf(stderr,"Error, unknown polynomial type!!! exit(2)");
           exit(1);
           break;
     }
   }

};

//A polynomial is zero if it is a constant polynomial and its value is 0
int isZeroExp(struct polynomial *p)
{
    if(p->eType!=T_CONSTANT)
        return 0;
    if(p->value==0.0)
        return 1;
    else
        return 0;
};

//A polynomial is one if it is a constant polynomial and its value is 1
int isOneExp(struct polynomial *p)
{
//    fprintf(stderr,"In isOneoExp,p->eType=%d p->value=%f index=%d\n",
//                    p->eType,p->value,p->indexd);
    if(p->eType!=T_CONSTANT)
        return 0;
    if(p->value==1.0)
        return 1;
    else
        return 0;
};

////A polynomial is minus one if it is a constant polynomial and its value is -1
int isMinusOneExp(struct polynomial *p)
{
    if(p->eType!=T_CONSTANT)
        return 0;
    if(p->value==-1.0)
        return 1;
    else
        return 0;
};

//All the constant polynomials
void printConstantList()
{
  int i;
  for(i=0;i<constantCount;i++)
   fprintf(stderr,"Constant No. %d value=%f\n",i,constantList[i]->value);
}

//Construct a constant polynomial
struct polynomial *constantExp(double con)
{
  int i;
  struct polynomial *p;
  int key;
  int hIndex,cIndex;
  int binaryStart, binaryMiddle, binaryEnd;
  int tempI1,tempI2;
  double tempD1, tempD2;


  tempI1 = tempI2 = -1;

  
  //Compute a key for the constant  
  tempD1=frexp(con,&tempI1);

  //The key is a number between 0 and MAX_POLYNOMIAL_KEY
  if(tempD1>0)
     key=floor((tempD1-0.5)*2*MAX_POLYNOMIAL_KEY)+tempI1;
  else
     key=floor((tempD1+0.5)*2*MAX_POLYNOMIAL_KEY)+tempI1;

  //Make sure the key is positive and between 0 and MAX_POLYNOMIAL_KEY
  if(key>=0)
       key=key % MAX_POLYNOMIAL_KEY;
  else
       key=(key + MAX_POLYNOMIAL_KEY) % MAX_POLYNOMIAL_KEY;

  if(key<0)
     key=0-key;


  //Compute the index of this constant polynomial in the hash table of constant polynomials
  hIndex=key % CONSTANT_HASH_SIZE;

//   fprintf(stderr,"       key=%d hIndex=%d con=%f tempD1=%f tempI1=%d Hash num=%d\n",
//                   key,hIndex,con,tempD1,tempI1,constantHash[hIndex].num);

//   if(debugInfo==1)
//   {
//      for(i=0;i<constantHash[hIndex].num;i++)
//         fprintf(stderr,"       %5d  %10.5f  %8d\n",
//                                constantHash[hIndex].index[i],
//                                constantList[constantHash[hIndex].index[i]]->value,
//                                constantHash[hIndex].key[i]);
//   }



      //If the hash table item is not empty, determine if the constant has been in the 
      //constant polynomial list.  If it is not, determine a position in the hash table to
      //save the key and the index in the constant list of this constant polynomial and 
      //save this constant polynomial in the constant polynomial list
      if(constantHash[hIndex].num>0)
      {
          //Check the location of the key of a new constant polynomial in a list of keys
          //using binary search
          binaryStart = 0;
          binaryEnd   = constantHash[hIndex].num-1;
          while(binaryStart<=binaryEnd)
          {
             binaryMiddle=floor((binaryStart+binaryEnd)/2.0);
             if(key==constantHash[hIndex].key[binaryMiddle])
             {
                 binaryStart = binaryMiddle;
                 binaryEnd   = binaryMiddle;
                 break;
             }
             else if(key>constantHash[hIndex].key[binaryMiddle])
                 binaryStart = binaryMiddle+1;
             else
                 binaryEnd   = binaryMiddle-1;
          } //end of while


          
//          if(debugInfo==1)
//          fprintf(stderr,"binaryStart=%d binaryMiddle=%d binaryEnd=%d constantHash[hIndex].key[binaryStart]=%d\n",
//                     binaryStart,binaryMiddle,binaryEnd,constantHash[hIndex].key[binaryStart]);

          //If the newly generated key is equal to the key of a polynomial that stays in
          //the same indexed position at the hash table
          if(binaryStart>=0 && binaryStart<=constantHash[hIndex].num-1 
                            && key==constantHash[hIndex].key[binaryStart])
          {
              binaryMiddle = binaryStart;
              binaryEnd    = binaryStart+1;
              binaryStart  = binaryStart-1;
              //Search for more keys that are equal to the newly generated key
              while(binaryStart>=0 && 
                   constantHash[hIndex].key[binaryStart]==constantHash[hIndex].key[binaryMiddle])
                  binaryStart--;
              //Search for more keys that are equal to the newly generated key
              while(binaryEnd<constantHash[hIndex].num && 
                   constantHash[hIndex].key[binaryEnd]==constantHash[hIndex].key[binaryMiddle])
                  binaryEnd++;
              if(binaryStart<0)
                  binaryStart++;
              if(binaryEnd>=constantHash[hIndex].num || key!=constantHash[hIndex].key[binaryEnd])
                  binaryEnd--;

              //if the key of this constant is equal to the keys of some polynomials in the 
              //constant list, compare if the value of the new constant is equal to
              //the value of an existing constant
             for(i=binaryStart;i<=binaryEnd;i++)
             {
	         cIndex=constantHash[hIndex].index[i];
                 tempD2=frexp(constantList[cIndex]->value,&tempI2);
                 //Compare if the two constants are the same
	         if(tempI2==tempI1 && (int)(tempD2*100000000)==(int)(tempD1*100000000))
	         {
	                 //if the two constants are the same, return the constant in the constant list 
	                 return constantList[cIndex];
	         } 
	      } 
	  } 
          //If the newly generated key is unique, we are sure that this constant appears
          //the first time.
          else
          {
              binaryEnd = binaryStart;
          }
       }
       //If the hash table item is empty, we are sure that this constant appears
       //the first time
       else
       {
              binaryEnd = 0;
       }


//  fprintf(stderr,"binaryEnd=%d\n",binaryEnd);

  //next, insert it into the constant list


  //Generate a constant polynomial
  p=(struct polynomial *)malloc(sizeof(struct polynomial));
  if(p==NULL)
  {
    fprintf(stderr,"Memory allocation error: malloc returned NULL!");
    exit(1);
  }
  p->eType = T_CONSTANT;
  p->value = con; 

  //check if the constant polynomial list is full.  Apply for more items if it is full
  if(constantCount>=constantListLength)
  {
         constantListLength+=1000;
         constantList=realloc(constantList,constantListLength*sizeof(struct polynomial *));
  }



  //save the constant in the constant polynomial list
  constantList[constantCount]=p;
  p->index=constantCount;
  p->id=nodeId;
  constantCount++;
  nodeId++;
  p->key=key;
  p->valid=0;

      //increase the size of the hash table item to record the newlt generated constant
      constantHash[hIndex].num++;
      constantHash[hIndex].key=realloc(constantHash[hIndex].key,
                             sizeof(int)*constantHash[hIndex].num);
      constantHash[hIndex].index=realloc(constantHash[hIndex].index,
                             sizeof(int)*constantHash[hIndex].num);
      if( constantHash[hIndex].key==NULL || constantHash[hIndex].index==NULL)
      {
           fprintf(stderr,"Memory allocation for constants' hash table failed!\n");
           exit(1);
      }

//      fprintf(stderr,"constantHash[hIndex].num=%d\n",constantHash[hIndex].num);




      //Prepare a place in hash table for the new polynomial
      if(binaryEnd<=constantHash[hIndex].num-2)
      for(i=constantHash[hIndex].num-2;i>=binaryEnd;i--)
      {
         constantHash[hIndex].key[i+1]=constantHash[hIndex].key[i];
         constantHash[hIndex].index[i+1]=constantHash[hIndex].index[i];
      }

      //Insert the new polynomiual in the hash table
      constantHash[hIndex].key[binaryEnd]=key;
      constantHash[hIndex].index[binaryEnd]=constantCount-1;


  return p;
};



struct polynomial *variableExp(double *vD, int *vI, char vType, char name[10])
{
  int i;
  struct polynomial *p;
  struct variablePoly *vPoly;
  int key;
  int hIndex,vIndex;
  int binaryStart, binaryMiddle, binaryEnd;

  
  //compuete a key for this variable polynomial
  if(vType=='I')
     key= (long int)vI;
  else if(vType=='D')
     key= (long int)vD+1;
  else
  {
      fprintf(stderr,"UNKNOWN variable type !");
      exit(1);
  }

  //Make sure the key is positive and between 0 and MAX_POLYNOMIAL_KEY
  if(key>=0)
       key=key % MAX_POLYNOMIAL_KEY;
  else
       key=(key + MAX_POLYNOMIAL_KEY) % MAX_POLYNOMIAL_KEY;

  if(key<0)
     key=0-key;
  //Compute the index of this variable polynomial in the hash table of variable polynomials
 

  hIndex=key % VARIABLE_HASH_SIZE;

      //If the hash table item is not empty, determine if the variable has been in the 
      //variable polynomial list.  If it is not, determine a position in the hash table to
      //save the key and the index in the variable list of this variable polynomial and 
      //save this variable polynomial in the variable polynomial list
      if(variableHash[hIndex].num>0)
      {
          //Check the location of the key of a new variable polynomial in a list of keys
          //using binary search
          binaryStart = 0;
          binaryEnd   = variableHash[hIndex].num-1;
          while(binaryStart<=binaryEnd)
          {
             binaryMiddle=floor((binaryStart+binaryEnd)/2.0);
             if(key==variableHash[hIndex].key[binaryMiddle])
             {
                 binaryStart = binaryMiddle;
                 binaryEnd   = binaryMiddle;
                 break;
             }
             else if(key>variableHash[hIndex].key[binaryMiddle])
                 binaryStart = binaryMiddle+1;
             else
                 binaryEnd   = binaryMiddle-1;
          } //end of while

          //If the newly generated key is equal to the key of a polynomial that stays in
          //the same indexed position at the hash table
          if(binaryStart>=0 && binaryStart<=variableHash[hIndex].num-1 
                            && key==variableHash[hIndex].key[binaryStart])
          {
              binaryMiddle = binaryStart;
              binaryEnd    = binaryStart+1;
              binaryStart  = binaryStart-1;
              //Search for more keys that are equal to the newly generated key
              while(binaryStart>=0 && 
                   variableHash[hIndex].key[binaryStart]==variableHash[hIndex].key[binaryMiddle])
                  binaryStart--;
              //Search for more keys that are equal to the newly generated key
              while(binaryEnd<variableHash[hIndex].num && 
                   variableHash[hIndex].key[binaryEnd]==variableHash[hIndex].key[binaryMiddle])
                  binaryEnd++;
              if(binaryStart<0)
                  binaryStart++;
              if(binaryEnd>=variableHash[hIndex].num || key!=variableHash[hIndex].key[binaryEnd])
                  binaryEnd--;

              //if the key of this variable is equal to the keys of some polynomials in the 
              //variable list, compare if the variable has already been there
             for(i=binaryStart;i<=binaryEnd;i++)
             {
	         vIndex=variableHash[hIndex].index[i];
                 //Compare if the two variables are the same
	         if((vType=='D' && variableList[vIndex]->e.v->vAddrD==vD) || (vType=='I' && variableList[vIndex]->e.v->vAddrI==vI))
	         {
	                 //if the two variables are the same, return the variable in the variable list 
	                 return variableList[vIndex];
	         } 
	      } 
	  } 
          //If the newly generated key is unique, we are sure that this variable appears
          //the first time.
          else
          {
              binaryEnd = binaryStart;
          }
       }
       //If the hash table item is empty, we are sure that this variable appears
       //the first time
       else
       {
              binaryEnd = 0;
       }

//  fprintf(stderr,"binaryEnd=%d\n",binaryEnd);

  //next, insert it into the variable list

  p=(struct polynomial *)malloc(sizeof(struct polynomial));
  vPoly=(struct variablePoly *)malloc(sizeof(struct variablePoly));
  if(p==NULL || vPoly==NULL)
    fprintf(stderr,"Memory allocation error: malloc returned NULL!");
  p->eType = T_VARIABLE;
  strcpy(vPoly->vName,name);
  if(vType=='D')
     vPoly->vAddrD = vD;
  else
     vPoly->vAddrI = vI;
  vPoly->vType = vType;
  p->e.v   = vPoly;

  if(variableCount>=variableListLength)
  {
         variableListLength+=50;
         variableList=realloc(variableList,variableListLength*sizeof(struct polynomial *));
  }
  if(variableList==NULL)
  {
         fprintf(stderr,"Memory allocation error in variableExp(), exit!");
         exit(1);
  }

  p->index=variableCount;
  p->id=nodeId;
  p->key=key;
  p->valid=0;

  variableList[variableCount]=p;
  variableCount++;
  nodeId++;

      //increase the size of the hash table item to record the newlt generated variable
      variableHash[hIndex].num++;
      variableHash[hIndex].key=realloc(variableHash[hIndex].key,
                             sizeof(int)*variableHash[hIndex].num);
      variableHash[hIndex].index=realloc(variableHash[hIndex].index,
                             sizeof(int)*variableHash[hIndex].num);
      if( variableHash[hIndex].key==NULL || variableHash[hIndex].index==NULL)
      {
           fprintf(stderr,"Memory allocation for variables' hash table failed!\n");
           exit(1);
      }

      //Prepare a place in hash table for the new polynomial
      if(binaryEnd<=variableHash[hIndex].num-2)
      for(i=variableHash[hIndex].num-2;i>=binaryEnd;i--)
      {
         variableHash[hIndex].key[i+1]=variableHash[hIndex].key[i];
         variableHash[hIndex].index[i+1]=variableHash[hIndex].index[i];
      }

      //Insert the new polynomiual in the hash table
      variableHash[hIndex].key[binaryEnd]=key;
      variableHash[hIndex].index[binaryEnd]=p->index;

  return p;
};

///////////////////////////////////////////////////////////////////////////////////////////
//This function generates a sum polynomial.  It accepts a group of <factor, poly> pairs. //
//The sum is represented as factor_1*poly_1+factor_2*poly_2+...+factor_n*poly_n.         //
//The sub polynomials are checked to see if it can be combined with other polynomials for//
//simplification.  Also, if the sum is a constant, the result will be a constant.        //
///////////////////////////////////////////////////////////////////////////////////////////
struct polynomial *plusExp(int num, ...)
{
   int i,j,k,l;
   va_list args;
   struct sumPoly *sP;
   struct polynomial *rp;
   double            *factor;
   struct polynomial     **p;
   int counter=0;
   int flag;
   double f1,f0;
   struct polynomial     *p1=0, *p0=0;
   double con=0;
   int key=0;
   int tempI,tempI2;
   double tempD,tempD2;
   int hIndex, sIndex;
   int binaryStart,binaryEnd, binaryMiddle;
   int p0SubHIndex=0, p0HIndex=0, p0Index=0, p0Id=0, p0Key, p0Count;
   enum expressionType p0EType;

   counter_v=0;
   counter_p=0;
   counter_f=0;

   //get the number of items for this sum
   va_start(args, num);
   
   //iterate through all the items in the parameter list
   for(i=0;i<num;i++)
   {  
      //get the coefficient
      f1= va_arg( args, double );
      //get the polynomial
      p1= va_arg( args, struct polynomial *);

      if(i==0)
      {
        f0=f1;
        p0=p1;

//        fprintf(stderr,"f0=%f p0=",f0);
//        expPrinting(p0);
//        fprintf(stderr,"\n");
      }

//        fprintf(stderr,"In Sum factor=%f  item No. %d of %d type=%d\n",f1,i,num,p1->eType);
//        expPrinting(p1);
//        fprintf(stderr,"\n");

      //If an item is contant 0, it has no effect on a sum, therefore we do nothing
      if(f1==0.0 )
      {
           continue;
      }
      switch(p1->eType)
      {
         //If an item is a non-zero constant, we collect it
         case T_CONSTANT:
              con+=p1->value*f1;         
              break;
         //The item is either a variable
         case T_VARIABLE:         
              //search for the position where the new item should be inserted
              for(j=0;j<counter_v;j++)
              {
                 if(p_v[j]->index>=p1->index)
                    break;
              }
              //this is a new item in the sum, insert it at the start or end of the sum
              if(counter_v==0 || j>=counter_v)
              {
	           if(counter_v>=containerLength_v-1)
	           {
	               containerLength_v+=50;
	               factor_v = (double *)realloc(factor_v,containerLength_v*sizeof(double));
	               p_v      = (struct polynomial **)realloc(p_v,     containerLength_v*sizeof(struct polynomial *));
	               if(factor_v==NULL || p_v==NULL)
	               {
	                   fprintf(stderr,"Momery allocation error!\n");
	                   exit(1);
	               }
	           }
                   p_v[counter_v]=p1;
                   factor_v[counter_v]=f1;
                   counter_v++;
              }
              //this item is currently in the sum, just merge their coefficients
              else if(p_v[j]==p1)
              {
                   factor_v[j]+=f1;
              }
              //insert the item in the middle of the sum
              else
              {
	           if(counter_v>=containerLength_v-1)
	           {
	                containerLength_v+=50;
	                factor_v = (double *)realloc(factor_v,containerLength_v*sizeof(double));
	                p_v      = (struct polynomial **)realloc(p_v,     containerLength_v*sizeof(struct polynomial *));
	                if(factor_v==NULL || p_v==NULL)
	                {
	                    fprintf(stderr,"Momery allocation error!\n");
	                    exit(1);
	                }
	           }
                   //Move the items backward
                   for(k=counter_v-1;k>=j;k--)
                   {
                        p_v[k+1]=p_v[k];
                        factor_v[k+1]=factor_v[k];
                   }

                   //insert the new item
                   p_v[j]=p1;
                   factor_v[j]=f1;
                   counter_v++;
              }       
              break;
         case T_PRODUCT:
              //search for the position where the new item should be inserted
              for(j=0;j<counter_p;j++)
              {
                 if(p_p[j]->index>=p1->index)
                    break;
              }
              //this is a new item in the sum, insert it at the start or end of the sum
              if(counter_p==0 || j>=counter_p)
              {
                   if(counter_p>=containerLength_p-1)
                   {
                       containerLength_p+=50;
                       factor_p = (double *)realloc(factor_p,containerLength_p*sizeof(double));
                       p_p      = (struct polynomial **)realloc(p_p,     containerLength_p*sizeof(struct polynomial *));
                       if(factor_p==NULL || p_p==NULL)
                       {
                           fprintf(stderr,"Momery allocation error!\n");
                           exit(1);
                       }
                   }
                   p_p[counter_p]=p1;
                   factor_p[counter_p]=f1;
                   counter_p++;
              }
              //this item is currently in the sum, just merge their coefficients
              else if(p_p[j]==p1)
              {
                   factor_p[j]+=f1;
              }
              //insert the item in the middle of the sum
              else
              {
                   if(counter_p>=containerLength_p-1)
                   {
                        containerLength_p+=50;
                        factor_p = (double *)realloc(factor_p,containerLength_p*sizeof(double));
                        p_p      = (struct polynomial **)realloc(p_p,     containerLength_p*sizeof(struct polynomial *));
                        if(factor_p==NULL || p_p==NULL)
                        {
                            fprintf(stderr,"Momery allocation error!\n");
                            exit(1);
                        }
                   }
                   //Move the items backward
                   for(k=counter_p-1;k>=j;k--)
                   {
                        p_p[k+1]=p_p[k];
                        factor_p[k+1]=factor_p[k];
                   }
                   //insert the new item
                   p_p[j]=p1;
                   factor_p[j]=f1;
                   counter_p++;
              }
              break;

         case T_FUNCTIONCALL:
              //search for the position where the new item should be inserted
              for(j=0;j<counter_f;j++)
              {
                 if(p_f[j]->index>=p1->index)
                    break;
              }
              //this is a new item in the sum, insert it at the start or end of the sum
              if(counter_f==0 || j>=counter_f)
              {
                   if(counter_f>=containerLength_f-1)
                   {
                       containerLength_f+=50;
                       factor_f = (double *)realloc(factor_f,containerLength_f*sizeof(double));
                       p_f      = (struct polynomial **)realloc(p_f,     containerLength_f*sizeof(struct polynomial *));
                       if(factor_f==NULL || p_f==NULL)
                       {
                           fprintf(stderr,"Momery allocation error!\n");
                           exit(1);
                       }
                   }
                   p_f[counter_p]=p1;
                   factor_f[counter_f]=f1;
                   counter_f++;
              }
              //this item is currently in the sum, just merge their coefficients
              else if(p_f[j]==p1)
              {
                   factor_f[j]+=f1;
              }
              //insert the item in the middle of the sum
              else
              {
                   if(counter_f>=containerLength_f-1)
                   {
                        containerLength_f+=50;
                        factor_f = (double *)realloc(factor_f,containerLength_f*sizeof(double));
                        p_f      = (struct polynomial **)realloc(p_f,     containerLength_f*sizeof(struct polynomial *));
                        if(factor_f==NULL || p_f==NULL)
                        {
                            fprintf(stderr,"Momery allocation error!\n");
                            exit(1);
                        }
                   }
                   //Move the items backward
                   for(k=counter_f-1;k>=j;k--)
                   {
                        p_f[k+1]=p_f[k];
                        factor_f[k+1]=factor_f[k];
                   }
                   //insert the new item
                   p_f[j]=p1;
                   factor_f[j]=f1;
                   counter_f++;
              }
              break;

         case T_SUM:
              //This item is a sum, we go through the items of this sum item
              //each item in this sum is called sum item
              for(l=0;l<p1->e.s->num;l++)
              {
                  switch(p1->e.s->sum[l]->eType)
                  {
                      //If this sum item is a constant, collected it into con
                      case T_CONSTANT:
                           con+=f1*p1->e.s->sum[l]->value*p1->e.s->factor[l];
                           break;
                      //Variable
                      case T_VARIABLE:              
                           //search for a position to insert this sum item
                           for(j=0;j<counter_v;j++)
              		   {
                               if(p_v[j]->index>=p1->e.s->sum[l]->index)
                                  break;
                           }
                           //insert this sum item in the start or end of the result sum
                           if(counter_v==0 || j>=counter_v)
                           {
	                       if(counter_v>=containerLength_v-1)
	                       {
	                           containerLength_v+=50;
	                           factor_v = (double *)realloc(factor_v,containerLength_v*sizeof(double));
	                           p_v      = (struct polynomial **)realloc(p_v,  containerLength_v*sizeof(struct polynomial *));
	                           if(factor_v==NULL || p_v==NULL)
	                           {
	                               fprintf(stderr,"Momery allocation error!\n");
	                               exit(1);
	                           }
	                       }
                               p_v[counter_v]=p1->e.s->sum[l];
                               factor_v[counter_v]=f1*p1->e.s->factor[l];
                               counter_v++;
                           }
                           //this sum item is in the result sum, merge the coefficients
                           else if(p_v[j]==p1->e.s->sum[l])
                           {
                               factor_v[j]+=f1*p1->e.s->factor[l];
                           }
                           //insert this sum item in the middle of the result sum
                           else
                           {
 	                       if(counter_v>=containerLength_v-1)
	                       {
	                           containerLength_v+=50;
	                           factor_v = (double *)realloc(factor_v,containerLength_v*sizeof(double));
	                           p_v      = (struct polynomial **)realloc(p_v,     containerLength_v*sizeof(struct polynomial *));
	                           if(factor_v==NULL || p_v==NULL)
	                           {
	                                fprintf(stderr,"Momery allocation error!\n");
	                                exit(1);
	                           }
                               }
                               //search for a position
                               for(k=counter_v-1;k>=j;k--)
                               {
                                   p_v[k+1]=p_v[k];
                                   factor_v[k+1]=factor_v[k];
                               }  
                               //insert the item
                               p_v[j]=p1->e.s->sum[l];
                               factor_v[j]=f1*p1->e.s->factor[l];
                               counter_v++;
                           }
                           break;
                      //product
                      case T_PRODUCT:
                           //search for a position to insert this sum item
                           for(j=0;j<counter_p;j++)
                           {
                               if(p_p[j]->index>=p1->e.s->sum[l]->index)
                                  break;
                           }
                           //insert this sum item in the start or end of the result sum
                           if(counter_p==0 || j>=counter_p)
                           {
                               if(counter_p>=containerLength_p-1)
                               {
                                   containerLength_p+=50;
                                   factor_p = (double *)realloc(factor_p,containerLength_p*sizeof(double));
                                   p_p      = (struct polynomial **)realloc(p_p,  containerLength_p*sizeof(struct polynomial *));
                                   if(factor_p==NULL || p_p==NULL)
                                   {
                                       fprintf(stderr,"Momery allocation error!\n");
                                       exit(1);
                                   }
                               }
                               p_p[counter_p]=p1->e.s->sum[l];
                               factor_p[counter_p]=f1*p1->e.s->factor[l];
                               counter_p++;
                           }
                           //this sum item is in the result sum, merge the coefficients
                           else if(p_p[j]==p1->e.s->sum[l])
                           {  
                               factor_p[j]+=f1*p1->e.s->factor[l];
                           }
                           //insert this sum item in the middle of the result sum
                           else
                           {
                               if(counter_p>=containerLength_p-1)
                               {
                                   containerLength_p+=50;
                                   factor_p = (double *)realloc(factor_p,containerLength_p*sizeof(double));
                                   p_p      = (struct polynomial **)realloc(p_p,     containerLength_p*sizeof(struct polynomial *));
                                   if(factor_p==NULL || p_p==NULL)
                                   {
                                        fprintf(stderr,"Momery allocation error!\n");
                                        exit(1);
                                   }
                               }
                               //search for a position
                               for(k=counter_p-1;k>=j;k--)
                               {
                                   p_p[k+1]=p_p[k];
                                   factor_p[k+1]=factor_p[k];
                               }
                               //insert the item
                               p_p[j]=p1->e.s->sum[l];
                               factor_p[j]=f1*p1->e.s->factor[l];
                               counter_p++;
                           }
                           break;
                     case T_FUNCTIONCALL:
                           //search for a position to insert this sum item
                           for(j=0;j<counter_f;j++)
                           {
                               if(p_f[j]->index>=p1->e.s->sum[l]->index)
                                  break;
                           }
                           //insert this sum item in the start or end of the result sum
                           if(counter_f==0 || j>=counter_f)
                           {
                               if(counter_f>=containerLength_f-1)
                               {
                                   containerLength_f+=50;
                                   factor_f = (double *)realloc(factor_f,containerLength_f*sizeof(double));
                                   p_f      = (struct polynomial **)realloc(p_f,  containerLength_f*sizeof(struct polynomial *));
                                   if(factor_f==NULL || p_f==NULL)
                                   {
                                       fprintf(stderr,"Momery allocation error!\n");
                                       exit(1);
                                   }
                               }
                               p_f[counter_f]=p1->e.s->sum[l];
                               factor_f[counter_f]=f1*p1->e.s->factor[l];
                               counter_f++;
                           }
                           //this sum item is in the result sum, merge the coefficients
                           else if(p_f[j]==p1->e.s->sum[l])
                           {
                               factor_f[j]+=f1*p1->e.s->factor[l];
                           }
                           //insert this sum item in the middle of the result sum
                           else
                           {
                               if(counter_f>=containerLength_f-1)
                               {
                                   containerLength_f+=50;
                                   factor_f = (double *)realloc(factor_f,containerLength_f*sizeof(double));
                                   p_f      = (struct polynomial **)realloc(p_f,     containerLength_f*sizeof(struct polynomial *));
                                   if(factor_f==NULL || p_f==NULL)
                                   {
                                        fprintf(stderr,"Momery allocation error!\n");
                                        exit(1);
                                   }
                               }
                               //search for a position
                               for(k=counter_f-1;k>=j;k--)
                               {
                                   p_f[k+1]=p_f[k];
                                   factor_f[k+1]=factor_f[k];
                               }
                               //insert the item
                               p_f[j]=p1->e.s->sum[l];
                               factor_f[j]=f1*p1->e.s->factor[l];
                               counter_f++;
                           }
                           break;
                      default:
                          fprintf(stderr,"UNKNOWN polynomial type, exit!!\n");
   			  exit(1);
                   }
              }
              break;
         default:
              fprintf(stderr,"UNKNOWN polynomial type, exit!!\n");
              exit(1);
     }

//              for(k=0;k<counter;k++)
//              {
//                expPrinting(p[k]);
//                fprintf(stderr," k=%d index=%d factor=%f\n",k,p[k]->index,factor[k]);
//              }
   }
   flag= va_arg( args, int );
   va_end(args);





//flag=0;









   if(flag==0)
       sum0++;
   else
       sum1++;


//  fprintf(stderr,"111111111\n");
 
   counter=counter_v+counter_p+counter_f;
   factor = (double *)malloc((counter+1)*sizeof(double));
   p      = (struct polynomial **)malloc((counter+1)*sizeof(struct polynomial *));
   j=0;
   for(i=0;i<counter_v;i++)
   {
      factor[j] = factor_v[i];
      p[j]      = p_v[i];
      j++;
   }
   for(i=0;i<counter_p;i++)
   {
      factor[j] = factor_p[i];
      p[j]      = p_p[i];
      j++;
   }
   for(i=0;i<counter_f;i++)
   {
      factor[j] = factor_f[i];
      p[j]      = p_f[i];
      j++;
   }

   //After we go through all the items in the sum,
   //we get only a constant
   if(counter==0)
   {
      rp=constantExp(con);
      free(factor);
      free(p);
      sum2++;
      return rp;
   }
   else if(con==0.0 && counter==1 && factor[0]==1.0)
   {
      rp=p[0];
      free(factor);
      free(p);
      sum3++;
      return rp;
   }
   else
   {
      if(con!=0.0)          
      {
         p1=constantExp(con);
         factor[counter]=1.0;
         p[counter]=p1;
         counter++;
      }
   }

    //compute the key for this polynomial
      key=0;

      double f=0;
      for(j=0;j<counter;j++)
      {
          f+=fabs(factor[j]);
      }
      for(j=0;j<counter;j++)
      {
          key+=p[j]->id+p[j]->key+(int)(factor[j]*12345*(int)(p[j]->eType))%10000;
          if(key>=0)
              key=key % MAX_POLYNOMIAL_KEY;
          else
              key=key % MAX_POLYNOMIAL_KEY+MAX_POLYNOMIAL_KEY;
      }


      if(key<0)
        key=0-key;

//   fprintf(stderr,"22222222222 sum key=%d\n",key);


      hIndex=key % SUM_HASH_SIZE;

 //   fprintf(stderr,"First poly statistics\n");
 //    polyStatistics(NULL);


//     fprintf(stderr,"second poly statistics\n");
//     polyStatistics(NULL);

//   fprintf(stderr,"3333333\n");
//   fprintf(stderr,"hIndex=%d sumHash[hIndex].num=%d\n",hIndex,sumHash[hIndex].num);
      //compare the key of this polynomial with the keys of polynomials in the
      //sum list
      if(sumHash[hIndex].num>0)
      {

          //Binary search is used to locate the a sum in ths hash table
          //that share the same hash table with the new sum
          binaryStart = 0;
          binaryEnd   = sumHash[hIndex].num-1;
          while(binaryStart<=binaryEnd)
          {
             binaryMiddle=floor((binaryStart+binaryEnd)/2.0);
             if(key==sumHash[hIndex].key[binaryMiddle])
             {
                 binaryStart = binaryMiddle;
                 binaryEnd   = binaryMiddle;
                 break;
             }
             else if(key>sumHash[hIndex].key[binaryMiddle])
                 binaryStart = binaryMiddle+1;
             else
                 binaryEnd   = binaryMiddle-1;
          } //end of while

 //         fprintf(stderr,"3333333AAAAA\n");

          //If we have found the key in the hash table,
          //we compare the new poly with the old polys recorded in the hash table
          //to see if they are the same
          if(binaryStart>=0 && binaryStart<=sumHash[hIndex].num-1 
                            && key==sumHash[hIndex].key[binaryStart])
          {
              binaryMiddle = binaryStart;
              binaryEnd    = binaryStart+1;
              binaryStart  = binaryStart-1;
              while(binaryStart>=0 && 
                   sumHash[hIndex].key[binaryStart]==sumHash[hIndex].key[binaryMiddle])
                  binaryStart--;
              while(binaryEnd<sumHash[hIndex].num && 
                   sumHash[hIndex].key[binaryEnd]==sumHash[hIndex].key[binaryMiddle])
                  binaryEnd++;
              if(binaryStart<0)
                  binaryStart++;
              if(binaryEnd>=sumHash[hIndex].num || key!=sumHash[hIndex].key[binaryEnd])
                  binaryEnd--;

//          fprintf(stderr,"3333333BBBBBBB\n");
              //if the key of this sum is equal to the key of a polynomial in the 
              //sum list, compare the number of items in the two polynomials
//          fprintf(stderr,"binaryStart=%d binaryEnd=%d\n",binaryStart, binaryEnd);
             for(i=binaryStart;i<=binaryEnd;i++)
             {
	             sIndex=sumHash[hIndex].index[i];

//                     expPrinting(sumList[sIndex]);
//                     fprintf(stderr,"\n");
//                     fprintf(stderr,"hIndex=%d sIndex=%d sumCount=%d counter=%d i=%d start=%d end=%d\n",
//                              hIndex, sIndex, sumCount, counter, i, binaryStart, binaryEnd);
//                     fprintf(stderr,"num=%d\n",sumList[sIndex]->e.s->num);
	             if(counter==sumList[sIndex]->e.s->num)
	             {
	                 //compare the two sums item by item
	                 for(k=0;k<counter;k++)
                         {                           
	                      if(p[k]==sumList[sIndex]->e.s->sum[k])
                              {                
                                    tempD  = frexp(factor[k],&tempI);
                                    tempD2 = frexp(sumList[sIndex]->e.s->factor[k],&tempI2);
                                    if(tempI2!=tempI || (int)(tempD2*100000000)!=(int)(tempD*100000000))
               	                      break;                              
                              }
                              else
                                    break;
                         }
//                    fprintf(stderr,"k=%d\n",k);
	                 //if the two sums are the same, return the sum in the sum list 
	                 if(k>=counter)
	                 {
//     fprintf(stderr,"third poly statistics\n");
//     polyStatistics(NULL);

	                     free(factor);
	                     free(p);
                             sumList[sIndex]->count=2;
//                             fprintf(stderr,"Between 3333---44444 sIndex=%d\n",sIndex);
//     fprintf(stderr,"fourth poly statistics\n");
//     polyStatistics(NULL);

                             sum4++;
	                     return sumList[sIndex];
	                 }
	             } //end of if(counter==sumList[sIndex]->e.s->num)
	      } //end of for(i=binaryStart;i<=binaryEnd;i++)
	  } //end of if(binaryStart>=0 && binarySt
          else
          {
              binaryEnd = binaryStart;
          }
       }//end of if(sumHash[hIndex].num>0)
       else
       {
              binaryEnd = 0;
       }

  //fprintf(stderr,"444444444444\n");

      if(flag!=0 && p0->eType==T_SUM && p0->count==1)
      {
           p0Index    = p0->index;
           p0Id       = p0->id;
           p0Key      = p0->key;
           p0HIndex    = p0->key%SUM_HASH_SIZE;
           p0SubHIndex = p0->subHIndex;
           p0EType     = p0->eType;
           p0Count     = p0->count;                   
//           fprintf(stderr,"p0Index=%d p0Id=%d hIndex0=%d subHIndex0=%d SUM_HASH_SIZE=%d key=%d\n",
//                           p0Index,p0Id, p0HIndex,p0SubHIndex,SUM_HASH_SIZE,p0->key);
//           fprintf(stderr,"sumHash[hIndex0].num=%d subHIndex0=%d\n",sumHash[p0HIndex].num,p0SubHIndex);

           free(p0->e.s->sum);
           free(p0->e.s->factor);
           free(p0->e.s);
           free(p0);

      }
      else
      {
           p0EType=999;
           p0Count=999;
      }

      //If the sum is not found in the sum list, insert it in the sum list
      //Build a new polynomial
      rp=(struct polynomial *)malloc(sizeof(struct polynomial));
      if(rp==NULL)
      {
        fprintf(stderr,"Memory allocation error: malloc returned NULL!");
        exit(1);
      }
      rp->eType=T_SUM;
      sP=(struct sumPoly *)malloc(sizeof(struct sumPoly));
      if(sP==NULL)
      {
         fprintf(stderr,"Memory allocation error: malloc returned NULL!");
         exit(1);
      }
      sP->num    = counter;
      sP->sum    = (struct polynomial **)malloc(counter*sizeof(struct polynomial *));
      sP->factor = (double *)malloc(counter*sizeof(double));
      if(sP->sum==NULL || sP->factor==NULL)
      {
         fprintf(stderr,"Memory allocation error: malloc returned NULL!");
         exit(1);
      }
      for(i=0;i<sP->num;i++)
      {
         sP->factor[i] = factor[i];
         sP->sum[i]    = p[i];
      }
      rp->e.s       = sP;
      rp->key       = key;
      rp->valid     = 0;
      rp->subHIndex = binaryEnd;
      rp->count     = 1;

      //Insert the new built polynomial in sum list
      if(sumCount>=sumListLength){
         sumListLength+=10000;
         sumList=realloc(sumList,sumListLength*sizeof(struct polynomial *));
      }
      if(sumList==NULL)
      {
         fprintf(stderr,"Memory allocation error in plusExp, exit!");
         exit(1);
      }

//      fprintf(stderr,"sumCount=%d sumListLength=%d key=%d hIndex=%d\n",sumCount, sumListLength, key, hIndex);

//   fprintf(stderr,"5555555555\n");

      if(flag!=0 && p0EType==T_SUM && p0Count==1)
      {
              sum11++;
              sumList[p0Index]=rp;
              sumList[p0Index]->index = p0Index;
              sumList[p0Index]->id   = p0Id;
      }
      else
      {
              sum00++;
	      sumList[sumCount]=rp;
	      sumList[sumCount]->index=sumCount;           
	      sumList[sumCount]->id=nodeId;
	      sumCount++;
	      nodeId++;
      }
     
//   fprintf(stderr,"6666666666\n");

//         fprintf(stderr,"Before insert into sumHash\n");
//         for(i=0;i<sumHash[hIndex].num;i++)
//           fprintf(stderr,"hIndex=%d i=%d this key=%d key=%d index=%d\n",
//                           hIndex,i,key,sumHash[hIndex].key[i],
//                           sumHash[hIndex].index[i]);
//         fprintf(stderr,"\n");



      //Insert the newly built polynomial into the Hash table
      if(flag!=0 && p0EType==T_SUM && p0Count==1 )
      {
//         fprintf(stderr,"AAAAAA\n");
         if(p0HIndex!=hIndex || p0SubHIndex!=binaryEnd)
         {
//fprintf(stderr,"BBBBBB\n");
	      sumHash[hIndex].num++;
	      sumHash[hIndex].key=realloc(sumHash[hIndex].key,
	                             sizeof(int)*sumHash[hIndex].num);
	      sumHash[hIndex].index=realloc(sumHash[hIndex].index,
	                             sizeof(int)*sumHash[hIndex].num);
	      if( sumHash[hIndex].key==NULL || sumHash[hIndex].index==NULL)
	      {
	           fprintf(stderr,"Memory allocation for sums' hash table failed!(1) sumHash[%d].num=%d\n",hIndex,sumHash[hIndex].num);
	           exit(1);
	      }      
              //else
              //      fprintf(stderr,"Memory allocation for sums' hash table succeeded!(1) sumHash[%d].num=%d\n",hIndex,sumHash[hIndex].num);
	      if(binaryEnd<=sumHash[hIndex].num-2)
	      for(i=sumHash[hIndex].num-2;i>=binaryEnd;i--)
	      {
	         sumHash[hIndex].key[i+1]=sumHash[hIndex].key[i];
	         sumHash[hIndex].index[i+1]=sumHash[hIndex].index[i];	
                 sumList[sumHash[hIndex].index[i+1]]->subHIndex++;
	      }
	      sumHash[hIndex].key[binaryEnd]=key;
	      sumHash[hIndex].index[binaryEnd]=rp->index;
//fprintf(stderr,"CCCC\n");
              if(p0HIndex!=hIndex)
              {
//fprintf(stderr,"DDDD, p0HIndex=%d hIndex=%d p0SubHIndex=%d num=%d\n",p0HIndex,hIndex,p0SubHIndex,sumHash[p0HIndex].num);
                  for(i=p0SubHIndex;i<sumHash[p0HIndex].num-1;i++)
                  {
                      sumHash[p0HIndex].key[i]=sumHash[p0HIndex].key[i+1];
                      sumHash[p0HIndex].index[i]=sumHash[p0HIndex].index[i+1];
                      sumList[sumHash[p0HIndex].index[i]]->subHIndex--;                    
                  }
                  sumHash[p0HIndex].num--;
                  sumHash[p0HIndex].key=realloc(sumHash[p0HIndex].key,
                                     sizeof(int)*sumHash[p0HIndex].num);
                  sumHash[p0HIndex].index=realloc(sumHash[p0HIndex].index,
                                     sizeof(int)*sumHash[p0HIndex].num);
                  if( sumHash[p0HIndex].num>0 && (sumHash[p0HIndex].key==NULL || sumHash[p0HIndex].index==NULL))
                  {
                       fprintf(stderr,"Memory allocation for sums' hash table failed!(2) sumHash[%d].num=%d\n",p0HIndex,sumHash[p0HIndex].num);
                       exit(1);
                  }
                  //else
                  //     fprintf(stderr,"Memory allocation for sums' hash table succeeded(2) sumHash[%d].num=%d\n",p0HIndex,sumHash[p0HIndex].num);

                  
//fprintf(stderr,"EEEE\n");
              }
              else //p0SubHIndex!=binaryEnd
              {
//fprintf(stderr,"FFFF\n");
                 if(p0SubHIndex<binaryEnd)
                 {
//fprintf(stderr,"GGGGGG\n");
                     for(i=p0SubHIndex;i<sumHash[p0HIndex].num-1;i++)
                     {
                          sumHash[p0HIndex].key[i]=sumHash[p0HIndex].key[i+1];
                          sumHash[p0HIndex].index[i]=sumHash[p0HIndex].index[i+1];
                          sumList[sumHash[p0HIndex].index[i]]->subHIndex--;
                     }
                     sumHash[p0HIndex].num--;
                     sumHash[p0HIndex].key=realloc(sumHash[p0HIndex].key,
                                     sizeof(int)*sumHash[p0HIndex].num);
                     sumHash[p0HIndex].index=realloc(sumHash[p0HIndex].index,
                                     sizeof(int)*sumHash[p0HIndex].num);
                    if(sumHash[p0HIndex].num>0 && (sumHash[p0HIndex].key==NULL || sumHash[p0HIndex].index==NULL))
                    {
                       fprintf(stderr,"Memory allocation for sums' hash table failed! (3) sumHash[%d].num=%d\n",p0HIndex,sumHash[p0HIndex].num);
                       exit(1);
                    }
                    //else
                    //   fprintf(stderr,"Memory allocation for sums' hash table succeeded(3) sumHash[%d].num=%d\n",p0HIndex,sumHash[p0HIndex].num);



//fprintf(stderr,"HHHHHH\n");
                 }
                 else
                 {
//fprintf(stderr,"IIIII\n");
                     for(i=p0SubHIndex+1;i<sumHash[p0HIndex].num-1;i++)
                     {
                          sumHash[p0HIndex].key[i]=sumHash[p0HIndex].key[i+1];
                          sumHash[p0HIndex].index[i]=sumHash[p0HIndex].index[i+1];
                          sumList[sumHash[p0HIndex].index[i]]->subHIndex--;
                     }
                     sumHash[p0HIndex].num--;
                     sumHash[p0HIndex].key=realloc(sumHash[p0HIndex].key,
                                     sizeof(int)*sumHash[p0HIndex].num);
                     sumHash[p0HIndex].index=realloc(sumHash[p0HIndex].index,
                                     sizeof(int)*sumHash[p0HIndex].num);

                     if( sumHash[p0HIndex].num>0 && (sumHash[p0HIndex].key==NULL || sumHash[p0HIndex].index==NULL))
                     {
                       fprintf(stderr,"Memory allocation for sums' hash table failed! (4) sumHash[%d].num=%d\n",p0HIndex,sumHash[p0HIndex].num);
                       exit(1);
                     }
                     //else
                     //  fprintf(stderr,"Memory allocation for sums' hash table succeeded(4) sumHash[%d].num=%d\n",p0HIndex,sumHash[p0HIndex].num);

//fprintf(stderr,"JJJJJJ\n");
                 }
              }
         }
         else //p0HIndex==hIndex && p0SubHIndex==binaryEnd
         {
//fprintf(stderr,"KKKKK\n");
              sumHash[hIndex].key[binaryEnd]=key;
         }
     }
     else 
     {
//fprintf(stderr,"LLLLLL\n");
              sumHash[hIndex].num++;
              sumHash[hIndex].key=realloc(sumHash[hIndex].key,
                                     sizeof(int)*sumHash[hIndex].num);
              sumHash[hIndex].index=realloc(sumHash[hIndex].index,
                                     sizeof(int)*sumHash[hIndex].num);
              if( sumHash[hIndex].key==NULL || sumHash[hIndex].index==NULL)
              {
                   fprintf(stderr,"Memory allocation for hash table failed!(5) sumHash[%d].sum=%d\n",hIndex,sumHash[hIndex].num);
                   exit(1);
              }
              //else
              //     fprintf(stderr,"Memory allocation for hash table succeeded(5) sumHash[%d].sum=%d\n",hIndex,sumHash[hIndex].num);
              if(binaryEnd<=sumHash[hIndex].num-2)
              for(i=sumHash[hIndex].num-2;i>=binaryEnd;i--)
              {
                 sumHash[hIndex].key[i+1]=sumHash[hIndex].key[i];
                 sumHash[hIndex].index[i+1]=sumHash[hIndex].index[i];
                 sumList[sumHash[hIndex].index[i+1]]->subHIndex++;
              }
              sumHash[hIndex].key[binaryEnd]=key;
              sumHash[hIndex].index[binaryEnd]=rp->index;
//fprintf(stderr,"mmmmmm\n");
     }
//         for(i=0;i<sumHash[hIndex].num;i++)
//           fprintf(stderr,"hIndex=%d i=%d this key=%d key=%d index=%d\n",
//                           hIndex,i,key,sumHash[hIndex].key[i],
//                           sumHash[hIndex].index[i]);
//         fprintf(stderr,"\n");

     free(p);
     free(factor);

//     expPrinting(rp);
//     fprintf(stderr,"\n");
//     fprintf(stderr,"888888888888\n");

//     fprintf(stderr,"last poly statistics\n");
//     polyStatistics(NULL);

     sum5++;
     return rp;
};
////////////////////////////////////////////////////////////////////////////////////////
//This function create a product polynomial.  It accepts a group of <exponent, poly>  //
//pairs and the product is in the form of                                             //
//poly_1^exponet1*poly_2^exponent_2*...*poly_n^exponent_n.                            //
//It checks the sub polynomials to see if they can be combined for simplification     //
//If the exponents of the sub polynomials can be cancelled out, then the result will  //
//be a constant                                                                       //
////////////////////////////////////////////////////////////////////////////////////////

struct polynomial *timesExp(int num,...)
{
    int i,j,k,l,counter;
   va_list args;
   struct productPoly *pP;
   struct polynomial *rp;
   struct polynomial     **p;
   int                   *exponent;

   struct polynomial     *p1=0, *p0=0;
   int                    e1, e0;
   int isZero=0;
   double factor=1;
   int key=0;


   int pIndex,hIndex;
   int binaryStart, binaryMiddle, binaryEnd;
   int flag;
   int p0SubHIndex=0, p0HIndex=0, p0Index=0, p0Id=0, p0Key, p0Count;
   enum expressionType p0EType;

   int    *exponent_v;
   struct polynomial **p_v;
   int    containerLength_v;
   int    counter_v=0;

   int    *exponent_s;
   struct polynomial **p_s;
   int    containerLength_s;
   int    counter_s=0;

   int    *exponent_f;
   struct polynomial **p_f;
   int    containerLength_f;
   int    counter_f=0;

   va_start(args, num);

   //apply memory for container to hold product items
   //Apply memory for containers to hold the sum items
   containerLength_v=50;
   exponent_v = (int *)malloc(containerLength_v*sizeof(int));
   p_v      = (struct polynomial **)malloc(containerLength_v*sizeof(struct polynomial *));
   if(exponent_v==NULL || p_v==NULL)
   {
        fprintf(stderr,"Momery allocation error!\n");
        exit(1);
   }

   containerLength_s=50;
   exponent_s = (int *)malloc(containerLength_s*sizeof(int));
   p_s      = (struct polynomial **)malloc(containerLength_s*sizeof(struct polynomial *));
   if(exponent_s==NULL || p_s==NULL)
   {     
        fprintf(stderr,"Momery allocation error!\n");
        exit(1);
   }

   containerLength_f=50;
   exponent_f = (int *)malloc(containerLength_f*sizeof(double));
   p_f      = (struct polynomial **)malloc(containerLength_f*sizeof(struct polynomial *));
   if(exponent_f==NULL || p_f==NULL)
   {
        fprintf(stderr,"Momery allocation error!\n");
        exit(1);
   }

   //go through each component and its exponent of the product
   for(i=0;i<num;i++)
   {
      //get the polynomial
      p1= va_arg( args, struct polynomial *);
      //get the exponent
      e1= va_arg( args, int);
      //fprintf(stderr," TimesExp  exponent=%d item No. %d of %d type=%d\n",
      //                e1,i,num,p1->eType);
      //expPrinting(p1);
      //fprintf(stderr,"\n");
      

      if(i==0)
      {
        e0=e1;
        p0=p1;
      }


      //If any of the component is zero, then the product is zero
      if(isZeroExp(p1))
      {
         isZero=1;
         break;
      }
      //If this component is a constant, this constant and its exponent
      //is accumulated into factor
      else if(p1->eType==T_CONSTANT)
      {
         factor*=pow(p1->value,e1);
         continue;
      }
      //If the component is variable/sum/product
      else
      {
         //If this component is a sum polynomial and this sum polynomial
         //has only one item
         if(p1->eType==T_SUM && p1->e.s->num==1)
         {
            factor*=pow(p1->e.s->factor[0],e1);
            p1=p1->e.s->sum[0];
         }

         //If this component is a variable or a sum that has
         //more than one items
         switch(p1->eType)
         {
             case T_VARIABLE:
	             //Search for a position for this item in the container
	             for(j=0;j<counter_v;j++)
	             {
	                if(p_v[j]->index>=p1->index)
	                   break;              
	             }
	             //If this item should be in the start of end of the product
	             if(counter_v==0 || j>=counter_v)
	             {
	
	                if(counter_v>=containerLength_v-1)
	                {
	                   containerLength_v+=50;
	                   exponent_v   = (int *)realloc(exponent_v,containerLength_v*sizeof(int));
	                   p_v          = (struct polynomial **)realloc(p_v,     containerLength_v*sizeof(struct polynomial *));
	                   if(exponent_v==NULL || p_v==NULL)
	                   {
	                      fprintf(stderr,"Momery allocation error in timesExp()!\n");
	                      exit(1);
	                   }
	                }
	                p_v[counter_v]=p1;
	                exponent_v[counter_v]=e1;
	                counter_v++;
	             }
	             //If this component has appeared in the container, accumulate the exponent
	             else if(p_v[j]==p1)
	             {
	                exponent_v[j]+=e1;
	             }
	             //Otherwise, insert it in the container
	             else
	             {
	                if(counter_v>=containerLength_v-1)
	                {
	                   containerLength_v += 50;
	                   exponent_v         = (int *)realloc(exponent_v,containerLength_v*sizeof(int));
	                   p_v   = (struct polynomial **)realloc(p_v,     containerLength_v*sizeof(struct polynomial *));
	                   if(exponent_v==NULL || p_v==NULL)
	                   {
	                      fprintf(stderr,"Momery allocation error in timesExp()!\n");
	                      exit(1);
	                   }
	                }
	                //search for a position for this component in the container
	                for(k=counter_v-1;k>=j;k--)
	                {
	                  p_v[k+1]=p_v[k];
	                  exponent_v[k+1]=exponent_v[k];
	                }
	                //put this component in the container
	                p_v[j]=p1;
	                exponent_v[j]=e1;
	                counter_v++;
	             }
	             break;

             case T_SUM:

	             //Search for a position for this item in the container
	             for(j=0;j<counter_s;j++)
	             {
	                if(p_s[j]->index>=p1->index)
	                   break;
	             }
	             //If this item should be in the start of end of the product
	             if(counter_s==0 || j>=counter_s)
	             {
	                if(counter_s>=containerLength_s-1)
	                {
	                   containerLength_s+=50;
	                   exponent_s   = (int *)realloc(exponent_s,containerLength_s*sizeof(int));
	                   p_s   = (struct polynomial **)realloc(p_s,     containerLength_s*sizeof(struct polynomial *));
	                   if(exponent_s==NULL || p_s==NULL)
	                   {
	                      fprintf(stderr,"Momery allocation error in timesExp()!\n");
	                      exit(1);
	                   }
	                }
	                p_s[counter_s]=p1;
	                exponent_s[counter_s]=e1;
	                counter_s++;
	             }
	             //If this component has appeared in the container, accumulate the exponent
	             else if(p_s[j]==p1)
	             {
	                exponent_s[j]+=e1;
	             }
	             //Otherwise, insert it in the container
	             else
	             {
	                if(counter_s>=containerLength_s-1)
	                {
	                   containerLength_s+=50;
	                   exponent_s   = (int *)realloc(exponent_s,containerLength_s*sizeof(int));
	                   p_s   = (struct polynomial **)realloc(p_s,     containerLength_s*sizeof(struct polynomial *));
	                   if(exponent_s==NULL || p_s==NULL)
	                   {
	                      fprintf(stderr,"Momery allocation error in timesExp()!\n");
	                      exit(1);
	                   }
	                }
	                //search for a position for this component in the container
	                for(k=counter_s-1;k>=j;k--)
	                {
	                  p_s[k+1]=p_s[k];
	                  exponent_s[k+1]=exponent_s[k];
	                }   
	                //put this component in the container
	                p_s[j]=p1;
	                exponent_s[j]=e1;
	                counter_s++;
	             }
	             break;

             case T_FUNCTIONCALL:
	             //Search for a position for this item in the container
	             for(j=0;j<counter_f;j++)
	             {
	                if(p_f[j]->index>=p1->index)
	                   break;
	             }
	             //If this item should be in the start of end of the product
	             if(counter_f==0 || j>=counter_f)
	             {
	
	                if(counter_f>=containerLength_f-1)
	                {
	                   containerLength_f+=50;
	                   exponent_f   = (int *)realloc(exponent_f,containerLength_f*sizeof(int));
	                   p_f   = (struct polynomial **)realloc(p_f,     containerLength_f*sizeof(struct polynomial *));
	                   if(exponent_f==NULL || p_f==NULL)
	                   {
	                      fprintf(stderr,"Momery allocation error in timesExp()!\n");
	                      exit(1);
	                   }
	                }
	                p_f[counter_f]=p1;
	                exponent_f[counter_f]=e1;
	                counter_f++;
	             }
	             //If this component has appeared in the container, accumulate the exponent
	             else if(p_f[j]==p1)
	             {
	                exponent_f[j]+=e1;
	             }
	             //Otherwise, insert it in the container
	             else
	             {
	                if(counter_f>=containerLength_f-1)
	                {
	       	            containerLength_f+=50;
	                   exponent_f   = (int *)realloc(exponent_f,containerLength_f*sizeof(int));
	                   p_f   = (struct polynomial **)realloc(p_f,     containerLength_f*sizeof(struct polynomial *));
	                   if(exponent_f==NULL || p_f==NULL)
	                   {
	                      fprintf(stderr,"Momery allocation error in timesExp()!\n");
	                      exit(1);
	                   }
	                }
	                //search for a position for this component in the container
	                for(k=counter_f-1;k>=j;k--)
	                {
	                  p_f[k+1]=p_f[k];
	                  exponent_f[k+1]=exponent_f[k];
	                }   
	                //put this component in the container
	                p_f[j]=p1;
	                exponent_f[j]=e1;
	                counter_f++;
	             }
	             break;

             case T_PRODUCT:

	             for(l=0;l<p1->e.p->num;l++)
	             {	
                         switch(p1->e.p->product[l]->eType)
                         {
                             case T_VARIABLE:	
		                 for(j=0;j<counter_v;j++)
		                 {
		                    if(p_v[j]->index>=p1->e.p->product[l]->index)
		                       break;	
		                 }
		                 //If this item is new, save it in the container
		                 if(counter_v==0 || j>=counter_v)
		                 {
		                    //container is full, then apply for more memory
			            if(counter_v>=containerLength_v-1)
			            {
			                   containerLength_v+=50;
			                   exponent_v   = (int *)realloc(exponent_v,containerLength_v*sizeof(int));
			                   p_v   = (struct polynomial **)realloc(p_v,     containerLength_v*sizeof(struct polynomial *));
			                   if(exponent_v==NULL || p_v==NULL)
			                   {
			                      fprintf(stderr,"Momery allocation error in timesExp()!\n");
			                      exit(1);
			                   }
			            }
		                    p_v[counter_v]=p1->e.p->product[l];
		                    exponent_v[counter_v]=e1*p1->e.p->exponent[l];
		                    counter_v++;
		                 }
		                 //merge this item with another item that has existed in the container
		                 //for example, y^2+y^2=2*y^2
		                 else if(p_v[j]==p1->e.p->product[l])
		                 {
		                    exponent_v[j]+=e1*p1->e.p->exponent[l];
		                 }
		                 else
		                 {
		                    //container is full, then apply for more memory
		                    if(counter_v>=containerLength_v-1)
		                    {
			                   containerLength_v+=50;
			                   exponent_v   = (int *)realloc(exponent_v,containerLength_v*sizeof(int));
			                   p_v   = (struct polynomial **)realloc(p_v,     containerLength_v*sizeof(struct polynomial *));
			                   if(exponent_v==NULL || p_v==NULL)
			                   {
			                      fprintf(stderr,"Momery allocation error in timesExp()!\n");
			                      exit(1);
			                   }
		                    }
		                    //create an empty space for the new item
		                    for(k=counter_v-1;k>=j;k--)
		                    {
		                      p_v[k+1]=p_v[k];
		                      exponent_v[k+1]=exponent_v[k];                    
		                    }
		                    p_v[j]=p1->e.p->product[l];
		                    exponent_v[j]=e1*p1->e.p->exponent[l];
		                    counter_v++;
                                 }
                                 break;

                             case T_SUM:   
                                 for(j=0;j<counter_s;j++)
                                 {
                                    if(p_s[j]->index>=p1->e.p->product[l]->index)
                                       break;
                                 }
                                 //If this item is new, save it in the container
                                 if(counter_s==0 || j>=counter_s)
                                 {
                                    //container is full, then apply for more memory
                                    if(counter_s>=containerLength_s-1)
                                    {
                                           containerLength_s+=50;
                                           exponent_s   = (int *)realloc(exponent_s,containerLength_s*sizeof(int));
                                           p_s   = (struct polynomial **)realloc(p_s,     containerLength_s*sizeof(struct polynomial *));
                                           if(exponent_s==NULL || p_s==NULL)
                                           {
                                              fprintf(stderr,"Momery allocation error in timesExp()!\n");
                                              exit(1);
                                           }
                                    }
                                    p_s[counter_s]=p1->e.p->product[l];
                                    exponent_s[counter_s]=e1*p1->e.p->exponent[l];
                                    counter_s++;
                                 }
                                 //merge this item with another item that has existed in the container
                                 //for example, y^2+y^2=2*y^2
                                 else if(p_s[j]==p1->e.p->product[l])
                                 {
                                    exponent_s[j]+=e1*p1->e.p->exponent[l];
                                 }
                                 else
                                 {
                                    //container is full, then apply for more memory
                                    if(counter_s>=containerLength_s-1)
                                    {
                                           containerLength_s+=50;
                                           exponent_s   = (int *)realloc(exponent_s,containerLength_s*sizeof(int));
                                           p_s   = (struct polynomial **)realloc(p_s,     containerLength_s*sizeof(struct polynomial *));
                                           if(exponent_s==NULL || p_s==NULL)
                                           {
                                              fprintf(stderr,"Momery allocation error in timesExp()!\n");
                                              exit(1);
                                           }
                                    }
                                    //create an empty space for the new item
                                    for(k=counter_s-1;k>=j;k--)
                                    {
                                      p_s[k+1]=p_s[k];
                                      exponent_s[k+1]=exponent_s[k];
                                    }
                                    p_s[j]=p1->e.p->product[l];
                                    exponent_s[j]=e1*p1->e.p->exponent[l];
                                    counter_s++;
                                 }
                                 break;

                             case T_FUNCTIONCALL:   
                                 for(j=0;j<counter_f;j++)
                                 {
                                    if(p_f[j]->index>=p1->e.p->product[l]->index)
                                       break;
                                 }
                                 //If this item is new, save it in the container
                                 if(counter_f==0 || j>=counter_f)
                                 {
                                    //container is full, then apply for more memory
                                    if(counter_f>=containerLength_f-1)
                                    {
                                           containerLength_f+=50;
                                           exponent_f   = (int *)realloc(exponent_f,containerLength_f*sizeof(int));
                                           p_f          = (struct polynomial **)realloc(p_f,     containerLength_f*sizeof(struct polynomial *));
                                           if(exponent_f==NULL || p_f==NULL)
                                           {
                                              fprintf(stderr,"Momery allocation error in timesExp()!\n");
                                              exit(1);
                                           }
                                    }
                                    p_f[counter_f]=p1->e.p->product[l];
                                    exponent_f[counter_f]=e1*p1->e.p->exponent[l];
                                    counter_f++;
                                 }
                                 //merge this item with another item that has existed in the container
                                 //for example, y^2+y^2=2*y^2
                                 else if(p_f[j]==p1->e.p->product[l])
                                 {
                                    exponent_f[j]+=e1*p1->e.p->exponent[l];
                                 }
                                 else
                                 {
                                    //container is full, then apply for more memory
                                    if(counter_f>=containerLength_f-1)
                                    {
                                           containerLength_f+=50;
                                           exponent_f   = (int *)realloc(exponent_f,containerLength_f*sizeof(int));
                                           p_f   = (struct polynomial **)realloc(p_f,     containerLength_f*sizeof(struct polynomial *));
                                           if(exponent_f==NULL || p_f==NULL)
                                           {
                                              fprintf(stderr,"Momery allocation error in timesExp()!\n");
                                              exit(1);
                                           }
                                    }
                                    //create an empty space for the new item
                                    for(k=counter_f-1;k>=j;k--)
                                    {
                                      p_f[k+1]=p_f[k];
                                      exponent_f[k+1]=exponent_f[k];
                                    }
                                    p_f[j]=p1->e.p->product[l];
                                    exponent_f[j]=e1*p1->e.p->exponent[l];
                                    counter_f++;
                                 }
                                 break;

                             default:
                                 fprintf(stderr,"unknown expression type\n");
 				 exit(1);
                           } //end of switch
                    } //end of for
                    break;

             default:
                    fprintf(stderr,"unknown expression type\n");
		    break;
            } //end of switch
        } //end of else
//      for(k=0;k<counter;k++)
//      {
//        expPrinting(p[k]);
//        fprintf(stderr," list member %d exponent %d\n",k,exponent[k]);
//      }
   } //end of for

   flag= va_arg( args, int );
   va_end(args);
   



//flag=0;


                                                           
   if(flag==0)
       product0++;
   else
       product1++;

   //fprintf(stderr,"isZero=%d counter=%d factor=%f\n",isZero,counter,factor);

   //The product is zero
   if(isZero)
   {
      rp=constantExp(0.0);
      free(exponent_v);
      free(p_v);
      free(exponent_s);
      free(p_s);
      free(exponent_f);
      free(p_f);
      product2++;
      return rp;
   }

   //combine all the items into one polynomial
   counter=counter_v+counter_s+counter_f;
   exponent = (int *)malloc((counter+1)*sizeof(int));
   p        = (struct polynomial **)malloc((counter+1)*sizeof(struct polynomial *));

   j=0;
   for(i=0;i<counter_v;i++)
   {
      exponent[j]=exponent_v[i];
      p[j]=p_v[i];
      j++;
   }
   for(i=0;i<counter_s;i++)
   {
      exponent[j]=exponent_s[i];
      p[j]=p_s[i];
      j++;
   }
   for(i=0;i<counter_f;i++)
   {
      exponent[j]=exponent_f[i];
      p[j]=p_f[i];
      j++;
   }

   //The product has 0 items
   if(counter==0)
   {
      rp=constantExp(factor);
      free(exponent_v);
      free(p_v);
      free(exponent_s);
      free(p_s);
      free(exponent_f);
      free(p_f);
      free(exponent);
      free(p);
      product3++;
      return rp;
   }
   else if(counter==1 && exponent[0]==1)
   {
      if(factor==1.0)
      {
         rp=p[0];
         free(exponent_v);
         free(p_v);
         free(exponent_s);
         free(p_s);
         free(exponent_f);
         free(p_f);
         free(exponent);
         free(p);
         product4++;
         return rp;
      }
      else
      {
         rp=plusExp(1,factor,p[0],0);
         free(exponent_v);      
         free(p_v);      
         free(exponent_s);      
         free(p_s);      
         free(exponent_f);      
         free(p_f);      
         free(exponent);
         free(p);
         product5++;
         return rp;
      }
   }
   else
   {

      //compute the key for this polynomial
      key=0;
      double e=0;
      for(j=0;j<counter;j++)
      {
         e+=abs(exponent[j]);
      }
      for(j=0;j<counter;j++)
      {
         key+=p[j]->id+p[j]->key+(exponent[j]*12345*(int)(p[j]->eType))%10000;
//         fprintf(stderr,"  product key: key=%d  p->key=%d exponent=%d e=%f\n",key,p[j]->key,exponent[j],e);
         if(key>=0)
            key=key % MAX_POLYNOMIAL_KEY;
         else
            key=key % MAX_POLYNOMIAL_KEY+MAX_POLYNOMIAL_KEY;

      }      

     if(key<0)
        key=0-key;

//      fprintf(stderr,"product Key=%d\n",key);


      hIndex=key % PRODUCT_HASH_SIZE;

      //compare the key of this polynomial with the keys of polynomials in the same hash item
      if(productHash[hIndex].num>0)
      {
          binaryStart = 0;
          binaryEnd   = productHash[hIndex].num-1;
          while(binaryStart<=binaryEnd)
          {
             binaryMiddle=floor((binaryStart+binaryEnd)/2.0);
             if(key==productHash[hIndex].key[binaryMiddle])
             {
                 binaryStart = binaryMiddle;
                 binaryEnd   = binaryMiddle;
                 break;
             }
             else if(key>productHash[hIndex].key[binaryMiddle])
                 binaryStart = binaryMiddle+1;
             else
                 binaryEnd   = binaryMiddle-1;
          } //end of while



          //search in the hash table to see if the new polynomial has been existing
          if(binaryStart>=0 && binaryStart<=productHash[hIndex].num-1 
                            && key==productHash[hIndex].key[binaryStart])
          {
              binaryMiddle = binaryStart;
              binaryEnd    = binaryStart+1;
              binaryStart  = binaryStart-1;
              while(binaryStart>=0 && 
                   productHash[hIndex].key[binaryStart]==productHash[hIndex].key[binaryMiddle])
                  binaryStart--;
              while(binaryEnd<productHash[hIndex].num && 
                   productHash[hIndex].key[binaryEnd]==productHash[hIndex].key[binaryMiddle])
                  binaryEnd++;
              if(binaryStart<0)
                  binaryStart++;
              if(binaryEnd>=productHash[hIndex].num ||  key!=productHash[hIndex].key[binaryEnd])
                  binaryEnd--;

              //if the key of this sum is equal to the key of a polynomial in the 
              //sum list, compare the number of items in the two polynomials
             for(i=binaryStart;i<=binaryEnd;i++)
             {

                 pIndex=productHash[hIndex].index[i];
//            fprintf(stderr,"counter=%d num=%d\n",counter,productList[pIndex]->e.p->num);
                 //If the numer of items of the new polynomial is equal to the number of 
                 //items of another polynomial whose key is equal to that of the new polynomial
                 if(counter==productList[pIndex]->e.p->num)
                 {
	               for(k=0;k<counter;k++)
	                   if(p[k]!=productList[pIndex]->e.p->product[k] 
	                       || exponent[k]!=productList[pIndex]->e.p->exponent[k])
	                     break;
	               //If the polynomials are the same
	               if(k>=counter)
                       {
                          //this item has been existing, return it
	                  if(factor==1.0)
	                  {

			      free(exponent_v);
			      free(p_v);
			      free(exponent_s);
			      free(p_s);
			      free(exponent_f);
			      free(p_f);
	                     free(exponent);
	                     free(p);
                             productList[pIndex]->count=2;
                             product6++;
	                     return productList[pIndex];
	                  }
                          //factor is different, return a sum
	                  else
	                  {
                              free(exponent_v);
                              free(p_v);
                              free(exponent_s);
                              free(p_s);
                              free(exponent_f);
                              free(p_f);
	                     free(exponent);
	                     free(p);
                             product7++;
                             productList[pIndex]->count=2;
	                     return plusExp(1,factor, productList[pIndex],0);
	                  }
                       }//end of if
	         }//end of if
             }//end of for
          }//end of if
          else
          {
              binaryEnd = binaryStart;
          }
      }//end of if
      else
      {
              binaryEnd = 0;
      }


      if(flag!=0 && p0->eType==T_PRODUCT && p0->count==1)
      {
           p0Index     = p0->index;
           p0Id        = p0->id;
           p0Key       = p0->key;
           p0HIndex    = p0->key%PRODUCT_HASH_SIZE;
           p0SubHIndex = p0->subHIndex;
           p0EType     = p0->eType;
           p0Count     = p0->count;           

                          
//           fprintf(stderr,"p0Index=%d p0Id=%d hIndex0=%d subHIndex0=%d SUM_HASH_SIZE=%d key=%d\n",
//                           p0Index,p0Id, p0HIndex,p0SubHIndex,SUM_HASH_SIZE,p0->key);
//           fprintf(stderr,"sumHash[hIndex0].num=%d subHIndex0=%d\n",sumHash[p0HIndex].num,p0SubHIndex);
                          
           free(p0->e.p->exponent);
           free(p0->e.p->product);
           free(p0->e.p);
           free(p0);                                                      

      }
      else
      {
           p0EType=999;
           p0Count=999;

      }


      //fprintf(stderr,"flag=%d p0EType=%d p0Count=%d p0Key=%d PRODUCT_HASH_SIZE=%d p0HIndex=[%d %d]\n",
      //        flag,p0EType,p0Count,p0Key,PRODUCT_HASH_SIZE,p0->key%PRODUCT_HASH_SIZE,p0HIndex);


      //if the new polynomial is not found, we build this polynomial from the
      //components we save in the container
      rp=(struct polynomial *)malloc(sizeof(struct polynomial));
      if(rp==NULL)
          fprintf(stderr,"Memory allocation error: malloc returned NULL!");
      rp->eType=T_PRODUCT;
      pP=(struct productPoly *)malloc(sizeof(struct productPoly));
      if(pP==NULL)
         fprintf(stderr,"Memory allocation error: malloc returned NULL!");
      pP->num        = counter;
      pP->product    = (struct polynomial **)malloc(counter*sizeof(struct polynomial *));
      pP->exponent   = (int *)malloc(counter*sizeof(int));
      if(pP->product==NULL || pP->exponent==NULL)
         fprintf(stderr,"Memory allocation error: malloc returned NULL!");
      for(i=0;i<counter;i++)
      {
         pP->product[i]    = p[i];
         pP->exponent[i]   = exponent[i];
      }
      rp->e.p=pP;
      rp->index=productCount;
      rp->id=nodeId;
      rp->key=key;
      rp->valid=0;
      rp->subHIndex = binaryEnd;
      rp->count=1;


//      rp->count=(int *)malloc(10*sizeof(int));
//      rp->values=(double *)malloc(sizeof(double)*10);
//      memset(rp->count,0,10*sizeof(int));

 
      //After the new polynomial is built, it is recorded in the product polynomial list
      if(productCount>=productListLength){
         productListLength+=10000;
         productList=realloc(productList,productListLength*sizeof(struct polynomial *));
      }
      if(productList==NULL){
         fprintf(stderr,"Memory allocation error in timesExp, exit!");
         exit(1);
      }      


      if(flag!=0 && p0EType==T_PRODUCT && p0Count==1)
      {
              product11++;
              productList[p0Index]=rp;
              productList[p0Index]->index = p0Index;
              productList[p0Index]->id   = p0Id;
      }
      else
      {
              product00++;
              productList[productCount]=rp;
              productList[productCount]->index=productCount;
              productList[productCount]->id=nodeId;
              productCount++;
              nodeId++;
      }

      //the new polynomial is also recorded in the hash table
      if(flag!=0 && p0EType==T_PRODUCT && p0Count==1 )
      {

         if(p0HIndex!=hIndex || p0SubHIndex!=binaryEnd)
         {

              //fprintf(stderr,"productHash[%d]=%d\n",hIndex,productHash[hIndex].num);
              productHash[hIndex].num++;
              productHash[hIndex].key=realloc(productHash[hIndex].key,
                                  sizeof(int)*productHash[hIndex].num);
              productHash[hIndex].index=realloc(productHash[hIndex].index,
                                    sizeof(int)*productHash[hIndex].num);
              if( productHash[hIndex].key==NULL || productHash[hIndex].index==NULL)
              {   
                   fprintf(stderr,"Memory allocation for products' hash table failed!(1) productHash[%d].num=%d\n",hIndex,productHash[hIndex].num);
                   exit(1);
              }
              //else
              //     fprintf(stderr,"Memory allocation for products' hash table succeeded!(1) productHash[%d].num=%d\n",hIndex,productHash[hIndex].num);

              if(binaryEnd<=productHash[hIndex].num-2)
              for(i=productHash[hIndex].num-2;i>=binaryEnd;i--)
              {
                 productHash[hIndex].key[i+1]   = productHash[hIndex].key[i];
                 productHash[hIndex].index[i+1] = productHash[hIndex].index[i];
                 productList[productHash[hIndex].index[i+1]]->subHIndex++;
              }   
              productHash[hIndex].key[binaryEnd]=key;
              productHash[hIndex].index[binaryEnd]=rp->index;

              if(p0HIndex!=hIndex)
              {
                  //fprintf(stderr,"DDDD, p0HIndex=%d hIndex=%d p0SubHIndex=%d\n",p0HIndex,hIndex,p0SubHIndex);
                  //fprintf(stderr,"productHash[%d].num=%d\n",p0HIndex,productHash[p0HIndex].num);
                  for(i=p0SubHIndex;i<productHash[p0HIndex].num-1;i++)
                  {
                      productHash[p0HIndex].key[i]   = productHash[p0HIndex].key[i+1];
                      productHash[p0HIndex].index[i] = productHash[p0HIndex].index[i+1];
                      productList[productHash[p0HIndex].index[i]]->subHIndex--;
                  }
                  productHash[p0HIndex].num--;
                  productHash[p0HIndex].key=realloc(productHash[p0HIndex].key,
                                  sizeof(int)*productHash[p0HIndex].num);
                  productHash[p0HIndex].index=realloc(productHash[p0HIndex].index,
                                    sizeof(int)*productHash[p0HIndex].num);
                 if( productHash[p0HIndex].num>0 &&( productHash[p0HIndex].key==NULL || productHash[p0HIndex].index==NULL))
                 {
                         fprintf(stderr,"Memory allocation for products' hash table failed!(2) productHash[%d].num=%d\n",p0HIndex,productHash[p0HIndex].num);
                         exit(1);
                 }
                 //else
                 //        fprintf(stderr,"Memory allocation for products' hash table succeeded!(2)  productHash[%d].num=%d\n",p0HIndex,productHash[p0HIndex].num);



              }
              else //p0SubHIndex!=binaryEnd
              {

                 if(p0SubHIndex<binaryEnd)
                 {

                     for(i=p0SubHIndex;i<productHash[p0HIndex].num-1;i++)
                     {
                          productHash[p0HIndex].key[i]   = productHash[p0HIndex].key[i+1];
                          productHash[p0HIndex].index[i] = productHash[p0HIndex].index[i+1];
                          productList[productHash[p0HIndex].index[i]]->subHIndex--;
                     }
                     productHash[p0HIndex].num--;
                     productHash[p0HIndex].key=realloc(productHash[p0HIndex].key,
                                  sizeof(int)*productHash[p0HIndex].num);
                     productHash[p0HIndex].index=realloc(productHash[p0HIndex].index,
                                    sizeof(int)*productHash[p0HIndex].num);

                     if( productHash[p0HIndex].num>0 && (productHash[p0HIndex].key==NULL || productHash[p0HIndex].index==NULL))
                     {
                         fprintf(stderr,"Memory allocation for products' hash table failed!(3) productHash[%d].num=%d\n",p0HIndex,productHash[p0HIndex].num);
                         exit(1);
                     }
                     //else
                     //    fprintf(stderr,"Memory allocation for products' hash table succeeded!(3)  productHash[%d].num=%d\n",p0HIndex,productHash[p0HIndex].num);

                 }
                 else
                 {

                     for(i=p0SubHIndex+1;i<productHash[p0HIndex].num-1;i++)
                     {
                          productHash[p0HIndex].key[i]   = productHash[p0HIndex].key[i+1];
                          productHash[p0HIndex].index[i] = productHash[p0HIndex].index[i+1];
                          productList[productHash[p0HIndex].index[i]]->subHIndex--;
                     }
                     productHash[p0HIndex].num--;
                     productHash[p0HIndex].key=realloc(productHash[p0HIndex].key,
                                  sizeof(int)*productHash[p0HIndex].num);
                     productHash[p0HIndex].index=realloc(productHash[p0HIndex].index,
                                    sizeof(int)*productHash[p0HIndex].num);
                     if( productHash[p0HIndex].num>0 && (productHash[p0HIndex].key==NULL || productHash[p0HIndex].index==NULL))
                     {
                         fprintf(stderr,"Memory allocation for products' hash table failed!(3) productHash[%d].num=%d\n",p0HIndex,productHash[p0HIndex].num);
                         exit(1);
                     }
                     //else
                     //    fprintf(stderr,"Memory allocation for products' hash table succeeded!(3)  productHash[%d].num=%d\n",p0HIndex,productHash[p0HIndex].num);



                 }
              }
         }
         else //p0HIndex==hIndex && p0SubHIndex==binaryEnd
         {

              productHash[hIndex].key[binaryEnd]=key;
              
         }
     }
     else
     {
              productHash[hIndex].num++;
              productHash[hIndex].key=realloc(productHash[hIndex].key,
                                     sizeof(int)*productHash[hIndex].num);
              productHash[hIndex].index=realloc(productHash[hIndex].index,
                                     sizeof(int)*productHash[hIndex].num);
              if( productHash[hIndex].key==NULL || productHash[hIndex].index==NULL)
              {
                    fprintf(stderr,"Memory allocation for products' hash table failed!(5) productHash[%d].num=%d\n",hIndex,productHash[hIndex].num);
                    exit(1);
              }
              //else
              //      fprintf(stderr,"Memory allocation for products' hash table succeeded!(5) productHash[%d]num=%d\n",hIndex,productHash[hIndex].num);

              if(binaryEnd<=productHash[hIndex].num-2)
              for(i=productHash[hIndex].num-2;i>=binaryEnd;i--)
              {
                 productHash[hIndex].key[i+1]   = productHash[hIndex].key[i];
                 productHash[hIndex].index[i+1] = productHash[hIndex].index[i];
                 productList[productHash[hIndex].index[i+1]]->subHIndex++;
              }
              productHash[hIndex].key[binaryEnd]=key;
              productHash[hIndex].index[binaryEnd]=rp->index;

     }

      //free all the memory before return
      if(factor==1.0)
      {
	      free(exponent_v);
	      free(p_v);
	      free(exponent_s);
	      free(p_s);
	      free(exponent_f);
	      free(p_f);
        free(exponent);
        free(p);
        product8++;
        return rp;
      }
      else
      {
	      free(exponent_v);
	      free(p_v);
	      free(exponent_s);
	      free(p_s);
	      free(exponent_f);
	      free(p_f);
        free(exponent);
        free(p);
        product9++;
        return plusExp(1,factor,rp,0);
      }
   }
};

////////////////////////////////////////////////////////////////////////////////////////
//This function create a function call.  It accepts the function name and parameters of
//the function
//
////////////////////////////////////////////////////////////////////////////////////////
struct polynomial *functionCallExp(int num, ...)
{
   int i,k,hIndex,fIndex;
   char *fName;
   int key=0;
   struct polynomial     **p;
   int binaryStart, binaryMiddle, binaryEnd;
   va_list args;
   struct polynomial *   rp;
   struct functionPoly  *fP;

   va_start(args, num);
   fName = va_arg( args, char *);

   p=(struct polynomial **)malloc((num-1)*sizeof(struct polynomial *));

//      fprintf(stderr,"Name=%s num=%d\n",fName,num);
   

   for(i=0;i<strlen(fName);i++)
       key+=(int)fName[i];   
   for(i=0;i<num-1;i++)
   {
      p[i]= va_arg( args, struct polynomial *);
      key = (key+p[i]->key) % FUNCTIONCALL_HASH_SIZE;
   }  

   if(key<0)
      key=0-key;

   va_end(args);

//   fprintf(stderr,"new key=%d\n",key);

   hIndex=key % FUNCTIONCALL_HASH_SIZE;

      //compare the key of this polynomial with the keys of polynomials in the
      //function list
      if(functionCallHash[hIndex].num>0)
      {
          binaryStart = 0;
          binaryEnd   = functionCallHash[hIndex].num-1;
          while(binaryStart<=binaryEnd)
          {
             binaryMiddle=floor((binaryStart+binaryEnd)/2.0);
             if(key==functionCallHash[hIndex].key[binaryMiddle])
             {
                 binaryStart = binaryMiddle;
                 binaryEnd   = binaryMiddle;
                 break;
             }
             else if(key>functionCallHash[hIndex].key[binaryMiddle])
                 binaryStart = binaryMiddle+1;
             else
                 binaryEnd   = binaryMiddle-1;
          } //end of while

          if(binaryStart>=0 && binaryStart<=functionCallHash[hIndex].num-1 
                            && key==functionCallHash[hIndex].key[binaryStart])
          {
              binaryMiddle = binaryStart;
              binaryEnd    = binaryStart+1;
              binaryStart  = binaryStart-1;
              while(binaryStart>=0 && 
                   functionCallHash[hIndex].key[binaryStart]==functionCallHash[hIndex].key[binaryMiddle])
                  binaryStart--;
              while(binaryEnd<functionCallHash[hIndex].num && 
                   functionCallHash[hIndex].key[binaryEnd]==functionCallHash[hIndex].key[binaryMiddle])
                  binaryEnd++;
              if(binaryStart<0)
                  binaryStart++;
              if(binaryEnd>=functionCallHash[hIndex].num)
                  binaryEnd--;
              //if the key of this function is equal to the key of a polynomial in the 
              //function call list, compare the number of items in the two polynomials

             for(i=binaryStart;i<=binaryEnd;i++)
             {
	             fIndex=functionCallHash[hIndex].index[i];


	             if(num-1==functionCallList[fIndex]->e.f->paraNum)
	             {
	                 //compare the two function calls item by item
	                 for(k=0;k<num-1;k++)
	                   if(p[k]!=functionCallList[fIndex]->e.f->para[k])
	                     break;


	                 //if the two function calls are the same, 
                         //return the function call in the function call list 
	                 if(k>=num-1)
	                 {
	                     free(p);
	                     return functionCallList[fIndex];
	                 }
	             } //end of if(num-1==functionCallList[fIndex]->e.f->paraNum)
	      } //end of for(i=binaryStart;i<=binaryEnd;i++)
	  } //end of if(binaryStart>=0 && binarySt
          else
          {
              binaryEnd = binaryStart;
          }
       }//end of if(functionCallHash[hIndex].num>0)
       else
       {
              binaryEnd = 0;
       }


      //If the function call is not found in the list, insert it in the list
      //Build a new polynomial
      rp=(struct polynomial *)malloc(sizeof(struct polynomial));
      rp->eType=T_FUNCTIONCALL;
      fP=(struct functionPoly *)malloc(sizeof(struct functionPoly));
      if( rp==NULL || fP==NULL )
         fprintf(stderr,"Memory allocation error: malloc returned NULL!");
      fP->paraNum    = num-1;
      fP->para    = (struct polynomial **)malloc((num-1)*sizeof(struct polynomial *));
      fP->name    = (char *)malloc(strlen(fName)+1);
      strcpy(fP->name,fName);

//      rp->count=(int *)malloc(10*sizeof(int));
//      rp->values=(double *)malloc(sizeof(double)*10);
//      memset(rp->count,0,10*sizeof(int));


      if(fP->para==NULL)
      {
         fprintf(stderr,"Memory allocation error: malloc returned NULL!");
         exit(1);
      }
      for(i=0;i<fP->paraNum;i++)
      {
         fP->para[i]    = p[i];
      }
      rp->e.f=fP;
      rp->index=functionCallCount;
      rp->id=nodeId;
      rp->key=key;
      rp->valid=0;
      //Insert the new built polynomial in function call list
      if(functionCallCount>=functionCallListLength){
         functionCallListLength+=10000;
         functionCallList=realloc(functionCallList,functionCallListLength*sizeof(struct polynomial *));
      }
      if(functionCallList==NULL){
         fprintf(stderr,"Memory allocation error in functionCallExp, exit!");
         exit(1);
      }      

//      fprintf(stderr,"functionCount=%d functionListLength=%d key=%d hIndex=%d\n",
//              functionCallCount, functionCallListLength, key, hIndex);

      functionCallList[functionCallCount]=rp;
      functionCallCount++;
      nodeId++;


//         fprintf(stderr,"Before insert into functionCallHash\n");
//         for(i=0;i<functionCallHash[hIndex].num;i++)
//           fprintf(stderr,"hIndex=%d i=%d this key=%d key=%d index=%d\n",
//                           hIndex,i,key,functionCallHash[hIndex].key[i],
//                           functionCallHash[hIndex].index[i]);
//         fprintf(stderr,"\n");


      //Insert the newly built polynomial into the Hash table
      functionCallHash[hIndex].num++;
      functionCallHash[hIndex].key=realloc(functionCallHash[hIndex].key,
                             sizeof(int)*functionCallHash[hIndex].num);
      functionCallHash[hIndex].index=realloc(functionCallHash[hIndex].index,
                             sizeof(int)*functionCallHash[hIndex].num);
      if( functionCallHash[hIndex].key==NULL || functionCallHash[hIndex].index==NULL)
      {
           fprintf(stderr,"Memory allocation for functionCalls' hash table failed!\n");
           exit(1);
      }      
      if(binaryEnd<=functionCallHash[hIndex].num-2)
      for(i=functionCallHash[hIndex].num-2;i>=binaryEnd;i--)
      {
         functionCallHash[hIndex].key[i+1]=functionCallHash[hIndex].key[i];
         functionCallHash[hIndex].index[i+1]=functionCallHash[hIndex].index[i];
      }

      functionCallHash[hIndex].key[binaryEnd]=key;
      functionCallHash[hIndex].index[binaryEnd]=functionCallCount-1;


/*
         fprintf(stderr,"After insert into functionCallHash binaryEnd=%d functionCallCount-1=%d\n",
                  binaryEnd,functionCallCount-1); 
         for(i=0;i<functionCallHash[hIndex].num;i++)
           fprintf(stderr,"hIndex=%d i=%d this key=%d key=%d index=%d name=%s\n",
                           hIndex,i,key,functionCallHash[hIndex].key[i],
                           functionCallHash[hIndex].index[i],functionCallList[functionCallCount-1]->e.f->name);
         fprintf(stderr,"\n");
*/
     free(p);

     return rp;

};

void expPrinting(struct polynomial *p)
{
   int i;

//   fprintf(stderr,"Inside expPrinting p->eType=%d\n",p->eType);

   switch(p->eType){
     case T_CONSTANT:
        fprintf(stderr,"%f",p->value);
        break;

     case T_VARIABLE:
        fprintf(stderr,"%s",p->e.v->vName);
        break;
  
     case T_SUM:
        if(p->e.s->num>1)
          fprintf(stderr,"(");

        if(fabs(p->e.s->factor[0])!=1.0)
           fprintf(stderr,"%f*",p->e.s->factor[0]);
        else if(p->e.s->factor[0]==-1.0)
           fprintf(stderr,"-");
        expPrinting(p->e.s->sum[0]);
        
        for(i=1;i<p->e.s->num;i++)
        {
          if(p->e.s->factor[i]==1.0)
             fprintf(stderr,"+");
          else if(p->e.s->factor[i]==-1.0)
             fprintf(stderr,"-");
          else if(p->e.s->factor[i]>=0.0)
             fprintf(stderr,"+%f*",p->e.s->factor[i]);
          else
             fprintf(stderr,"%f*",p->e.s->factor[i]);

          
          expPrinting(p->e.s->sum[i]);
        }
        if(p->e.s->num>1)
          fprintf(stderr,")");
        break;

     case T_PRODUCT:

//        if(p->e.p->num>1)
//          fprintf(stderr,"(");

        expPrinting(p->e.p->product[0]);
        if(p->e.p->exponent[0]!=1)
           fprintf(stderr,"^%d",p->e.p->exponent[0]);
        for(i=1;i<p->e.s->num;i++)
        {
          fprintf(stderr,"*");
          expPrinting(p->e.p->product[i]);
          if(p->e.p->exponent[i]!=1)
             fprintf(stderr,"^%d",p->e.p->exponent[i]);
        }

//        if(p->e.p->num>1)
//          fprintf(stderr,")");
        break;
      case T_FUNCTIONCALL:
         fprintf(stderr,"%s(",p->e.f->name);
         for(i=0;i<p->e.f->paraNum-1;i++)
         {
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


struct polyList *buildPolyList()
{
   struct polyList *l;
   l=(struct polyList *)malloc(sizeof(struct polyList));
   l->listSize=100;
   l->listNext=0;
   l->pList=(struct polynomial **)malloc(sizeof(struct polynomial *)*l->listSize);
   return l;
};

void polyListAppend(struct polyList *l, struct polynomial *p)
{
    //valid is a mark showing that this polynomial appears on a sorting list
    p->valid=1;


    if(l->listNext>=l->listSize)  
    {
        l->pList=realloc(l->pList,sizeof(struct polynomial *)*(l->listSize+100));
        l->listSize=l->listSize+100;
//        fprintf(stderr,"Length of polyList=%d size of polyList=%d\n",
//                l->listNext, sizeof(struct polynomial *)*l->listSize);
    } 
    l->pList[l->listNext]=p;
    l->listNext++;
};

/////////////////////////////////////////////////////////
// Sort the polynomial list for preparation for evaluation.
// Polynomails are evaluated in order.  This function is
// to determine this order
/////////////////////////////////////////////////////////
void polyListSorting(struct polynomial *p, struct polyList *l)
{
   int i;

   switch(p->eType){
           //If the polynomial is a constant, put it in the evaluation list
           case T_CONSTANT:
                polyListAppend(l,p);
                break;

           //If the polynomial is a variable, put it in the evaluation list
           case T_VARIABLE:
                if(p->valid!=1)
                {
                   polyListAppend(l,p);
                }
                break;

           //If the polynomial is a sum, put all the items of the sum in the evaluation list
           //except constants
           case T_SUM:

              if(p->valid==1)
                break;

              for(i=0;i<p->e.s->num;i++)
                  if(p->e.s->sum[i]->eType!=T_CONSTANT && p->e.s->sum[i]->valid!=1)
                  {
                         polyListSorting(p->e.s->sum[i],l);
                  }
              polyListAppend(l,p);
              break;

           //If the polynomial is a product, put all the items of the product in the 
           //evaluation list except constants
           case T_PRODUCT:

              if(p->valid==1)
                break;

              for(i=0;i<p->e.p->num;i++)
                  if(p->e.p->product[i]->eType!=T_CONSTANT && p->e.p->product[i]->valid!=1)
                  {
                         polyListSorting(p->e.p->product[i],l);
                  }
              polyListAppend(l,p);
              break;
 
           //If the polynomial is a functionCall, put the parameters in the evaluation list
           case T_FUNCTIONCALL:

                if(p->valid==1)
                break;


                for(i=0;i<p->e.f->paraNum;i++)
                  if(p->e.f->para[i]->eType!=T_CONSTANT && p->e.f->para[i]->valid!=1)
                  {
                         polyListSorting(p->e.f->para[i],l);
                  }
                polyListAppend(l,p);
                break;

           default:
              break;
    }  
//  for(i=0;i<l->listNext;i++)
//     if(l->pList[i]->eType==T_VARIABLE)
//       fprintf(stderr,"List element i=%d of %d l->pList[i].eType=%d  index=%d name=%s\n",
//                                  i,l->listNext,l->pList[i]->eType,l->pList[i]->index,l->pList[i]->vName);
//     else
//       fprintf(stderr,"List element i=%d of %d l->pList[i].eType=%d  index=%d \n",
//                                 i,l->listNext,l->pList[i]->eType,l->pList[i]->index);

};

void polyListSorting2(struct polynomial *p, struct polyList *l)
{
   int i;

   switch(p->eType){
           case T_CONSTANT:
                polyListAppend(l,p);
                break;

           case T_VARIABLE:
                polyListAppend(l,p);
                break;

           case T_SUM:
              for(i=0;i<p->e.s->num;i++)
                  if(p->e.s->sum[i]->eType!=T_CONSTANT)
                  {
                         polyListSorting2(p->e.s->sum[i],l);
                  }
              polyListAppend(l,p);
              break;

           case T_PRODUCT:
              for(i=0;i<p->e.p->num;i++)
                  if(p->e.p->product[i]->eType!=T_CONSTANT)
                  {
                         polyListSorting2(p->e.p->product[i],l);
                  }
              polyListAppend(l,p);
              break;

           case T_FUNCTIONCALL:
                for(i=0;i<p->e.f->paraNum;i++)
                  if(p->e.f->para[i]->eType!=T_CONSTANT )
                  {
                         polyListSorting(p->e.f->para[i],l);
                  }
                polyListAppend(l,p);
                break;

           default:
              break;
    }

};

/////////////////////////////////////////////////////////////
//This function compute the value of a polynomial.  It evaluate
//all the polynomials in the polynomial list for this polynomial.
//////////////////////////////////////////////////////////////
double evaluatePoly(struct polynomial *pp, struct polyList *l)
{
   struct polynomial *p;
   register int i,j;
   int num;
   double pV,re;
   int pE;
   //   double tempD1,tempD2;
   //   int tempI1,tempI2;
//  fprintf(stderr,"l->listNext=%d\n",l->listNext);
//   for(i=0;i<10;i++)
//   {
//      sum_Count[i]=0;
//      product_Count[i]=0;
//   }

   if(l->listNext==0)
   {
       return pp->value;
   }

   for(j=0; j<=l->listNext-1;j++)
   {
      p=l->pList[j];
//      expPrinting(p);
//      fprintf(stderr,"\n");
      switch(p->eType){
//        case T_CONSTANT:
//             break;
        //Read the value of the variable
        case T_VARIABLE:
             if(p->e.v->vType=='D')
                p->value=*(p->e.v->vAddrD);
             else if(p->e.v->vType=='I')
                p->value=*(p->e.v->vAddrI);
             else
             {
                fprintf(stderr,"Wrong variable type, exit!\n");
                exit(1);
             }
             break;

        //Sum up all the items in a sum
        case T_SUM:
           p->value=0;
           num=p->e.s->num;
//           sum_Count[num]++;
           for(i=0;i<num;i++)
               p->value += p->e.s->sum[i]->value*p->e.s->factor[i];
           break;

        //Product of all the items
        case T_PRODUCT:
           p->value=1;
           num=p->e.p->num;
           for(i=0;i<num;i++)
//             p->value*=pow(p->e.p->product[i]->value,p->e.p->exponent[i]);
           {

                 pV=p->e.p->product[i]->value;
                 pE=p->e.p->exponent[i];


/*
                 if(pE<=11 && pE>=2)
                 {
                    if(p->e.p->product[i]->count[pE-2]!=0)
                    {
                        p->value*=p->e.p->product[i]->values[pE-2];
                        continue;
                    }
                 }
*/


/*
                 if(pE>=1 && pE<=9)
                    p->e.p->product[i]->count[pE-1]++;
                 else if(pE>9)
                    p->e.p->product[i]->count[9]++;
*/
                 switch(pE){
                  case 1:
                        re=pV;
//                        fprintf(stderr,"Here,1, re=%10.8f \n",re);
                        break;
                  case 2:
                        re=pV*pV;
//                        fprintf(stderr,"Here,2, re=%10.8f \n",re);
                        break;
                  case 3:
                        re=pV*pV*pV;
//                        fprintf(stderr,"Here,3, re=%10.8f \n",re);
                        break;
                  case 4:
                        re=pV*pV;
                        re=re*re;
//                        fprintf(stderr,"Here,4, re=%10.8f \n",re);
                        break;
                  case 5:
                        re=pV*pV;
                        re=re*re;
                        re=re*pV;
//                        fprintf(stderr,"Here,5, re=%10.8f \n",re);
                        break;
                  case 6:
                        re=pV*pV*pV;
                        re=re*re;
//                        fprintf(stderr,"Here,6, re=%10.8f \n",re);
                        break;
                  default:
                        re = pow(pV,pE);
//                        fprintf(stderr,"Here, default re=%10.8f \n",re);
                        break;
                  }
                  p->value*=re;

/*
                  if(pE<=11 && pE>=2)
                  {
                    if(p->e.p->product[i]->count[pE-2]==0)
                    {
                        p->e.p->product[i]->count[pE-2]=1;
                        p->e.p->product[i]->values[pE-2]=re;
                    }
                  }
*/
           }
           break;
        
        //Function calls are evaluated by the library functions.
        case T_FUNCTIONCALL:
           if(strcmp(p->e.f->name,"log10")==0)
           {
              p->value=log10(p->e.f->para[0]->value);
           }
           else if(strcmp(p->e.f->name,"gsl_ran_tdist_pdf")==0)
           {
              p->value=gsl_ran_tdist_pdf(p->e.f->para[0]->value, p->e.f->para[1]->value);
           }
           else if(strcmp(p->e.f->name,"gsl_cdf_tdist_Q")==0)
           {
              p->value=gsl_cdf_tdist_Q(p->e.f->para[0]->value,p->e.f->para[1]->value);
           }   
           else if(strcmp(p->e.f->name,"gsl_cdf_tdist_P")==0)
           {
              p->value=gsl_cdf_tdist_P(p->e.f->para[0]->value,p->e.f->para[1]->value);
           }
           else if(strcmp(p->e.f->name,"gsl_ran_ugaussian_pdf")==0)
           {
              p->value=gsl_ran_ugaussian_pdf(p->e.f->para[0]->value);
           }
           else if(strcmp(p->e.f->name,"gsl_cdf_ugaussian_Q")==0)
           {
              p->value=gsl_cdf_ugaussian_Q(p->e.f->para[0]->value);
           }
           else if(strcmp(p->e.f->name,"gsl_cdf_ugaussian_P")==0)
           {
              p->value=gsl_cdf_ugaussian_P(p->e.f->para[0]->value);
           }
           else if(strcmp(p->e.f->name,"gsl_cdf_chisq_P")==0)
           {
              p->value=gsl_cdf_chisq_P(p->e.f->para[0]->value, p->e.f->para[1]->value);
           }
           else if(strcmp(p->e.f->name,"gsl_cdf_chisq_Q")==0)
           {
              p->value=gsl_cdf_chisq_Q(p->e.f->para[0]->value, p->e.f->para[1]->value);
           }
           else if(strcmp(p->e.f->name,"gsl_ran_chisq_pdf")==0)
           {
              p->value=gsl_ran_chisq_pdf(p->e.f->para[0]->value, p->e.f->para[1]->value);
           }
           else if(strcmp(p->e.f->name,"pow")==0)
           {
              p->value=pow(p->e.f->para[0]->value,p->e.f->para[1]->value);

           }
           else if(strcmp(p->e.f->name,"exp")==0)
           {
//              fprintf(stderr,"exp parameter:\n");
//              expPrinting(p->e.f->para[0]);
//              fprintf(stderr,"\n");
              p->value=exp(evaluateValue(p->e.f->para[0]));              
           }
           else if(strcmp(p->e.f->name,"sqrt")==0)
           {
              p->value=sqrt(p->e.f->para[0]->value);
           }
           else
           {
              fprintf(stderr,"unknown function name %s in polynomials\n",p->e.f->name);
              exit(1);
           }
           break;

        default:
           fprintf(stderr,"Error, unknown polynomial type!!!! exit(7)");
           break;
     }
   }
   return pp->value;
//   return l->pList[l->listNext-1]->value;
}
/////////////////////////////////////////////
//Print the polynomial with variales, constants,
//and symbols such as +,-,*,^,(, and ).
/////////////////////////////////////////////
void printPolyList(struct polyList *l)
{
   struct polynomial *p;
   register int i,j;

   for(j=l->listNext-1;j>=0;j--)
   {
      p=l->pList[j];
      
      expPrinting(p);
      fprintf(stderr,"\n");
      switch(p->eType){
        case T_CONSTANT:
             fprintf(stderr,"list No. %d type:CONSTANT   value=%f\n",
                            j,p->value);
             break;

        case T_VARIABLE:
             fprintf(stderr,"list No. %d type:VAR  name=%s  value=%f\n",
                            j,p->e.v->vName,p->value);
             break;

        case T_SUM:
           for(i=0;i<p->e.s->num;i++)
           {
             fprintf(stderr,"list No. %d Type:SUM  item %d factor %f item value%f  sum value %f \n",
                             j,i,p->e.s->factor[i],p->e.s->sum[i]->value, p->value);
           }
           break;

        case T_PRODUCT:
           for(i=0;i<p->e.p->num;i++)
           {
             fprintf(stderr,"list No. %d Type:Product  item %d item value %f exponent %d product value %f \n",
                             j,i,p->e.p->product[i]->value,p->e.p->exponent[i], p->value);
           }
           break;

        case T_FUNCTIONCALL:
             fprintf(stderr,"list No. %d type:CONSTANT CALL  name=%s  value=%f\n",
                            j,p->e.f->name,p->value);
             break;

        default:
           fprintf(stderr,"Error, unknown polynomial type!!!! exit(8)");
           break;
     }
   }
}

////////////////////////////////////////////////////////////////
// We need to allocate some memory and initiate some variables 
// to start polynomials.
//
////////////////////////////////////////////////////////////////
void polynomialInitialization()
{
   int i;

  

   sum0=0;
   sum1=0;
   sum2=0;
   sum3=0;
   sum4=0;
   sum5=0;
   sum00=0;
   sum11=0;

   product0=0;
   product1=0;
   product2=0;
   product3=0;
   product4=0;
   product5=0;
   product6=0;
   product7=0;
   product8=0;
   product9=0;

   product00=0;
   product11=0;



   nodeId=0;
 
   //allocate memory for each type of polynomials and set the counter of each
   //type of polynomials to be 0
   constantListLength = 1000;
   constantCount      = 0;
   constantList = (struct polynomial **)malloc(constantListLength*sizeof(struct polynomial *));

   variableListLength = 50;
   variableCount    = 0;
   variableList = (struct polynomial **)malloc(variableListLength*sizeof(struct polynomial *));

   sumListLength   =  10000;
   sumCount        =  0;
   sumList      = (struct polynomial **)malloc(sumListLength*sizeof(struct polynomial *));

   productListLength = 10000;
   productCount      = 0;
   productList  = (struct polynomial **)malloc(productListLength*sizeof(struct polynomial *));

   functionCallListLength = 10000;
   functionCallCount      = 0;
   functionCallList = (struct polynomial **)malloc(functionCallListLength*sizeof(struct polynomial *));

   //allocate memory for hash tables
   constantHash    = (struct hashStruct *)malloc(CONSTANT_HASH_SIZE     * sizeof(struct hashStruct));
   variableHash    = (struct hashStruct *)malloc(VARIABLE_HASH_SIZE     * sizeof(struct hashStruct));
   sumHash         = (struct hashStruct *)malloc(SUM_HASH_SIZE          * sizeof(struct hashStruct));
   productHash     = (struct hashStruct *)malloc(PRODUCT_HASH_SIZE      * sizeof(struct hashStruct));
   functionCallHash= (struct hashStruct *)malloc(FUNCTIONCALL_HASH_SIZE * sizeof(struct hashStruct));

   fprintf(stderr,"Sum hash size: %d product hash size: %d\n",(int)SUM_HASH_SIZE*(int)sizeof(struct hashStruct),(int)PRODUCT_HASH_SIZE*(int)sizeof(struct hashStruct));

   //Initialize the constant hash table
   for(i=0;i<CONSTANT_HASH_SIZE;i++)
   {
     constantHash[i].num=0;
     constantHash[i].index=NULL;
     constantHash[i].key=NULL;
   }

   //Initialize the variable hash table 
   for(i=0;i<VARIABLE_HASH_SIZE;i++)
   {
     variableHash[i].num=0;
     variableHash[i].index=NULL;
     variableHash[i].key=NULL;
   }

   //Initialize the sum hash table 
   for(i=0;i<SUM_HASH_SIZE;i++)
   {
     sumHash[i].num=0;
     sumHash[i].index=NULL;
     sumHash[i].key=NULL;
   }

   //Initialize the product hash table 
   for(i=0;i<PRODUCT_HASH_SIZE;i++)
   {
     productHash[i].num=0;
     productHash[i].index=NULL;
     productHash[i].key=NULL;
   }

   //Initialize the function call hash table 
   for(i=0;i<FUNCTIONCALL_HASH_SIZE;i++)
   {
     functionCallHash[i].num=0;
     functionCallHash[i].index=NULL;
     functionCallHash[i].key=NULL;
   }

   //Apply memory for containers to hold the sum items
   containerLength_v=100;
   factor_v = (double *)malloc(containerLength_v*sizeof(double));
   p_v      = (struct polynomial **)malloc(containerLength_v*sizeof(struct polynomial *));
   if(factor_v==NULL || p_v==NULL)
   {
        fprintf(stderr,"Momery allocation error!\n");
        exit(1);
   }
   containerLength_p=100;
   factor_p = (double *)malloc(containerLength_p*sizeof(double));
   p_p      = (struct polynomial **)malloc(containerLength_p*sizeof(struct polynomial *));
   if(factor_p==NULL || p_p==NULL)
   {
        fprintf(stderr,"Momery allocation error!\n");
        exit(1);
   }
   containerLength_f=100;
   factor_f = (double *)malloc(containerLength_f*sizeof(double));
   p_f      = (struct polynomial **)malloc(containerLength_f*sizeof(struct polynomial *));
   if(factor_f==NULL || p_f==NULL)
   {
        fprintf(stderr,"Momery allocation error!\n");
        exit(1);
   }



}

///////////////////////////////////////////////
//record the current status of the polynomials
///////////////////////////////////////////////
void makePolynomialStamp()
{
   
//           nodeIdStamp            = nodeId;
           constantCountStamp     = constantCount;
           variableCountStamp     = variableCount;
           sumCountStamp          = sumCount;
           productCountStamp      = productCount;
           functionCallCountStamp = functionCallCount;
                 
//           printAllPolynomials();  
}

//record the current status of the polynomials
void makePolynomialStamp2()
{

//           nodeIdStamp2            = nodeId;
           constantCountStamp2     = constantCount;
           variableCountStamp2     = variableCount;
           sumCountStamp2          = sumCount;
           productCountStamp2      = productCount;
           functionCallCountStamp2 = functionCallCount;
//           printAllPolynomials();
}


void partialPolynomialClearance()
{
   int i,j,k;
    
//   nodeId=nodeIdStamp;
                
  //partially clear constants
   for(j=constantCountStamp;j<constantCount;j++)
              free(constantList[j]);
   constantCount=constantCountStamp;

   //partially clear variables
   for(j=variableCountStamp;j<variableCount;j++)
   {
              free(variableList[j]->e.v);
              free(variableList[j]);
   }

   variableCount=variableCountStamp;
   //partially clear sums
   for(j=sumCountStamp;j<sumCount;j++)
   {
              free(sumList[j]->e.s->sum);
              free(sumList[j]->e.s->factor);
              free(sumList[j]->e.s);
              free(sumList[j]);
   }
   sumCount=sumCountStamp;

   //partially clear products
   for(j=productCountStamp;j<productCount;j++)
   {
              free(productList[j]->e.p->product);
              free(productList[j]->e.p->exponent);
              free(productList[j]->e.p);
              free(productList[j]);
   }
   productCount=productCountStamp;
                 
   //partially clear function-calls
   for(j=functionCallCountStamp;j<functionCallCount;j++)
   {
              free(functionCallList[j]->e.f->name);
              free(functionCallList[j]->e.f->para);
              free(functionCallList[j]->e.f);
              free(functionCallList[j]);
   }
   functionCallCount=functionCallCountStamp;

   //adjust the memory for constants
   constantListLength = constantCount+10;
   constantList = realloc(constantList, constantListLength*sizeof(struct polynomial *));

   //adjust the memory for variables
   variableListLength = variableCount+10;
   variableList = realloc(variableList, variableListLength*sizeof(struct polynomial *));
              
   //adjust the memory for sums
   sumListLength   =  sumCount+10;
   sumList      = realloc(sumList, sumListLength*sizeof(struct polynomial *));

   //adjust the memory for prodcuts
   productListLength = productCount+10;
   productList       = realloc(productList,productListLength*sizeof(struct polynomial *));
    
   //adjust the memory for functioncalls
   functionCallListLength = functionCallCount+10; 
   functionCallList = realloc(functionCallList,functionCallListLength*sizeof(struct polynomial *));

   //rearrange the hash table of constants
   for(j=0;j<CONSTANT_HASH_SIZE;j++)
   {
      if(constantHash[j].num>0)
      {
                k=0;
                for(i=0;i<constantHash[j].num;i++)
                {
                   if(constantHash[j].index[i]<constantCountStamp)
                   {
                        constantHash[j].index[k] = constantHash[j].index[i];
                        constantHash[j].key[k]   = constantHash[j].key[i];
                        k++;
                   }
                }
                constantHash[j].num=k;
                constantHash[j].key=realloc(constantHash[j].key,
                             sizeof(int)*constantHash[j].num);
                constantHash[j].index=realloc(constantHash[j].index,
                             sizeof(int)*constantHash[j].num);

      }
   }


   //rearrange the hash table of variables
   for(j=0;j<VARIABLE_HASH_SIZE;j++)
   {
      if(variableHash[j].num>0)
      {
                k=0;
                for(i=0;i<variableHash[j].num;i++)
                {
                   if(variableHash[j].index[i]<variableCountStamp)
                   {
                        variableHash[j].index[k] = variableHash[j].index[i];
                        variableHash[j].key[k]   = variableHash[j].key[i];
                        k++;
                   }
                }
                variableHash[j].num=k;
                variableHash[j].key=realloc(variableHash[j].key,
                             sizeof(int)*variableHash[j].num);
                variableHash[j].index=realloc(variableHash[j].index,
                             sizeof(int)*variableHash[j].num);

      }
   }

   //rearrange hash table of sums
   for(j=0;j<SUM_HASH_SIZE;j++)
   {
      if(sumHash[j].num>0)
      {
                k=0;
                for(i=0;i<sumHash[j].num;i++)
                {
                   if(sumHash[j].index[i]<sumCountStamp)
                   {
                        sumHash[j].index[k] = sumHash[j].index[i];
                        sumHash[j].key[k]   = sumHash[j].key[i];
                        k++;
                   }
                }
                sumHash[j].num=k;
                sumHash[j].key=realloc(sumHash[j].key,
                             sizeof(int)*sumHash[j].num);
                sumHash[j].index=realloc(sumHash[j].index,
                             sizeof(int)*sumHash[j].num);
      }
   }
   //rearrange hash table of products
   for(j=0;j<PRODUCT_HASH_SIZE;j++)
   {
       if(productHash[j].num>0)
       {
                k=0;
                for(i=0;i<productHash[j].num;i++)
                {
                   if(productHash[j].index[i]<productCountStamp)
                   {
                        productHash[j].index[k] = productHash[j].index[i];
                        productHash[j].key[k]   = productHash[j].key[i];
                        k++;
                   }
                }
                productHash[j].num=k;
                productHash[j].key=realloc(productHash[j].key,
                             sizeof(int)*productHash[j].num);
                productHash[j].index=realloc(productHash[j].index,
                             sizeof(int)*productHash[j].num);
       }
   }
   //rearrange hash table of function calls
   for(j=0;j<FUNCTIONCALL_HASH_SIZE;j++)
   {
       if(functionCallHash[j].num>0)
       {
                k=0;
                for(i=0;i<functionCallHash[j].num;i++)   
                {
                   if(functionCallHash[j].index[i]<functionCallCountStamp)
                   {
                        functionCallHash[j].index[k] = functionCallHash[j].index[i];
                        functionCallHash[j].key[k]   = functionCallHash[j].key[i];
                        k++;   
                   }
                }   
                functionCallHash[j].num=k;
                functionCallHash[j].key=realloc(functionCallHash[j].key,
                             sizeof(int)*functionCallHash[j].num);
                functionCallHash[j].index=realloc(functionCallHash[j].index,
                             sizeof(int)*functionCallHash[j].num);
       }
   }
   for(j=0;j<constantCount;j++)
       constantList[j]->valid=0;
   for(j=0;j<variableCount;j++)
       variableList[j]->valid=0;
   for(j=0;j<sumCount;j++)
       sumList[j]->valid=0;
   for(j=0;j<productCount;j++)
       productList[j]->valid=0;
   for(j=0;j<functionCallCount;j++)
       functionCallList[j]->valid=0;
                    
};

void partialPolynomialClearance2()
{
   int i,j,k;
   int *sumIndex;
   int *productIndex;
   int *functionCallIndex;
   int sumRecount;   
   int productRecount;
   int functionCallRecount;   

//   fprintf(stderr,"partialPolynomialClearance2 sumIndex          Size=%d\n",sizeof(int)*(sumCount-sumCountStamp2+1));
//   fprintf(stderr,"partialPolynomialClearance2 productIndex      size=%d\n",sizeof(int)*(productCount-productCountStamp2+1));
//   fprintf(stderr,"partialPolynomialClearance2 functionCallIndex size=%d\n",sizeof(int)*(functionCallCount-functionCallCountStamp2+1));

   //scanf("%d",&i);

   sumIndex          = (int *)malloc(sizeof(int)*(sumCount-sumCountStamp2+1));

//   fprintf(stderr,"sumIndex application done\n");
   if(sumIndex==NULL)
   {
       fprintf(stderr,"sumIndex memory application failed!\n");
       exit(1);
   }

      //scanf("%d",&i);

   productIndex      = (int *)malloc(sizeof(int)*(productCount-productCountStamp2+1));

//   fprintf(stderr,"productIndex application done\n");

   if(productIndex==NULL)
   {
       fprintf(stderr,"productIndex memory application failed!\n");
       exit(1);
   }

     //scanf("%d",&i);

   functionCallIndex = (int *)malloc(sizeof(int)*(functionCallCount-functionCallCountStamp2+1));
   
//   fprintf(stderr,"functionCallIndex application done\n");
   if(functionCallIndex==NULL)
   {
       fprintf(stderr,"functionCallIndex memory application failed!\n");
       exit(1);
   }

   //partially clear sums
   sumRecount=sumCountStamp2;
   for(j=sumCountStamp2;j<sumCount;j++)
   {
       if(sumList[j]->valid>=1)
       {
              sumList[sumRecount]=sumList[j];
              sumList[sumRecount]->index=sumRecount;
              sumIndex[j-sumCountStamp2]=sumRecount;
              sumRecount++;
       }
       else
       {
              sumIndex[j-sumCountStamp2]=-1;
              free(sumList[j]->e.s->sum);
              free(sumList[j]->e.s->factor);
              free(sumList[j]->e.s);
              free(sumList[j]);
       }
   }
   sumCount=sumRecount;

//   fprintf(stderr,"sumList clearance completed!\n");

   //partially clear products 
   productRecount=productCountStamp2;
   for(j=productCountStamp2;j<productCount;j++)    
   {
       if(productList[j]->valid>=1)
       {
              productList[productRecount]=productList[j];
              productList[productRecount]->index=productRecount;
              productIndex[j-productCountStamp2]=productRecount;
              productRecount++;
       }
       else 
       {
              productIndex[j-productCountStamp2]=-1;
              free(productList[j]->e.p->product);
              free(productList[j]->e.p->exponent);
              free(productList[j]->e.p);
              free(productList[j]);
        }
   }
   productCount=productRecount;

//   fprintf(stderr,"productList clearance completed!\n");

   //partially clear function-calls
   functionCallRecount=functionCallCountStamp2;
   for(j=functionCallCountStamp2;j<functionCallCount;j++)
   {
       if(functionCallList[j]->valid>=1)
       {
              functionCallList[functionCallRecount]=functionCallList[j];
              functionCallList[functionCallRecount]->index=functionCallRecount;
              functionCallIndex[j-functionCallCountStamp2]=functionCallRecount;
              functionCallRecount++;
       }
       else 
       {
              functionCallIndex[j-functionCallCountStamp2]=-1;
              free(functionCallList[j]->e.f->name);
              free(functionCallList[j]->e.f->para);
              free(functionCallList[j]->e.f);
              free(functionCallList[j]);
       }
   }       
   functionCallCount=functionCallRecount;

//   fprintf(stderr,"function List clearance completed!\n");

/*
   //adjust the memory for constants
   constantListLength = constantCount+10;
   constantList = realloc(constantList, constantListLength*sizeof(struct polynomial *));
   //adjust the memory for variables
   variableListLength = variableCount+10;
   variableList = realloc(variableList, variableListLength*sizeof(struct polynomial *));                 
   //adjust the memory for sums
   sumListLength   =  sumCount+10;
   sumList      = realloc(sumList, sumListLength*sizeof(struct polynomial *));                        
   //adjust the memory for prodcuts
   productListLength = productCount+10;
   productList       = realloc(productList,productListLength*sizeof(struct polynomial *));
   //adjust the memory for functioncalls
   functionCallListLength = functionCallCount+10;
   functionCallList = realloc(functionCallList,functionCallListLength*sizeof(struct polynomial *));
*/

   //rearrange hash table of sums
   for(j=0;j<SUM_HASH_SIZE;j++)
   {
      if(sumHash[j].num>0)
      {
                k=0;
                for(i=0;i<sumHash[j].num;i++)
                {
                   if(sumHash[j].index[i]<sumCountStamp2)
                   {
                        sumHash[j].index[k] = sumHash[j].index[i];
                        sumHash[j].key[k]   = sumHash[j].key[i];
                        sumList[sumHash[j].index[k]]->subHIndex=k;
                        k++;
                   }
                   else if(sumIndex[sumHash[j].index[i]-sumCountStamp2]!=-1)
                   {
                        sumHash[j].index[k] = sumIndex[sumHash[j].index[i]-sumCountStamp2];
                        sumHash[j].key[k]   = sumHash[j].key[i];
                        sumList[sumHash[j].index[k]]->subHIndex=k;
                        k++;                                           
                   }
                }
                sumHash[j].num=k;
                sumHash[j].key=realloc(sumHash[j].key,
                             sizeof(int)*sumHash[j].num);
                sumHash[j].index=realloc(sumHash[j].index,
                             sizeof(int)*sumHash[j].num);
      }
   }

 //  fprintf(stderr,"sumHash rearrangement completed!\n");

   //rearrange hash table of products
   for(j=0;j<PRODUCT_HASH_SIZE;j++)
   {
       if(productHash[j].num>0)
       {
                k=0;
                for(i=0;i<productHash[j].num;i++)
                {
                   if(productHash[j].index[i]<productCountStamp2)
                   {
                        productHash[j].index[k] = productHash[j].index[i];
                        productHash[j].key[k]   = productHash[j].key[i];
                        productList[productHash[j].index[k]]->subHIndex=k;
                        k++;
                   }
                   else if(productIndex[productHash[j].index[i]-productCountStamp2]!=-1)
                   {
                        productHash[j].index[k] = productIndex[productHash[j].index[i]-productCountStamp2];
                        productHash[j].key[k]   = productHash[j].key[i];
                        productList[productHash[j].index[k]]->subHIndex=k;
                        k++;
                   }
                }
                productHash[j].num=k;
                productHash[j].key=realloc(productHash[j].key,
                             sizeof(int)*productHash[j].num);
                productHash[j].index=realloc(productHash[j].index,
                             sizeof(int)*productHash[j].num);
       }
   }

//   fprintf(stderr,"productHash rearrangement completed!\n");


   //rearrange hash table of function calls
   for(j=0;j<FUNCTIONCALL_HASH_SIZE;j++)
   {
       if(functionCallHash[j].num>0)
       {
                k=0;
                for(i=0;i<functionCallHash[j].num;i++)
                {
                   if(functionCallHash[j].index[i]<functionCallCountStamp2)
                   {
                        functionCallHash[j].index[k] = functionCallHash[j].index[i];
                        functionCallHash[j].key[k]   = functionCallHash[j].key[i];
                        functionCallList[functionCallHash[j].index[k]]->subHIndex=k;
                        k++;
                   }
                   else if(functionCallIndex[functionCallHash[j].index[i]-functionCallCountStamp2]!=-1)
                   {
                        functionCallHash[j].index[k] = functionCallIndex[functionCallHash[j].index[i]-functionCallCountStamp2];
                        functionCallHash[j].key[k]   = functionCallHash[j].key[i];
                        functionCallList[functionCallHash[j].index[k]]->subHIndex=k;
                        k++;
                   }
                }
                functionCallHash[j].num=k;
                functionCallHash[j].key=realloc(functionCallHash[j].key,
                             sizeof(int)*functionCallHash[j].num);
                functionCallHash[j].index=realloc(functionCallHash[j].index,
                             sizeof(int)*functionCallHash[j].num);
       }
   }

//   fprintf(stderr,"functionCallHash rearrangement completed!\n");

   free(sumIndex);
   free(productIndex);
   free(functionCallIndex);
	
};


void polynomialClearance()
{
   int j;
   nodeId=0;
   for(j=0;j<constantCount;j++)
              free(constantList[j]);
   constantCount=0;
   for(j=0;j<variableCount;j++)
   {
              free(variableList[j]->e.v);
              free(variableList[j]);
   }
   variableCount=0;

   for(j=0;j<sumCount;j++)
   {
              free(sumList[j]->e.s->sum);
              free(sumList[j]->e.s->factor);
              free(sumList[j]->e.s);
              free(sumList[j]);
   }
   sumCount=0;
   for(j=0;j<productCount;j++)
   {
              free(productList[j]->e.p->product);
              free(productList[j]->e.p->exponent);
              free(productList[j]->e.p); 
              free(productList[j]);
   }
   productCount=0;

   for(j=0;j<functionCallCount;j++)
   {
              free(functionCallList[j]->e.f->name);
              free(functionCallList[j]->e.f->para);
              free(functionCallList[j]->e.f);
              free(functionCallList[j]);
   }
   functionCallCount=0;

   for(j=0;j<CONSTANT_HASH_SIZE;j++)
   {
      if(constantHash[j].num>0)
      {
                free(constantHash[j].index);
                free(constantHash[j].key);
      }
   }


   for(j=0;j<VARIABLE_HASH_SIZE;j++)
   {
      if(variableHash[j].num>0)
      {
                free(variableHash[j].index);
                free(variableHash[j].key);
      }          
   }


   for(j=0;j<SUM_HASH_SIZE;j++)
   {
      if(sumHash[j].num>0)
      {
                free(sumHash[j].index);
                free(sumHash[j].key);
      }
   }

   for(j=0;j<PRODUCT_HASH_SIZE;j++)
   {
       if(productHash[j].num>0)
       {
                free(productHash[j].index);
                free(productHash[j].key);
       }
   }

   for(j=0;j<FUNCTIONCALL_HASH_SIZE;j++)
   {
       if(functionCallHash[j].num>0)
       {
                free(functionCallHash[j].index);
                free(functionCallHash[j].key);
       }
   }

   free(constantList);
   free(variableList);
   free(sumList);
   free(productList);
   free(functionCallList);

   free(constantHash);
   free(variableHash);
   free(sumHash);
   free(productHash);
   free(functionCallHash);

   free(factor_v);
   free(p_v);
   free(factor_p);
   free(p_p);
   free(factor_f);
   free(p_f);
   fprintf(stderr,"containerLength=[%d %d %d]\n",containerLength_v,containerLength_p,containerLength_f);

}



void dismantle()
{ 
   int i;
   for(i=0;i<constantCount;i++)
     fprintf(stderr,"Constant index=%d  value=%10.8f key=%d\n",constantList[i]->index,constantList[i]->value,constantList[i]->key);
   for(i=0;i<variableCount;i++)
     fprintf(stderr,"Variable index=%d  value=%10.8f key=%d\n",variableList[i]->index,variableList[i]->value,variableList[i]->key);
   for(i=0;i<sumCount;i++)
   {
     fprintf(stderr,"Sum      index=%d  value=%10.8f key=%d\n",sumList[i]->index,sumList[i]->value,sumList[i]->key);
/*
     for(j=0;j<sumList[i]->e.s->num;j++)
     {
        fprintf(stderr,"   %d  %ld  %10.8f  %10.8f\n",
              j,sumList[i]->e.s->sum[j]->index,sumList[i]->e.s->sum[j]->value, sumList[i]->e.s->factor[j]);
     }
*/
   }
   for(i=0;i<productCount;i++)
   {
     fprintf(stderr,"Product  index=%d  value=%10.8f key=%d\n",productList[i]->index,productList[i]->value,productList[i]->key);
/*
     for(j=0;j<productList[i]->e.p->num;j++)
     {
        fprintf(stderr,"   %d  %ld  %10.8f  %d\n",
                j,productList[i]->e.p->product[j]->index,
                productList[i]->e.p->product[j]->value, 
                productList[i]->e.p->exponent[j]);
     }
*/
   }
   for(i=0;i<functionCallCount;i++)
   {
     fprintf(stderr,"Function Call  id=%d  name: %s  value=%10.8f key=%d\n",
             functionCallList[i]->index,functionCallList[i]->e.f->name,functionCallList[i]->value,functionCallList[i]->key);
/*
     for(j=0;j<functionCallList[i]->e.f->paraNum;j++)
     {
        fprintf(stderr,"   %d  %ld  %10.8f \n",
                j,functionCallList[i]->e.f->para[j]->index,
                functionCallList[i]->e.f->para[j]->value);
     }
*/
   }
};

void polyStatistics(FILE *fp)
{ 
   long constantSize,variableSize,sumSize,productSize,functionCallSize;
   int  constantHashCount = 0, variableHashCount = 0, sumHashCount = 0, productHashCount = 0, functionHashCount = 0;
   int  constantHashSize  = 0, variableHashSize  = 0, sumHashSize  = 0, productHashSize  = 0, functionHashSize  = 0;
   int i;

   fprintf(stderr,"Number of constants:       %d\n",constantCount);
   fprintf(stderr,"Number of variables:       %d\n",variableCount);
   fprintf(stderr,"Number of sums:            %d\n",sumCount);
   fprintf(stderr,"Number of products:        %d\n",productCount);
   fprintf(stderr,"Number of function Calls:  %d\n",functionCallCount);

   fprintf(stderr,"sum0=%d sum1=%d sum2=%d sum3=%d sum4=%d sum5=%d sum00=%d sum11=%d\n",
                   sum0,sum1,sum2,sum3,sum4,sum5,sum00,sum11);

   fprintf(stderr,"product0=%d product1=%d product2=%d product3=%d product4=%d product5=%d product6=%d product7=%d product8=%d product9=%d product00=%d product11=%d \n", 
                    product0, product1, product2, product3,product4,product5,product6,product7,product8,product9, product00,product11);
   if(fp!=NULL)
   fprintf(fp," %5d   %5d   %5d   %5d  %5d",
	   constantCount,variableCount,sumCount,productCount,functionCallCount);


   constantSize=constantCount*sizeof(Polynomial);
   variableSize=variableCount*sizeof(Polynomial);
   sumSize=sumCount*sizeof(Polynomial)+sumCount*sizeof(struct sumPoly);
   for(i=0;i<sumCount;i++)
   {     
     sumSize+=sumList[i]->e.s->num*(sizeof(Polynomial *)+sizeof(double));
   }
   productSize=productCount*sizeof(Polynomial)+productCount*sizeof(struct productPoly);
   for(i=0;i<productCount;i++)
   {
      productSize+=productList[i]->e.p->num*(sizeof(Polynomial *)+sizeof(int));
   }
   functionCallSize=functionCallCount*sizeof(Polynomial)+functionCallCount*sizeof(struct functionPoly);
   for(i=0;i<functionCallCount;i++)
   {  
      functionCallSize+=strlen(functionCallList[i]->e.f->name)+sizeof(int)+sizeof(Polynomial *);
   }
   fprintf(stderr,"size of constants:       %ld\n",constantSize);
   fprintf(stderr,"size of variables:       %ld\n",variableSize);
   fprintf(stderr,"size of sums:            %ld\n",sumSize);
   fprintf(stderr,"size of products:        %ld\n",productSize);
   fprintf(stderr,"size of function Calls:  %ld\n",functionCallSize);


   if(fp!=NULL)
   fprintf(fp," %10ld   %10ld   %10ld   %10ld   %10ld   %10ld\n",constantSize,variableSize,sumSize,productSize,functionCallSize,
                       constantSize+variableSize+sumSize+productSize+functionCallSize);

/*
   fprintf(stderr,"Hash table usage for constants:\n");
   for(i=0;i<CONSTANT_HASH_SIZE;i++)
   {
      if(constantHash[i].num>0)
      {
         fprintf(stderr,"%d  ",constantHash[i].num);
         constantHashCount++;
         constantHashSize+=constantHash[i].num*sizeof(int)*2;
      }
   }
   fprintf(stderr,"\n");

   fprintf(stderr,"Hash table usage for variables:\n");
   for(i=0;i<VARIABLE_HASH_SIZE;i++)
   {
      if(variableHash[i].num>0)
      {
         fprintf(stderr,"%d  ",variableHash[i].num);
         variableHashCount++;
         variableHashSize+=variableHash[i].num*sizeof(int)*2;
      }
   }
   fprintf(stderr,"\n");


   fprintf(stderr,"Hash table usage for sums:\n");
   for(i=0;i<SUM_HASH_SIZE;i++)
   {
      if(sumHash[i].num>0)
      {
        fprintf(stderr,"%d  ",sumHash[i].num);
        sumHashCount++;
        sumHashSize+=sumHash[i].num*sizeof(int)*2;
      }
   }
   fprintf(stderr,"\n");

   fprintf(stderr,"Hash table usage for products:\n");
   for(i=0;i<PRODUCT_HASH_SIZE;i++)
   {
      if(productHash[i].num>0)
      {
        fprintf(stderr,"%d  ",productHash[i].num);
        productHashCount++;
        productHashSize+=productHash[i].num*sizeof(int)*2;
      }
   }
   fprintf(stderr,"\n");

   fprintf(stderr,"Hash table usage for function calls:\n");
   for(i=0;i<FUNCTIONCALL_HASH_SIZE;i++)
   {
      if(functionCallHash[i].num>0)
      {
        fprintf(stderr,"%d  ",functionCallHash[i].num);
        functionHashCount++;
        functionHashSize+=functionCallHash[i].num*sizeof(int)*2;
      }
   }
   fprintf(stderr,"\n");
*/

   fprintf(stderr,"constantHashSize: %d   variableHashSize: %d   sumHashSize: %d    productHashSize: %d    functionHashSize: %d\n",
                   constantHashSize, variableHashSize, sumHashSize, productHashSize,functionHashSize);

   fprintf(stderr,"constantHashCount:[%d  %f]  variableHashCount:[%d  %f] sumHashCount:[%d  %f]  productHashCount:[%d  %f]  functionHashCount:[%d  %f]\n",
                   constantHashCount, constantHashCount*1.0/CONSTANT_HASH_SIZE,
                   variableHashCount, variableHashCount*1.0/VARIABLE_HASH_SIZE,
                   sumHashCount,  sumHashCount*1.0/SUM_HASH_SIZE,
                   productHashCount,productHashCount*1.0/PRODUCT_HASH_SIZE,
                   functionHashCount,functionHashCount*1.0/FUNCTIONCALL_HASH_SIZE);
   if(fp!=NULL)
   {
       fprintf(fp,"constantHashSize: %d   variableHashSize: %d   sumHashSize: %d    productHashSize: %d    functionHashSize: %d\n",
                   constantHashSize, variableHashSize, sumHashSize, productHashSize,functionHashSize);

       fprintf(fp,"constantHashCount:[%d  %f]  variableHashCount:[%d  %f] sumHashCount:[%d  %f]  productHashCount:[%d  %f]  functionHashCount:[%d  %f]\n",
                   constantHashCount, constantHashCount*1.0/CONSTANT_HASH_SIZE,
                   variableHashCount, variableHashCount*1.0/VARIABLE_HASH_SIZE,
                   sumHashCount,  sumHashCount*1.0/SUM_HASH_SIZE,
                   productHashCount,productHashCount*1.0/PRODUCT_HASH_SIZE,
                   functionHashCount,functionHashCount*1.0/FUNCTIONCALL_HASH_SIZE);
   }

/*
   fprintf(stderr,"variable's exponents' refer times");
   for(i=0;i<variableCount;i++)
   {
      if(variableList[i]->valid==1)
      {
         fprintf(stderr,"variable %d: name: %s  ",i,variableList[i]->e.v->vName);
         for(j=0;j<10;j++)
            if(variableList[i]->count[j]!=0)
              fprintf(stderr,"(%d %d) ", j+1,(int)variableList[i]->count[j]);
         fprintf(stderr,"\n");
      }
   }
   for(i=0;i<sumCount;i++)
   {
      if(sumList[i]->valid==1)
      {
         fprintf(stderr,"sum %d: ",i);
         for(j=0;j<10;j++)
           if(sumList[i]->count[j]!=0)
              fprintf(stderr,"(%d %d) ", j+1,sumList[i]->count[j]);
         fprintf(stderr,"\n");
      }
   }
*/

/*
   for(i=0;i<productCount;i++)
   {    
      if(productList[i]->count[19]!=0) 
      {
         fprintf(stderr,"product %d: ",i);
         for(j=0;j<20;j++)
            fprintf(stderr,"(%d %d) ", j,productList[i]->count[j]);
         fprintf(stderr,"\n");
      }
   }
 */


}

//print out all variables and their values
void printAllVariables()
{
   int i;

   fprintf(stderr,"All the variables:\n");
   for(i=0;i<variableCount;i++)
   {
        fprintf(stderr,"index=%d variable: ",i);
        expPrinting(variableList[i]);

         if(variableList[i]->e.v->vType=='D')
         {
            fprintf(stderr,"  value: %f  %f\n", *(variableList[i]->e.v->vAddrD),variableList[i]->value);
         }
         else if(variableList[i]->e.v->vType=='I')
         {
            fprintf(stderr,"  value: %d  %f\n", *(variableList[i]->e.v->vAddrI),variableList[i]->value);
         }
   }
   fprintf(stderr,"\n");


}



//Print all the polynomials
void printAllPolynomials()
{
   int i,j;
   fprintf(stderr,"All the constants:\n");
   for(i=0;i<CONSTANT_HASH_SIZE;i++)
   {  
      if(constantHash[i].num<=0)
         continue;
      for(j=0;j<constantHash[i].num;j++)
      {
        fprintf(stderr,"Key=%d index=%d constant: ",constantHash[i].key[j],constantHash[i].index[j]);
        expPrinting(constantList[constantHash[i].index[j]]);
        fprintf(stderr,"\n");
      }
      fprintf(stderr,"\n");
   }
   fprintf(stderr,"\n");
   fprintf(stderr,"All the variables:\n");
   for(i=0;i<variableCount;i++)
   {
        fprintf(stderr,"index=%d variable: ",i);                                                
        expPrinting(variableList[i]);
        fprintf(stderr,"\n");
   }
   fprintf(stderr,"\n");      
   fprintf(stderr,"All the sums:\n");
   for(i=0;i<SUM_HASH_SIZE;i++)
   {  
      if(sumHash[i].num<=0)
         continue;
      for(j=0;j<sumHash[i].num;j++)
      {
        fprintf(stderr,"Key=%d index=%d sum: ",sumHash[i].key[j],sumHash[i].index[j]);
        expPrinting(sumList[sumHash[i].index[j]]);
        fprintf(stderr,"\n");
      }
      fprintf(stderr,"\n");
   }
   fprintf(stderr,"\n");
   fprintf(stderr,"All the products:\n");
   for(i=0;i<PRODUCT_HASH_SIZE;i++)
   {  
      if(productHash[i].num<=0)
         continue;
      for(j=0;j<productHash[i].num;j++)
      {
        fprintf(stderr,"Key=%d index=%d product: ",productHash[i].key[j],productHash[i].index[j]);
        expPrinting(productList[productHash[i].index[j]]);
        fprintf(stderr,"\n");
      }
      fprintf(stderr,"\n");
   }
   fprintf(stderr,"\n");

   
   fprintf(stderr,"All function calls:\n");
   for(i=0;i<FUNCTIONCALL_HASH_SIZE;i++)
   {
      if(functionCallHash[i].num<=0)
         continue;
      for(j=0;j<functionCallHash[i].num;j++)
      {
        fprintf(stderr,"Key=%d index=%d functionCall: ",functionCallHash[i].key[j],functionCallHash[i].index[j]);
        expPrinting(functionCallList[functionCallHash[i].index[j]]);
        fprintf(stderr,"\n");
      }
      fprintf(stderr,"\n");
   }
   fprintf(stderr,"\n");   
}


void polyEvaluationInitialization()
{
/*
   int i;
   for(i=0;i<sumCount;i++)
   {
     if(sumList[i]->valid==1)
     {
        memset(sumList[i]->count,0,10*sizeof(int));
     }
   }
   for(i=0;i<variableCount;i++)
   {    
     if(variableList[i]->valid==1)
     {  
        memset(variableList[i]->count,0,10*sizeof(int));
     }
   }
*/

}


// print out a polynomial and its sorting list
void dismantlePolynomialAndSortingList(struct polynomial *p, struct polyList *l)
{   
 int j;
 fprintf(stderr,"Polynomial Value:  %e\n",p->value);
 expPrinting(p);
 fprintf(stderr,"\n");
              
 for(j=0; j<=l->listNext-1;j++)
 {
      fprintf(stderr,"%4d  value: %e ",j, l->pList[j]->value);
      expPrinting(l->pList[j]);
      fprintf(stderr,"\n");
 }
};


void polyListAppend3(struct polyList *l, struct polynomial *p, int signature)
{
    //valid is a mark showing that this polynomial appears on a sorting list
    p->valid=signature;

//    expPrinting(p);
//    fprintf(stderr,"\n");

    /*******************By Honglong on Aug 8, 2006
    if(p->eType==T_VARIABLE)
    {         
      if(p->id==0)
         return;
      else
         p->id=0;  
    }
    **************************************************/
                
    if(l->listNext>=l->listSize)
    {         
        l->pList=realloc(l->pList,sizeof(struct polynomial *)*(l->listSize+100));
        l->listSize=l->listSize+100;
    }
    l->pList[l->listNext]=p;
    l->listNext++;
};


void polyListStatistics(struct polyList *l)
{
    int j;
    int numConstant=0,numVariable=0,numSum=0, numProduct=0,numFunctionCall=0;
    for(j=0; j<=l->listNext-1;j++)
    {
       switch(l->pList[j]->eType){
         case T_CONSTANT:
              numConstant++;
              break;
         case T_VARIABLE:
              numVariable++;
              break;
         case T_SUM:
              numSum++;
              break;
         case T_PRODUCT:
              numProduct++;
              break;
         case T_FUNCTIONCALL:
              numFunctionCall++;
              break;
         default:
              fprintf(stderr,"Unknown expression type, exit!\n");
              exit(1);
       }
   }
   fprintf(stderr,"%d  %d  %d  %d  %d\n",numConstant,numVariable,numSum, numProduct,numFunctionCall);
};


void polyListSorting3(struct polynomial *p, struct polyList *l, int signature)
{    
   int i;
 
//   fprintf(stderr,"sorting3 valid=%d signature=%d \n",p->valid,signature);
         
   switch(p->eType){
           case T_CONSTANT:
                polyListAppend3(l,p,signature);
                break;
     
           case T_VARIABLE:
                if(p->valid!=signature)
                {
                   polyListAppend3(l,p,signature);
                }
                break;
          
           case T_SUM:
              if(p->valid==signature)
                break;
              
              for(i=0;i<p->e.s->num;i++)
                  if(p->e.s->sum[i]->eType!=T_CONSTANT && p->e.s->sum[i]->valid!=signature)
                  {
                         polyListSorting3(p->e.s->sum[i],l,signature);
                  }
              polyListAppend3(l,p,signature);
              break;
                   
           case T_PRODUCT:

              if(p->valid==signature)
                break;
                
              for(i=0;i<p->e.p->num;i++)
                  if(p->e.p->product[i]->eType!=T_CONSTANT && p->e.p->product[i]->valid!=signature)
                  {
                         polyListSorting3(p->e.p->product[i],l,signature);
                  }
              polyListAppend3(l,p,signature);
              break;
           
           case T_FUNCTIONCALL:
     
                if(p->valid==signature)
                break;

                for(i=0;i<p->e.f->paraNum;i++)
                  if(p->e.f->para[i]->eType!=T_CONSTANT && p->e.f->para[i]->valid!=signature)
                  {
                         polyListSorting3(p->e.f->para[i],l,signature);
                  }
                polyListAppend3(l,p,signature);
                break;
                   
           default: 
              break;
    }
//  for(i=0;i<l->listNext;i++)
//     if(l->pList[i]->eType==T_VARIABLE)
//       fprintf(stderr,"List element i=%d of %d l->pList[i].eType=%d  id=%d name=%s\n",
//                                  i,l->listNext,l->pList[i]->eType,l->pList[i]->id,l->pList[i]->vName);
//     else
//       fprintf(stderr,"List element i=%d of %d l->pList[i].eType=%d  id=%d \n",
//                                 i,l->listNext,l->pList[i]->eType,l->pList[i]->id);
                
};                       



