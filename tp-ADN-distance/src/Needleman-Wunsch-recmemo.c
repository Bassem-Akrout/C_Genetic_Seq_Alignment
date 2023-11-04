/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 *
 * Documentation: see Needleman-Wunsch-recmemo.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in Needleman-Wunsch-recmemo.h
 */


#include "Needleman-Wunsch-recmemo.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <stdint.h>
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */
//#include "DequeNode.h" /* double-ended queue struct*/

/*****************************************************************************/
   
/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L 

/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm 
*/
struct NW_MemoContext 
{
    char *X ; /*!< the longest genetic sequences */
    char *Y ; /*!< the shortest genetic sequences */
    size_t M; /*!< length of X */
    size_t N; /*!< length of Y,  N <= M */
    long **memo; /*!< memoization table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
} ;


// Fonction pour créer une nouvelle deque vide
queue* createFile() {
    queue *deque = (queue *)malloc(sizeof(queue)); // Allocation de mémoire pour la deque
    if (deque == NULL) {
        fprintf(stderr, "Erreur d'allocation de mémoire\n"); // En cas d'erreur, affiche un message d'erreur et quitte le programme
        exit(1);
    }
    deque->front = NULL; // Initialise l'extrémité avant de la deque à NULL (vide)
    deque->rear = NULL;  // Initialise l'extrémité arrière de la deque à NULL (vide)
    return deque;
}

queue* copyQueue(queue *originalDeque) {
    queue *copiedDeque = createFile(); // Crée une nouvelle deque vide pour la copie
    queueNode *currentNode = originalDeque->front;

    while (currentNode != NULL) {
        pushBack(copiedDeque, currentNode->value); // Ajoute chaque élément de la deque originale à la copie
        currentNode = currentNode->next;
    }

    return copiedDeque;
}
// Fonction pour ajouter un élément à l'extrémité arrière de la deque
void pushBack(queue *deque, long value) {
    queueNode *newNode = (queueNode *)malloc(sizeof(queueNode)); // Alloue de la mémoire pour un nouveau nœud
    if (newNode == NULL) {
        fprintf(stderr, "Erreur d'allocation de mémoire\n");
        exit(1);
    }
    newNode->value = value; // Attribue la valeur à ce nouveau nœud
    newNode->next = NULL;   // Initialise le pointeur vers le nœud suivant à NULL
    newNode->prev = deque->rear; // Le nœud précédent est l'actuelle extrémité arrière de la deque

    // Met à jour les pointeurs de l'extrémité arrière de la deque
    if (deque->rear != NULL) {
        deque->rear->next = newNode;
    }
    deque->rear = newNode;

    if (deque->front == NULL) {
        deque->front = newNode;
    }
}

// Fonction pour retirer un élément de l'extrémité avant de la deque
long popFront(queue *deque) {
    if (deque->front == NULL) {
        fprintf(stderr, "Le deque est vide\n");
        exit(1);
    }
    queueNode *node = deque->front;
    long value = node->value; // Récupère la valeur du nœud à l'extrémité avant
    deque->front = node->next; // Met à jour l'extrémité avant de la deque

    if (deque->front != NULL) {
        deque->front->prev = NULL; // Met à jour le pointeur vers le nœud précédent du nouveau nœud à l'extrémité avant
    } else {
        deque->rear = NULL; // Si la deque est vide après le retrait, l'extrémité arrière doit également être NULL
    }

    free(node); // Libère la mémoire du nœud retiré
    return value;
}
/*
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoization \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoization array 
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ] 
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ] 
 */ 
static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
/* compute and returns phi(i,j) using data in c -allocated and initialized by EditDistance_NW_Rec */
{
   if (c->memo[i][j] == NOT_YET_COMPUTED)
   {  
      long res ;
      char Xi = c->X[i] ;
      char Yj = c->Y[j] ;
      if (i == c->M) /* Reach end of X */
      {  if (j == c->N) res = 0;  /* Reach end of Y too */
         else res = (isBase(Yj) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else if (j == c->N) /* Reach end of Y but not end of X */
      {  res = (isBase(Xi) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
      {  ManageBaseError( Xi ) ;
         res = EditDistance_NW_RecMemo(c, i+1, j) ;
      }
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res = EditDistance_NW_RecMemo(c, i, j+1) ;
      }
      else  
      {  /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */ 
         long min = /* initialization  with cas 1*/
                   ( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                   )
                   + EditDistance_NW_RecMemo(c, i+1, j+1) ; 
         { long cas2 = INSERTION_COST + EditDistance_NW_RecMemo(c, i+1, j) ;      
           if (cas2 < min) min = cas2 ;
         }
         { long cas3 = INSERTION_COST + EditDistance_NW_RecMemo(c, i, j+1) ;      
           if (cas3 < min) min = cas3 ; 
         }
         res = min ;
      }
       c->memo[i][j] = res ;
   }
   return c->memo[i][j] ;
}

/* EditDistance_NW_Rec :  is the main function to call, cf .h for specification 
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the 
 * recursivefunction EditDistance_NW_RecMemo 
 * See .h file for documentation
 */
long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB)
{
   _init_base_match() ;
   struct NW_MemoContext ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;
   {  /* Allocation and initialization of ctx.memo to NOT_YET_COMPUTED*/
      /* Note: memo is of size (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements 
       * It would have been possible to allocate only one big array memezone of (M+1)*(N+1) elements 
       * and then memo as an array of (M+1) pointers, the memo[i] being the address of memzone[i*(N+1)].
       */ 
      ctx.memo = (long **) malloc ( (M+1) * sizeof(long *)) ;
      if (ctx.memo == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo" ); exit(EXIT_FAILURE); }
      for (int i=0; i <= M; ++i) 
      {  ctx.memo[i] = (long*) malloc( (N+1) * sizeof(long));
         if (ctx.memo[i] == NULL) { perror("EditDistance_NW_Rec: malloc of ctx_memo[i]" ); exit(EXIT_FAILURE); }
         for (int j=0; j<=N; ++j) ctx.memo[i][j] = NOT_YET_COMPUTED ;
      }   }    
   /* Compute phi(0,0) = ctx.memo[0][0] by calling the recursive function EditDistance_NW_RecMemo */
   long res = EditDistance_NW_RecMemo( &ctx, 0, 0 ) ;
   { /* Deallocation of ctx.memo */
      for (int i=0; i <= M; ++i) free( ctx.memo[i] ) ;
      free( ctx.memo ) ;
   }
   return res ;
}

long iteratif(char* A, size_t lengthA, char* B, size_t lengthB) {

   _init_base_match() ;
   struct NW_MemoContext ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;
   queue *deque = createFile(); // Crée une nouvelle deque vide

   // Initialisation de la derniere ligne de Ø avec les valeurs appropriées
   pushBack(deque, (long)0);

   for (int64_t j = N-1; j >= 0; j--) {
      char Yj = ctx.Y[j] ;      
      
      pushBack(deque, (isBase(Yj) ? INSERTION_COST : 0) + deque->rear->value);
      
   }
   //remplissage de Ø en ne gardant en mémoire que la derniere ligne 
   for (int64_t i = M - 1; i >= 0; i--) {
      char Xi = ctx.X[i] ;
      pushBack(deque, (isBase(Xi) ? INSERTION_COST : 0) + deque->front->value);
         for (int64_t j = N-1; j >= 0; j--) {
            char Yj = ctx.Y[j] ;
            long res ;
            if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
               {  ManageBaseError( Xi ) ;
                  res = deque->front->next->value ;
               }
            else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
            {  ManageBaseError( Yj ) ;
               res =  deque->rear->value ;
            }
            else{  
               long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                   ) + deque->front->value ;                             //i+1, j+1
               long cas2 = INSERTION_COST + deque->front->next->value ; //i+1, j
               if (cas2 < min) min = cas2 ;
               long cas3 = INSERTION_COST + deque->rear->value ;        //i  , j+1
               if (cas3 < min) min = cas3 ;

               res=min;
               } 
            pushBack(deque, res);
            popFront(deque);
         }
      popFront(deque);
   }
   long result=deque->rear->value ;
   // Libération de la mémoire de la deque
   while (deque->front != NULL) {
      popFront(deque); // Retire les éléments restants de la deque et libère leur mémoire
   }
   free(deque); // Libère la mémoire de la structure de la deque

   return result;
}
int64_t max(int64_t a, int64_t b) {
    return (a > b) ? a : b;
}

// Fonction pour créer une nouvelle deque vide
queueBlock* createFileBlock() {
    queueBlock * deque = (queueBlock *)malloc(sizeof(queueBlock)); // Allocation de mémoire pour la deque
    if (deque == NULL) {
        fprintf(stderr, "Erreur d'allocation de mémoire\n"); // En cas d'erreur, affiche un message d'erreur et quitte le programme
        exit(1);
    }
    deque->front = NULL; // Initialise l'extrémité avant de la deque à NULL (vide)
    deque->rear = NULL;  // Initialise l'extrémité arrière de la deque à NULL (vide)
    return deque;
}

// Fonction pour ajouter un élément à l'extrémité arrière de la deque
void pushBackBlock(queueBlock *deque, queue value) {
    Block *newNode = (Block *)malloc(sizeof(Block)); // Alloue de la mémoire pour un nouveau nœud
    if (newNode == NULL) {
        fprintf(stderr, "Erreur d'allocation de mémoire\n");
        exit(1);
    }
    newNode->value = value; // Attribue la valeur à ce nouveau nœud
    newNode->next = NULL;   // Initialise le pointeur vers le nœud suivant à NULL
    newNode->prev = deque->rear; // Le nœud précédent est l'actuelle extrémité arrière de la deque

    // Met à jour les pointeurs de l'extrémité arrière de la deque
    if (deque->rear != NULL) {
        deque->rear->next = newNode;
    }
    deque->rear = newNode;

    if (deque->front == NULL) {
        deque->front = newNode;
    }
}
// Fonction pour retirer un élément de l'extrémité avant de la deque
queue popFrontBlock(queueBlock *deque) {
    if (deque->front == NULL) {
        fprintf(stderr, "Le deque est vide\n");
        exit(1);
    }
    Block *node = deque->front;
    queue value = node->value; // Récupère la valeur du nœud à l'extrémité avant
    deque->front = node->next; // Met à jour l'extrémité avant de la deque

    if (deque->front != NULL) {
        deque->front->prev = NULL; // Met à jour le pointeur vers le nœud précédent du nouveau nœud à l'extrémité avant
    } else {
        deque->rear = NULL; // Si la deque est vide après le retrait, l'extrémité arrière doit également être NULL
    }

    free(node); // Libère la mémoire du nœud retiré
    return value;
}

queue*  block_down_right(size_t N,size_t M,struct NW_MemoContext ctx,queueBlock* listBlock){
   queue *deque = createFile(); // Crée une nouvelle deque vide
   queue *tempo_deque = createFile();//utilisé par le block d'à gauche de celui là
   // Initialisation de la derniere ligne de Ø avec les valeurs appropriées
   pushBack(deque, (long)0);
   for (int64_t j = N-1; j >  max(N-K,-1); j--) {
      char Yj = ctx.Y[j] ;      
      pushBack(deque, (isBase(Yj) ? INSERTION_COST : 0) + deque->rear->value);
      
   }
   //remplissage de Ø en ne gardant en mémoire que la derniere ligne 
   for (int64_t i = M - 1; i >  max(M-K,-1); i--) {
      char Xi = ctx.X[i] ;
      pushBack(deque, (isBase(Xi) ? INSERTION_COST : 0) + deque->front->value);
         for (int64_t j = N-1; j > max(N-K,-1); j--) {
            char Yj = ctx.Y[j] ;
            long res ;
            if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
               {  ManageBaseError( Xi ) ;
                  res = deque->front->next->value ;
               }
            else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
            {  ManageBaseError( Yj ) ;
               res =  deque->rear->value ;
            }
            else{  
               long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                          : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                   ) + deque->front->value ;                             //i+1, j+1
               long cas2 = INSERTION_COST + deque->front->next->value ; //i+1, j
               if (cas2 < min) min = cas2 ;
               long cas3 = INSERTION_COST + deque->rear->value ;        //i  , j+1
               if (cas3 < min) min = cas3 ;

               res=min;
               } 
            pushBack(deque, res);
            popFront(deque);
         }
      long elt= popFront(deque);
      pushBack(tempo_deque,elt);
   }

   pushBack(tempo_deque,deque->rear->value);
   pushBackBlock(listBlock,*deque);
   return tempo_deque;
}

queue*  block_down(int J,queue* tempo_deque,size_t N,size_t M,struct NW_MemoContext ctx,queueBlock* listBlock){
   
   queue *deque = createFile(); // Crée une nouvelle deque vide
   // Initialisation de la derniere ligne de Ø avec les valeurs appropriées
   char YJ = ctx.Y[J] ;
   pushBack(deque,(isBase(YJ) ? INSERTION_COST : 0) + tempo_deque->front->value );
   for (int64_t j = J-1; j > max(J-K,-1); j--) {
      char Yj = ctx.Y[j] ;      
      pushBack(deque, (isBase(Yj) ? INSERTION_COST : 0) + deque->rear->value);
   }

   for (int64_t i = M - 1; i > max(M-K,-1); i--) {

      char Yj = ctx.Y[J] ;
      long res ;
      char Xi = ctx.X[i] ;
      if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
         {  ManageBaseError( Xi ) ;
            res = deque->front->value ;
         }
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res =  tempo_deque->front->next->value ;
      }
      else{  
         long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                     : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
               ) + tempo_deque->front->value ;                             //i+1, j+1
         long cas2 = INSERTION_COST + deque->front->value ;             //i+1, j
         if (cas2 < min) min = cas2 ;
         long cas3 = INSERTION_COST + tempo_deque->front->next->value ;        //i  , j+1
         if (cas3 < min) min = cas3 ;

         res=min;
         } 
      pushBack(deque, res);
      for (int64_t j = J-1; j > max(J-K,-1); j--) {
         char Yj = ctx.Y[j] ;
         long res ;
         if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
            {  ManageBaseError( Xi ) ;
               res = deque->front->next->value ;
            }
         else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
         {  ManageBaseError( Yj ) ;
            res =  deque->rear->value ;
         }
         else{  
            long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                        : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                  ) + deque->front->value ;                             //i+1, j+1
            long cas2 = INSERTION_COST + deque->front->next->value ; //i+1, j
            if (cas2 < min) min = cas2 ;
            long cas3 = INSERTION_COST + deque->rear->value ;        //i  , j+1
            if (cas3 < min) min = cas3 ;

            res=min;
            } 
         pushBack(deque, res);
         popFront(deque);
      }
      long elt= popFront(deque);
      pushBack(tempo_deque,elt);
      popFront(tempo_deque);
   }
   popFront(tempo_deque);
   pushBack(tempo_deque,deque->rear->value);
   pushBackBlock(listBlock,*deque);

   return tempo_deque;
}

queue*  block_right(int I,size_t N,size_t M,struct NW_MemoContext ctx,queueBlock* listBlock){
   queue* originalDeque=&(listBlock->front->value);
   queue *deque =copyQueue(originalDeque) ; // Crée une nouvelle deque vide
   
   queue *tempo_deque = createFile();//utilisé par le block d'à gauche de celui là
   // Initialisation de la derniere ligne de Ø avec les valeurs appropriées
   char Xi = ctx.X[I] ;
   pushBack(deque, (isBase(Xi) ? INSERTION_COST : 0) + deque->front->value);

   for (int64_t j = N-1; j > max(N-K,-1); j--) {
      char Yj = ctx.Y[j] ;
      long res ;
      if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
         {  ManageBaseError( Xi ) ;
            res = deque->front->next->value ;
         }
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res =  deque->rear->value ;
      }
      else{  
         long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                  : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
            ) + deque->front->value ;                             //i+1, j+1
         long cas2 = INSERTION_COST + deque->front->next->value ; //i+1, j
         if (cas2 < min) min = cas2 ;
         long cas3 = INSERTION_COST + deque->rear->value ;        //i  , j+1
         if (cas3 < min) min = cas3 ;

         res=min;
         } 
      pushBack(deque, res);
      popFront(deque);
   }
   popFront(deque);


   for (int64_t i = I-1; i >  max(I-K,-1); i--) {
      char Xi = ctx.X[i] ;
      pushBack(deque, (isBase(Xi) ? INSERTION_COST : 0) + deque->front->value);
      for (int64_t j = N-1; j > max(N-K,-1); j--) {
               char Yj = ctx.Y[j] ;
               long res ;
               if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
                  {  ManageBaseError( Xi ) ;
                     res = deque->front->next->value ;
                  }
               else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
               {  ManageBaseError( Yj ) ;
                  res =  deque->rear->value ;
               }
               else{  
                  long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                           : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
                     ) + deque->front->value ;                             //i+1, j+1
                  long cas2 = INSERTION_COST + deque->front->next->value ; //i+1, j
                  if (cas2 < min) min = cas2 ;
                  long cas3 = INSERTION_COST + deque->rear->value ;        //i  , j+1
                  if (cas3 < min) min = cas3 ;

                  res=min;
                  } 
               pushBack(deque, res);
               popFront(deque);
            }
      long elt= popFront(deque);
      pushBack(tempo_deque,elt);
   }
   pushBack(tempo_deque,deque->rear->value);
   pushBackBlock(listBlock,*deque);
   return tempo_deque;  
}

queue*  block(int I,int J,queue* tempo_deque,size_t N,size_t M,struct NW_MemoContext ctx,queueBlock* listBlock){
   long elt =(popFrontBlock(listBlock)).rear->value ; // Crée une nouvelle deque vide
   
   queue *originalDeque =&(listBlock->front->value);
   queue * deque=copyQueue(originalDeque);

   char Xi = ctx.X[I] ;
   char Yj = ctx.Y[J] ;
   long res ;
   if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
      {  ManageBaseError( Xi ) ;
         res = deque->front->value;
      }
   else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
   {  ManageBaseError( Yj ) ;
      res =  tempo_deque->front->value ;
   }
   else{  
      long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
               : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
         ) + elt ;                             //i+1, j+1
      long cas2 = INSERTION_COST + deque->front->value ; //i+1, j
      if (cas2 < min) min = cas2 ;
      long cas3 = INSERTION_COST + tempo_deque->front->value ;        //i  , j+1
      if (cas3 < min) min = cas3 ;

      res=min;
      } 
   pushBack(deque, res);

   for (int64_t j = J-1; j >max(J-K,-1); j--) {
      char Yj = ctx.Y[j] ;
      long res ;
      if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
         {  ManageBaseError( Xi ) ;
            res = deque->front->next->value ;
         }
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res =  deque->rear->value ;
      }
      else{  
         long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                  : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
            ) + deque->front->value ;                             //i+1, j+1
         long cas2 = INSERTION_COST + deque->front->next->value ; //i+1, j
         if (cas2 < min) min = cas2 ;
         long cas3 = INSERTION_COST + deque->rear->value ;        //i  , j+1
         if (cas3 < min) min = cas3 ;

         res=min;
         } 
      pushBack(deque, res);
      popFront(deque);
   }
   popFront(deque);

   for (int64_t i = I-1; i >  max(I-K,-1); i--) {

      char Yj = ctx.Y[J] ;
      long res ;
      char Xi = ctx.X[i] ;
      if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
         {  ManageBaseError( Xi ) ;
            res = deque->front->value ;
         }
      else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
      {  ManageBaseError( Yj ) ;
         res =  tempo_deque->front->next->value ;
      }
      else{  
         long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                     : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
               ) + tempo_deque->front->value ;                             //i+1, j+1
         long cas2 = INSERTION_COST + deque->front->value ;             //i+1, j
         if (cas2 < min) min = cas2 ;
         long cas3 = INSERTION_COST + tempo_deque->front->next->value ;        //i  , j+1
         if (cas3 < min) min = cas3 ;

         res=min;
         } 
      pushBack(deque, res);

      for (int64_t j = J-1; j > max(J-K,-1); j--) {
         char Yj = ctx.Y[j] ;
         long res ;
         if (! isBase(Xi))  /* skip ccharacter in Xi that is not a base */
            {  ManageBaseError( Xi ) ;
               res = deque->front->next->value ;
            }
         else if (! isBase(Yj))  /* skip ccharacter in Yj that is not a base */
         {  ManageBaseError( Yj ) ;
            res =  deque->rear->value ;
         }
         else{  
            long min =( isUnknownBase(Xi) ?  SUBSTITUTION_UNKNOWN_COST 
                     : ( isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST ) 
               ) + deque->front->value ;                             //i+1, j+1
            long cas2 = INSERTION_COST + deque->front->next->value ; //i+1, j
            if (cas2 < min) min = cas2 ;
            long cas3 = INSERTION_COST + deque->rear->value ;        //i  , j+1
            if (cas3 < min) min = cas3 ;

            res=min;
            } 
         pushBack(deque, res);
         popFront(deque);
      }
      long elt= popFront(deque);
      pushBack(tempo_deque,elt);
      popFront(tempo_deque);

      }
   pushBack(tempo_deque,deque->rear->value);
   popFront(tempo_deque);
   pushBackBlock(listBlock,*deque);

   return tempo_deque;
}

long iteratif_cache_aware(char* A, size_t lengthA, char* B, size_t lengthB) {
   _init_base_match() ;
   struct NW_MemoContext ctx;
   if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
   {  ctx.X = A ;
      ctx.M = lengthA ;
      ctx.Y = B ;
      ctx.N = lengthB ;
   }
   else
   {  ctx.X = B ;
      ctx.M = lengthB ;
      ctx.Y = A ;
      ctx.N = lengthA ;
   }
   size_t M = ctx.M ;
   size_t N = ctx.N ;
   queueBlock *listBlock =createFileBlock();
   queue * tempo_queue =block_down_right(N,M,ctx,listBlock);
   for (int J =N-K; J >= 0; J -= K){
      tempo_queue=block_down(J,tempo_queue,N,M,ctx,listBlock);
      }

   for (int I =M-K; I >= 0; I -= K){
         tempo_queue=block_right(I,N,M,ctx,listBlock);

         for (int J =N-K; J >= 0; J -= K){
            tempo_queue=block(I,J,tempo_queue,N,M,ctx,listBlock);
         }
         popFrontBlock(listBlock);
         
      }
   return (listBlock->rear->value.rear->value);

}

