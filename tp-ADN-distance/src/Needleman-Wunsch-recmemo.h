/**
 * \file Needleman-Wunsch-recmemo.h
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 */

#include <stdlib.h> /* for size_t */

/*
 * Costs for operations on canonical bases
 * Three  operations: insertion and sustitution of one base by an another 
 * Note= substitution of an unknown base N by another one (known or unknown) as the same cost than substitution between 2 different known bases
 */
/** \def SUBSTITUTION_COST
 *  \brief Cost of substitution of one canonical base by another
 */
#define SUBSTITUTION_COST	1

/** \def SUBSTITUTION_UNKNOWN_COST
 *  \brief Cost of substitution of an unknown base (N) by another one (canonical or unknown)
 */
#define SUBSTITUTION_UNKNOWN_COST	1  /* Cost for sustitition of an Unknown bas N by another on -known or unkown- */ 

/** \def INSERTION_COST
 *  \brief Cost of insertion of a canonical base 
 */
#define INSERTION_COST		2


#define D1_CACHE_SIZE (int64_t) 262144
#define CACHE_ALIGNEMENT (int64_t) 64
#define K  (int64_t) 64 //sqrt(262144/64) on suppose que les données sont non contigues

#define S  (int64_t) 64

// Définition of the structure for a deque node
typedef struct queueNode {
    long value;
    struct queueNode *prev;
    struct queueNode *next;
} queueNode;

// Définition of the structure for the deque
typedef struct queue {
    queueNode *front;
    queueNode *rear;
} queue;

/********************************************************************************
 * Recursive implementation of NeedlemanWunsch with memoization
 */
/**
 * \fn long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB);
 * \brief computes the edit distance between A[0 .. lengthA-1] and B[0 .. lengthB-1]
 * \param A  : array of char represneting a genetic sequence A 
 * \param lengthA :  number of elements in A 
 * \param B  : array of char represneting a genetic sequence B
 * \param lengthB :  number of elements in B 
 * \return :  edit distance between A and B }
 *
 * editDistance_RecMemo is a memoized recursive immplementatioin of Needleman-Wunsch algorithm.
 * It allocates the data structure for memoization table and calls the internal recursive function _editDistance_memo
 * that fills in the memoization table.
 * 
 * If lengthA < lengthB, the sequences A and B are swapped.
 *
 */
long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB);

long iteratif(char* A, size_t lengthA, char* B, size_t lengthB);

// Function to create a new empty deque
queue* createFile();

// Function to add an element to the rear end of the deque
void pushBack(queue *deque, long value);

// Function to remove an element from the front end of the deque
long popFront(queue *deque);
typedef struct Block {
    queue value;
    struct Block *prev;
    struct Block *next;
} Block;

// Définition of the structure for the deque
typedef struct queueBlock {
    Block *front;
    Block *rear;
} queueBlock;

// Structure pour représenter une paire ((i, j), long valeur)
typedef struct {
    int i;
    int j;
    long valeur;
} Pair;

// Structure pour représenter un élément de la table de hachage
typedef struct HashNode {
    Pair pair;
    struct HashNode* next;
} HashNode;

// Structure pour représenter la table de hachage
typedef struct HashTable {
    HashNode** table;
    size_t size;
} HashTable;


queueBlock* createFileBlock();
void pushBackBlock(queueBlock *deque, queue value) ;
queue popFrontBlock(queueBlock *deque) ;

long iteratif_cache_aware(char* A, size_t lengthA, char* B, size_t lengthB) ;
long rec_cache_oblivious(char* A, size_t lengthA, char* B, size_t lengthB) ;
