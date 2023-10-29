#include "DequeNode.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <stdint.h>

// Fonction pour créer une nouvelle deque vide
Deque* createDeque() {
    Deque *deque = (Deque *)malloc(sizeof(Deque)); // Allocation de mémoire pour la deque
    if (deque == NULL) {
        fprintf(stderr, "Erreur d'allocation de mémoire\n"); // En cas d'erreur, affiche un message d'erreur et quitte le programme
        exit(1);
    }
    deque->front = NULL; // Initialise l'extrémité avant de la deque à NULL (vide)
    deque->rear = NULL;  // Initialise l'extrémité arrière de la deque à NULL (vide)
    return deque;
}

// Fonction pour ajouter un élément à l'extrémité arrière de la deque
void pushBack(Deque *deque, char value) {
    DequeNode *newNode = (DequeNode *)malloc(sizeof(DequeNode)); // Alloue de la mémoire pour un nouveau nœud
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
char popFront(Deque *deque) {
    if (deque->front == NULL) {
        fprintf(stderr, "Le deque est vide\n");
        exit(1);
    }
    DequeNode *node = deque->front;
    char value = node->value; // Récupère la valeur du nœud à l'extrémité avant
    deque->front = node->next; // Met à jour l'extrémité avant de la deque

    if (deque->front != NULL) {
        deque->front->prev = NULL; // Met à jour le pointeur vers le nœud précédent du nouveau nœud à l'extrémité avant
    } else {
        deque->rear = NULL; // Si la deque est vide après le retrait, l'extrémité arrière doit également être NULL
    }

    free(node); // Libère la mémoire du nœud retiré
    return value;
}