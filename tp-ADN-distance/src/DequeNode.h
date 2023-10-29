#ifndef DEQUENODE_H
#define DEQUENODE_H

// Définition of the structure for a deque node
typedef struct DequeNode {
    char value;
    struct DequeNode *prev;
    struct DequeNode *next;
} DequeNode;

// Définition of the structure for the deque
typedef struct Deque {
    DequeNode *front;
    DequeNode *rear;
} Deque;

// Function to create a new empty deque
Deque* createDeque();

// Function to add an element to the rear end of the deque
void pushBack(Deque *deque, char value);

// Function to remove an element from the front end of the deque
char popFront(Deque *deque);

#endif // DEQUENODE_H
