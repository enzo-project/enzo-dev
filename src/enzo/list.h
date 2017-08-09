/********************************************************
c GENERIC LINKED-LIST CLASS
c  
c Written by: Iryna Butsky 11/2016
c
c PURPOSE: 
c    Makes a linked list of arbitrary type. 
c    Current use case: SuperNovaList in Grid.h
c 
c FUNCTIONS: 
c    append() - appends to the head of the list
c    clear() - clears the list
c    size() - returns the size of the list (int)
c    
c EXAMPLES: 
c    List<SuperNova> SuperNovaList; // creates a new list of type SuperNova
c    SuperNovaList.append(Sn1); // appends a SuperNova object
c    SuperNovaList.append(Sn2); // appends a second SuperNova object
c    
c    // iterating through the list
c    List<SuperNova>::Iterator *It = SuperNovaList.begin()
c    while(It != SuperNovaList.end()){
c       Current_Supernova = It->get(); 
c       //....do something ....
c       It = It->next();  // It 
c 
c   SuperNovaList.clear(); //clears list and deletes from memory
***************************************************/


#ifndef ____List__
#define ____List__

template <class T>
class List
{
 public:

  class Iterator
  {
  public:
    Iterator *next()
    {
      return nextObject;
    }
    T *get()
    {
      return wrappedObject;
    }
  private:
    T *wrappedObject;
    Iterator *nextObject;
    friend class List;
  } ;


  List()
    {
      theList = new Iterator;
      theList->wrappedObject = NULL;
      theList->nextObject = NULL;
      theEnd = theList;
      listSize = 0;
    }

  ~List()
    {
      clear();
      delete theList;
    }

  void append(T *newThing)
  {
    Iterator *head = new Iterator;
    head->wrappedObject = newThing;
    head->nextObject = theList;

    theList = head;

    listSize += 1;
  }

  Iterator *begin()
  {
    return theList;
  }

  T *front() {
    return theList->wrappedObject;
  }

  void clear()
  {
    Iterator *it = theList;

    while (it != NULL) {
      delete it->wrappedObject;

      Iterator *next = it->nextObject;      
      delete it;
      it = next;
    }

    theList = new Iterator;
    theList->wrappedObject = NULL;
    theList->nextObject = NULL;
    theEnd = theList;
    listSize = 0;
  }

  void remove(T *thingToRemove);

  int size()
  {
    return listSize;
  }

  Iterator *end()
  {
    return theEnd;
  }
 private:
  Iterator *theList;
  Iterator *theEnd;
  int listSize;
  
} ;

#endif
