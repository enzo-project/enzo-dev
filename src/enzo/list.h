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
