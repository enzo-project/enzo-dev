#ifndef EAL_INT_H
#define EAL_INT_H

class ealInt{

 public:
  ealInt(int size=1);

  ~ealInt();

  ealInt( const ealInt & );

  void Zero();
  int Sum();
  int Min();
  int Max();
  void ReduceSum();
  void ReduceMin();
  void ReduceMax();
  int GlobalSum();
  int GlobalMin();
  int GlobalMax();

  void Bcast(int FromProcessor);
  void Print(const char *fmt="%6"ISYM, FILE *stream=NULL);
  void PrintWithIndex(const char *fmt="%"ISYM"  %6"ISYM, FILE *stream=NULL);
  int *Array;
 
  int &operator[](int);
  int operator[](int) const;

  /******** Assignment ********/
  const ealInt &operator=(const ealInt &);
  const ealInt &operator=(const int &);

  const ealInt &operator+=(const ealInt &);
  const ealInt &operator+=(const int &);

  const ealInt &operator-=(const ealInt &);
  const ealInt &operator-=(const int &);

  const ealInt &operator*=(const ealInt &);
  const ealInt &operator*=(const int &);

  const ealInt &operator/=(const ealInt &);
  const ealInt &operator/=(const int &);


  /******** Mathematical ********/
  const ealInt operator+(const ealInt &);
  const ealInt operator+(const int &);

  const ealInt operator-(const ealInt &);
  const ealInt operator-(const int &);

  const ealInt operator*(const ealInt &);
  const ealInt operator*(const int &);

  const ealInt operator/(const ealInt &);
  const ealInt operator/(const int &);

  inline int ReturnSize(){ return Size; };

 private:
  int Size;
};

#endif 
