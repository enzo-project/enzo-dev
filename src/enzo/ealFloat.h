#ifndef EAL_FLOAT_H
#define EAL_FLOAT_H

class ealFloat{

 public:
  ealFloat(int size=1,
			 float start=0.,
			 float end=1.,
			 const char *axis_type="none");

  ~ealFloat();

  ealFloat( const ealFloat & );

  void Zero();
  void Log(float start, float end);
  void Linear(float start, float end);
  float Sum();
  float Min();
  float Max();
  void ReduceSum();
  void ReduceMin();
  void ReduceMax();
  float GlobalSum();
  float GlobalMin();
  float GlobalMax();

  void Bcast(int FromProcessor);
  void Print(const char *fmt="%6.3"FSYM, FILE *stream=NULL);
  void PrintWithIndex(const char *fmt="%"ISYM"  %6.3"FSYM, FILE *stream=NULL);
  float *Array;
 
  float &operator[](int);
  float operator[](int) const;

  /******** Assignment ********/
  const ealFloat &operator=(const ealFloat &);
  const ealFloat &operator=(const float &);

  const ealFloat &operator+=(const ealFloat &);
  const ealFloat &operator+=(const float &);

  const ealFloat &operator-=(const ealFloat &);
  const ealFloat &operator-=(const float &);

  const ealFloat &operator*=(const ealFloat &);
  const ealFloat &operator*=(const float &);

  const ealFloat &operator/=(const ealFloat &);
  const ealFloat &operator/=(const float &);


  /******** Mathematical ********/
  const ealFloat operator+(const ealFloat &);
  const ealFloat operator+(const float &);

  const ealFloat operator-(const ealFloat &);
  const ealFloat operator-(const float &);

  const ealFloat operator*(const ealFloat &);
  const ealFloat operator*(const float &);

  const ealFloat operator/(const ealFloat &);
  const ealFloat operator/(const float &);

  inline int ReturnSize(){ return Size; };

 private:
  int Size;
};

#endif /* ENZO_ANALYSIS_FLOAT_ARRAY_H */
