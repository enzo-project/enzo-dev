// test setup hack
#define START_VECTOR_TEST(name) \
  START_TEST(name) \
  ealFloat scalar; \
  ealFloat A(5); \
  ealFloat B(5); \
  ealFloat C(5); \
  ealFloat D(5); \
  ealFloat twelve; \
  twelve = 12; \
  for( int i=0; i < 5; i++){ \
    A[i] = 1+i; \
    B[i] = 2+i; \
    C[i] = 2*(1+i); \
  } \


#include "enzo_unit_tests.h"
#include "ealFloat.h"

// Start the overall test suite
START_TESTS()

  CommunicationInitialize(&argc, &argv);
  // Initialize Communications   

START_VECTOR_TEST("print test")

  /* Initialize Communications package. */

  if( MyProcessorNumber == ROOT_PROCESSOR ){
    printf("Scalar   = ");  
    scalar.Print();
    printf("\n");

    printf("Twelve   = ");  
    twelve.Print();
    printf("\n");

    printf("Vector A = ");  
    A.Print();
    printf("\n");

    printf("Vector B = ");  
    B.Print();
    printf("\n");

    printf("Vector C = ");  
    C.Print();
    printf("\n");

    printf("Vector D = ");  
    D.Print();
    printf("\n");

  }

ASSERT(TRUE);
END_TEST()

START_VECTOR_TEST("default scalar creation")
   ASSERT_EQUALS_INT(scalar.ReturnSize(), 1);
   ASSERT_EQ(scalar.Array[0], 0);
END_TEST()

START_VECTOR_TEST("create then assign scalar")
   ASSERT_EQUALS_INT(twelve.ReturnSize(), 1);
   ASSERT_EQUALS_FLOAT(twelve.Array[0], 12);
END_TEST()

START_VECTOR_TEST("default array creation")
   ASSERT_EQUALS_INT(D.ReturnSize(), 5);
   ASSERT_EQUALS_FLOAT(D.Array[4], 0);
END_TEST()

START_VECTOR_TEST("index operator")
   ASSERT_EQUALS_FLOAT(A[4], 5);
END_TEST()

START_VECTOR_TEST("vector = vector")
  A = B;
  ASSERT_EQUALS_FLOAT(A[4], 6);
END_TEST()

START_VECTOR_TEST("scalar = float")
  scalar = 2;
  ASSERT_EQUALS_FLOAT(scalar[0], 2);
END_TEST()

START_VECTOR_TEST("vector = vector + vector")
    A = B+C;
    ASSERT_EQUALS_FLOAT(A[3], 13);
END_TEST()

START_VECTOR_TEST("vector = vector - vector")
    A = B-C;
    ASSERT_EQUALS_FLOAT(A[3], -3);
END_TEST()

START_VECTOR_TEST("vector = vector * vector")
    A = B*C;
    ASSERT_EQUALS_FLOAT(A[3], 40);
END_TEST()

START_VECTOR_TEST("vector = vector / vector")
    B[3] = 17;
    float checkval = 17.0/8.0;
    A = B/C;
    ASSERT_EQUALS_FLOAT(A[3], checkval);
    ASSERT_EQUALS_FLOAT(A[4], 0.6);
END_TEST()

START_VECTOR_TEST("vector = float")
  A = 2;
  ASSERT_EQUALS_FLOAT(A[2], 2);
END_TEST()

START_VECTOR_TEST("vector = vector * float")
    B = A*5;
    ASSERT_EQUALS_FLOAT(A[4], 5);
    ASSERT_EQUALS_FLOAT(B[4], 25);
END_TEST()

START_VECTOR_TEST("vector = vector * scalar")
    B = A*scalar;
    ASSERT_EQUALS_FLOAT(A[4], 5);
    ASSERT_EQUALS_FLOAT(B[4], 0);
END_TEST()

START_VECTOR_TEST("vector -= scalar")
    C-=twelve;
    ASSERT_EQUALS_FLOAT(C[4], -2);
END_TEST()

START_VECTOR_TEST("vector += scalar")
    C+=twelve;
    ASSERT_EQUALS_FLOAT(C[4], 22);
END_TEST()

START_VECTOR_TEST("vector *= scalar")
    C*=twelve;
    ASSERT_EQUALS_FLOAT(C[4], 120);
END_TEST()

START_VECTOR_TEST("vector /= scalar")
    C[4] = 25;
    C[3] = 12;
    C[2] = 11;
    C/=twelve;
    ASSERT_EQUALS_FLOAT(C[4], 25.0/12.0);
    ASSERT_EQUALS_FLOAT(C[3], 1.0);
    ASSERT_EQUALS_FLOAT(C[2], 11.0/12.0);
END_TEST()

START_VECTOR_TEST("vector -= float")
    C-=12;
    ASSERT_EQUALS_FLOAT(C[4], -2);
END_TEST()

START_VECTOR_TEST("vector += float")
    C+=12;
    ASSERT_EQUALS_FLOAT(C[4], 22);
END_TEST()

START_VECTOR_TEST("vector *= float")
    C*=12;
    ASSERT_EQUALS_FLOAT(C[4], 120);
END_TEST()

START_VECTOR_TEST("vector /= float")
    C[4] = 25;
    C[3] = 12;
    C[2] = 11;
    C/=12;
    ASSERT_EQUALS_FLOAT(C[4], 25.0/12.0);
    ASSERT_EQUALS_FLOAT(C[3], 1.0);
    ASSERT_EQUALS_FLOAT(C[2], 11.0/12.0);
END_TEST()


START_VECTOR_TEST("vector -= vector")
    C-=B;
    ASSERT_EQUALS_FLOAT(C[4], 4);
END_TEST()

START_VECTOR_TEST("vector += vector")
    C+=B;
    ASSERT_EQUALS_FLOAT(C[4], 16);
END_TEST()

START_VECTOR_TEST("vector *= vector")
    C*=B;
    ASSERT_EQUALS_FLOAT(C[4], 60);
END_TEST()

START_VECTOR_TEST("vector /= vector")
    C[4] = 18;
    C[2] = 3;
    C/=B;
    ASSERT_EQUALS_FLOAT(C[4], 3.0);
    ASSERT_EQUALS_FLOAT(C[3], 8.0/5.0);
    ASSERT_EQUALS_FLOAT(C[2], 0.75);
END_TEST()


START_VECTOR_TEST("vector[i] = vector[j]")
    B[3] = A[2];
    ASSERT_EQUALS_FLOAT(B[3], 3);
END_TEST()

START_VECTOR_TEST("vector.Min()")
    A[4] = -2;
   float minval = A.Min();
    ASSERT_EQUALS_FLOAT(minval, -2);
END_TEST()

START_VECTOR_TEST("vector.Max()")
    A[1] = 99;
    float maxval = A.Max();
    ASSERT_EQUALS_FLOAT(maxval, 99);
END_TEST()

START_VECTOR_TEST("vector.Sum()")
    float sumval = A.Sum();
    ASSERT_EQUALS_FLOAT(sumval, 15.0);
END_TEST()

START_VECTOR_TEST("vector.Zero()")
    A.Zero();
    ASSERT_EQ(A[2], 0);
END_TEST()


START_VECTOR_TEST("vector.Bcast()")
    if(MyProcessorNumber != ROOT_PROCESSOR){
       A.Zero();
   }
   A.Bcast(ROOT_PROCESSOR);
//   A.Print();
//   printf("\n");
   ASSERT_EQUALS_FLOAT(A[4], 5);
END_TEST()



START_VECTOR_TEST("vector.ReduceSum()")
    A.ReduceSum();
    float tot = NumberOfProcessors * 5;
    ASSERT_EQUALS_FLOAT(A[4], tot);
END_TEST()


START_VECTOR_TEST("vector.ReduceMin()")
    A -= MyProcessorNumber;
    A.ReduceMin();
ASSERT_EQUALS_FLOAT(A[3], 5 - NumberOfProcessors);
END_TEST()

START_VECTOR_TEST("vector.ReduceMax()")
    A += MyProcessorNumber;
    A.ReduceMax();
    ASSERT_EQUALS_FLOAT(A[4], 4+NumberOfProcessors);
END_TEST()

START_VECTOR_TEST("vector.GlobalSum()")
    float sumval = A.GlobalSum();
    ASSERT_EQUALS_FLOAT(sumval, NumberOfProcessors*15);    
END_TEST()

START_VECTOR_TEST("vector.GlobalMin()")
    A -= MyProcessorNumber;
    float minval = A.GlobalMin();
    ASSERT_EQUALS_FLOAT(minval, 2 - NumberOfProcessors);
END_TEST()

START_VECTOR_TEST("vector.GlobalMax()")
    A += MyProcessorNumber;
    float maxval = A.GlobalMax();
    A.GlobalMax();
ASSERT_EQUALS_FLOAT(maxval, 4+NumberOfProcessors);
END_TEST()

  CommunicationFinalize();
// End the overall test suite
END_TESTS()
