// test setup hack
#define START_VECTOR_TEST(name) \
  START_TEST(name) \
  ealInt scalar; \
  ealInt A(5); \
  ealInt B(5); \
  ealInt C(5); \
  ealInt D(5); \
  ealInt twelve; \
  twelve = 12; \
  for( int i=0; i < 5; i++){ \
    A[i] = 1+i; \
    B[i] = 2+i; \
    C[i] = 2*(1+i); \
  } \


#include "enzo_unit_tests.h"
#include "ealInt.h"

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
   ASSERT_EQ(scalar.ReturnSize(), 1);
   ASSERT_EQ(scalar.Array[0], 0);
END_TEST()

START_VECTOR_TEST("create then assign scalar")
   ASSERT_EQ(twelve.ReturnSize(), 1);
   ASSERT_EQ(twelve.Array[0], 12);
END_TEST()

START_VECTOR_TEST("default array creation")
   ASSERT_EQ(D.ReturnSize(), 5);
   ASSERT_EQ(D.Array[4], 0);
END_TEST()

START_VECTOR_TEST("index operator")
   ASSERT_EQ(A[4], 5);
END_TEST()

START_VECTOR_TEST("vector = vector")
  A = B;
  ASSERT_EQ(A[4], 6);
END_TEST()

START_VECTOR_TEST("scalar = int")
  scalar = 2;
  ASSERT_EQ(scalar[0], 2);
END_TEST()

START_VECTOR_TEST("vector = vector + vector")
    A = B+C;
    ASSERT_EQ(A[3], 13);
END_TEST()

START_VECTOR_TEST("vector = vector - vector")
    A = B-C;
    ASSERT_EQ(A[3], -3);
END_TEST()

START_VECTOR_TEST("vector = vector * vector")
    A = B*C;
    ASSERT_EQ(A[3], 40);
END_TEST()

START_VECTOR_TEST("vector = vector / vector")
    B[3] = 17;
    A = B/C;
    ASSERT_EQ(A[3], 2);
    ASSERT_EQ(A[4], 0);
END_TEST()

START_VECTOR_TEST("vector = int")
  A = 2;
  ASSERT_EQ(A[2], 2);
END_TEST()

START_VECTOR_TEST("vector = vector * int")
    B = A*5;
    ASSERT_EQ(A[4], 5);
    ASSERT_EQ(B[4], 25);
END_TEST()

START_VECTOR_TEST("vector = vector * scalar")
    B = A*scalar;
    ASSERT_EQ(A[4], 5);
    ASSERT_EQ(B[4], 0);
END_TEST()

START_VECTOR_TEST("vector -= scalar")
    C-=twelve;
    ASSERT_EQ(C[4], -2);
END_TEST()

START_VECTOR_TEST("vector += scalar")
    C+=twelve;
    ASSERT_EQ(C[4], 22);
END_TEST()

START_VECTOR_TEST("vector *= scalar")
    C*=twelve;
    ASSERT_EQ(C[4], 120);
END_TEST()

START_VECTOR_TEST("vector /= scalar")
    C[4] = 25;
    C[3] = 12;
    C[2] = 11;
    C/=twelve;
    ASSERT_EQ(C[4], 2);
    ASSERT_EQ(C[3], 1);
    ASSERT_EQ(C[2], 0);
END_TEST()

START_VECTOR_TEST("vector -= int")
    C-=12;
    ASSERT_EQ(C[4], -2);
END_TEST()

START_VECTOR_TEST("vector += int")
    C+=12;
    ASSERT_EQ(C[4], 22);
END_TEST()

START_VECTOR_TEST("vector *= int")
    C*=12;
    ASSERT_EQ(C[4], 120);
END_TEST()

START_VECTOR_TEST("vector /= int")
    C[4] = 25;
    C[3] = 12;
    C[2] = 11;
    C/=12;
    ASSERT_EQ(C[4], 2);
    ASSERT_EQ(C[3], 1);
    ASSERT_EQ(C[2], 0);
END_TEST()


START_VECTOR_TEST("vector -= vector")
    C-=B;
    ASSERT_EQ(C[4], 4);
END_TEST()

START_VECTOR_TEST("vector += vector")
    C+=B;
    ASSERT_EQ(C[4], 16);
END_TEST()

START_VECTOR_TEST("vector *= vector")
    C*=B;
    ASSERT_EQ(C[4], 60);
END_TEST()

START_VECTOR_TEST("vector /= vector")
    C[4] = 18;
    C[2] = 3;
    C/=B;
    ASSERT_EQ(C[4], 3);
    ASSERT_EQ(C[3], 1);
    ASSERT_EQ(C[2], 0);
END_TEST()


START_VECTOR_TEST("vector[i] = vector[j]")
    B[3] = A[2];
    ASSERT_EQ(B[3], 3);
END_TEST()

START_VECTOR_TEST("vector.Min()")
    A[4] = -2;
int minval = A.Min();
    ASSERT_EQ(minval, -2);
END_TEST()

START_VECTOR_TEST("vector.Max()")
    A[1] = 99;
    int maxval = A.Max();
    ASSERT_EQ(maxval, 99);
END_TEST()

START_VECTOR_TEST("vector.Sum()")
    int sumval = A.Sum();
    ASSERT_EQ(sumval, 15);
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
   ASSERT_EQ(A[4], 5);
END_TEST()



START_VECTOR_TEST("vector.ReduceSum()")
    A.ReduceSum();
    int tot = NumberOfProcessors * 5;
    ASSERT_EQ(A[4], tot);
END_TEST()


START_VECTOR_TEST("vector.ReduceMin()")
    A -= MyProcessorNumber;
    A.ReduceMin();
ASSERT_EQ(A[3], 5 - NumberOfProcessors);
END_TEST()

START_VECTOR_TEST("vector.ReduceMax()")
    A += MyProcessorNumber;
    A.ReduceMax();
    ASSERT_EQ(A[4], 4+NumberOfProcessors);
END_TEST()

START_VECTOR_TEST("vector.GlobalSum()")
    int sumval = A.GlobalSum();
    ASSERT_EQ(sumval, NumberOfProcessors*15);    
END_TEST()

START_VECTOR_TEST("vector.GlobalMin()")
    A -= MyProcessorNumber;
    int minval = A.GlobalMin();
    ASSERT_EQ(minval, 2 - NumberOfProcessors);
END_TEST()

START_VECTOR_TEST("vector.GlobalMax()")
    A += MyProcessorNumber;
    int maxval = A.GlobalMax();
    A.GlobalMax();
ASSERT_EQ(maxval, 4+NumberOfProcessors);
END_TEST()

  CommunicationFinalize();
// End the overall test suite
END_TESTS()


