#include <stdio.h>
#include <stdlib.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#define SWAP(a,b)  temp =(a);(a)=(b);(b)=temp;
#define SWAPI(a,b) tempi=(a);(a)=(b);(b)=tempi;

float selectb(unsigned long k, unsigned long n, float arr[], int ind[])
{
	unsigned long i,ir,j,l,mid;
	float a,temp;
	int    ai,tempi;


	l=1;
	ir=n;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir])
				SWAPI(ind[l],ind[ir])  
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;

			SWAP(arr[mid],arr[l+1])
			SWAPI(ind[mid],ind[l+1])

			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
				SWAPI(ind[l],ind[ir])
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
				SWAPI(ind[l+1],ind[ir])
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1])
				SWAPI(ind[l],ind[l+1])
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			ai=ind[l+1];

			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
				SWAPI(ind[i],ind[j])
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			ind[l+1]=ind[j];
			ind[j]=ai;

			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}
#undef SWAP
#undef SWAPI
