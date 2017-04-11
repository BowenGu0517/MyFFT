/****************************************************

Project Name: MyFFT function
Author: BOWEN GU
Function: Using FFT algorithm to do DFT for discrete signals

****************************************************/

#include <iostream>
using namespace std;

#define BIT_NUM (32)
#ifndef PI
#define PI (3.1415926536)
#endif 

#define TEST_ARR_NUM (17)

struct Complex
{
	float real;
	float imag;
};

struct DataInfo
{
	int dataLength;
	Complex *pArr;
};

Complex complexAdd(Complex A, Complex B);
Complex complexSub(Complex A, Complex B);
Complex complexMultiply(Complex A, Complex B);
Complex getW(int N, int k);
DataInfo* doFFT(DataInfo *pDataIn);
DataInfo* reverse(DataInfo *pDataIn);
void display(DataInfo *pDataIn);
bool isPowerOfTwo(int dataIn);
DataInfo* doZeroTapping(DataInfo* pDataIn);

int main()
{
	DataInfo *pAll = new DataInfo;
	pAll->pArr = new Complex[TEST_ARR_NUM];
	for (int i = 0; i < TEST_ARR_NUM; i++)
	{
		(pAll->pArr + i)->real = (float)i;
		(pAll->pArr + i)->imag = (float)i;
	}
	pAll->dataLength = TEST_ARR_NUM;

	pAll = doFFT(pAll);
	//pAll = doZeroTapping(pAll);
	display(pAll);

	system("pause");
	return 0;
}

Complex complexAdd(Complex A, Complex B)
{
	Complex result;
	result.real = A.real + B.real;
	result.imag = A.imag + B.imag;
	return result;
}

Complex complexSub(Complex A, Complex B)
{
	Complex result;
	result.real = A.real - B.real;
	result.imag = A.imag - B.imag;
	return result;
}

Complex complexMultiply(Complex A, Complex B)
{
	Complex result;
	result.real = A.real*B.real - A.imag*B.imag;
	result.imag = A.real*B.imag + A.imag*B.real;
	return result;
}

Complex getW(int N, int k)
{
	Complex W;
	W.real = cosf((float)(2 * PI*k) / N);
	W.imag = sinf((float)(-2 * PI*k) / N);
	return W;
}

DataInfo* doFFT(DataInfo *pDataIn)
{
	pDataIn = doZeroTapping(pDataIn);
	int totalLength = pDataIn->dataLength;
	int curLength = 1;

	Complex *WstoredArr = new Complex[totalLength / 2];
	for (int i = 0; i < totalLength / 2; i++)
	{
		*(WstoredArr + i) = getW(totalLength, i);
	}

	DataInfo *pFinalResult = reverse(pDataIn);
	Complex* pProcessedArr = pFinalResult->pArr;
	Complex *pFinalAns = new Complex[totalLength];

	while (curLength != totalLength)
	{
		curLength *= 2;
		int curArrNum = totalLength/curLength;
		
		for (int i = 0; i < curArrNum; i++)
		{

			for (int j = 0; j < (curLength / 2); j++)
			{	
				Complex X1 = *(pProcessedArr + i*curLength + j);
				Complex X2 = *(pProcessedArr + i*curLength + j + curLength / 2);
				X2 = complexMultiply(*(WstoredArr+j*totalLength/curLength),X2);
				*(pFinalAns + i*curLength + j) = complexAdd(X1,X2);
				*(pFinalAns + i*curLength + j +curLength/2) = complexSub(X1,X2);
			}
		}

		for (int i = 0; i < totalLength; i++)
		{
			*(pProcessedArr + i) = *(pFinalAns + i);
		}
	}

	for (int i = 0; i < totalLength; i++)
	{
		*(pFinalResult->pArr + i) = *(pFinalAns + i);
	}

	return pFinalResult;
}

DataInfo* reverse(DataInfo *pDataIn)
{
	DataInfo *pResult = pDataIn;
	int arrLength = pDataIn->dataLength;
	Complex *pArr = new Complex[arrLength];
	int N = log2(arrLength);
	for (int i = 0; i < arrLength; ++i)
	{
		int *ptempArr=new int[N];
		for (int j = 0; j < N; j++)
		{
			*(ptempArr + N - j - 1) = (i >> j) & 1;
		}

		int tempNum = 0;
		for (int j = 0; j < N; j++)
		{
			tempNum += *(ptempArr + j)*pow(2, j);
		}
		*(pArr + tempNum) = *(pDataIn->pArr + i);
	}

	for (int i = 0; i < arrLength; i++)
	{
		*(pResult->pArr + i) = *(pArr + i);
	}
	return pResult;
}

void display(DataInfo *pDataIn)
{
	for (int i = 0; i < pDataIn->dataLength; ++i)
	{
		cout << (pDataIn->pArr+i)->real<<" + j"<< (pDataIn->pArr + i)->imag << endl;
	}
}

bool isPowerOfTwo(int dataIn)
{
	if (dataIn == 0)
	{
		printf("There is no data input! Program exit!\n");
		getchar();
		exit(0);
	}

	if (dataIn == 1)
	{
		return true;
	}

	else
	{
		if (dataIn % 2 != 0)
		{
			return false;
		}
		else
		{
			int j = dataIn;
			int k = 0;
			while (j != 1 && k == 0)
			{
				j = j / 2;
				k = j % 2;
			}
			if (j == 1)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
}

DataInfo* doZeroTapping(DataInfo* pDataIn)
{
	DataInfo *pResult = pDataIn;
	if (!isPowerOfTwo(pDataIn->dataLength))
	{
		int curLength = pDataIn->dataLength;
		int intArr[BIT_NUM] = { 0 };
		for (int i = 0; i < BIT_NUM; i++)
		{
			intArr[i] = (curLength >> i & 1);
		}

		int maxBit = 0;
		for (int i = BIT_NUM; i >= 0; i--)
		{
			if (intArr[i] == 1)
			{
				maxBit = i + 1;
				break;
			}
		}
		int totalLength = pow(2, maxBit);

		DataInfo *pNew=new DataInfo;
		pNew->dataLength = totalLength;
		pNew->pArr = new Complex[totalLength];
		for (int i =0;i<curLength;i++)
		{
			*(pNew->pArr + i) = *(pDataIn->pArr + i);
		}
		for (int i = curLength; i < totalLength; i++)
		{
			(pNew->pArr + i)->real = 0.0f;
			(pNew->pArr + i)->imag = 0.0f;
		}
		return pNew;
	}
	return pResult;
}
