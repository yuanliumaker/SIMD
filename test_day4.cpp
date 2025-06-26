#include <iostream>
#include <x86intrin.h>
using namespace std;
/*AVX512f(foundation)
avx-512 提供了32个512位的ZMM寄存器，可以一次性出来更多数据，比之前的avx 和avx2 指令集所支持的
256位数据宽度更大*/ 
/*1 __mm512_2intersect_epi32(__m512i a,__m512i b,_mmask16* k1,_mmask16* k2)
作用：计算两个向量a 和b 的交集，k1 存储交集中属于a 的元素，k2存储属于b的元素
类似的函数有 __mm512_2intersect_epi64*/
void test()
{
    int a[16],b[16];
    for(int i=0;i<16;i++)
    {
        a[i]=i;
        b[i]=i+10;
    }
    __m512i a1=_mm512_loadu_si512(a);
    __m512i b1=_mm512_loadu_si512(b);
    unsigned short k1,k2;
    _mm512_2intersect_epi32(a1,b1,&k1,&k2);
    
    
} 
int main()
{   
    test();
    return 0;
}