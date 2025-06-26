#include <x86intrin.h>
#include <iostream>
using namespace std;
/*sse4.1*/ 
/*__m128i _mm_blend_epi16(__m128i a,__m128i b,const int imm8)
作用：根据控制掩码M从a 或b 中合并整数值
如果掩码位是0，则从a中复制相应的元素到结果
掩码位是1，则从b中复制相应的元素到结果
类似的函数有double 类型：__mm128d _mm_blend_pd(__m128d a,__m128d b,const int imm8)
float 类型：__mm128 _mm_blend__ps(__m128 a,__m128 b,const int imm8)
*/ 
void test(){
    __m128i a=_mm_set_epi16(1,2,3,4,5,6,7,8);
    __m128i b=_mm_set_epi16(9,10,11,12,13,14,15,16);
    __m128d c=_mm_set_pd(1.0,2.0);
    __m128d d=_mm_set_pd(3.0,4.0);
    
    __m128i result=_mm_blend_epi16(a,b,0b01010101);
    __m128d result1=_mm_blend_pd(c,d,0b01);
    int16_t* res=(int16_t*)&result;
    double* res1=(double*)&result1;
    for (int i=0;i<8;i++)
    {
        cout<<res[i]<<" "; //16 7 14 5 12 3 10 1
        
    }
    cout<<res1[0]<<res1[1]<<endl; //4 1 
}
/*__m128i _mm_blendv_epi8 (__m128i a,__m128ib,__m128i mask)
和上述作用一样，区别在于第三个参数掩码的不同
类似的有 _mm_blendv_pd
_mm_blendv_ps*/

void test1()
{
    __m128i a=_mm_set_epi8(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
    __m128i b=_mm_set_epi8(17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32);
    __m128i m=_mm_set_epi8(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1);
    __m128i result=_mm_blendv_epi8(a,b,m);
    int8_t* res=(int8_t*)&result;
    for (int i=0;i<16;i++)
    {
        cout<<res[i]<<" ";
    }    
}

/*_mm_ceil_pd(__m128d a)
作用:对寄存器中的每个数向上取整
类似的float：_mm_ceil_ps()
_mm_ceil_sd(__m128d a,__m128d b)
作用；scalar double :将b的低64位（即第一个双精度浮点值）向上取整为整数值，
并将结果复制到返回值的低64位，高64位保持a的原值
类似的 __mm_ceil_ss()*/ 
void test2()
{
    __m128d a=_mm_set_pd(1.3,2.5); //a=[2.5,1.3]，a[0]=2.5 低64位 a[1]=1.3 高64位
    __m128d result=_mm_ceil_pd(a);
    double* res=(double*)&result;
    // cout<<res[0]<<" "<<res[1]<<endl; //3 2
    __m128d b=_mm_set_pd(3.2,4.7);
    __m128d re=_mm_ceil_sd(a,b);
    double* res1=(double*)&re;
    cout<<re[0]<<re[1]<<endl;

}
/*_mm_cmpeq_epi64 比较两个寄存器中的对应元素是否相等，相等则1，否则0 
easy 略*/ 


/*_mm_cvtepi16_epi32(__128i a)
作用：将128位向量a 中的每个16位整数元素拓展为32 位整数（显然从低位开始取），并返回
类似的函数还有 _mm_cvtepi16_epi64
_mm_cvtepi32_epi64
_mm_cvtepi8_epi16等
_mm_cvtepu32_epi64 32位无符号转成64位整数*/
void test3()
{
    __m128i a=_mm_set_epi16(1,2,3,4,5,6,7,8);//a=[8,7,6,5,4,3,2,1]
    __m128i result=_mm_cvtepi16_epi32(a);
    int32_t* re=(int32_t*)&result;
    cout<<re[0]<<re[1]<<re[2]<<re[3]<<endl;//8 7 6 5
} 
/*提取
int _mm_extract_epi32(__m128i a,const int N)
N:一个立即数 根据N从向量a 中提取指定的32位整数元素
类似的还有：_mm_extract_epi64 _mm_extract_epi8
float _mm_extract_ps
*/
void test4() 
{
    __m128i a=_mm_set_epi32(1,2,3,4);//a=[4,3,2,1]
    int result=_mm_extract_epi32(a,1);//提取a[1]
    cout<<result<<endl;

}
/*__mm_insert_epi32(__mm128i a,int x,const int N)
作用：该函数首先通过复制a来构造一个向量，然后由立即数N指定的偏移处插入整数x
类似的 __mm_insert_epi64
__mm_insert_epi8
__mm_insert_ps*/
void test5()
{
    __m128i a=_mm_set_epi32(1,2,3,4);//a=[4,3,2,1]
    __m128i result=_mm_insert_epi32(a,5,2);//a[2]处插5
    int* res=(int*)&result;
    cout<<res[0]<<res[1]<<res[2]<<res[3]<<endl;
} 
/*_mm_minpos_epu16 
作用：128位 8个无符号整数 查找最小值及其索引 
*/
void test6()
{
    __m128i a=_mm_set_epi16(10,20,4,13,27,32,40,17);
    __m128i result=_mm_minpos_epu16(a);
    // 寄存器的前两个分别是最小值和索引，其余部分为0
    
    int16_t* res=(int16_t*)&result;
    cout<<"最小值"<<res[0]<<endl;
    cout<<"最小值索引"<<res[1]<<endl;
} 
/*_mm_mpsadbw_epu8(__m128i a,__m128i b,const int imm8)
作用: 用于计算a 和b 向量的绝对差值之和,参数imm8 用于指定如何计算这些差值的具体规则
通常用于定义哪一对元素将参与计算*/
void test7()
{
    __m128i a=_mm_set_epi8(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
    __m128i b=_mm_set_epi8(-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-16);
    __m128i result=_mm_mpsadbw_epu8(a,b,0x0);
    int16_t* res=(int16_t*)&result;
    for (int i=0;i<8;i++)
    {
        cout<<res[i]<<" ";
    }
}
/*_mm_mullo_epi32(__m128i a,__m128i b)
两个向量对应元素相乘，返回低32位(low)的结果*/
void test8()
{
    __m128i a=_mm_set_epi32(1,2,3,4);
    __m128i b=_mm_set_epi32(5,6,7,8);
    __m128i result=_mm_mullo_epi32(a,b);
    int res[4];
    _mm_store_si128((__m128i*)res,result);
   for (int i=0;i<4;i++){
    cout<<res[i]<<" ";//32 21, 12 5
   }
}  
/*_mm_testc_si128(__m128i a,__m128ib)
a和b 进行按位与操作返回结果
*/ 
int main()
{
    // test();
    // test1();
    // test2();
    // test3();
    // test4();
    // test5();
    // test6();
    // test7();
    test8();
    return 0;
}