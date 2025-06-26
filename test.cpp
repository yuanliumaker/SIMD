#include <iostream>
#include <xmmintrin.h>
#include <smmintrin.h> /* SSE4.1 which includes __mm_fmadd_ps*/
using namespace std;
typedef struct vec4{
    float x;
    float y;
    float z;
    float w;
}vec4;
float mul( vec4* v1, vec4* v2)
{
    float v3=v1->x*v2->x+v1->y*v2->y+v1->z*v2->z+v1->w*v2->w;
    return v3;
}
void simd_dot(vec4* v1,vec4* v2)
{
    // 寄存器m1,m2 
    // 注意 _mm_load_ps 需要的是指向float 的指针
    __m128 m1=_mm_load_ps(reinterpret_cast<float*>(&v1->x));
    __m128 m2=_mm_load_ps(reinterpret_cast<float*>(&v2->x));
    __m128 m3=_mm_mul_ps(m1, m2);
    float* res=(float*)&m3;
    cout<<res[0]<<res[1]<<res[2]<<res[3]<<endl;
    
}
/*
_mm_fmadd_ps(Fused multiply-add)浮点乘加操作，能够在一次操作中同时完成乘法和加法
dst[i]=a[i]*b[i]+c[i]
*/ 
// void fmadd(__m128 a,__m128 b,__m128 c,float* res)
// {
//     __m128 result=_mm_fmadd_ps(a,b,c);
//     // res= (float*)&result;
//     // cout<<res[0]<<res[1]<<res[2]<<res[3]<<endl;    
//     // // 将结果存储到res
//     _mm_store_ps(res,result);
// }
/*按位指令集 如果需要按位NOT 最快的方法可能是与全1 进行xor 示例如下*/
__m128i bitwiseNot(__m128i x){
    // si (signed integer)
    const __m128i zero=_mm_setzero_si128();
    // 得到全1
    const __m128i one=_mm_cmpeq_epi32(zero,zero);
    return _mm_xor_si128(x,one);

} 
/*shuffle
_mm_movehl_ps 将寄存器中的a 中的高64位(即a的后两个浮点数)移动到寄存器b的低64位，而寄存器b的低64位(即前两个浮点数)将被丢弃*/ 

void test(){
    __m128 a=_mm_set_ps(1.0,2.0,3.0,4.0); //[4,3,2,1]
    float* r=(float*)&a;
    cout<<r[0]<<r[1]<<r[2]<<r[3]<<endl;
    __m128 b=_mm_set_ps(8.0,5.0,6.0,7.0);//[7,6,5,8]
    // move high to low
    __m128 res=_mm_movehl_ps(a,b);
    float result[4];
    _mm_store_ps(result,res);
    cout<<result[0]<<result[1]<<result[2]<<result[3]<<endl;
}
int main(){
    vec4 v1,v2;
    // /float v3;
    v1.x=1.0; v1.y=2.0;v1.z=3.0; v1.w=4.0;
    v2.x=2.0;v2.y=2.0;v2.z=2.0;v2.w=2.0;
    // v3=mul(&v1,&v2);
    // cout<<v3<<endl;
    
    // simd_dot(&v1,&v2);


    __m128 a=_mm_set_ps(1.0,2.0,3.0,4.0);
    __m128 b=_mm_set_ps(2.0,2.0,2.0,2.0);
    __m128 c=_mm_set_ps(3.0,3.0,3.0,3.0);
    float result;
    // fmadd(a,b,c,&result);
    test();

    return 0;
}