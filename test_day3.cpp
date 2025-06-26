#include <iostream>
#include <x86intrin.h>
using namespace std;
/*avx2*/

/*1 __256i _mm256_abs_epi16(__m256i a)
返回一个类型为__m256i的256位寄存器，其中包含16位整数的绝对值
类似的 _mm256_abs_epi32
_mm256_abs_epi8
easy 略*/
/*2. _mm256_add_epi16(__m256i a,__m256i b)
     _mm256_add_epi32
     _mm256_add_epi64
     _mm256_add_epi8
     两个寄存器相加返回结果*/ 


/*3. _mm256_adds_epi16(__256i a,__m256i b)
作用：用于对256位宽AVX2寄存器中打包的16位整数进行饱和加法运算
饱和加法确保结果不会溢出，而是被截断到可表示的最小值或最大值
类似的：_mm256_adds_epi8 _mm256_adds_epu16 _mm256_adds_epu8*/

void test1()
{
    // 初始化两个包含打包的 16 位整数的 256 位寄存器
    __m256i a = _mm256_set_epi16(30000, -30000, 20000, -20000, 10000, -10000, 5000, -5000, 3000, -3000, 2000, -2000, 1000, -1000, 500, -500);
    __m256i b = _mm256_set_epi16(10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000);

    // 执行饱和加法运算
    __m256i result = _mm256_adds_epi16(a, b);
    int16_t* res=(int16_t*)&result;
    for(int i=0;i<16;i++)
    {
        cout<<res[i]<<endl; //30000+10000 溢出 保留最大值32767
    }

}
/*__256i _mm256_alignr_epi8(__m256i a,__m256i b,int count)
作用：返回一个新的256位整数向量，该向量由从b和a 中提取的数字组成，count 表示从b开始提取的位置

alignr 即align right*/ 
void test2()
{
    int8_t a[32];
    int8_t b[32];
    for(int i=0;i<32;i++)
    {
        a[i]=i;
        b[i]=i+1;
    }
    // reinterpret_cast :指针类型之间的转换，将一种指针类型转换为另一种指针类型
    __m256i a1=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(a));
    __m256i b1=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(b));
    __m256i result=_mm256_alignr_epi8(a1,b1,8);//从b[8]开始提取
    int8_t res[32];
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(res), result);

    for (int i = 0; i < 32; ++i) {
        // 这里int8_t 在编译器和标准库视线中被视为char 类型，为确保能看到数值，
        // 将其转换为一个更大的整数类型
        std::cout << static_cast<int>(res[i]) << " ";
    }
    std::cout << std::endl;
}
/*_mm256_broadcastb_epi8(__m128i a)
作用:将一个8位整数(一个字节)broadcast byte的值广播到256位向量，
其中 a 只使用最低8位的值
类似的：_mm256_broadcastd_epi32
_mm256_broadcastq_epi64
注意上述 epi8 对应b （byte）
epi32 对应d(doubleword)四字节 即32bit
epi64 对应quadword 八字节 64bit
同理128位寄存器的有
__m128d _mm_broadcastsd_pd(__m128d a)
sd即scalar double 
__m128i _mm_broadcastq_epi64等*/
void test3()
{
    int8_t a[16];
    for (int i=0;i<16;i++)
    {
        a[i]=i;
    }
    __m128i a1=_mm_load_si128(reinterpret_cast<__m128i*>(a));
    __m256i result=_mm256_broadcastb_epi8(a1);
    int8_t res[32];
    // res 未对齐 因此storeu
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(res),result);
    for(int i=0;i<32;i++){
        cout<<static_cast<int>(res[i])<<" ";
    }
    
} 
/*__256i _mm256_avg_epu16(__m256i a,__m256i b)
计算a和b 的均值并向上取整*/
void test4()
{
    int16_t a[16],b[16];
    for(int i=0;i<16;i++)
    {
        a[i]=i;
        b[i]=i+1;
    }
    __m256i a1=_mm256_loadu_si256(reinterpret_cast<__m256i*>(a));
    __m256i b1=_mm256_loadu_si256(reinterpret_cast<__m256i*>(b));
    __m256i result=_mm256_avg_epu16(a1,b1);
    int16_t* res=(int16_t*)&result;
    for(int i=0;i<16;i++)
    {
        cout<<static_cast<int>(res[i])<<" ";
    }
} 

int main()
{
    // test1();
    // test2();
    // test3();
    test4();
    return 0;
}