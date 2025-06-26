#include <iostream>
#include <immintrin.h>
using namespace std;
/*SVML (short vector math library)
其中的指令主要用于加速向量化数学运算，如三角函数、指数和对数等*/
/*1 _mm_acos_pd(__m128d a)
计算向量的反余弦值
编译器不支持*/
// void test()
// {
//     __m128d a=_mm_set_pd(-0.5,0.5);
//     __m128d result=_mm_acos_pd(a);


// } 

// 辅助函数：打印 __m128d 向量
void print_m128d(__m128d vec) {
    double* elements = (double*)&vec;
    for (int i = 0; i < 2; ++i) {
        std::cout << elements[i] << " ";
    }
    std::cout << std::endl;
}
void test1()
{
    // 初始化向量
    __m128d a = _mm_set_pd(1.0, 2.0);
    __m128d b = _mm_set_pd(3.0, 4.0);
    __m128d c = _mm_set_pd(5.0, 6.0);

    // 使用 _mm_fmadd_pd 函数计算融合乘加运算 注意：此函数需要在执行命令中添加 -mfma
    __m128d result = _mm_fmadd_pd(a, b, c);

    // 打印结果
    std::cout << "结果: ";
    print_m128d(result);
}
int main() {
    
    test1();
    return 0;
}