#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// 标量版本的 t8_vec_axpy 函数
static inline void
t8_vec_axpy_scalar(const double vec_x[3], double vec_y[3], const double alpha)
{
    for (int i = 0; i < 3; i++) {
        vec_y[i] += alpha * vec_x[i];
    }
}

// SIMD 版本的 t8_vec_axpy 函数
static inline void 
t8_vec_axpy_simd(const double vec_x[3], double vec_y[3], const double alpha)
{
    __m256d alpha_ = _mm256_set1_pd(alpha);
    __m256d x_ = _mm256_maskload_pd(vec_x, _mm256_set_epi64x(0, -1, -1, -1));
    __m256d y_ = _mm256_maskload_pd(vec_y, _mm256_set_epi64x(0, -1, -1, -1));
    __m256d result = _mm256_fmadd_pd(x_, alpha_, y_);
    _mm256_maskstore_pd(vec_y, _mm256_set_epi64x(0, -1, -1, -1), result);
}

// 高分辨率计时器
double get_time() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main() {
    const int iterations = 1000000;
    
    double vec_x[3] = {12345.6789012345, 67890.1234567890, 23456.7890123456};
    double vec_y_scalar[3] = {10000.0000123456, 20000.0000234567, 30000.0000345678};
    double vec_y_simd[3] = {10000.0000123456, 20000.0000234567, 30000.0000345678};
    double alpha = 2.0;

    // 测试标量版本的执行时间
    double start_time = get_time();
    for (int iter = 0; iter < iterations; iter++) {
        t8_vec_axpy_scalar(vec_x, vec_y_scalar, alpha);
    }
    double scalar_time = get_time() - start_time;

    // 测试 SIMD 版本的执行时间
    start_time = get_time();
    for (int iter = 0; iter < iterations; iter++) {
        t8_vec_axpy_simd(vec_x, vec_y_simd, alpha);
    }
    double simd_time = get_time() - start_time;

    // 打印结果
    printf("Scalar version time: %f seconds\n", scalar_time);
    printf("SIMD version time: %f seconds\n", simd_time);

    // 打印最终结果以确保正确性
    for (int i = 0; i < 3; i++) {
        if (vec_y_scalar[i] != vec_y_simd[i]) {
            printf("Mismatch at index %d: scalar = %f, simd = %f\n", i, vec_y_scalar[i], vec_y_simd[i]);
            return -1;
        }
    }

    printf("Results match!\n");

    return 0;
}