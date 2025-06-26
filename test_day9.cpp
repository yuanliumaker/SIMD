#include <cblas.h>
#include <iostream>
#include <cblas.h>
#include <stddef.h>
#include <stdint.h>
using namespace std;
/*坐标计算*/
typedef struct t8_dtri
{
  int8_t level;
  int8_t type;
  int32_t x,y,z; //anchor node
}t8_dtri_t;


void t8_dtri_compute_reference_coords_openblas(const t8_dtri_t *elem, const double* ref_coords, const size_t num_coords,
#ifndef t8_tet
        const size_t skip_coords,
#endif
        double* out_coords)
{
    int8_t type;
    int32_t h;
    type = elem->type;
    h = (1 << (5 - (elem->level)));  // 根据元素的级别计算 h

#ifndef t8_tet
    const int tri_orientation = type;  // 三角形方向
#else
    const int tet_orientation0 = type / 2;
    const int tet_orientation1 = (tet_orientation0 + ((type % 2 == 0) ? 1 : 2)) % 3;
    const int tet_orientation2 = (tet_orientation0 + ((type % 2 == 0) ? 2 : 1)) % 3;
#endif

    // 遍历每个坐标并计算
    for (size_t coord = 0; coord < num_coords; ++coord) {
        // 计算偏移量
#ifndef t8_tet
        const size_t offset = (2 + skip_coords) * coord;
        const size_t offset_3d = 3 * coord;  // 对于三维坐标使用的偏移量
#else
        const size_t offset = 3 * coord;  // 对于四面体使用的偏移量
#endif

        // 初始化输出坐标（x, y, z）
        out_coords[offset + 0] = elem->x;
        out_coords[offset + 1] = elem->y;
#ifdef t8_tet
        out_coords[offset + 2] = elem->z;  // 四面体有 3 个分量
#endif

#ifndef t8_tet
        // 使用 OpenBLAS 进行加速计算，并根据方向控制选择性更新
        // 用 OpenBLAS 进行加速计算 double precision AX+Y
        // 即y[i*incy]=y[i*incy]+alpha*x[i*incx] incy 和incx 分别表示x 和y 的步长
        /*void cblas_daxpy(int n,double alpha,const double* x,int incx,double *y,int incy)
        n:表示向量的元素数量，即x和y 向量的大小
        alpha：即A
        incx: 向量x 的步长，通常是 1，表示每次访问相邻的元素。如果你希望从向量 X 中跳过元素（比如读取 X[0], X[2], X[4] 等），可以设置 incx 为 2
        incy:向量y的步长*/
         cblas_daxpy(1, h, &ref_coords[offset_3d + 0], 1, &out_coords[offset + tri_orientation], 1);
         cblas_daxpy(1, h, &ref_coords[offset_3d + 1], 1, &out_coords[offset + 1 - tri_orientation], 1);
#else
        // 四面体的方向控制并使用 OpenBLAS 加速计算
        cblas_daxpy(1, h, &ref_coords[offset + 0], 1, &out_coords[offset + tet_orientation0], 1);
        cblas_daxpy(1, h, &ref_coords[offset + 1], 1, &out_coords[offset + tet_orientation1], 1);
        cblas_daxpy(1, h, &ref_coords[offset + 2], 1, &out_coords[offset + tet_orientation2], 1);
#endif

        // 归一化操作：除以 T8_DTRI_ROOT_LEN (1<<2)
        out_coords[offset + 0] /= (double)(1 << 2);
        out_coords[offset + 1] /= (double)(1 << 2);
#ifdef t8_tet
        out_coords[offset + 2] /= (double)(1 << 2);
#endif
    }
     for(int i=0;i<9;i++)
        {
        cout<<out_coords[i]<<" ";
        }
        cout<<endl;
}

void
t8_dtri_compute_reference_coords(const t8_dtri_t *elem,const double* ref_coords,const size_t num_coords,
#ifndef t8_tet
        const size_t skip_coords,
#endif
        double* out_coords)
{
  int8_t type;
  int32_t h;
  type=elem->type;
  h=(1<<(5-(elem->level)));
#ifndef t8_tet
  const int tri_orientation=type;
#else
  const int tet_orientation0 = type / 2;
  const int tet_orientation1 = (tet_orientation0 + ((type % 2 == 0) ? 1 : 2)) % 3;
  const int tet_orientation2 = (tet_orientation0 + ((type % 2 == 0) ? 2 : 1)) % 3;
#endif
for (size_t coord = 0; coord < num_coords; ++coord) {
    /* offset defines, how many coordinates to skip in an iteration. */
#ifndef t8_tet
    const size_t offset = (2 + skip_coords) * coord;
    const size_t offset_3d = 3 * coord;
#else
    const size_t offset = 3 * coord;
#endif
    out_coords[offset + 0] = elem->x;
    out_coords[offset + 1] = elem->y;
#ifdef t8_tet
    out_coords[offset + 2] = elem->z;
#endif
#ifndef t8_tet
    out_coords[offset + tri_orientation] += h * ref_coords[offset_3d + 0];
    out_coords[offset + 1 - tri_orientation] += h * ref_coords[offset_3d + 1];
#else
    out_coords[offset + tet_orientation0] += h * ref_coords[offset + 0];
    out_coords[offset + tet_orientation1] += h * ref_coords[offset + 1];
    out_coords[offset + tet_orientation2] += h * ref_coords[offset + 2];

    /* done 3D */
#endif
    /* Since the integer coordinates are coordinates w.r.t to
     * the embedding into [0,T8_DTRI_ROOT_LEN]^d, we just need
     * to divide them by the root length. */
    out_coords[offset + 0] /= (double) (1<<2);
    out_coords[offset + 1] /= (double) (1<<2);
#ifdef t8_tet
    out_coords[offset + 2] /= (double) (1<<2);
#endif

}
    for(int i=0;i<9;i++)
        {
        cout<<out_coords[i]<<" ";
        }
        cout<<endl;
}



void t8_dtri_compute_reference_coords1(const t8_dtri_t *elem, const double* ref_coords, const size_t num_coords,
#ifndef t8_tet
        const size_t skip_coords,
#endif
        double* out_coords)
{
    int8_t type;
    int32_t h;
    type = elem->type;
    h = (1 << (5 - (elem->level)));
#ifndef t8_tet
    const int tri_orientation = type;
    const int dim = 2;
#else
    const int tet_orientation0 = type / 2;
    const int tet_orientation1 = (tet_orientation0 + ((type % 2 == 0) ? 1 : 2)) % 3;
    const int tet_orientation2 = (tet_orientation0 + ((type % 2 == 0) ? 2 : 1)) % 3;
    const int dim = 3;
#endif

    // Initialize out_coords with elem coordinates
    for (size_t coord = 0; coord < num_coords; ++coord) {
        size_t offset = dim * coord;
        out_coords[offset + 0] = elem->x;
        out_coords[offset + 1] = elem->y;
#ifdef t8_tet
        out_coords[offset + 2] = elem->z;
#endif
    }

    // Create a temporary matrix to store scaled ref_coords
    double* scaled_ref_coords = (double*)malloc(num_coords * dim * sizeof(double));
    for (size_t i = 0; i < num_coords * dim; ++i) {
        scaled_ref_coords[i] = h * ref_coords[i];
    }

    // Create orientation matrix
    double orientation_matrix[dim * dim];
#ifndef t8_tet
    orientation_matrix[0] = (tri_orientation == 0) ? 1 : 0;
    orientation_matrix[1] = (tri_orientation == 1) ? 1 : 0;
    orientation_matrix[2] = (tri_orientation == 1) ? 1 : 0;
    orientation_matrix[3] = (tri_orientation == 0) ? 1 : 0;
#else
    orientation_matrix[0] = (tet_orientation0 == 0) ? 1 : 0;
    orientation_matrix[1] = (tet_orientation1 == 0) ? 1 : 0;
    orientation_matrix[2] = (tet_orientation2 == 0) ? 1 : 0;
    orientation_matrix[3] = (tet_orientation0 == 1) ? 1 : 0;
    orientation_matrix[4] = (tet_orientation1 == 1) ? 1 : 0;
    orientation_matrix[5] = (tet_orientation2 == 1) ? 1 : 0;
    orientation_matrix[6] = (tet_orientation0 == 2) ? 1 : 0;
    orientation_matrix[7] = (tet_orientation1 == 2) ? 1 : 0;
    orientation_matrix[8] = (tet_orientation2 == 2) ? 1 : 0;
#endif

    // Use cblas_dgemm to perform matrix multiplication and addition
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, num_coords, dim, dim, 1.0,
                scaled_ref_coords, dim, orientation_matrix, dim, 1.0, out_coords, dim);

    // Normalize out_coords
    for (size_t coord = 0; coord < num_coords * dim; ++coord) {
        out_coords[coord] /= (double)(1 << 2);
    }

    free(scaled_ref_coords);
    for(int i=0;i<9;i++)
        {
        cout<<out_coords[i]<<" ";
        }
        cout<<endl;
}
int main()
{
    t8_dtri tri;
    tri.level=1;
    tri.type=0;
    tri.x=1; tri.y=0,tri.z=0;
    double ref_coords[9]={0.0,0.0,0.0,1.0,0.0,0.0,1.0,1.0,0.0};
    double out_coords[9]={0.0};
    t8_dtri_compute_reference_coords(&tri,ref_coords,3,1,out_coords);
    t8_dtri_compute_reference_coords_openblas(&tri,ref_coords,3,1, out_coords);
    t8_dtri_compute_reference_coords1(&tri,ref_coords,3,1,out_coords);
    
    return 0;
}
