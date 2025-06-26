#include <iostream>
#include <x86intrin.h>
#include <cmath>
#include <cassert>
#include <chrono>
using namespace std;
static inline double t8_vec_dist(const double vec_x[3],const double vec_y[3])
{
    __m256d vec_x_part=_mm256_set_pd(0.0,vec_x[2],vec_x[1],vec_x[0]);

 
    __m256d vec_y_part=_mm256_set_pd(0.0,vec_y[2],vec_y[1],vec_y[0]);
  // 计算diff
  __m256d vec_diff_part=_mm256_sub_pd(vec_x_part,vec_y_part);
  // 计算平方和
  __m256d dist=_mm256_mul_pd(vec_diff_part,vec_diff_part);
  double res[4];
  _mm256_storeu_pd(res,dist);
  double result=sqrt(res[0]+res[1]+res[2]+res[3]);
  return result;
}
static inline void
t8_vec_copy(const double vec_in[3],double vec_out[3])
{
  // 创建32字节对齐的临时数组
  alignas(32) double temp_in[4]={vec_in[0],vec_in[1],vec_in[2],0.0};
  __m256d vec=_mm256_load_pd(temp_in);
  // 将寄存器的前3个值存储到vec_out
  // _mm256_set_epi64x(long long  a,long long b,long long c,long long d)
  _mm256_maskstore_pd(vec_out,_mm256_set_epi64x(0,-1,-1,-1),vec);
}
/*double _mm256_cvtsd_f64(__m256d a)
作用：用于从256位双精度浮点数向量中提取最低有效的64位双精度浮点数，
并将其转换为标量双精度浮点数
cvtsd convert scalar double */ 
/*_mm256d _mm256_permute_pd(__256d a,int imm8)
作用：用于对256位双精度浮点数向量中的元素进行重新排列，根据掩码来指定如何重排元素
例如 _256d a=_mm256_set_pd(4.0,3.0,2.0,1.0);
_mm256_permute_pd(a,0x5)
则0x5=0101 
即第0个元素的位置设置为原第1个元素
第1个元素的位置设置为原第0个元素
第2个元素的位置设置为原第1个元素
第3个元素的位置设置为原第0个元素
类似的 _mm256_permutex64_pd
 */ 
double test1()
{
  double vec[3]={1.0,2.0,3.0};
  __m256d vec_part=_mm256_set_pd(0.0,vec[2],vec[1],vec[0]);
  __m256d vec_norm=_mm256_mul_pd(vec_part,vec_part);
  // 提取并累加平方值
  double norm = _mm256_cvtsd_f64(vec_norm) + _mm256_cvtsd_f64(_mm256_permute_pd(vec_norm, 0x5)) +
                  _mm256_cvtsd_f64(_mm256_permute4x64_pd(vec_norm, 0x4E));
  cout<<norm<<endl;
  return sqrt(norm);
}
/*_mm256_maskload_pd()
作用：根据掩码中的位来加载相应的元素
在计算机系统中，整数通常使用二进制补码表示法
负数的表示方法是对该数的绝对值取反，然后加1
如-1的二进制表示如下
首先将1 的二进制表示（00000001）取反得到11111110
然后加1 得到11111111
拓展到64位整数，-1 的所有位都是1 即0xFFFFFFFFFFFFFFFF
*/ 
static inline void 
t8_vec_ax (double vec_x[3],const double alpha)
{
  __m256d vec=_mm256_maskload_pd(vec_x,_mm256_set_epi64x(0,-1,-1,-1));
  __m256d scalar=_mm256_set1_pd(alpha);
  vec=_mm256_mul_pd(vec,scalar);
  _mm256_maskstore_pd(vec_x,_mm256_set_epi64x(0,-1,-1,-1),vec);
  cout<<vec_x[0]<<" "<<vec_x[1]<<" "<<vec_x[2]<<endl;
}
static inline void
t8_vec_axy(const double vec_x[3],double vec_y[3],const double alpha)
{
  cout<<vec_y[0]<<" "<<vec_y[1]<<" "<<vec_y[2]<<endl;
  __m256d vec=_mm256_maskload_pd(vec_x,_mm256_set_epi64x(0,-1,-1,-1));
  __m256d scalar=_mm256_set1_pd(alpha);
  __m256d vec1=_mm256_mul_pd(vec,scalar);
  _mm256_maskstore_pd(vec_y,_mm256_set_epi64x(0,-1,-1,-1),vec1);
  cout<<vec_y[0]<<" "<<vec_y[1]<<" "<<vec_y[2]<<endl;
}
static inline void t8_vec_axb(const double vec_x[3],double vec_y[3],const double alpha,const double b)
{
  cout<<vec_x[0]<<" "<<vec_x[1]<<" "<<vec_x[2]<<endl;
  __m256d vec=_mm256_maskload_pd(vec_x,_mm256_set_epi64x(0,-1,-1,-1));
  __m256d alpha_=_mm256_set1_pd(alpha);
  __m256d b_=_mm256_set1_pd(b);
  vec=_mm256_fmadd_pd(vec,alpha_,b_);
  _mm256_maskstore_pd(vec_y,_mm256_set_epi64x(0,-1,-1,-1),vec);
  cout<<vec_y[0]<<" "<<vec_y[1]<<" "<<vec_y[2]<<endl;

}
static inline void 
t8_vec_axpy(const double vec_x[3],double vec_y[3],const double alpha)
{
  __m256d vec=_mm256_maskload_pd(vec_x,_mm256_set_epi64x(0,-1,-1,-1));
  __m256d alpha_=_mm256_set1_pd(alpha);
  __m256d y_=_mm256_maskload_pd(vec_y,_mm256_set_epi64x(0,-1,-1,-1));
  vec=_mm256_fmadd_pd(vec,alpha_,y_);
  _mm256_maskstore_pd(vec_y,_mm256_set_epi64x(0,-1,-1,-1),vec);  
  cout<<vec_y[0]<<" "<<vec_y[1]<<" "<<vec_y[2]<<endl;


}
static inline void t8_vec_axpyz(const double vec_x[3],double vec_y[3],double vec_z[3],const double alpha)
{
  __m256d x_=_mm256_maskload_pd(vec_x,_mm256_set_epi64x(0,-1,-1,-1));
  __m256d alpha_=_mm256_set1_pd(alpha);
  __m256d y_=_mm256_maskload_pd(vec_y,_mm256_set_epi64x(0,-1,-1,-1));
  x_=_mm256_fmadd_pd(x_,alpha_,y_);
  _mm256_maskstore_pd(vec_z,_mm256_set_epi64x(0,-1,-1,-1),x_);  
  cout<<vec_z[0]<<" "<<vec_z[1]<<" "<<vec_z[2]<<endl;


}
static inline double 
t8_vec_dot(const double vec_x[3],const double vec_y[3])
{
  __m256d x_=_mm256_maskload_pd(vec_x,_mm256_set_epi64x(0,-1,-1,-1));
  __m256d y_=_mm256_maskload_pd(vec_y,_mm256_set_epi64x(0,-1,-1,-1));
  __m256d dot_=_mm256_mul_pd(x_,y_);
  double dot = _mm256_cvtsd_f64(dot_) + _mm256_cvtsd_f64(_mm256_permute_pd(dot_, 0x5)) +
                  _mm256_cvtsd_f64(_mm256_permute4x64_pd(dot_, 0x4E));
  cout<<dot<<endl;
  return dot;

}
/*计算叉乘，permute 和shuffle 搭配 置换向量顺序*/
static inline void 
t8_vec_cross(const double vec_x[3],const double vec_y[3],double cross[3])
{
  // __m256d x_=_mm256_maskload_pd(vec_x,_mm256_set_epi64x(0,-1,-1,-1));
  // __m256d y_=_mm256_maskload_pd(vec_y,_mm256_set_epi64x(0,-1,-1,-1));
  __m256d x_=_mm256_setr_pd(vec_x[0],vec_x[1],vec_x[2],0.0);
  __m256d y_=_mm256_setr_pd(vec_y[0],vec_y[1],vec_y[2],0.0);
  __m256d x1=_mm256_permute4x64_pd(x_,_MM_SHUFFLE(1,2,0,3));
  double* x1_=(double*)&x1;
  cout<<x1_[0]<<" "<<x1_[1]<<" "<<x1_[2]<<" "<<x1_[3]<<endl;
  __m256d x2=_mm256_permute4x64_pd(x_,_MM_SHUFFLE(2,0,1,3));
  __m256d y1=_mm256_permute4x64_pd(y_,_MM_SHUFFLE(1,2,0,3));
  __m256d y2=_mm256_permute4x64_pd(y_,_MM_SHUFFLE(2,0,1,3));
  double* y2_=(double*)&y2;
  cout<<y2_[0]<<" "<<y2_[1]<<" "<<y2_[2]<<" "<<y2_[3]<<endl;
  
  __m256d a=_mm256_mul_pd(x1,y2);
  __m256d b=_mm256_mul_pd(x2,y1);
  __m256d c=_mm256_sub_pd(a,b);
  double temp[4];
  _mm256_storeu_pd(temp,c);
  cross[0]=temp[3];
  cross[1]=temp[2];
  cross[2]=temp[1];
  
  cout<<cross[0]<<cross[1]<<cross[2]<<endl;
} 
static  inline void
t8_vec_diff(const double vec_x[3],const double vec_y[3], double diff[3])
{
  __m256d x_=_mm256_set_pd(0.0,vec_x[2],vec_x[1],vec_x[0]);
  __m256d y_=_mm256_set_pd(0.0,vec_y[2],vec_y[1],vec_y[0]);
  __m256d diff_=_mm256_sub_pd(x_,y_);
  _mm256_maskstore_pd(diff,_mm256_set_epi64x(0,-1,-1,-1),diff_);
  cout<<diff[0]<<" "<<diff[1]<<" "<<diff[2]<<endl;
}
static inline void
t8_vec_swap(double p1[3], double p2[3])
{
    
    __m256d vec1 = _mm256_maskload_pd(p1, _mm256_set_epi64x(0, -1, -1, -1));
    __m256d vec2 = _mm256_maskload_pd(p2, _mm256_set_epi64x(0, -1, -1, -1));
    __m256d tmp = vec1;
    vec1 = vec2;
    vec2 = tmp;
    _mm256_maskstore_pd(p1, _mm256_set_epi64x(0, -1, -1, -1), vec1);
    _mm256_maskstore_pd(p2, _mm256_set_epi64x(0, -1, -1, -1), vec2);
    cout<<p1[0]<<" "<<p1[1]<<" "<<p1[2]<<endl;
    cout<<p2[0]<<" "<<p2[1]<<" "<<p2[2]<<endl;
}

/*坐标计算*/
typedef struct t8_dtri
{
  int8_t level;
  int8_t type;
  int32_t x,y,z; //anchor node
}t8_dtri_t;
void t8_dtri_compute_interger_coords(const t8_dtri_t *elem,const int vertex,int32_t coordinates[3]) 
{
  int8_t type;
  int ei;
#ifdef t8_tet
  int ej;
#endif
  int32_t h;
  type=elem->type;
  h=(1<<(5-(elem->level)));
#ifndef t8_tet
  ei=type;
#else
  ei=type/2;
  ej=(ei+((type%2==0)? 2:1))%3;
#endif
  coordinates[0]=elem->x;
  coordinates[1]=elem->y;
#ifdef t8_tet 
  coordinates[2]=elem->z;
#endif
  if(vertex==0){
    cout<<"vertex0: "<<coordinates[0]<<" "<<coordinates[1]<<endl;
    return ;
  }
  coordinates[ei]+=h;
  cout<<"vertex1: "<<coordinates[0]<<" "<<coordinates[1]<<" "<<endl;
#ifndef t8_tet
  if(vertex==2)
  {
    coordinates[1-ei]+=h;
    cout<<"vertex2: "<<coordinates[0]<<" "<<coordinates[1]<<" "<<endl;
    return ;
  }
#else
  if (vertex == 2) {
    coordinates[ej] += h;
    cout<<"vertex2: "<<coordinates[0]<<" "<<coordinates[1]<<" "<<coordinates[2]<<endl;
    return;
  }
  if (vertex == 3) {
    coordinates[(ei + 1) % 3] += h;
    coordinates[(ei + 2) % 3] += h;
    cout<<"vertex3 "<<coordinates[0]<<" "<<coordinates[1]<<" "<<coordinates[2]<<endl;
  }
#endif

}
void t8_dtri_compute_integer_coords1(const t8_dtri_t* elem, const int vertex, int32_t coordinates[3]) {
  int8_t type = elem->type;
  int ei;
#ifdef t8_tet
  int ej;
#endif
  int32_t h=(1<<(5-(elem->level)));

#ifndef t8_tet
  ei = type;
  #else
  ei = type / 2;
  ej = (ei + ((type % 2 == 0) ? 2 : 1)) % 3;
#endif
  coordinates[0]=elem->x;
  coordinates[1]=elem->y;
#ifdef t8_tet
  coordinates[2]=elem->z;
#endif
  // Load coordinates into SIMD registers
  __m128i coords = _mm_set_epi32(0, coordinates[1], coordinates[0], 0);//coords[0,coordinates[0],coordinates[1],0]
#ifdef t8_tet
   coords=_mm_insert_epi32(coords,coordinates[2],3);//coords[0,coordinates[0],coordinates[1],coordinates[2]]
#endif
  int32_t a[4];
  if (vertex == 0) {
    _mm_storeu_si128((__m128i*)a, coords);
    coordinates[0]=a[1];
    coordinates[1]=a[2];
    cout<<"vertex0: "<<coordinates[0]<<" "<<coordinates[1]<<endl;                    
    return;
  }

  // Update coordinates based on vertex
#ifndef t8_tet
  coords = _mm_add_epi32(coords, _mm_set_epi32(0, (ei == 1) ? h : 0, (ei == 0) ? h : 0, 0));
#else
  coords=_mm_add_epi32(coords,_mm_set_epi32((ei==2)? h:0,(ei==1)? h:0,(ei==0)? h:0,0));
#endif
  if (vertex == 2) {
#ifndef t8_tet
    coords = _mm_add_epi32(coords, _mm_set_epi32(0, (ei == 0) ? h : 0, (ei == 1) ? h : 0, 0));
    _mm_store_si128((__m128i*)a,coords);
    coordinates[0]=a[1];
    coordinates[1]=a[2];
#else
    coords = _mm_add_epi32(coords, _mm_set_epi32((ej == 2) ? h : 0, (ej == 1) ? h : 0, (ej==0)? h:0, 0));
    _mm_store_si128((__m128i*)a,coords);
    coordinates[0]=a[1];
    coordinates[1]=a[2];
    coordinates[2]=a[3];
#endif
    cout<<"vertex2: "<<coordinates[0]<<" "<<coordinates[1]<<" "<<coordinates[2]<<" "<<endl;
    return;
  }

#ifdef t8_tet
  if (vertex == 3) {
    coords = _mm_add_epi32(coords, _mm_set_epi32(h, h, 0, 0)); // Add h to both coordinates
  }
#endif
  _mm_storeu_si128((__m128i*)a, coords);
  coordinates[0]=static_cast<int32_t>(a[1]);
  coordinates[1]=static_cast<int32_t>(a[2]);
  cout<<coordinates[0]<<" "<<coordinates[1]<<endl;
#ifdef t8_tet
  coordinates[2]=static_cast<int32_t>(a[3]);
#endif
  cout<<"vertex3 "<<coordinates[0]<<" "<<coordinates[1]<<" "<<coordinates[2]<<endl;
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
}

#define t8_dtri_maxlevel 29
#define T8_DTRI_LEN(l) (1 << (t8_dtri_maxlevel - (l)))
#ifndef t8_tet
  int t8_dtri_dim=2;
#else
  int t8_dtri_dim=3;
  #define t8_dtri_type_cid_to_Iloc t8_dtet_type_cid_to_Iloc
  #define t8_dtri_cid_type_to_parenttype t8_dtet_cid_type_to_parenttype
#endif
const int t8_dtri_cid_type_to_parenttype[4][2] = { 
  { 0, 1 }, 
  { 0, 0 }, 
  { 1, 1 }, 
  { 0, 1 } };
const int t8_dtri_type_cid_to_Iloc[2][4] = { 
  { 0, 1, 1, 3 }, 
  { 0, 2, 2, 3 } };
const int t8_dtet_cid_type_to_parenttype[8][6] = { 
  { 0, 1, 2, 3, 4, 5 }, 
  { 0, 1, 1, 1, 0, 0 }, 
  { 2, 2, 2, 3, 3, 3 }, 
  { 1, 1, 2, 2, 2, 1 },
  { 5, 5, 4, 4, 4, 5 }, 
  { 0, 0, 0, 5, 5, 5 }, 
  { 4, 3, 3, 3, 4, 4 }, 
  { 0, 1, 2, 3, 4, 5 } };
const int t8_dtet_type_cid_to_Iloc[6][8]= { 
  { 0, 1, 1, 4, 1, 4, 4, 7 }, 
  { 0, 1, 2, 5, 2, 5, 4, 7 }, 
  { 0, 2, 3, 4, 1, 6, 5, 7 },
  { 0, 3, 1, 5, 2, 4, 6, 7 }, 
  { 0, 2, 2, 6, 3, 5, 5, 7 }, 
  { 0, 3, 3, 6, 3, 6, 6, 7 } 
};
static int8_t
compute_cubeid (const t8_dtri_t *t, int level)
{
  int8_t id = 0;
  int32_t h;

  /* TODO: assert that 0 < level? This may simplify code elsewhere */

  assert(0 <= level && level <= t8_dtri_maxlevel);
  h = T8_DTRI_LEN (level);

  if (level == 0) {
    return 0;
  }

  id |= ((t->x & h) ? 0x01 : 0);
  id |= ((t->y & h) ? 0x02 : 0);
#ifdef t8_tet
  id |= ((t->z & h) ? 0x04 : 0);
#endif

  return id;
}
 
static int8_t
compute_type_ext (const t8_dtri_t *t, int level, int8_t known_type, int known_level)
{
  int8_t type = known_type;
  int8_t cid;
  int i;

  assert (0 <= level && level <= known_level);
  assert (known_level <= t->level);
  if (level == known_level) {
    return known_type;
  }
  if (level == 0) {
    /* TODO: the type of the root tet is hardcoded to 0
     *       maybe once we want to allow the root tet to have different types */
    return 0;
  }
  for (i = known_level; i > level; i--) {
    cid = compute_cubeid (t, i);
    /* compute type as the type of T^{i+1}, that is T's ancestor of level i+1 */
    type = t8_dtri_cid_type_to_parenttype[cid][type];
  }
  return type;
}

static int8_t
compute_type (const t8_dtri_t *t, int level)
{
  return compute_type_ext (t, level, t->type, t->level);
}

uint64_t
t8_dtri_linear_id (const t8_dtri_t *t, int level)
{
  uint64_t id = 0;
  int8_t type_temp = 0;
  int8_t cid;
  int i;
  int exponent;
  int my_level;

  assert (0 <= level && level <= t8_dtri_maxlevel);
  my_level = t->level;
  exponent = 0;
  /* If the given level is bigger than t's level
   * we first fill up with the ids of t's descendants at t's
   * origin with the same type as t */
  if (level > my_level) {
    exponent = (level - my_level) * t8_dtri_dim;
    type_temp = t->type;
    level = my_level;
  }
  else {
    type_temp = compute_type (t, level);
  }
  for (i = level; i > 0; i--) {
    cid = compute_cubeid (t, i);
    id |= ((uint64_t) t8_dtet_type_cid_to_Iloc[type_temp][cid]) << exponent;
    exponent += t8_dtri_dim; /* multiply with 4 (2d) resp. 8  (3d) */
    type_temp = t8_dtri_cid_type_to_parenttype[cid][type_temp];
  }
  cout<<static_cast<unsigned long long>(id)<<" ";
  return id;
}

int main()
{

    double a[3]={2.0,3.0,4.0};
    double b[3]={1.0,2.0,5.0};
    // double res=t8_vec_dist(a,b);
    // t8_vec_copy(a,b);
    // cout<<b[0]<<b[1]<<b[2];
    // cout<<res<<endl;
    // test1();
    // t8_vec_ax(a,3);
    // t8_vec_axy(a,b,3);
    // t8_vec_axb(a,b,3,3);
    // t8_vec_axpy(a,b,3);
    double c[3];
    // t8_vec_axpyz(a,b,c,3);
    // t8_vec_dot(a,b);
    // t8_vec_cross(a,b,c);
    // t8_vec_diff(a,b,c);
    // t8_vec_swap(a,b);
    t8_dtri tri;
    tri.level=24;
    tri.type=1;
    tri.x=200; tri.y=128,tri.z=3;
    const int vertex=1;
    int32_t coordinates[3];
    // t8_dtri_compute_interger_coords(&tri,vertex,coordinates);
    // t8_dtri_compute_integer_coords1(&tri,vertex,coordinates);
    double ref_coords[9]={0.0,0.0,0.0,1.0,0.0,0.0,1.0,1.0,0.0};
    double out_coords[9]={0.0};
    // t8_dtri_compute_reference_coords(&tri,ref_coords,3,0,out_coords);
    // if_tet 四面体
    // t8_dtri_compute_reference_coords(&tri,ref_coords,3,out_coords);

    
    // for(int i=0;i<9;i++)
    // {
    //   cout<<out_coords[i]<<" ";
    // }
    // cout<<endl;
    int level=24;
    auto start=std::chrono::steady_clock::now();
    for(int i=1;i<level;i++){
      uint64_t id=t8_dtri_linear_id(&tri,i);
    }
    auto end=std::chrono::steady_clock::now();
    auto diff=end-start;
    double time=std::chrono::duration<double, std::milli>(diff).count()/60000;
    cout<<"time "<<time<<endl;

    return 0;
    
}