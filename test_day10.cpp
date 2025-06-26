#include <iostream>
#include <cassert>
#include <thread>
#include <mutex>
#include <omp.h>
#include <chrono>
#define t8_dtri_dim 3
typedef struct t8_dtri
{
  int8_t level;
  int8_t type;
  int32_t x,y,z; //anchor node
}t8_dtri_t;
#ifdef t8_tet
#define t8_dtri_type_cid_to_Iloc t8_dtet_type_cid_to_Iloc
#define t8_dtri_cid_type_to_parenttype t8_dtet_cid_type_to_parenttype
#endif
using namespace std;
/* 预计算数据结构定义 */
typedef uint8_t t8_dtri_type_t;
const int T8_DTRI_MAX_PRECOMP_LEVEL = 29; // 支持最大32层预计算

// 类型传递预计算表 [level][cid][current_type] -> next_type
alignas(64) static t8_dtri_type_t pre_parent_type[T8_DTRI_MAX_PRECOMP_LEVEL][4][2];

// Iloc值预计算表 [level][type][cid] -> Iloc 
alignas(64) static uint8_t pre_Iloc[T8_DTRI_MAX_PRECOMP_LEVEL][2][4];

// CubeID掩码预计算 [level] -> h值
alignas(64) static uint32_t pre_cubeid_h[T8_DTRI_MAX_PRECOMP_LEVEL];

/* 初始化预计算表 */
void t8_dtri_precompute() {
    for (int level = 0; level < T8_DTRI_MAX_PRECOMP_LEVEL; ++level) {
        // 预计算h值（论文Algorithm 4.1）
        pre_cubeid_h[level] = (level == 0) ? 0 : (1 << (29 - level));

        // 严格匹配原始类型传递表（论文式46）
        const int8_t parent_type_table[4][2] = {{0,1}, {0,0}, {1,1}, {0,1}};
        for (int cid = 0; cid < 4; ++cid) {
            for (int type = 0; type < 2; ++type) {
                pre_parent_type[level][cid][type] = parent_type_table[cid][type];
            }
        }

        // 严格匹配原始Iloc表（论文表6）
        const uint8_t iloc_table[2][4] = {{0,1,1,3}, {0,2,2,3}};
        for (int type = 0; type < 2; ++type) {
            for (int cid = 0; cid < 4; ++cid) {
                pre_Iloc[level][type][cid] = iloc_table[type][cid];
            }
        }
    }
}
#define T8_DTRI_LEN(l) (1 << (29 - (l)))
const int t8_dtri_cid_type_to_parenttype[4][2] = { 
  { 0, 1 }, 
  { 0, 0 }, 
  { 1, 1 }, 
  { 0, 1 } };
const int t8_dtri_type_cid_to_Iloc[2][4] = { 
  { 0, 1, 1, 3 }, 
  { 0, 2, 2, 3 } };
const int t8_dtet_type_cid_to_Iloc[6][8]= { 
  { 0, 1, 1, 4, 1, 4, 4, 7 }, 
  { 0, 1, 2, 5, 2, 5, 4, 7 }, 
  { 0, 2, 3, 4, 1, 6, 5, 7 },
  { 0, 3, 1, 5, 2, 4, 6, 7 }, 
  { 0, 2, 2, 6, 3, 5, 5, 7 }, 
  { 0, 3, 3, 6, 3, 6, 6, 7 } 
};
const int t8_dtet_cid_type_to_parenttype[8][6] = { 
  { 0, 1, 2, 3, 4, 5 }, 
  { 0, 1, 1, 1, 0, 0 }, 
  { 2, 2, 2, 3, 3, 3 }, 
  { 1, 1, 2, 2, 2, 1 },
  { 5, 5, 4, 4, 4, 5 }, 
  { 0, 0, 0, 5, 5, 5 }, 
  { 4, 3, 3, 3, 4, 4 }, 
  { 0, 1, 2, 3, 4, 5 } };
static int8_t
compute_cubeid (const t8_dtri_t *t, int level)
{
  int8_t id = 0;
  int32_t h;

  /* TODO: assert that 0 < level? This may simplify code elsewhere */

  assert(0 <= level && level <= 29);
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

/* 优化后的线性ID计算函数 */
uint64_t t8_dtri_linear_id_opt(const t8_dtri_t *t, int level) {
    static std::once_flag precompute_flag;
    std::call_once(precompute_flag, t8_dtri_precompute);

    uint64_t id = 0;
    t8_dtri_type_t type_temp = 0;
    int exponent = 0;
    const int my_level = t->level;

    assert(0 <= level && level < T8_DTRI_MAX_PRECOMP_LEVEL);

    /* 处理超层级情况（论文引理19） */
    if (level > my_level) {
        const int level_diff = level - my_level;
        // 修正后的填充模式（论文式37）
        const uint64_t fill_pattern = (t->type == 0) ? 0x5555555555555555ULL : 
                                    (t8_dtri_dim == 3 ? 0x9249249249249249ULL : 
                                                      0xAAAAAAAAAAAAAAAALL);
        id = fill_pattern >> (64 - level_diff * t8_dtri_dim);
        exponent = level_diff * t8_dtri_dim;
        type_temp = t->type;
        level = my_level;
    } else {
        type_temp = compute_type(t, level); // 保留原始类型计算路径
    }

    /* 修正层级索引处理 */
    int current_level = level;
    while (current_level > 0) {
        uint8_t cid = ((t->x & pre_cubeid_h[current_level]) ? 0x1 : 0) |
                     ((t->y & pre_cubeid_h[current_level]) ? 0x2 : 0);
        #ifdef t8_tet
            cid |= ((t->z & pre_cubeid_h[current_level]) ? 0x4 : 0);
        #endif

        // 查表操作严格对齐原始数据
        id |= (uint64_t)pre_Iloc[current_level][type_temp][cid] << exponent;
        exponent += 3;
        type_temp = pre_parent_type[current_level][cid][type_temp];
        
        current_level--;
    }
    cout<<id<<" ";
    return id;
}
int main()
{
    t8_dtri tri;
    tri.level=24;
    tri.type=1;
    tri.x=200; tri.y=128,tri.z=0;
    int level=24;
    auto start=std::chrono::steady_clock::now();
    for(int i=1;i<level;i++){
      uint64_t id=t8_dtri_linear_id_opt(&tri,i);
    
    }
    auto end=std::chrono::steady_clock::now();
    auto diff=end-start;
    double time=std::chrono::duration<double, std::milli>(diff).count()/60000;
    cout<<"time "<<time<<endl;
    return 0;
}