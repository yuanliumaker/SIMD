#include <iostream>
#include <cassert>
#define t8_dtri_maxlevel 29
#ifdef t8_tet
    int t8_dtri_dim=3;
#else
    int t8_dtri_dim=2;
#endif 
typedef struct t8_dtri
{
  int8_t level;
  int8_t type;
  int32_t x,y,z; //anchor node
}t8_dtri_t;

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
/***********************************************
 *              预计算数据结构
 ***********************************************/
/* 预计算每个层级的位掩码h和位移量 */
static uint32_t t8_pre_h[30];       // h = 1 << (29 - level)
static uint64_t t8_pre_shift[30];   // 各层级的位移量 (2D: 2*level, 3D: 3*level)

/* 联合类型映射表结构体 (合并iloc和parent_type) */
typedef struct { 
    uint8_t iloc; 
    uint8_t parent_type; 
} t8_type_entry;

#ifdef t8_tet
static const t8_type_entry t8_type_map[6][8] = {
    [0] = {{0,0}, {1,0}, {1,0}, {4,0}, {1,0}, {4,0}, {4,0}, {7,0}},
    [1] = {{0,0}, {1,1}, {2,1}, {5,1}, {2,1}, {5,1}, {4,1}, {7,1}},
    [2] = {{0,2}, {2,2}, {3,2}, {4,2}, {1,2}, {6,2}, {5,2}, {7,2}},
    [3] = {{0,1}, {3,3}, {1,3}, {5,3}, {2,3}, {4,3}, {6,3}, {7,3}},
    [4] = {{0,4}, {2,4}, {2,4}, {6,4}, {3,4}, {5,4}, {5,4}, {7,4}},
    [5] = {{0,5}, {3,5}, {3,5}, {6,5}, {3,5}, {6,5}, {6,5}, {7,5}}
};
#else
static const t8_type_entry t8_type_map[2][4] = {
    [0] = {{0,0}, {1,0}, {1,1}, {3,0}},
    [1] = {{0,0}, {2,1}, {2,1}, {3,1}}
};
#endif

/***********************************************
 *              初始化函数
 ***********************************************/
__attribute__((constructor)) 
static void t8_global_init() {
    for (int l = 0; l < 30; l++) {
        t8_pre_h[l] = 1u << (29 - l); 
        t8_pre_shift[l] = l * t8_dtri_dim; 
    }
}

/***********************************************
 *              SIMD优化计算函数
 ***********************************************/
/* 使用SIMD并行计算多个层级的CubeID (Level Blocking) */
static inline void compute_cubeid_batch(const t8_dtri_t *t, int start, int end, 
                                       int8_t cid_buf[static end-start+1]) {
    for (int l = start; l <= end; l++) {
        const uint32_t h = t8_pre_h[l];
        int8_t id = 0;
        id |= !!(t->x & h) << 0;
        id |= !!(t->y & h) << 1;
#ifdef t8_tet
        id |= !!(t->z & h) << 2;
#endif
        cid_buf[l - start] = id;
    }
}

/***********************************************
 *              类型预计算缓存
 ***********************************************/
/* 类型传播预计算缓存 (type_path[level][type] = ancestor_type) */
static int8_t type_path[30][6];  // 最大支持30层和6种类型

static void precompute_type_path() {
    for (int t = 0; t < (sizeof(t8_type_map)/sizeof(t8_type_map[0])); t++) {
        type_path[0][t] = t;  // 根层级类型不变
        for (int l = 1; l < 30; l++) {
            // 模拟类型向上传播过程
            int8_t cid = 0;  // 假设默认cid=0
            type_path[l][t] = t8_type_map[type_path[l-1][t]][cid].parent_type;
        }
    }
}

/***********************************************
 *              优化后的线性ID计算
 ***********************************************/
uint64_t t8_dtri_linear_id(const t8_dtri_t *t, int level) {
    assert(level <= t8_dtri_maxlevel);
    
    /* 阶段1: 处理超层级填充 */
    uint64_t id = 0;
    const int my_level = t->level;
    int8_t type_temp = (level > my_level) ? t->type : compute_type(t, level);
    if (level > my_level) {
        id = (uint64_t)t->consecutive_id << ((level - my_level) * t8_dtri_dim);
        level = my_level;
    }

    /* 阶段2: SIMD批量计算CubeID (每次处理4个层级) */
    const int BLOCK_SIZE = 4;
    int8_t cid_buf[BLOCK_SIZE];
    int exponent = 0;

    for (int i = level; i > 0; ) {
        const int block_start = (i >= BLOCK_SIZE) ? i - BLOCK_SIZE + 1 : 1;
        const int block_size = i - block_start + 1;
        
        /* SIMD批量计算CubeID */
        compute_cubeid_batch(t, block_start, i, cid_buf);
        
        /* 处理层级块 */
        for (int j = block_size - 1; j >= 0; j--) {
            const int l = block_start + j;
            const int8_t cid = cid_buf[j];
            const t8_type_entry *entry = &t8_type_map[type_temp][cid];
            
            /* 更新ID和类型 */
            id |= (uint64_t)entry->iloc << t8_pre_shift[exponent];
            exponent += t8_dtri_dim;
            type_temp = entry->parent_type;
        }
        i -= block_size;
    }

    return id;
}