#include <iostream>
#include <x86intrin.h>
#include <omp.h>
using namespace std;

typedef struct t8_dtri
{
  int8_t level;
  int8_t type;
  int32_t x,y,z; //anchor node
}t8_dtri_t;
const int t8_dtri_parenttype_Iloc_to_cid[2][4] = { 
  { 0, 1, 1, 3 }, 
  { 0, 2, 2, 3 } };
const int t8_dtri_parenttype_Iloc_to_type[2][4] = { 
  { 0, 0, 1, 0 }, 
  { 1, 0, 1, 1 } };

void 
t8_dtri_init_linear_id(t8_dtri_t *t, uint64_t id,int level)
{
    int i;
    int offset_coords,offset_index;
    const int children_m1=3;
    uint64_t local_index;
    int8_t cid;
    int8_t type;
    t->level=level;
    t->x=0;
    t->y=0;
#ifdef t8_tet
    t->z=0;
#endif
    type=0;
    int local_x=0;
    int local_y=0;
    int local_z=0;
    /*#pragma omp parallel for: 创建一个并行区域，并将接下来的for循环分配给多个线程执行
    private(list)指定在并行区域中每个线程都有自己的独立的变量副本，
    shared(list) 指定在并行区域中多个线程共享同一个变量，这些变量在线程之间共享数据，确保一致性*/ 
    #pragma omp parallel for private(offset_coords,offset_index,local_index,cid) shared(t,id,level,type) reduction(|:local_x,local_y,local_z)
    for(i=1;i<=level;i++){
        offset_coords=5-i;
        offset_index=level-i;
        local_index=(id>>(2*offset_index))&children_m1;
        cid = t8_dtri_parenttype_Iloc_to_cid[type][local_index];
        type = t8_dtri_parenttype_Iloc_to_type[type][local_index];
        // 原子操作，确保特定的内存操作以原子方式执行，原子操作开销大
        // 从而避免多线程环境下的竞态条件
        // #pragma omp atomic
        local_x |= (cid & 1) ? 1 << offset_coords : 0;
        // #pragma omp atomic
        local_y |= (cid & 2) ? 1 << offset_coords : 0;
#ifdef t8_tet
        // #pragma omp atomic
        local_z |=(cid & 4)? 1<<offset_coords:0;
#endif
    }
    t->x = local_x;
    t->y = local_y;
#ifdef t8_tet
    t->z = local_z;
#endif 
    t->type=type;
    cout<<"t x y z "<<t->x<<" "<<t->y<<" "<<t->z<<endl;
    cout<<"t type "<<static_cast<int>(t->type)<<endl;

}
static int8_t 
compute_cubeid(const t8_dtri_t *t,int level)
{
    int8_t id=0;
    int32_t h;
    h=(1<<(5-(t->level)));
    if(level==0)
    {
        return 0;
    }
    id |= ((t->x & h) ? 0x01 : 0);
    
    id |= ((t->y & h) ? 0x02 : 0);
#ifdef t8_tet
    id |= ((t->z & h) ? 0x04 : 0);
#endif
    cout<<"id "<<+id<<endl;
  return id;
}
static int8_t compute_cubeid1(const t8_dtri_t *t, int level) {
    int8_t id = 0;
    int32_t h = (1 << (5 - (t->level)));
    if (level == 0) {
        return 0;
    }

    __m128i t_xyz = _mm_set_epi32(t->z, t->y, t->x, 0);
    __m128i h_vec = _mm_set1_epi32(h);
    __m128i mask = _mm_and_si128(t_xyz, h_vec);
    mask = _mm_cmpeq_epi32(mask, h_vec);

    id |= (_mm_extract_epi16(mask, 2) ? 0x01 : 0);
    id |= (_mm_extract_epi16(mask, 3) ? 0x02 : 0);
#ifdef t8_tet
    id |= (_mm_extract_epi16(mask, 4) ? 0x04 : 0);
#endif
    cout<<"id "<<+id<<endl;
    return id;
    
}


int main()
{   t8_dtri tri;
    tri.level=2;
    tri.type=1;
    tri.x=3; tri.y=4,tri.z=5;
    uint64_t id=3;
    // t8_dtri_init_linear_id(&tri,id,2);
    int8_t id1=compute_cubeid(&tri,1);
    int8_t id2=compute_cubeid1(&tri,1);
    return 0;
}