#include <iostream>
#include <string.h>
#include <cassert>
using namespace std;
#define YH_DTRI_CHILDREN 4
#define YH_DTRI_MAXLEVEL 29
#define YH_DTRI_LEN(l) (1 << (YH_DTRI_MAXLEVEL - (l)))
typedef  int8_t yh_dtri_type_t;
typedef int32_t yh_dtri_coord_t;
typedef int8_t yh_dtri_cube_id_t;
const int yh_dtri_cid_type_to_parenttype[4][2] = { 
  { 0, 1 }, 
  { 0, 0 }, 
  { 1, 1 }, 
  { 0, 1 } };
const int yh_dtri_type_cid_to_Iloc[2][4] = { 
  { 0, 1, 1, 3 }, 
  { 0, 2, 2, 3 } };
const int yh_dtri_parenttype_Iloc_to_type[2][4] = { 
  { 0, 0, 1, 0 }, 
  { 1, 0, 1, 1 } };
const int yh_dtri_parenttype_Iloc_to_cid[2][4] = { 
  { 0, 1, 1, 3 }, 
  { 0, 2, 2, 3 } };
typedef struct yh_dtri
{
  int8_t level;
  yh_dtri_type_t type;
  yh_dtri_coord_t x, y,z;
} yh_dtri_t;

void yh_dtri_copy(const yh_dtri_t *t,yh_dtri_t *dest)
{
    if(t==dest){
        return;
    }
    memcpy(dest,t,sizeof(yh_dtri_t));
}
static yh_dtri_cube_id_t
compute_cubeid (const yh_dtri_t *t, int level)
{
  yh_dtri_cube_id_t id = 0;
  yh_dtri_coord_t h;

  /* TODO: assert that 0 < level? This may simplify code elsewhere */

  assert (0 <= level && level <= YH_DTRI_MAXLEVEL);
  h = YH_DTRI_LEN (level);

  if (level == 0) {
    return 0;
  }

  id |= ((t->x & h) ? 0x01 : 0);
  id |= ((t->y & h) ? 0x02 : 0);
#ifdef YH_DTRI_TO_DTET
  id |= ((t->z & h) ? 0x04 : 0);
#endif
  return id;
}
static yh_dtri_type_t
compute_type_ext (const yh_dtri_t *t, int level, yh_dtri_type_t known_type, int known_level)
{
  int8_t type = known_type;
  yh_dtri_cube_id_t cid;
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
    type = yh_dtri_cid_type_to_parenttype[cid][type];
  }
  return type;
}
static yh_dtri_type_t
compute_type (const yh_dtri_t *t, int level)
{
  return compute_type_ext (t, level, t->type, t->level);
}
static void
yh_dtri_succ_pred_recursion (const yh_dtri_t *t, yh_dtri_t *s, int level, int increment)
{
  yh_dtri_type_t type_level, type_level_p1;
  yh_dtri_cube_id_t cid;
  int local_index;
  int sign;

  /* We exclude the case level = 0, because the root triangle does
   * not have a successor. */
  assert(1 <= level && level <= t->level);
  assert(-YH_DTRI_CHILDREN < increment && increment < YH_DTRI_CHILDREN);

  if (increment == 0) {
    yh_dtri_copy (t, s);
    return;
  }
  cid = compute_cubeid (t, level);
  type_level = compute_type (t, level);
  // 当前t的局部索引
  local_index = yh_dtri_type_cid_to_Iloc[type_level][cid];
  // t的前驱或后继s的局部索引
  local_index = (local_index + YH_DTRI_CHILDREN + increment) % YH_DTRI_CHILDREN;
  if (local_index == 0) {
    sign = increment < 0 ? -1 : increment > 0;
    // 此时s为t的前驱且s在t 的上一层级，因此到上一层级上递归生成
    yh_dtri_succ_pred_recursion (t, s, level - 1, sign);
    type_level_p1 = s->type; /* We stored the type of s at level-1 in s->type */
  }
  else {
    type_level_p1 = yh_dtri_cid_type_to_parenttype[cid][type_level];
  }
  type_level = yh_dtri_parenttype_Iloc_to_type[type_level_p1][local_index];
  cid = yh_dtri_parenttype_Iloc_to_cid[type_level_p1][local_index];
  s->type = type_level;
  s->level = level;
  /* Set the x,y(,z) coordinates at level to the cube-id. */
  /* TODO: check if we set the correct bits here! */
  s->x = (cid & 1 ? s->x | 1 << (YH_DTRI_MAXLEVEL - level) : s->x & ~(1 << (YH_DTRI_MAXLEVEL - level)));
  s->y = (cid & 2 ? s->y | 1 << (YH_DTRI_MAXLEVEL - level) : s->y & ~(1 << (YH_DTRI_MAXLEVEL - level)));
#ifdef YH_DTRI_TO_DTET
  s->z = (cid & 4 ? s->z | 1 << (YH_DTRI_MAXLEVEL - level) : s->z & ~(1 << (YH_DTRI_MAXLEVEL - level)));
#endif

}
int main()
{   yh_dtri tri,s;
    tri.level=5;
    tri.type=1;
    tri.x=0; tri.y=0,tri.z=0;
    int level=5;
    yh_dtri_succ_pred_recursion(&tri,&s,level,1);
    cout<<static_cast<int>(s.level)<<endl;
    cout<<static_cast<int>(s.type)<<endl;
    cout<<static_cast<int>(s.x)<<" "<<static_cast<int>(s.y)<<endl;
    return 0;
}