#ifndef GASTONANY
#define GASTONANY
inline bool any(const std::vector<bool> x) {
  for(int i = 0; i < x.size(); i++)
    if(x[i]) return true;
  return false;
}
#endif
