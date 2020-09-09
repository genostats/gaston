#include<string>

#ifndef FLIP_STRAND
#define FLIP_STRAND

inline char flip_strand(char x) {
  if(x == 'A') return 'T';
  if(x == 'C') return 'G';
  if(x == 'G') return 'C';
  if(x == 'T') return 'A';
  return x;
}

inline std::string flip_strand(const char * str) {
  std::string s;
  while(*str) {
    s += flip_strand(*str++);
  }
  return s;
}

inline std::string flip_strand(const SEXP a) {
  return flip_strand(CHAR(a));
}

inline std::string flip_strand(const std::string & a) {
  return flip_strand(a.c_str());
}
#endif
