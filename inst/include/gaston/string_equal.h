#include <string>

#ifndef STRING_EQUAL
#define STRING_EQUAL
inline bool string_equal(SEXP a, SEXP b) {
  return !(std::strcmp(CHAR(a), CHAR(b)));
}

inline bool string_equal(char * a, SEXP b) {
  return !(std::strcmp(a, CHAR(b)));
}

inline bool string_equal(std::string & a, SEXP b) {
  return (a == CHAR(b));
}
#endif

/*
bool string_equal(std::string & a, std::string & b) {
  return (a.compare(b) == 0);
}

bool string_equal(SEXP a, std::string & b) {
  return (CHAR(a) == b);
}

bool string_equal(char * a, std::string & b) {
  return (b == a);
}


bool string_equal(std::string & a, char * b) {
  return (a == b);
}

bool string_equal(SEXP a, char * b) {
  return !(std::strcmp(CHAR(a), b));
}

bool string_equal(char * a, char * b) {
  return !std::strcmp(a, b);
}
*/
