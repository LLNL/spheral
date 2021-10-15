#include <fstream>

void main() {
  fstream* f;
  ios::openmode mode;
  mode = ios::out;
  f = new fstream("teststream.txt", mode);
  *f << mode << endl;
  for (int i = 0; i < 10; ++i) {
    *f << i << " ";
  }
  *f << endl;
  delete f;
}
