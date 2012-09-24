#include <iostream>
#include "aligner.h"
#include "aligner_local.h"

using namespace std;

int main(int argc, char *argv[])
{
  LocalAligner a("ACCGATG", "CAACGTCGATGTCA", 1);
  a.print_matches();
  return 0;
}
