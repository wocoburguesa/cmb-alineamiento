#include <iostream>
#include "star.h"
#include "aligner_local.h"

using namespace std;

int main(int argc, char *argv[])
{
  //  GlobalAligner a("ATTGCCATT", "ATCCAATTTT", 1, 1);
  vector<string> seqs;
  seqs.push_back("ATTGCCATT");
  seqs.push_back("ATGGCCATT");
  seqs.push_back("ATCCAATTTT");
  seqs.push_back("ATCTTCTT");
  seqs.push_back("ACTGACC");
  StarAligner a(seqs);
  a.run();
  a.print_matches();
  return 0;
}
