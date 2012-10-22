#include <iostream>
#include "progressive.h"

using namespace std;

int main(int argc, char *argv[])
{
  vector<string> seqs;
  string in;
  bool good = 0, show_matrices = 0;
  int counter = 0;
  while(!good){
    seqs.clear();
    cout <<
      "Ingrese una serie de cadenas (por lo menos 3 y separadas por espacio o por línea). Termine la entrada escribiendo el número 0." << endl;
    while(cin >> in && in != "0"){
      seqs.push_back(in);
      ++counter;
    }
    if(counter > 2)
      good = 1;
  }

  
  cout << "¿Desea que se muestren las matrices de distancias en cada paso? (s/n):" << endl;
  cin >> in;
  if(in == "s" || in == "S")
    show_matrices = 1;
  else;
  ProgAligner a(seqs, show_matrices);
  a.run();
  a.print_matches();
  return 0;
}
