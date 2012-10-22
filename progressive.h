#include <vector>
#include <algorithm>
#include <fstream>
#include "aligner.h"

#define INF 0xffffffff

using namespace std;

class ProgAligner{
 private:
  vector<string> seqs;
  vector<string> matches;
  bool print_matrices;

 public:
  ProgAligner(vector<string> input, bool print){
    seqs = input;
    print_matrices = print;
  }

  void fix_lengths(){
    int max_length = 0;
    for(int i = 0; i < matches.size(); ++i)
      if(matches[i].length() > max_length)
	max_length = matches[i].length();
    
    for(int i = 0; i < matches.size(); ++i){
      while(matches[i].length() < max_length)
	matches[i] += "_";
    }
  }

  float get_distance(string a, string b){
    int diff_counter = 0.0f;
    for(int i = 0; i < a.length(); ++i)
      if(a[i] != b[i] && a[i] != '_' && b[i] != '_')
	diff_counter += 1.0f;
    return diff_counter/(float)a.length();
  }

  vector< vector<float> > get_distance_matrix(vector<string> in){
    vector< vector<float> > dm(in.size(), vector<float>(in.size(), 0.0f));

    for(int i = 0; i < in.size(); ++i)
      for(int j = i; j < in.size(); ++j)
	if(i != j){
	  GlobalAligner aligner(in[i], in[j], 1, 1);
	  dm[i][j] = dm[j][i] = get_distance(aligner.get_matches()[0].first,
					     aligner.get_matches()[0].second);
	}

    return dm;
  }

  vector< vector<float> > get_avg_distance_matrix(vector< vector<float> > dm){
    vector< vector<float> > adm(dm.size(), vector<float>(dm.size(), 0.0f));

    for(int i = 0; i < dm.size(); ++i)
      for(int j = i; j < dm.size(); ++j)
	if(i != j){
	  float sum1 = 0.0f, sum2 = 0.0f;
	  for(int k = 0; k < dm.size(); ++k){
	    if(i != k) 
	      sum1 += dm[i][k];
	    else;
	    if(j != k)
	      sum2 += dm[j][k];
	    else;
	  }
	  float q = dm[i][j] - (1.0f/((float)dm.size() - 2.0f))*(sum1 + sum2);
	  adm[i][j] = adm[j][i] = q;
	}
	else
	  adm[i][j] = INF;

    return adm;
  }

  void run(){
    vector< vector<float> > dm = get_distance_matrix(seqs);
    vector< vector<float> > adm = get_avg_distance_matrix(dm);

    vector< pair< vector<float>, string > > aux;
    for(int i = 0; i < dm.size(); ++i)
      aux.push_back(pair< vector<float>, string >(dm[i], seqs[i]));

    while(dm.size() > 2){
      if(print_matrices){
	//informative prints
	cout << "CURRENT MATRIX SIZE: " << dm.size() << endl << endl;
	cout << "CURRENT DISTANCE MATRIX:" << endl;
	for(int i = 0; i < dm.size(); ++i){
	  for(int j = 0; j < dm.size(); ++j)
	    cout << dm[i][j] << " ";
	  cout << endl;
	}
	cout << endl;
      }

      //getting min avg distance
      int idx_i, idx_j;
      float max = INF;
      for(int i = 0; i < adm.size(); ++i)
	for(int j = 0; j < adm.size(); ++j)
	  if(adm[i][j] < max){
	    max = adm[i][j];
	    idx_i = i;
	    idx_j = j;
	  }

      if(aux[idx_i].second != "" && aux[idx_j].second != ""){
	GlobalAligner aligner(aux[idx_i].second, aux[idx_j].second, 1, 1);
	matches.push_back(aligner.get_matches()[0].first);
	matches.push_back(aligner.get_matches()[0].second);
      }
      else if(aux[idx_i].second != "")
	matches.push_back(aux[idx_i].second);
      else if(aux[idx_j].second != "")
	matches.push_back(aux[idx_j].second);
	else;
      
      vector<float> new_;

      for(int i = 0; i < dm.size(); ++i){
	new_.push_back((1.0f/2.0f)*
		      (dm[idx_i][i] + dm[idx_j][i] + dm[idx_i][idx_j]));
      }
      new_.push_back(INF);

      aux.erase(aux.begin()+idx_i);
      aux.erase(aux.begin()+idx_j-1);
      dm.erase(dm.begin()+idx_i);
      dm.erase(dm.begin()+idx_j-1);
      for(int i = 0; i < aux.size(); ++i){
	aux[i].first.erase(aux[i].first.begin()+idx_i);
	aux[i].first.erase(aux[i].first.begin()+idx_j-1);
	aux[i].first.push_back(new_[i]);
	dm[i].erase(dm[i].begin()+idx_i);
	dm[i].erase(dm[i].begin()+idx_j-1);
	dm[i].push_back(new_[i]);
      }
      dm.push_back(new_);
      aux.push_back(pair< vector<float>, string >(new_, ""));
      
      adm = get_avg_distance_matrix(dm);
    }
    
    matches.push_back((aux[0].second != "") ? aux[0].second : aux[1].second);

    for(int i = 0; i < matches.size(); ++i)
      fix_lengths();
  }

  void print_matches(){
    for(int i = 0; i < matches.size(); ++i)
      cout << matches[i] << endl;
  }
 

  void make_dot_file(){
    ifstream in = open("grafo.dot");
    in << "hola"  << endl;
    in.close();
  }
};
