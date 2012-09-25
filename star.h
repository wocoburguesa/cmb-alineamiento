#include <vector>
#include <algorithm>
#include "aligner.h"

using namespace std;

class StarAligner{
 private:
  vector<string> seqs;
  vector<string> matches;

 public:
  StarAligner(vector<string> input){
    seqs = input;
  }

  int get_center_idx(){
    int scores[seqs.size()][seqs.size()];
    for(int i = 0; i < seqs.size(); ++i){
      for(int j = i; j < seqs.size(); ++j){
	if(i != j){
	  GlobalAligner aligner(seqs[i], seqs[j], 1, 1);
	  scores[i][j] = aligner.get_score();	
	  scores[j][i] = aligner.get_score();  
	}
	else
	  scores[i][j] = 0;
      }
    }
    
    int center_idx = 0;
    int max_scores = 0;
    for(int i = 0; i < seqs.size(); ++i)
      max_scores += scores[0][i];

    for(int i = 1; i < seqs.size(); ++i){
      int buff_scores = 0;
      for(int j = 0; j < seqs.size(); ++j){
	buff_scores += scores[i][j];
      }
      if(buff_scores > max_scores){
	max_scores = buff_scores;
	center_idx = i;
      }
    }
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

  void run(){
    int center_idx = get_center_idx();
    vector< pair<string, string> > buff_matches;
    for(int i = 0; i < seqs.size(); ++i){
      if(i != center_idx){
	GlobalAligner aligner(seqs[center_idx], seqs[i], 1, 1);
	buff_matches.push_back(aligner.get_matches()[0]);
      }
    }

    matches.push_back(seqs[center_idx]);

    for(int i = 0; i < buff_matches.size(); ++i){
      matches.push_back(buff_matches[i].first);
      fix_lengths();
    }
  }

  void print_matches(){
    for(int i = 0; i < matches.size(); ++i)
      cout << matches[i] << endl;
  }

};
