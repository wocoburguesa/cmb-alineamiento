#include <vector>
#include <algorithm>

using namespace std;

class LocalAligner{
 private:
  string a;
  string b;
  vector< vector<int> > dp;
  vector< vector<int*> > backtrack;
  vector< pair<string,string> > matches;
  vector< pair<int,int> > max_positions;

 public:
  LocalAligner(string a_, string b_, bool auto_run=0){
    a = a_;
    b = b_;

    //setting the dp table
    for(int i = 0; i < b.length()+1; ++i){
      vector<int> row;
      for(int j = 0; j < a.length()+1; ++j)
	row.push_back(0);
      dp.push_back(row);
    }

    //setting the backtrack table
    for(int i = 0; i < b.length()+1; ++i){
      vector<int*> row;
      for(int j = 0; j < a.length()+1; ++j){
	int * backrefs = new int [3];
	backrefs[0] = backrefs[1] = backrefs[2] = 0;
	row.push_back(backrefs);
      }
      backtrack.push_back(row);
    }

    if(auto_run)
      run();
  }

  void rate_match(int i, int j){
    int max = dp[i-1][j-1] + ((b[i-1] == a[j-1]) ? 1 : -1);
    max = (dp[i-1][j] - 2 > max) ? dp[i-1][j]-2 : max;
    max = (dp[i][j-1] - 2 > max) ? dp[i][j-1]-2 : max;
    dp[i][j] = (max > 0) ? max : 0;
    
    if(max == dp[i-1][j-1] + ((b[i-1] == a[j-1]) ? 1 : -1))
      backtrack[i][j][0] = 1;
    else;
    if(max == dp[i-1][j] - 2)
      backtrack[i][j][1] = 1;
    else;
    if(max == dp[i][j-1] - 2)
      backtrack[i][j][2] = 1;
  }

  void traceback(int i, int j, string match_u, string match_d){
    if(dp[i][j] == 0){
      matches.push_back(pair<string, string>(match_u, match_d));
      return;
    }
    else{
      if(backtrack[i][j][0]){
	string new_u = match_u + ((i<=0) ? "_" : b.substr(i-1, 1));
	string new_d = match_d + ((j<=0) ? "_" : a.substr(j-1, 1));
	traceback(i-1, j-1, new_u, new_d);
      }
      else;
      if(backtrack[i][j][1]){
	string new_u = match_u + ((i<=0) ? "_" : b.substr(i-1, 1));
	traceback(i-1, j, new_u, match_d+"_");
      }
      else;
      if(backtrack[i][j][2]){
	string new_d = match_u + ((j<=0) ? "_" : a.substr(j-1, 1));
	traceback(i, j-1, match_u+"_", new_d);
      }
      else;
    }
  }

  int get_max(){
    int max = 0;
    for(int i = 0; i < dp.size(); ++i)
      for(int j = 0; j < dp[i].size(); ++j)
	if(dp[i][j] > max)
	  max = dp[i][j];
    return max;
  }

  void run(){
    //filling the table
    for(int i = 1; i < b.length()+1; ++i)
      for(int j = 1; j < a.length()+1; ++j)
	rate_match(i,j);

    //getting max positions
    int max = get_max();
    for(int i = 1; i < b.length()+1; ++i)
      for(int j = 1; j < a.length()+1; ++j)
	if(dp[i][j] == max)
	  max_positions.push_back(pair<int,int>(i,j));

    //backtracking matches
    for(int i = 0; i < max_positions.size(); ++i)    
      traceback(max_positions[i].first, max_positions[i].second, "", "");
    
    for(int i = 0; i < matches.size(); ++i){
      reverse(matches[i].first.begin(), matches[i].first.end());
      reverse(matches[i].second.begin(), matches[i].second.end());
    }
  }

  void print_matches(){
    for(int i = 0; i < matches.size(); ++i)
      cout << matches[i].first << endl << matches[i].second << endl << endl;
  }
};
