#include <vector>
#include <algorithm>

using namespace std;

class Filter{
 public:
  Filter(){};

  int count_cuts(string s){
    int count = 0;
    int line = 0;
    for(int i = 0; i < s.length(); ++i){
      if(s[i] == '_'){
	if(i == 0){
	  line = 1;
	}
	else{
	  count += !line;
	  line = 1;
	  if(i == (s.length()-1))
	    count--;
	}
      }
      else
	line = 0;
    }
    
    return count;
  }

  void operator()(vector< pair<string, string> > &matches){
    vector<int> underscores;
    for(int i = 0; i < matches.size(); ++i)
      underscores.push_back(count_cuts(matches[i].first) +
			    count_cuts(matches[i].second));
    
    int min = underscores[0];
    for(int i = 1; i < underscores.size(); ++i)
      min = (underscores[i] < min) ? underscores[i] : min;

    vector< pair<string, string> > buff;
    for(int i = 0; i < underscores.size(); ++i)
      if(underscores[i] == min)
	buff.push_back(matches[i]);

    matches = buff;
  }
};
