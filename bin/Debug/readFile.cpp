#include <fstream>
#include <iostream>
#include <string>
using namespace std;

int main()
  {
  string line;
  ifstream myfile( "paramDistribution.txt" );
  if (myfile)  // same as: if (myfile.good())
    {
    while (getline( myfile, line ))  // same as: while (getline( myfile, line ).good())
      {
	      cout<<line<<endl;
      }
    myfile.close();
    }
  else cout << "fooey\n";

  return 0;
  }
