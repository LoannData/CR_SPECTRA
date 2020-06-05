#include <iostream> //For console output/input
#include <fstream>  //Allows to read/write files 
#include <math.h>   //Basic mathematic functions 
#include <vector>   //For dynamic arrays 
#include <string>   //Operations on strings 
#include <tuple>
#include <cmath>
#include <sstream>
using namespace std;


std::string parameters = "./parameters.dat";
std::string snx = search(parameters, "NX"); int nx = stoi(snx);
std::string sne = search(parameters, "NE"); int ne = stoi(sne);
std::string sni = search(parameters,"ni");  double ni = stod(sni);
std::string sX  = search(parameters,"X");   double Xi = stod(sX); 
double nt = ni/Xi;
std::string smn = search(parameters,"mn");  double m_neutral = stod(smn); 
std::string sT   = search(parameters,"T");  double T  = stod(sT);
std::string scenter = search(parameters, "center"); double x_center = stod(scenter);
std::string scenter_index = search(parameters, "center_index"); int x_center_index = stod(scenter_index);
std::string sBcenter = search(parameters, "B"); double Bcenter = stod(sBcenter);
