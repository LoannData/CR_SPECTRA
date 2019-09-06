#include <string>
#include <vector>
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream


void showLog_0(double time, double Tmax, int nstep, int time_id, double duration)
{
    if (time_id % nstep == 0)
    {
        cout<<"###########################################################"<<endl;
        cout<<"# LOG                                                     #"<<endl;
        cout<<"###########################################################"<<endl;
        cout<<"Physical time = "<<time/kyr<<" kyr over "<<Tmax/kyr<<" kyr"<<endl;
        cout<<"Time step number : "<<time_id<<". Time step duration : "<<duration<<" seconds"<<endl;
    }
}