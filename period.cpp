#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cfloat>
#include <vector>
#include <algorithm>
#include <iomanip>
using namespace std;

string read_file;
string write_file;

struct atom
{
    int atom_type;
    double x,y,z;

    atom (int t,double x_,double y_,double z_)
        :atom_type(t),x(x_),y(y_),z(z_) {}

    bool operator<(const atom& other) const
    {
        return atom_type<other.atom_type;
    }
};

void read_data(vector<atom>& R,int *atom_number,double *a);
void boundary(vector<atom>& R,int *atom_number,double *a);
void write_data(vector<atom>& R,int *atom_number,double *a);

int main(int argc,char *argv[])
{
    read_file=argv[1];
    write_file=argv[2];
    double a[3];
    int atom_number[5]={0,0,0,0,0};
    vector<atom> R;

    read_data(R,atom_number,a);
    cout << "read_data finished!!" << endl;
    boundary(R,atom_number,a);
    cout << "boundary finished!!" << endl;
    write_data(R,atom_number,a);
    cout << "wirte_data finished!" << endl;
}




void read_data(vector<atom>& R,int *atom_number,double *a)
{
    double x_b[6];
    string line,word;
    int lines;
    ifstream file;
    file.open(read_file);

    if (!file.is_open())
    {
        cout << "Can't open dump.atom file." << endl;
        exit(1); 
    }

    for (int i=0;i<2;i++)
    {
        getline(file,line);
    }

    {
        getline(file,line);
        istringstream WORDS(line);
        WORDS >> word;
        atom_number[0]=stoi(word);
    }

    for (int i=0;i<2;i++)
    {
        getline(file,line);
    }

    for (int i=0;i<3;i++)
    {
        getline(file,line);
        istringstream WORDS(line);
        for (int j=0;j<2;j++)
        {
            if (WORDS >> word)
            {
                x_b[i*2+j]=stod(word);
            }
        }
        a[i]=x_b[i*2+1]-x_b[i*2];
    }
    for (int i=0;i<10;i++)
    {
        getline(file,line);
    }

    for (int i=0;i<atom_number[0];i++)
    {
        getline(file,line);
        istringstream WORDS(line);
        int idx,type;
        double x,y,z;
        for (int j=0;j<5;j++)
        {
            if (WORDS >> word)
            {
                if (j==0) idx=stoi(word);
                else if (j==1) type=stoi(word)-1;
                else if (j==2) x=stod(word)-x_b[0];
                else if (j==3) y=stod(word)-x_b[2];
                else if (j==4) z=stod(word)-x_b[4];
            }
        }
        atom_number[type+1]+=1;
        R.push_back(atom(type,x,y,z));
    }
    sort(R.begin(),R.end());
    file.close();
}

void boundary(vector<atom>& R,int *atom_number,double *a)
{
    const double EPSILON=1e-9;
    const int number=atom_number[0];
    for (int i=0;i<number;i++)
    {
        double x,y,z;
        int atom_type;
        if (R[i].z<EPSILON)
        {
            atom_type=R[i].atom_type;
            x=R[i].x;
            y=R[i].y;
            z=a[2];
            atom_number[0]+=1;
            atom_number[atom_type+1]+=1;
            R.push_back(atom(atom_type,x,y,z));
        }
        else if (fabs(R[i].z-a[2])<EPSILON)
        {
            atom_type=R[i].atom_type;
            x=R[i].x;
            y=R[i].y;
            z=0.0;
            atom_number[0]+=1;
            atom_number[atom_type+1]+=1;
            R.push_back(atom(atom_type,x,y,z));
        }
    }
    sort(R.begin(),R.end());
}

void write_data(vector<atom>& R,int *atom_number,double *a)
{
    ofstream file;
    file.open(write_file);
    //file.precision(9);
    //file.setf(ios::fixed,ios::floatfield);

    file << "Lammps data file written by whitecrn" << endl;
    file << endl;
    file << atom_number[0] << "  atoms" << endl;
    file << "4 atom types" << endl;
    file << endl;
    file << 0.0 << "  " << a[0] << "  " << "xlo " << " xhi" << endl;
    file << 0.0 << "  " << a[1] << "  " << "ylo " << " yhi" << endl;
    file << 0.0 << "  " << a[2] << "  " << "zlo " << " zhi" << endl;
    file << endl;
    file << "Masses" << endl;
    file << endl;
    file << "1 6.941 # Li" << endl;
    file << "2 138.9055 # La" << endl;
    file << "3 91.224 # Zr" << endl;
    file << "4 15.9994  # O" << endl;
    file << endl;
    file << "Atoms # atomic" << endl;
    file << endl;

    for (int i=0;i<atom_number[0];i++)
    {
        file << i+1 << " " << (R[i].atom_type)+1 << " "
            << fixed << setprecision(9)
            << R[i].x << " " << R[i].y << " " << R[i].z << endl;
    }
    file.close();
}