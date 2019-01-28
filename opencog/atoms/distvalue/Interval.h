
#include <vector>


namespace opencog
{
namespace interval
{

typedef std::vector<double> DVec;

struct Interval
{
	double lower;
	double upper;
};

typedef Interval Bin;
//A N-dimensional Bin
typedef std::vector<Bin> NBin;
//A Sequence of N-dimensional Bins
typedef std::vector<NBin> NBinSeq;

#define EMPTY_INTERVAL Interval{1,-1}

bool is_empty(const Interval& i);

bool is_empty(const NBin& i);

double median(const Interval& i);

DVec median(const NBin& i);

double lower(const Interval& i);

DVec lower(const NBin& i);

double upper(const Interval& i);

DVec upper(const NBin& i);

double width(const Interval& i);

DVec width(const NBin& i);

//If one of the argumenst is empty the result should also be emtpy
Interval intersect(const Interval& i1,const Interval& i2);

NBin intersect(const NBin& b1,const NBin& b2);;

NBin createNBin(const DVec& v1,const DVec& v2);

bool operator<(const Interval& i1, const Interval& i2);
bool operator>(const Interval& i1, const Interval& i2);
bool operator==(const Interval& i1, const Interval& i2);

DVec operator+(const DVec& a, const DVec& b);
DVec operator-(const DVec& a, const DVec& b);
DVec operator/(const DVec& a, const DVec& b);
DVec operator*(const DVec& a, const double& b);
DVec operator/(const DVec& a, const double& b);

double conditional_probability(const NBin&, const NBin&);

std::string to_string(const DVec &v,const std::string& indent);
std::string to_string(const NBin &b,const std::string& indent);
} //namespace interval
} //namespace opencog
