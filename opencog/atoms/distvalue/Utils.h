#include <algorithm>

#include <opencog/util/exceptions.h>


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

bool is_empty(const Interval& i)
{
	return i.lower > i.upper;
}

double median(const Interval& i)
{
	return (i.lower + i.upper) / 2;
}

DVec median(const NBin& i)
{
	DVec res;
	transform(i.begin(),i.end(),std::back_inserter(res),
	          (double (*)(const Interval&))median);
	return res;
}

double lower(const Interval& i)
{
	return i.lower;
}

DVec lower(const NBin& i)
{
	DVec res;
	transform(i.begin(),i.end(),std::back_inserter(res),
	          (double (*)(const Interval&))lower);
	return res;
}

double upper(const Interval& i)
{
	return i.upper;
}

DVec upper(const NBin& i)
{
	DVec res;
	transform(i.begin(),i.end(),std::back_inserter(res),
	          (double (*)(const Interval&))upper);
	return res;
}

double width(const Interval& i)
{
	if (is_empty(i)) throw RuntimeException(TRACE_INFO,
											"Empty interval has no width.");
	return (i.upper - i.lower);
}

DVec width(const NBin& i)
{
	DVec res;
	transform(i.begin(),i.end(),std::back_inserter(res),
	          (double (*)(const Interval&))width);
	return res;
}

//If one of the argumenst is empty the result should also be emtpy
Interval intersect(const Interval& i1,const Interval& i2)
{
	double lower = std::max(i1.lower,i2.lower);
	double upper = std::min(i1.upper,i2.upper);
	return Interval{lower,upper};
}

NBin intersect(const NBin& b1,const NBin& b2)
{
	NBin res;
	transform(b1.begin(),b1.end(),b2.begin(),std::back_inserter(res),
	          (Interval (*)(const Interval&,const Interval&))intersect);
	return res;
}

NBin createNBin(const DVec& v1,const DVec& v2)
{
	size_t size = v1.size();
	if (size != v2.size())
		throw RuntimeException(TRACE_INFO,"Vectors must be same lenght.");
	NBin res;
	res.resever(size);
	for (size_t i; i < size;i++)
		res.push_back(Interval{v1[i],v2[i]});
	return res;
}


bool operator<(const Interval& i1, const Interval& i2)
{
	if (is_empty(i1)) return true;
	if (is_empty(i2)) return false;
	if (i1.upper <= i2.lower) return true;
	if (i2.upper <= i1.lower) return false;
	throw RuntimeException(TRACE_INFO,
	                       "Comparison not defined for overlapping intervals!");
}

bool operator>(const Interval& i1, const Interval& i2)
{
	if (is_empty(i1)) return false;
	if (is_empty(i2)) return true;
	if (i1.upper <= i2.lower) return false;
	if (i2.upper <= i1.lower) return true;
	throw RuntimeException(TRACE_INFO,
	                       "Comparison not defined for overlapping intervals!");
}

bool operator==(const Interval& i1, const Interval& i2)
{
	if (is_empty(i1) && is_empty(i2)) return true;
	//Check if either of them is empty can't be both if we reach this line.
	if (is_empty(i1) || is_empty(i2)) return false;
	if (i1.lower == i2.lower && i1.upper == i2.upper) return true;
	return false;
}

DVec operator+(const DVec& a, const DVec& b)
{
    DVec result;

    const std::size_t n = std::min(a.size(), b.size()) ;
    std::transform(std::begin(a), std::begin(a)+n, std::begin(b),
                   std::back_inserter(result), std::plus<double>{});
    return result;
}

DVec operator-(const DVec& a, const DVec& b)
{
    DVec result;

    const std::size_t n = std::min(a.size(), b.size()) ;
    std::transform(std::begin(a), std::begin(a)+n, std::begin(b),
                   std::back_inserter(result), std::minus<double>{});
    return result;
}

DVec operator/(const DVec& a, const DVec& b)
{
    DVec result;

    const std::size_t n = std::min(a.size(), b.size()) ;
    std::transform(std::begin(a), std::begin(a)+n, std::begin(b),
                   std::back_inserter(result), std::divides<double>{});
    return result;
}

DVec operator*(const DVec& a, const double& b)
{
    DVec result;

	for (double elem : a)
		result.push_back(elem * b);

    return result;
}

DVec operator/(const DVec& a, const double& b)
{
    DVec result;

	for (double elem : a)
		result.push_back(elem / b);

    return result;
}
