#include <algorithm>
#include <numeric>
#include <cmath>
#include<sstream>

#include <opencog/util/exceptions.h>

#include <opencog/atoms/distvalue/Interval.h>

namespace opencog
{
namespace interval
{

bool is_empty(const Interval& i)
{
	return i.lower > i.upper;
}

bool is_empty(const NBin& b)
{
	return std::any_of(b.begin(),b.end(),(bool (*)(const Interval&))is_empty);
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
	if (b1.size() != b2.size())
		throw RuntimeException(TRACE_INFO,"Same Number of Dimensions Required.");
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
	res.reserve(size);
	for (size_t i = 0; i < size;i++)
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

DVec all_lower(const NBin& b1, const NBin& b2)
{

	if (b1.size() != b2.size())
		throw RuntimeException(TRACE_INFO,"Vectors must be same lenght.");
	auto it1 = b1.begin();
	auto it2 = b2.begin();
	while (it1 != b1.end())
	{
	}

}

//Find out how much of Key1 is contained in Key2
//The Keys can be considered continuous uniform distributions
//This calculates the conditional_probabilty of Key2 given Key1
double conditional_probability(const NBin &ks1,const NBin &ks2)
{
	NBin intersection = intersect(ks1,ks2);
	if (is_empty(intersection))
		return 0;
	DVec overlaps = width(intersection) / width(ks1);
	//If ks1 is a point-interval the above will result in nan
	//Because we know the intersection is not empty it must be completely inside
	//so replace all nans with 1
	for (size_t i = 0; i < overlaps.size(); i++)
	{
		if (std::isnan(overlaps[i]))
			overlaps[i] = 1;
	}
	return std::accumulate(overlaps.begin(),overlaps.end()
	                       ,1.0,std::multiplies<double>());
}

std::string to_string(const DVec &v,const std::string& indent)
{
	std::stringstream ss;
	ss << indent << "[";
	for (double val : v)
	{
		ss << val
		   << ",";
	}
	ss.seekp(-1,std::ios_base::end);
	ss << "]";
	return ss.str();
}
std::string to_string(const NBin &b,const std::string& indent)
{
	std::stringstream ss;
	ss << indent << "{";
	for (auto bin : b)
	{
		ss << "["
		   << bin.lower
		   << ","
		   << bin.upper
		   << ");";
	}
	ss.seekp(-1,std::ios_base::end);
	ss << "}";
	return ss.str();
}
} //namespace interval
} //namespace opencog
