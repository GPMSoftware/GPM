/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * DataSet.h
 *
 *  Created on: Sep 27, 2010
 *      Author: gareth
 */

#ifndef DataSet_H_
#define DataSet_H_

#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <math.h>

class Histogram
{
public:
    Histogram(std::vector<double> const &lower_limits = std::vector<double>(1) );

    void addPoint(double data);

    int operator[](int i) const;

    int points() const { return m_points; }

    int bins() const;

    double percentile(double p) const; // a fraction (not a percentage)

    void clear();

private:

    int m_points;                   // total number of points added;
    std::vector<double> m_ll;       // lower limits for each bin (last one is upper limit to last bin)
    std::vector<int> m_counts;      // counts in each bin
};


class DataSet;

// A data point (measured value of something) resulting from a number of observations.
// Reports mean, variance and other statistics.
struct DataPoint
{
    DataPoint(); // must not be used

    DataPoint(Histogram const &h, double y);

    void setPoint(double y); // resets to single observation
    void addPoint(double y); // adds an observation

    // statistics
    double sum()        const { return m_sum; }
    double sumsq()      const { return m_sum_sq; }
    int    num()        const { return m_num; }

    double mean()       const { return m_sum/m_num; }
    double var()        const { return m_sum_sq/m_num - mean() * mean(); }
    double stdDev()     const { return sqrt(var()); }

    double meanLog()    const { return m_sum_log/m_num; }
    double varLog()     const { return m_sum_sq_log/m_num - meanLog() * meanLog(); }
    double stdDevLog()  const { return sqrt(varLog()); }

    double gmean()      const { return pow(10.0, meanLog()); }
    double gVar()       const { return pow(10.0, varLog()); }    // ?? not sure this makes sense
    double gStdDev()    const { return pow(10.0, stdDevLog()); }

    double max()        const { return m_max; }
    double min()        const { return m_min; }

    double percentile(double p) const;

    Histogram m_hist; // keep a histogram for each point

    int    m_num;
    double m_sum;
    double m_sum_sq;
    double m_sum_log;
    double m_sum_sq_log;
    double m_max;
    double m_min;
};

// A series of DataPoints, each with an x value
class DataSeries
{
public:
    DataSeries();
    DataSeries(Histogram const &hist);
    DataSeries(Histogram const &hist, std::map<double, double> const &series);

    void setSeries(std::map<double, double> const &series);
    void addSeries(std::map<double, double> const &series);

    void setPoint(double x, double y);
    void addPoint(double x, double y);

    bool has(double x)      const { return m_series.find(x) != m_series.end(); }

    // statistics (will throw if value not found)
    double sum(double x)        const { return m_series.at(x).sum(); }
    double sumsq(double x)      const { return m_series.at(x).sumsq(); }
    int    num(double x)        const { return m_series.at(x).num(); }

    double mean(double x)       const { return m_series.at(x).mean(); }
    double var(double x)        const { return m_series.at(x).var(); }
    double stdDev(double x)     const { return m_series.at(x).stdDev(); }

    double meanLog(double x)    const { return m_series.at(x).meanLog(); }
    double varLog(double x)     const { return m_series.at(x).varLog(); }
    double stdDevLog(double x)  const { return m_series.at(x).stdDevLog(); }

    double gmean(double x)      const { return m_series.at(x).gmean(); }
    double gVar(double x)       const { return m_series.at(x).gVar(); }
    double gStdDev(double x)    const { return m_series.at(x).gStdDev(); }

    double max(double x)        const { return m_series.at(x).max(); }
    double min(double x)        const { return m_series.at(x).min(); }

    double percentile(double x, double p) const { return m_series.at(x).percentile(p); }

private:
    std::string m_name;

    std::string m_comment;

    std::map<double,    // x value
             DataPoint
            > m_series;

    Histogram m_hist;

    friend class DataSet;
};

// forward declarations are needed
template <typename T>
class CSVOutputter;

template <typename T>
std::ostream&
operator<<(std::ostream &os, CSVOutputter<T> const &op);

template <typename T>
class Fobj // function object
{
public:
    Fobj() {}
    virtual ~Fobj() {}
    virtual double operator()(DataPoint const &dp) = 0;
};


template <typename T>
class Fobj_void : public Fobj<T>
{
public:
    typedef T (DataPoint::*DataPointPTMFT) () const; // pointer to member of DataPoint returning T and taking (void)
    Fobj_void(DataPointPTMFT fp) : m_fp(fp) {}
    virtual double operator()(DataPoint const &dp) { return (dp.*m_fp)(); }
private:
    DataPointPTMFT m_fp; // the DataPoint member function to use
};

template <typename T>
class Fobj_double : public Fobj<T>
{
public:
    typedef T (DataPoint::*DataPointPTMFT_d) (double) const; // pointer to member of DataPoint returning T and taking (double)
    Fobj_double(DataPointPTMFT_d fp, double val) : m_fp(fp), m_val(val) {}
    virtual double operator()(DataPoint const &dp) { return (dp.*m_fp)(m_val); }
private:
    DataPointPTMFT_d m_fp; // the DataPoint member function to use
    const double m_val;    // argument
};

template <typename T>
class CSVOutputter
{
public:
    typedef T (DataPoint::*DataPointPTMFT) () const; // pointer to member of DataPoint returning T and taking (void)
    typedef T (DataPoint::*DataPointPTMFT_d) (double) const; // pointer to member of DataPoint returning T and taking (double)

    CSVOutputter(DataSet const &ds, DataPointPTMFT fp)
    : m_ds(ds), m_fobj(new Fobj_void<T>(fp)) {}

    CSVOutputter(DataSet const &ds, DataPointPTMFT_d fp, double p)
    : m_ds(ds), m_fobj(new Fobj_double<T>(fp, p)) {}

    ~CSVOutputter() { delete m_fobj; }

private:
    DataSet const &m_ds;   // the DataSet to output
    Fobj<T> *m_fobj; // the DataPoint member function to use

    friend std::ostream& operator<< <>(std::ostream &os, CSVOutputter const &op);
};

// Multiple data series each with a label.
// Provides output as a CSV file.
// The output contains a column for every x value occurring in any of the individual series.
class DataSet
{
public:
    DataSet(std::string const &xlabel,
            std::string const &ylabel,
            std::string name = "",
            std::string comment = "",
            Histogram const &hist = Histogram())
    : m_xlabel(xlabel)
    , m_ylabel(ylabel)
    , m_name(name)
    , m_comment(comment)
    , m_hist(hist)
    {}
//    virtual ~DataSet();

    void setXLabel(std::string const &xlabel)
    { m_xlabel = xlabel; }

    void setYLabel(std::string const &ylabel)
    { m_ylabel = ylabel; }

    void setName(std::string const &name)
    { m_name = name; }

    void setComment(std::string const &comment)
    { m_comment = comment; }

    void setSeries(std::string label, std::map<double, double> const &series);

    void addSeries(std::string label, std::map<double, double> const &series);

    void setPoint(std::string label, double x, double y);

    void addPoint(std::string label, double x, double y);

    bool has(std::string label)                const { return m_data.find(label) != m_data.end(); }
    bool has(std::string label, double x)      const { return has(label) && m_data.at(label).has(x); }

    // statistics (will throw if value not found)
    double     sum(std::string label, double x) const { return m_data.at(label).sum(x); }
    double   sumsq(std::string label, double x) const { return m_data.at(label).sumsq(x); }
    int        num(std::string label, double x) const { return m_data.at(label).num(x); }

    double    mean(std::string label, double x) const { return m_data.at(label).mean(x); }
    double     var(std::string label, double x) const { return m_data.at(label).var(x); }
    double  stdDev(std::string label, double x) const { return m_data.at(label).stdDev(x); }

    double meanLog(std::string label, double x) const   { return m_data.at(label).meanLog(x); }
    double varLog(std::string label, double x) const    { return m_data.at(label).varLog(x); }
    double stdDevLog(std::string label, double x) const { return m_data.at(label).stdDevLog(x); }

    double   gmean(std::string label, double x) const   { return m_data.at(label).gmean(x); }
    double   gVar(std::string label, double x) const    { return m_data.at(label).gVar(x); }
    double   gStdDev(std::string label, double x) const { return m_data.at(label).gStdDev(x); }

    double   max(std::string label, double x) const { return m_data.at(label).max(x); }
    double   min(std::string label, double x) const { return m_data.at(label).min(x); }

    double percentile(std::string label, double x, double p) const
    { return m_data.at(label).percentile(x, p); }

    // output the DataSet for a given DataPoint statistic
    template <class T>
    std::ostream& output(std::ostream &os, Fobj<T> *fp) const;

    // outputters for each of the DataPoint statistics
    CSVOutputter<double>     sumOut() { return CSVOutputter<double>(*this, &DataPoint::sum); }
    CSVOutputter<double>   sumsqOut() { return CSVOutputter<double>(*this, &DataPoint::sumsq); }
    CSVOutputter<int>        numOut() { return CSVOutputter<int>   (*this, &DataPoint::num); }

    CSVOutputter<double>    meanOut() { return CSVOutputter<double>(*this, &DataPoint::mean); }
    CSVOutputter<double>     varOut() { return CSVOutputter<double>(*this, &DataPoint::var); }
    CSVOutputter<double>  stdDevOut() { return CSVOutputter<double>(*this, &DataPoint::stdDev); }

    CSVOutputter<double> meanLogOut()   { return CSVOutputter<double>(*this, &DataPoint::meanLog); }
    CSVOutputter<double> varLogOut()    { return CSVOutputter<double>(*this, &DataPoint::varLog); }
    CSVOutputter<double> stdDevLogOut() { return CSVOutputter<double>(*this, &DataPoint::stdDevLog); }

    CSVOutputter<double>   gmeanOut()   { return CSVOutputter<double>(*this, &DataPoint::gmean); }
    CSVOutputter<double>   gVarOut()    { return CSVOutputter<double>(*this, &DataPoint::gVar); }
    CSVOutputter<double>   gStdDevOut() { return CSVOutputter<double>(*this, &DataPoint::gStdDev); }

    CSVOutputter<double>   maxOut() { return CSVOutputter<double>   (*this, &DataPoint::max); }
    CSVOutputter<double>   minOut() { return CSVOutputter<double>   (*this, &DataPoint::min); }

    CSVOutputter<double>   percentileOut(double p) { return CSVOutputter<double>(*this, &DataPoint::percentile, p); }

private:
    std::string m_xlabel;    // label for the x axis
    std::string m_ylabel;    // label for the y axis
    std::string m_name;
    std::string m_comment;

    std::map<std::string,    // Series label
             DataSeries      // Series
            > m_data;

    Histogram m_hist;
};

//#define Curve DataSeries;
//#define Graph DataSet;

#endif /* DataSet_H_ */
