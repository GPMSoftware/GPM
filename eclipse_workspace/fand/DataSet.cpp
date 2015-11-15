/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * DataSet.cpp
 *
 *  Created on: Sep 27, 2010
 *      Author: gareth
 */

#include "DataSet.h"
#include "fand/MessageStream.h"
#include "Assert.h"
#include <UnitTest++/UnitTest++.h>

using namespace std;

DataPoint::DataPoint()
{
    Assert(0); // default constructor must not be used
}

DataPoint::DataPoint(Histogram const &h, double y)
: m_hist(h), m_num(0), m_sum(0), m_sum_sq(0), m_sum_log(0), m_sum_sq_log(0), m_max(0), m_min(0)
{
    setPoint(y);
}

void DataPoint::setPoint(double y)
{
    m_num = 1; m_sum = y; m_sum_sq = y*y;
    double logy(log10(y)); m_sum_log = logy; m_sum_sq_log = logy*logy;
    m_max = m_min = y;

    m_hist.clear();

    if (m_hist.bins() > 0)
    {
        m_hist.addPoint(y);
    }
}

void DataPoint::addPoint(double y)
{
    if (m_num == 0)
    {
        setPoint(y);
    }
    else
    {
        m_num++; m_sum += y; m_sum_sq += y*y;
        double logy(log10(y)); m_sum_log += logy; m_sum_sq_log += logy*logy;
        m_max = std::max(m_max, y);
        m_min = std::min(m_min, y);
    }

    if (m_hist.bins() > 0)
    {
        m_hist.addPoint(y);
    }
}

double
DataPoint::percentile(double p) const
{
    Assert (m_hist.bins() > 0);

    return m_hist.percentile(p);
}

DataSeries::DataSeries()
{
    Assert(0); // default constructor must not be used
}

DataSeries::DataSeries(Histogram const &hist)
: m_hist(hist)
{
}

DataSeries::DataSeries(Histogram const &hist, std::map<double, double> const &series)
: m_hist(hist)
{
    setSeries(series);
}

void
DataSeries::setPoint(double x, double y)
{
    if (has(x))
    {
        m_series[x].setPoint(y);
    }
    else
    {
        m_series.insert(make_pair(x, DataPoint(m_hist, y)));
    }
}

void
DataSeries::addPoint(double x, double y)
{
    if (has(x))
    {
        m_series[x].addPoint(y);
    }
    else
    {
        setPoint(x, y);
    }
}

void
DataSeries::setSeries(std::map<double, double> const &series)
{
    m_series.clear();

    std::map<double, double>::const_iterator it;
    for (it = series.begin(); it != series.end(); ++it)
    {
        setPoint(it->first, it->second);
    }
}

void
DataSeries::addSeries(std::map<double, double> const &series)
{
    std::map<double, double>::const_iterator it;
    for (it = series.begin(); it != series.end(); ++it)
    {
        addPoint(it->first, it->second);
    }
}

template <typename T>
std::ostream&
operator<<(std::ostream &os, CSVOutputter<T> const &op)
{
    return op.m_ds.output(os, op.m_fobj);
}

void
DataSet::setSeries(std::string label, std::map<double, double> const &series)
{
    if (!has(label))
    {
        m_data.insert(make_pair(label, DataSeries(m_hist)));
    }
    m_data[label].setSeries(series);
}

void
DataSet::addSeries(std::string label, std::map<double, double> const &series)
{
    if (has(label))
    {
        m_data[label].addSeries(series);
    }
    else
    {
        setSeries(label, series);
    }
}

void
DataSet::setPoint(std::string label, double x, double y)
{
    if (!has(label))
    {
        m_data.insert(make_pair(label, DataSeries(m_hist)));
    }
    m_data[label].setPoint(x, y);
}

void
DataSet::addPoint(std::string label, double x, double y)
{
    if (!has(label))
    {
        m_data.insert(make_pair(label, DataSeries(m_hist)));
    }
    m_data[label].addPoint(x, y);
}

template <class T>
std::ostream&
DataSet::output(std::ostream &os, Fobj<T> *fp) const
{
    // get all the x values used in all datasets
    std::set<double> x_ords;

    std::map<std::string, DataSeries>::const_iterator it;
    for (it = m_data.begin(); it != m_data.end(); ++it)
    {
        DataSeries const &series = it->second;

        std::map<double, DataPoint>::const_iterator it2;
        for (it2 = series.m_series.begin(); it2 != series.m_series.end(); ++it2)
        {
            x_ords.insert(it2->first);
        }
    }

    // output name and comment
    if (!m_name.empty())
    {
        os << m_name;

        if (!m_comment.empty())
        {
            os << "," << m_comment;
        }
        os << "\n";
    }

    // output legend
    os << m_xlabel << " vs " << m_ylabel << "\n\n";

    // output x values
    os.unsetf(ios_base::floatfield);
    os.precision(3);

    os << m_xlabel;

    std::set<double>::const_iterator it3;
    for (it3 = x_ords.begin(); it3 != x_ords.end(); ++it3)
    {
        os << "," << *it3;
    }
    os << "\n";

    // output data series
    for (it = m_data.begin(); it != m_data.end(); ++it)
    {
        os << it->first; // Series label
        DataSeries const &series = it->second;

        std::map<double, DataPoint>::const_iterator it2;
        for (it3 = x_ords.begin(); it3 != x_ords.end(); ++it3)
        {
            os << ",";
            if ((it2 = series.m_series.find(*it3)) != series.m_series.end())
            {
                Assert(it2->second.m_num != 0);
                os << (*fp)(it2->second);
            }
        }
        os << "\n";
    }

    os << endl;

    return os;
}

char csv1[] = {
        "DataSet1,a comment\n"
        "XVALS vs YVALS\n\n"
        "XVALS,0,1,2,3,4\n"
        "EVENS,0.1,,2.1,,4.1\n"
        "ODDS,,1.1,,3.1,\n\n"
        };

char csv2[] = {
        "DataSet1,a comment\n"
        "XVALS vs YVALS\n\n"
        "XVALS,0,1,2,3,4\n"
        "EVENS,0.1,,2.1,,4.1\n"
        "ODDS,0.1,1.1,2.1,3.1,4.1\n\n"
        };

TEST(DataSet1)
{
    DataSet csv("XVALS", "YVALS", "DataSet1", "a comment");

    map<double, double> EVENS;
    EVENS[0] = 0.1;
    EVENS[2] = 2.1;
    EVENS[4] = 4.1;
    csv.addSeries("EVENS", EVENS);

    map<double, double> ODDS;
    ODDS[1] = 1.1;
    ODDS[3] = 3.1;
    csv.setSeries("ODDS", ODDS);

    ostringstream oss;
    oss << csv.meanOut();

    CHECK_EQUAL(std::string(csv1), oss.str());

    csv.addSeries("EVENS", EVENS);  // EVENS unchanged (mean of two copies of itself)
    csv.addSeries("ODDS", EVENS);   // ODDS now has points from both EVENS and ODDS

    oss.str(""); // clear contents
    oss << csv.meanOut();

    CHECK_EQUAL(std::string(csv2), oss.str());

    CHECK_EQUAL(1, csv.num("ODDS", 1));
    CHECK_EQUAL(1.1, csv.mean("ODDS", 1));
    csv.addPoint("ODDS", 1, 0.9);          // now mean of (1.1, 0.9) = 1;
    CHECK_EQUAL(2, csv.num("ODDS", 1));
    CHECK_EQUAL(1, csv.mean("ODDS", 1));
    CHECK_CLOSE(0.01, csv.var("ODDS", 1), 1e-6); // ( (1.1-1)^2 + (0.9-1)^2 ) / 2

    CHECK_EQUAL(true, csv.has("ODDS"));
    CHECK_EQUAL(false, csv.has("SODS"));
    CHECK_EQUAL(true, csv.has("ODDS", 1));
    CHECK_EQUAL(true, csv.has("ODDS", 2));
    CHECK_EQUAL(false, csv.has("EVENS", 1));

//    csv.output(cout, &DataPoint::var);
}

Histogram::Histogram(std::vector<double> const &lower_limits)
: m_points(0)
, m_ll(lower_limits)
, m_counts(lower_limits.size()-1)
{
}

int
Histogram::bins() const
{
    Assert(m_ll.size() > 0);
    return m_ll.size() - 1;
}

void Histogram::addPoint(double data)
{
    Assert(data >= m_ll.front() && data < m_ll.back());

    int d = m_counts.size();

    int i = 0;
    while( i < (int)m_ll.size() && m_ll[i] < data)
    {
        ++i;
    }

    Assert(i > 0 && i <= bins());

    m_counts[i-1]++;
    m_points++;
}

int Histogram::operator[](int i) const
{
    Assert(i>=0 && i < bins());
    return m_counts[i];
}

void
Histogram::clear()
{
    m_counts = std::vector<int>(m_counts.size(), 0);
    m_points = 0;
}

double
Histogram::percentile(double p) const // a fraction (not a percentage)
{
    double count_p = p * m_points;

    int c = 0; // cumulative count
    double tmp = count_p;
    size_t lim_u = 0, lim_l = 0; // lower and upper limits between which p occurs (indices)
    while( (tmp>0) && lim_u<m_counts.size() )
    {
        c   += m_counts[lim_u];
        tmp -= m_counts[lim_u++];
    }

    if (lim_u == 0) return m_ll.front();

    lim_l = lim_u - 1;
    c -= m_counts[lim_l];

//    cout << "lim_l           = " << lim_l           << " lim_u           = " << lim_u << endl;
//    cout << "m_ll[lim_l]     = " << m_ll[lim_l]     << " m_ll[lim_u]     = " << m_ll[lim_u] << endl;
//    cout << "count_p = " << count_p << " c = " << c << " m_counts[lim_u] = " << m_counts[lim_u] << endl;

    // interpolate
    return m_ll[lim_l] + (m_ll[lim_u] - m_ll[lim_l]) * (count_p - c) / m_counts[lim_l];
}

TEST(Histogram)
{
    vector<double> lims;
    lims.push_back(0);
    lims.push_back(1e1);
    lims.push_back(1e2);
    lims.push_back(1e3);
    lims.push_back(1e4);
    lims.push_back(1e5);

    Histogram hist(lims);
    CHECK_EQUAL(5, hist.bins());
    CHECK_EQUAL(0,  hist.points());

    hist.addPoint(1);  // 1
    hist.addPoint(11); hist.addPoint(12); hist.addPoint(13); hist.addPoint(14); // 4
    hist.addPoint(200); hist.addPoint(300); hist.addPoint(400); hist.addPoint(500); hist.addPoint(600);  // 5
    hist.addPoint(1001); hist.addPoint(3000); hist.addPoint(4000);  // 3
    hist.addPoint(10001); hist.addPoint(10001);  // 2

    CHECK_EQUAL(15,  hist.points());

    CHECK_EQUAL(1, hist[0]);
    CHECK_EQUAL(4, hist[1]);
    CHECK_EQUAL(5, hist[2]);
    CHECK_EQUAL(3, hist[3]);
    CHECK_EQUAL(2, hist[4]);

    CHECK_EQUAL(0,      hist.percentile(0));
    CHECK_EQUAL(1.5,    hist.percentile(0.01));
    CHECK_EQUAL(71.875, hist.percentile(0.25));
    CHECK_EQUAL(4750,   hist.percentile(0.75));
    CHECK_EQUAL(100000, hist.percentile(1));

}
