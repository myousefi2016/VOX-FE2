
/*
  from http://www.johndcook.com/blog/standard_deviation/
  A simple class to compute running statistics

  Slightly modified by R Holbrey to track min/max too.
  Also added histogram.
*/

#include <math.h>
#include <limits>
using std::numeric_limits;

class RunningStat {
  public:
    RunningStat() {
      Clear();
    }

    void Clear()
    {
      m_n = 0;
      m_min = numeric_limits<double>::max();
      m_max = -( numeric_limits<double>::max() - 1 );
    }

    void Push(const double& x)
    {
      m_n++;

      // See Knuth TAOCP vol 2, 3rd edition, page 232
      if (m_n == 1)
      {
        m_oldM = x;
        m_newM = x;
        m_oldS = 0.0;
      }
      else
      {
        m_newM = m_oldM + (x - m_oldM)/(double(m_n));
        m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

        // set up for next iteration
        m_oldM = m_newM;
        m_oldS = m_newS;
      }

      if( x < m_min ) m_min = x;
      if( x > m_max ) m_max = x;
    }

    long int NumDataValues() const
    {
      return m_n;
    }

    double Mean() const
    {
      return (m_n > 0) ? m_newM : 0.0;
    }

    double Variance() const
    {
      return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
    }

    double StandardDeviation() const
    {
      return sqrt( Variance() );
    }

    void MinMax( double& min, double& max ) {
      min = m_min;
      max = m_max;
    }

    double Min() { return m_min; }
    double Max() { return m_max; }

  private:
    long int m_n;
    double m_oldM, m_newM, m_oldS, m_newS;

    double m_min, m_max;
};


class RunningHistogram {
  public:
    RunningHistogram(long nbins, double hmin, double hmax ):
      bin(0),
      hist_min(hmin),
      hist_max(hmax + numeric_limits<double>::epsilon()),  // we want to *just* include the max value
      hist_count(0), num_bins(0)  {

      if( nbins > 0 )
        num_bins = nbins;
      else return;

      hist_range = hist_max - hist_min;

      //initialize the bins
      bin = new (std::nothrow) long[ num_bins ];
      for(long i=0; i<num_bins; i++ ) bin[i] = 0;
      outside_bin[0] = 0;
      outside_bin[1] = 0;
    }

    ~RunningHistogram() {

      if( bin ) delete[] bin;
    }

    void AddValue( const double& val ) {

      if( val < hist_min ) //add to the outside min bin
        outside_bin[0]++;
      else if( val >= hist_max ) //add to the outside max bin
        outside_bin[1]++;
      else {

        long index = long( double(num_bins) * ((val - hist_min)/hist_range) );
        bin[index]++;
      }

      hist_count++;
    }


    double GetBinValueAtCumulativeDensity( const double& density, double& sum ) {

      sum = 0;
      long index = 0;

      for(long i=0; i<num_bins; i++ ) {

        sum += (double(bin[i])/double(hist_count));
        if( sum > density ) {
          index = i;
          break;
        }
      }

      //return the mid-point of the selected bin
      return hist_min + ( (index + 0.5) * (hist_range/num_bins) );

    }

    void Print( std::ostream& out ) {

      double bin_width = hist_range/num_bins;
      for(long i=0; i<num_bins; i++ ) {
        out << "Bin: " << i << " Mid-pt: "
            << (hist_min + (i + 0.5) * bin_width) << "  count: "
            << bin[i] << "\n";
      }

      out << "\n\n";
      out << "Left-bin (below hmin):  " << outside_bin[0] << "\n";
      out << "Right-bin (above hmax): " << outside_bin[1] << "\n";
      out << "Total-count: " << hist_count << "\n";
    }


  protected:
    RunningHistogram();  //not implemented

    long* bin;           //accumulate counts
    long outside_bin[2];

    double hist_min;
    double hist_max;
    double hist_range;
    long hist_count;
    long num_bins;



};
