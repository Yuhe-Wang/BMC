#pragma once

#include <vector>
#include <cmath>
using namespace std;

#include "zeusConfig.h"

namespace Zeus {    

  /*! A simple structure to hold data for linear interpolations */
  struct LIData {
    double  p1, p2;    //!< Interpolations coefficients for left and right grid points
    int     i1, i2;    //!< Indeces of left and right grid points
  };

/*! \brief A simple structure to hold information for the voxel/pixel location, size and number along one dimension */
struct ZEUS_EXPORT Positions {

    double  _xmin;        //!< The left (first) voxel position. Note: this is the voxel center, so left image edge is _xmin - _delta/2
    double  _delta;       //!< The voxel size
    double  _deltai;      //!< 1/_delta
    double  _xmax;        //!< The right (last) voxel center, _xmax = _xmin + _delta*(_nx-1). Right image edge is _xmax + _delta/2
    int     _nx;          //!< Number of voxels.

    /*! A dummy ctor */
    Positions() : _xmin(0), _delta(0), _deltai(1), _xmax(0), _nx(0) {};

    /*! Ctor from a minimum position \a xmin, voxel size \a delta and number of voxels \a nx */
    Positions(double xmin, double delta, int nx) : _xmin(xmin), _delta(delta), _deltai(1/delta), _xmax(_xmin + delta*(nx-1)), _nx(nx) {};

    /*! Comparison operator */
    bool operator==(const Positions &p) const {
        return (_nx == p._nx && fabs(_xmin - p._xmin) < 1e-4 && fabs(_delta - p._delta) < 1e-5);
    };

    /*! This functions stores in \a lidata the linear interpolation coefficients of the 
     * voxel centers of the input positions \a p. 
     * It's main purpose is to be used in linear interpolation resampling (a.k.a. point resampling) of images.
     */
    void setLinearInterpolation(const Positions &p, vector<LIData> &lidata, bool extend = true) const {
      int N = p._nx; lidata.resize(N); double deltai = 1/_delta;
      int ind = extend ? 0 : -1;
      double leftBoundary = _xmin - 0.5*_delta, rightBoundary = _xmax + 0.5*_delta;
      for(int i=0; i<N; ++i) {
	lidata[i].p1 = lidata[i].p2 = 0;
	lidata[i].i1 = lidata[i].i2 = ind;
	double x = p._xmin + p._delta*i;
	if( x >= _xmin && x < _xmax ) {
	  double xv = (x - _xmin)*deltai; int j = (int) xv; xv -= j;
	  lidata[i].p1 = 1 - xv; lidata[i].p2 = xv; 
	  lidata[i].i1 = j; lidata[i].i2 = j+1 <= _nx-1 ? j+1 : _nx-1;
	}
	else if( x < _xmin && x >= leftBoundary ) lidata[i].p1 = 1; 
	else if( x >= _xmax && x <= rightBoundary ) {
	  lidata[i].p1 = 1; lidata[i].i1 = _nx-1; lidata[i].i2 = _nx-1;
	}
      }
    };

    /*! This functions stores in \a indeces the indeces of the nearest neighbours of the voxel centers of the input positions \a p.
        It's main purpose is to be used in nearest neighbour resampling of images.
     */
    void setNearestNeighbours(const Positions &p, vector<int> &indeces) const {
        int N = p._nx; indeces.resize(N); double deltai = 1/_delta;
        double leftBoundary = _xmin - 0.5*_delta, rightBoundary = _xmax + 0.5*_delta;
        for(int i=0; i<N; ++i) {
            double x = p._xmin + p._delta*i;
            if( x >= _xmin && x < _xmax ) {
                double xv = (x - _xmin)*deltai; int j = (int) xv;
                if( xv - j > 0.5 ) { ++j; if( j > _nx-1 ) j = _nx-1; }
                indeces[i] = j;
            }
            else if( x < _xmin ) indeces[i] = x >= leftBoundary ? 0 : -1;
            else                 indeces[i] = x <= rightBoundary ? _nx-1 : -1;
        }
    };

 /*! This functions stores in \a adata the indeces and fractions of the voxels covered by the voxels of the input positions \a p.
        It's main purpose is to be used in averaging resampling of images.
     */
    void setAveragingInterpolation(const Positions &p, vector<vector<pair<int,double> > > &adata) const {
        int N = p._nx; adata.resize(N); double deltai = 1/_delta;
        double leftBoundary = _xmin - 0.5*_delta, rightBoundary = _xmax + 0.5*_delta;
        double scale = _delta/p._delta;
        for(int i=0; i<N; ++i) {
            double x1 = p._xmin + p._delta*(i - 0.5), x2 = x1 + p._delta;
            if( x1 >= rightBoundary || x2 <= leftBoundary ) continue;
            int i1 = 0;
            if( x1 >= leftBoundary ) {
                x1 = (x1 - leftBoundary)*deltai;
                int ix1 = (int) x1; x1 -= ix1; adata[i].push_back(pair<int,double>(ix1,scale*(1-x1)));
                i1 = ix1+1;
            }
            int i2 = _nx;
            if( x2 < rightBoundary ) {
                x2 = (x2 - leftBoundary)*deltai;
                int ix2 = (int) x2;
                if( ix2 == _nx ) ix2 = _nx-1; // avoid roundoff errors
                x2 -= ix2; i2 = ix2;
            }
            if( i2 == i1-1 ) adata[i][0].second = scale*(x2 - x1);
            else {
                for(int ii=i1; ii<i2; ++ii) adata[i].push_back(pair<int,double>(ii,scale));
                if( i2 < _nx ) adata[i].push_back(pair<int,double>(i2,scale*x2));
            }
        }
    };


};
  

}
