/****************************************************************************
* $Id$
*
* *************************************************************************/

#include "zeusRandomGenerator.h"
#include "zeusLogProviders.h"

#include <cmath>
#include <vector>
using namespace std;

namespace Zeus { namespace Random {

/**
* \brief MT199937 random number generator.
 
A C-program for MT19937, with initialization improved 2002/1/26.
Coded by Takuji Nishimura and Makoto Matsumoto.

Before using, initialize the state by using init_genrand(seed)  
or init_by_array(init_key, key_length).

Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
All rights reserved.                          
Copyright (C) 2005, Mutsuo Saito,
All rights reserved.                          

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. The names of its contributors may not be used to endorse or promote 
products derived from this software without specific prior written 
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Any feedback is very welcome.
http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/
class ZEUS_EXPORT SimpleMT19937AR : public RandomGenerator {
public:
	/**
	* \brief Constructor
	* 
	* When you call getGenerator thru the interface, this class is constructed.
	* \param aSeed The seed of the random number generator.
	* \param aBufferSize The size of the buffer.
	*/
    SimpleMT19937AR(unsigned long aSeed, int aBufferSize );
	/**
	* \brief Destructor.
	*
	* Empty the buffer of random numbers.
	*/
    ~SimpleMT19937AR();
	/**
	* Fill the buffer with random numbers.
	*/
    void fillBuffer( );
#ifdef _WINDOWS
    /**
    * To stop Microsoft Warning C4512
    */
    SimpleMT19937AR & operator=( const SimpleMT19937AR & );
#endif
private:

    /* Period parameters */  
    const int _Np;
    const int _M;

    const unsigned long _MATRIX_A; /* constant vector a */
    const unsigned long _UPPER_MASK; /* most significant w-r bits */
    const unsigned long _LOWER_MASK; /* least significant r bits */

    /* initializes mt with a seed */
    void init_genrand(unsigned long s);

    /* generates a random number on [0,1)-real-interval */
    inline double genrand_real2();

    unsigned long *mt;
    int mti;

    /* 1/(2^32-1) to convert 32 bit integer into float in [0,1.0] */
    const double _oneOver2ToPower32Minus1;

};

SimpleMT19937AR::SimpleMT19937AR(unsigned long aSeed, int aBufferSize ) : RandomGenerator( aBufferSize ), _Np(624), _M(397), 
    _MATRIX_A(0x9908b0dfUL), _UPPER_MASK(0x80000000UL), _LOWER_MASK(0x7fffffffUL), mti(0), _oneOver2ToPower32Minus1(1.0/4294967295.0) {
        mt = new unsigned long[_Np];
        init_genrand( aSeed );
}

SimpleMT19937AR::~SimpleMT19937AR() { 
    _randomBuffer.clear(); 
    delete [] mt;
}

void SimpleMT19937AR::fillBuffer( ) {
    for ( int i = 0 ; i < _bufferSize ; ++i ) _randomBuffer[i] = genrand_real2( );
    _index = 0;
}

/* initializes mt with a seed */
void SimpleMT19937AR::init_genrand(unsigned long s) {
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<_Np; mti++) {
        mt[ mti ] = ( 1812433253UL * ( mt[mti-1] ^ (mt[mti-1] >> 30) ) + mti ); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[ mti ] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
    mti = 0;
}

/* generates a random number on [0,0xffffffff]-interval */
double SimpleMT19937AR::genrand_real2() {
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, _MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if ( mti >= _Np) { /* generate N words at one time */
        int kk;

        for (kk=0;kk<_Np-_M;kk++) {
            y = (mt[kk]&_UPPER_MASK)|(mt[kk+1]&_LOWER_MASK);
            mt[kk] = mt[kk+_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<_Np-1;kk++) {
            y = (mt[kk]&_UPPER_MASK)|(mt[kk+1]&_LOWER_MASK);
            mt[kk] = mt[kk+(_M-_Np)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[_Np-1]&_UPPER_MASK)|(mt[0]&_LOWER_MASK);
        mt[_Np-1] = mt[_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[ mti++ ];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y * _oneOver2ToPower32Minus1;

}
#ifdef _WINDOWS
SimpleMT19937AR & SimpleMT19937AR::operator=( const SimpleMT19937AR & ) { return *this; }
#endif

/********************************************** end of SimpleMT19937AR ********************************************************/

/*! The RANMAR generator by Marsaglia and Zaman. This generator is the default random number generator in the 
    EGSnrc package.
*/
class Ranmar : public RandomGenerator {

public:

    Ranmar(int sequence, int bufferSize) : RandomGenerator(bufferSize) {
        initialize(sequence);
    };

    void fillBuffer();

private:

    void initialize(int sequence);

    vector<int>     _buffer;
    int             _i1, _i2;
    int             _c;

    static int      _cd;
    static int      _cm;
    static int      _imax;
    static double   _intToFloat;

};

int Ranmar::_cd   = 7654321;
int Ranmar::_cm   = 16777213;
int Ranmar::_imax = 16777216;
double Ranmar::_intToFloat = 1./16777216.;

void Ranmar::initialize(int sequence) {
    int ixx = 1802;
    int jxx = 9373 + sequence;
    if( jxx <= 0 ) jxx += 30080;
    else if( jxx >= 30081 ) jxx -= 30080;
    if( jxx <= 0 || jxx >= 30081 ) {
        Zeus::LogProviders::FatalLog("Ranmar::initialize: illegal sequence index %d\n",sequence);
    }

    int i = (ixx/177)%177 + 2;
    int j = ixx%177 + 2;
    int k = (jxx/169)%178 + 1;
    int l = jxx%169;

    _buffer.resize(97);
    for(int ii=0; ii<97; ++ii) {
        int s = 0; int t = 8388608;
        for(int jj=0; jj<24; ++jj) {
            int m = (((i*j)%179)*k)%179;
            i = j; j = k; k = m;
            l = (53*l+1)%169;
            if( (l*m)%64 >= 32 ) s += t;
            t /= 2;
        }
        _buffer[ii] = s;
    }
    _c = 362436;
    _i1 = 96; _i2 = 32;
}

void Ranmar::fillBuffer() {
    for(int i=0; i<_bufferSize; ++i) {
        int r = _buffer[_i1] - _buffer[_i2--];
#ifdef MSVC
        if( r > 16777219 ) Zeus::LogProviders::InformationLog("r = %d\n",r);
#endif
        if( r < 0 ) r += _imax;
        _buffer[_i1--] = r;
        if(_i1 < 0) _i1 = 96;
        if(_i2 < 0) _i2 = 96;
        _c -= _cd; if( _c < 0 ) _c += _cm;
        r -= _c; if( r < 0 ) r += _imax;
#ifdef MSVC
        if( r > 16777219 ) Zeus::LogProviders::InformationLog("r = %d\n",r);
#endif
        ++r;
        _randomBuffer[i] = _intToFloat*r;
    }
    _index = 0;
}

/********************************************* end of Ranmar *******************************************************************/


RandomGenerator::RandomGenerator( int aBufferSize ) {
    if( aBufferSize < 1 ) Zeus::LogProviders::FatalLog("RandomGenerator::RandomGenerator: attempt to construct with a zero buffer size");
    _bufferSize = aBufferSize;
    _randomBuffer.resize(aBufferSize);
    _index = aBufferSize;
}

RandomGenerator* RandomGenerator::getGenerator(int seed, int bufferSize, int type) { 
    RandomGenerator *result;
    if( type == 0 ) result = new SimpleMT19937AR(seed,bufferSize);
    else result = new Ranmar(seed,bufferSize);
    return result;
}

class RandomAzimuthHelper {

public:

    RandomAzimuthHelper(int nbin);

    inline void compute(double phi, double &cphi, double &sphi) const {
        int bin = (int) (phi*_invBin);
        cphi = _table[bin]._ac + _table[bin]._bc*phi;
        sphi = _table[bin]._as + _table[bin]._bs*phi;
    };

private:

    struct Data {
        double _ac, _bc, _as, _bs;
        Data() {};
    };

    vector<Data> _table;
    double       _invBin;

};

#ifndef M_PI
#include "zeusConstants.h"
#define M_PI Zeus::Constants::PI
#endif

RandomAzimuthHelper::RandomAzimuthHelper(int nbin) {
    _table.resize(nbin);
    double dphi = 1./(nbin-1); _invBin = 1/dphi;
    double cold = 1, sold = 0;
    for(int i=1; i<nbin; ++i) {
        double phi = dphi*i;
        double c = cos(2*M_PI*phi), s = sin(2*M_PI*phi);
        _table[i-1]._bc = (c - cold)*_invBin;
        _table[i-1]._ac = c - _table[i-1]._bc*phi;
        _table[i-1]._bs = (s - sold)*_invBin;
        _table[i-1]._as = s - _table[i-1]._bs*phi;
        cold = c; sold = s;
    }
    _table[nbin-1] = _table[nbin-2];
}


void RandomGenerator::getRandomAzimuth(double &cosphi, double &sinphi) {
//#ifdef HAVE_SINCOS
//    double phi = 2*M_PI*getUniform();
//    sincos(phi,&sinphi,&cosphi);
//#else
//    double x,x2,y,y2,r;
//    do {
//        x = 2*getUniform() - 1; x2 = x*x;
//        y = getUniform(); y2 = y*y;
//        r = x2 + y2;
//    } while ( r > 1 );
//    double ri = 1/r;
//    cosphi = (x2 - y2)*ri;
//    sinphi = 2*x*y*ri;
//#endif
    static RandomAzimuthHelper helper(2048);
    helper.compute(getUniform(),cosphi,sinphi);
}

} // namespace Random
} // namespace Zeus
