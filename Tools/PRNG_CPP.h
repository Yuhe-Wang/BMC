
class RandomAzimuthHelper {

public:

	RandomAzimuthHelper(int nbin);

	inline void compute(double phi, double &cphi, double &sphi) const {
		int bin = (int)(phi*_invBin);
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


RandomAzimuthHelper::RandomAzimuthHelper(int nbin) {
	_table.resize(nbin);
	double dphi = 1. / (nbin - 1); _invBin = 1 / dphi;
	double cold = 1, sold = 0;
	for (int i = 1; i < nbin; ++i) {
		double phi = dphi*i;
		double c = cos(2 * PI*phi), s = sin(2 * PI*phi);
		_table[i - 1]._bc = (c - cold)*_invBin;
		_table[i - 1]._ac = c - _table[i - 1]._bc*phi;
		_table[i - 1]._bs = (s - sold)*_invBin;
		_table[i - 1]._as = s - _table[i - 1]._bs*phi;
		cold = c; sold = s;
	}
	_table[nbin - 1] = _table[nbin - 2];
}

/*
void PRNG::init(int NR, int BufferSize) //init the generator by its index in an array or threads
{
	
	if (sizeof(long) != 4) exitApp("not 32bit platform, random number generator cannot work!");
	long In = 123456789 + NR;

	In = 16807 * (In % 127773) - 2836 * (In / 127773);
	if (In < 0) In += 2147483647;
	ISEED1 = In;
	In = 16807 * (In % 127773) - 2836 * (In / 127773);
	if (In < 0) In += 2147483647;
	ISEED2 = In;

	buffer.resize(BufferSize);
	fillBuffer();
}

inline double PRNG::operator () () //default return random number in (0,1)
{
	if (cur == PRNGBufferSize)
	{
		fillBuffer();
	}
	return buffer[cur++];
}

void PRNG::fillBuffer()
{
	const double USCALE = 1.0 / 2.147483563e9;
	long I1, I2, IZ;
	int BufferSize = buffer.size();
	for (int i = 0; i < BufferSize; ++i)
	{
		I1 = ISEED1 / 53668;
		ISEED1 = 40014 * (ISEED1 - I1 * 53668) - I1 * 12211;
		if (ISEED1 < 0) ISEED1 += 2147483563;

		I2 = ISEED2 / 52774;
		ISEED2 = 40692 * (ISEED2 - I2 * 52774) - I2 * 3791;
		if (ISEED2 < 0) ISEED2 += 2147483399;

		IZ = ISEED1 - ISEED2;
		if (IZ < 1) IZ += 2147483562;
		buffer[i] = IZ*USCALE;
	}
	cur = 0;
}

void PRNG::getRandomAzimuth(double &cosphi, double &sinphi) 
{
	static RandomAzimuthHelper helper(2048);
	helper.compute((*this)(), cosphi, sinphi);
}

*/

void PRNG::init(unsigned long aSeed, int aBufferSize) {
	_Np = 624;
	_M = 397;
	_MATRIX_A = 0x9908b0dfUL;
	_UPPER_MASK = 0x80000000UL;
	_LOWER_MASK = 0x7fffffffUL; 
	mti = 0;
	_oneOver2ToPower32Minus1 = 1.0 / 4294967295.0;
	_bufferSize = aBufferSize;
	_randomBuffer.resize(aBufferSize);
	_index = aBufferSize;

	mt = new unsigned long[_Np];
	init_genrand(aSeed);
}

PRNG::~PRNG() {
	_randomBuffer.clear();
	delete[] mt;
}

void PRNG::fillBuffer() {
	for (int i = 0; i < _bufferSize; ++i) _randomBuffer[i] = genrand_real2();
	_index = 0;
}

/* initializes mt with a seed */
void PRNG::init_genrand(unsigned long s) {
	mt[0] = s & 0xffffffffUL;
	for (mti = 1; mti < _Np; mti++) {
		mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
	mti = 0;
}

/* generates a random number on [0,0xffffffff]-interval */
double PRNG::genrand_real2() {
	unsigned long y;
	static unsigned long mag01[2] = { 0x0UL, _MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= _Np) { /* generate N words at one time */
		int kk;

		for (kk = 0; kk < _Np - _M; kk++) {
			y = (mt[kk] & _UPPER_MASK) | (mt[kk + 1] & _LOWER_MASK);
			mt[kk] = mt[kk + _M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (; kk < _Np - 1; kk++) {
			y = (mt[kk] & _UPPER_MASK) | (mt[kk + 1] & _LOWER_MASK);
			mt[kk] = mt[kk + (_M - _Np)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[_Np - 1] & _UPPER_MASK) | (mt[0] & _LOWER_MASK);
		mt[_Np - 1] = mt[_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y * _oneOver2ToPower32Minus1;

}

static RandomAzimuthHelper helper(2048);

void PRNG::getRandomAzimuth(double &cosphi, double &sinphi) {
	
	helper.compute((*this)(), cosphi, sinphi);
}