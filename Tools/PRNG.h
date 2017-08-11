
/*
class DECLSPECIFIER PRNG
{
public:
	void init(int NR, int BufferSize = 512);
	inline double operator () ();
	void getRandomAzimuth(double &cosphi, double &sinphi);
private:
	int ISEED1;
	int ISEED2;
	int cur;
	vector<double> buffer;

	void fillBuffer();
};
*/

class  DECLSPECIFIER PRNG {
public:

	PRNG(){};
	void init(unsigned long aSeed, int aBufferSize = 624);
	PRNG(unsigned long aSeed, int aBufferSize = 624){ init(aSeed, aBufferSize); }
	void getRandomAzimuth(double &cosphi, double &sinphi);
	/**
	* \brief Destructor.
	*
	* Empty the buffer of random numbers.
	*/
	~PRNG();
	/**
	* Fill the buffer with random numbers.
	*/
	void fillBuffer();
	inline double  operator () () {
		if (_index == _bufferSize) fillBuffer();
		return _randomBuffer[_index++];
	};
private:

	/* Period parameters */
	int _Np;
	int _M;

	unsigned long _MATRIX_A; /* constant vector a */
	unsigned long _UPPER_MASK; /* most significant w-r bits */
	unsigned long _LOWER_MASK; /* least significant r bits */

	/* initializes mt with a seed */
	void init_genrand(unsigned long s);

	/* generates a random number on [0,1)-real-interval */
	inline double genrand_real2();

	unsigned long *mt;
	int mti;

	/* 1/(2^32-1) to convert 32 bit integer into float in [0,1.0] */
	double _oneOver2ToPower32Minus1;

	vector<double>    _randomBuffer;
	/**
	* The next random number in the buffer to return to the caller.
	*/
	int                 _index;
	/**
	* The size of the buffer of random numbers.
	*/
	int                 _bufferSize;
};