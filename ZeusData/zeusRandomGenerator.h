#pragma once

#include "zeusConfig.h"

#include <vector>
using std::vector;

namespace Zeus { 
/**
* \brief Namespace for all related to random number generation.
*
* We are still using MT19997, just have packaged it differently, now each thread
* has its own random number generator, as to current product DPM, where a single random
* number generator class kept the state for all threads.
*/
namespace Random {

/** \brief Interface for random number generators.
*
* Simulation code will get a handle to a generator thru this interface.
* Use method getGenerator to get the handle, and then use method getUniform 
* to get the random numbers
* <PRE>
*   RandomGenerator rng = RandomGenerator::getGenerator( seed, 256 );
*   double vz = 2.0 * aRng->getUniform() - 1.0; 
* </PRE>
*/
class ZEUS_EXPORT RandomGenerator {
public:
	/**
	* \brief Constructor.  
	*
	* \param bufferSize The size of the buffer used to store random numbers. Instead of calling the 
       *                   random number generating function each time a random number is needed, 
       *                   the implementation generates at once \a bufferSize random numbers kept in the buffer 
       *                   and then returns the numbers in thee buffer until it is empty.
	*/
    RandomGenerator(int bufferSize = 256);
	/**
	* Destructor.
	*/
    virtual ~RandomGenerator() {};
	/**
	* \brief Get a uniform random number.
	*
	* If there are still random numbers left in the buffer, you get one.  If all the random
	* numbers in the buffer have been consumed, the buffer is filled again, and you
	* get the first one.
	*/
    inline ZFloat getUniform() {
        if( _index == _bufferSize ) fillBuffer();
        return _randomBuffer[_index++];
    };
	/**
	* \brief Returns a poiter to a random number generator object.
	*
	* Call thru this interface to get a handle to a random number generator.
	* \param seed The seed for the random number generator.
	* \param bufferSize Size of the buffer (see constructor above)
       * \param type The generator type. 0->SimpleMT19937AR, 1->Ranmar
       * It is the responsibility of the user to delete the generator when no longer needed.
	*/
    static RandomGenerator* getGenerator(int seed, int bufferSize = 256, int type = 0);

    /**
     * Get the cosine and the sine of a randomly generated azimuthal angle 
     */
    void getRandomAzimuth(double &cosphi, double &sinphi);

protected:
	/** 
	* \brief Generate a buffer full of random numbers.
	*
	* If you provide your own random number generator, you must implement this method to 
	* generate a buffer of random numbers.  
	*/
    virtual void fillBuffer() = 0;
	/**
	* This array stores all the random numbers.
	* The type of the random numbers is described here \sa ZFloat
	*/
    vector<ZFloat>    _randomBuffer;
	/**
	* The next random number in the buffer to return to the caller.
	*/
    int                 _index;
	/**
	* The size of the buffer of random numbers.
	*/
    int                 _bufferSize;

};

}}
