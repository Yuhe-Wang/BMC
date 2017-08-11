#pragma once

namespace Zeus {
	/**
	* \brief Namespace containing all the constants used in ZeusMC.
	*/
	namespace Constants {
		/**
		* Electron mass energy equivalent [eV].
		* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?mec2mev|search_for=electron+mass">
		* Click here for the NIST reference, but note the units are MeV.</a>
		*/
		const double ELECTRON_MASS = 510998.910;

		/*! 1/ELECTRON_MASS in eV^-1*/
		const double INV_ELECTRON_MASS = 1.956951337e-6;

		/**
		* Classical electron radius [centimeters].
		* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?re|search_for=electron+radius">
		* Click here for the NIST reference, but note the units are in meters.</a>
		*/
		const double ELECTRON_RADIUS = 2.8179402894e-13;
		/**
		* Avogadro constant [1/mol].
		* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=Avogadro">
		* Click here for the NIST reference.</a>
		*/
		const double AVOGADRO = 6.02214179e23;
		/**
		* Speed of light in a vacuum [centimeters/sec].
		* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=speed+of+light">
		* Click here for the NIST reference, but note the units are [meters/sec].</a>
		*/
		const double SPEED_OF_LIGHT = 2.99792458e10;
		/**
		* Energy of 1 eV [J].
		* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?tevj|search_for=electron+volt">
		* Click here for the NIST reference.</a>
		*/
		const double ELECTRON_VOLT = 1.602176487e-19;
		/**
		* The mathematical constant PI.
		* N[Pi,30] from Mathematica 7.
		*/
		const double PI = 3.14159265358979323846264338328;
		/**
		* The square root of two.
		* N[Sqrt[2],30] from Mathematica 7.
		*/
		const double SQRT2 = 1.41421356237309;
		/**
		* This constant, the max value of a double precision floating point number,
		* is returned by method LambdaMoller::lmbdamo when the input energy is zero or less,
		* and is used to initialize a variable in method LambdaPhoton::SetupWoodcock, before
		* using that variable to hold a minimum value.
		*/
		const double DBLMAX = 1.79769313486231e308;
		/**
		* Zero divides are avoided in DoseWorker::inters by setting quotient to this value.
		*/
		const double FLTMAX = 3.402823466e38;
		/**
		* From microsoft, smallest number such that 1.0+FLTEPSILON != 1.0.
		*/
		const double FLTEPSILON = 1.192092896e-07;
		/**
		* A very small number, used in LambdaPhoton::SetupWoodcock.
		*/
		const double EPS = 1.0e-10;
		/**
		* The length of the table in LambdaPhoton::SetupWoodcock.
		*/
		const int WOODCOCKARRAYSIZE = 4096;
		/**
		* Number of leaf pairs in the MLC. Only the middle 30 are real.
		*/
		const int NLEAF_PAIRS_REAL = 30;
		/**
		* Number of leaf pairs, on each side, which are not real.
		*/
		const int N_LEAF_PAIRS_FAKE = 4;
		/**
		* Number of leaf pairs in the MLC. Only the middle 30 are real.
		*/
		const int NLEAF_PAIRS_TOTAL = 38;
		/**
		* Distance in CM from source to isocenter.
		*/
		const double SOURCE_TO_ISO_DISTANCE = 105;
		/**
		* Leaf thickness in CM.
		*/
		const double LEAF_THICKNESS = 1.05;
		/**
		* Leaf length in CM.
		*/
		const double LEAF_LENGTH = 100.0;
		/**
		* MLC model, the radius of the cylinder representing the inner leaves.
		*/
		const double MLC_INNER_RADIUS = 41;
		/**
		* MLC model, the radius of the cylinder representing the inner leaves.
		*/
		const double MLC_OUTER_RADIUS = 50;
		/**
		* The middle of the MLC, in the Y dimension.
		*/
		const double MLC_MIDDLE = (MLC_INNER_RADIUS + MLC_OUTER_RADIUS) / 2.0;
		/**
		* Minimum X location of MLC.
		*/
		const double MLC_MINX = -15;
		/**
		* Maximum X location of MLC.
		*/
		const double MLC_MAXX = 15;
		/**
		* MLC design constant
		*/
		const double INNER_RIGHT_FOCUS = 0.03;
		/**
		* MLC design constant
		*/
		const double INNER_GAP_FOCUS = -0.035;
		/**
		* MLC design constant
		*/
		const double INNER_Z_SHIFT = -0.015;
		/**
		* MLC design constant
		*/
		const double OUTER_LEFT_FOCUS = -0.03;
		/**
		* MLC design constant
		*/
		const double OUTER_GAP_FOCUS = 0.015;
		/**
		* MLC design constant
		*/
		const double OUTER_Z_SHIFT = 0.035;

	}
} // namespace

