/// <summary>
/// Definition of a class to contain an image with double precision pixels.
/// Includes pixel dimensions, resolutions and location arrays.
/// </summary>

#ifndef TEMPLATE_IMAGE_3D_
#define TEMPLATE_IMAGE_3D_

#include "zeusPositions.h"
#include "zeusConfig.h"

namespace Zeus {

/*! \brief A 3D image class template. */
template<class PixelType> class ZEUS_EXPORT Image3D {

public:

    /// <summary>
    /// Which modality are these voxels.
    /// </summary>
    Modality getModality( ) const { return _modality; }

    /// <summary>
    /// Query number of pixels in the x dimension.
    /// </summary>
    int Xsize( ) const { return _xPositions._nx; }

    /// <summary>
    /// Query number of pixels in the y dimension.
    /// </summary>
    int Ysize( ) const { return _yPositions._nx; }

    /// <summary>
    /// Query number of pixels in the z dimension.
    /// </summary>
    int Zsize( ) const { return _zPositions._nx; }

    /// <summary>
    /// Query total number of pixels in the volume.
    /// </summary>
    int NumberOfPixels( ) const { return _nelements; }

    /// <summary>
    /// Query voxel resolution in the x dimension.
    /// </summary>
    double Xresolution( ) const { return _xPositions._delta; }

    /// <summary>
    /// Query voxel resolution in the y dimension.
    /// </summary>
    double Yresolution( ) const { return _yPositions._delta; }

    /// <summary>
    /// Query voxel resolution in the y dimension.
    /// </summary>
    double Zresolution( ) const { return _zPositions._delta; }

    /// <summary>
    /// Returns the first x-location
    /// </summary>
    double FirstXlocation() const { return _xPositions._xmin; };

    /// <summary>
    /// Returns the last x-location
    /// </summary>
    double LastXlocation() const { return _xPositions._xmax; };

    /// <summary>
    /// Returns the first y-location
    /// </summary>
    double FirstYlocation() const { return _yPositions._xmin; };

    /// <summary>
    /// Returns the last y-location
    /// </summary>
    double LastYlocation() const { return _yPositions._xmax; };

    /// <summary>
    /// Returns the first y-location
    /// </summary>
    double FirstZlocation() const { return _zPositions._xmin; };

    /// <summary>
    /// Returns the last y-location
    /// </summary>
    double LastZlocation() const { return _zPositions._xmax; };

    /// <summary>
    /// Query pointer to the pixels
    /// </summary>
    PixelType *Buffer() { return _buffer; }

    /// <summary>
    /// Query pointer to the pixels
    /// </summary>
    const PixelType *Buffer() const { return _buffer; }

    /// <summary>
    /// Set all pixels to zero
    /// </summary>
    void Clear( PixelType value = 0 ) { for(int j=0; j<_nelements; ++j) _buffer[j] = value; };

    /// <summary>
    /// Constructor.  Allocates the volume.  
    /// </summary>
    /// <param name="aNx"> Number of pixels in x dimension. </param>
    /// <param name="aNy"> Number of pixels in y dimension.  </param>
    /// <param name="aNz"> Number of pixels in z dimension.  </param>
    /// <param name="aDx"> Resolution in x dimension. </param>
    /// <param name="aDy"> Resolution in y dimension. </param>
    /// <param name="aDz"> Resolution in z dimension. </param>
    /// <param name="aOx"> Location of first column. </param>
    /// <param name="aOy"> Location of first row. </param>
    /// <param name="aOz"> Location of first plane. </param>
    /// <param name="aModality"> Which modality are the voxels. </param> 
    Image3D( int aNx, int aNy, int aNz, double aDx, double aDy, double aDz, double aOx, double aOy, double aOz, Modality aModality ) : 
    _xPositions(aOx,aDx,aNx), _yPositions(aOy,aDy,aNy), _zPositions(aOz,aDz,aNz), 
        _nelements(_xPositions._nx*_yPositions._nx*_zPositions._nx), _buffer(_nelements > 0 ? new PixelType [_nelements] : 0), 
        _modality( aModality ) {};

    /// <summary>
    /// Constructor.  Allocates the volume.  
    /// </summary>
    /// <param name="xPositions"> The voxel grid in x direction. </param>
    /// <param name="yPositions"> The voxel grid in y direction. </param>
    /// <param name="zPositions"> The voxel grid in z direction. </param>
    /// <param name="aModality"> Which modality are the voxels. </param>
    Image3D(const Positions &xPositions, const Positions &yPositions, const Positions &zPositions, Modality aModality ) : 
    _xPositions(xPositions), _yPositions(yPositions), _zPositions(zPositions), 
        _nelements(_xPositions._nx*_yPositions._nx*_zPositions._nx), _buffer(_nelements > 0 ? new PixelType [_nelements] : 0), 
        _modality( aModality ) {};

    /// <summary>
    /// Destructor, delete the memory of the class.
    /// </summary>
    virtual ~Image3D( ) { if( _buffer ) delete [] _buffer; }

    /// <summary>
    /// Clone the class
    /// </summary>
    Image3D *Clone() const {
      Image3D* ptr = new Image3D(_xPositions, _yPositions, _zPositions, _modality);
        for(int j=0; j<_nelements; ++j) ptr->_buffer[j] = _buffer[j]; 
        return ptr;
    };

    /// <summary>
    /// Get the x-positions of the image
    /// </summary>
    const Positions& getXpositions() const { return _xPositions; };

    /// <summary>
    /// Get the y-positions of the image
    /// </summary>
    const Positions& getYpositions() const { return _yPositions; };

    /// <summary>
    /// Get the z-positions of the image
    /// </summary>
    const Positions& getZpositions() const { return _zPositions; };

    //enum ResamplingMethod {
    //    NearestNeighbour,
    //    Interpolation,
    //    Averaging
    //};
    //void resample(ResamplingMethod method, int nthread, Image3D<char> &outImage) { resampleT(method,nthread,outImage); };
    //void resample(ResamplingMethod method, int nthread, Image3D<unsigned char> &outImage) { resampleT(method,nthread,outImage); };
    //void resample(ResamplingMethod method, int nthread, Image3D<short> &outImage) { resampleT(method,nthread,outImage); };
    //void resample(ResamplingMethod method, int nthread, Image3D<unsigned short> &outImage) { resampleT(method,nthread,outImage); };
    //void resample(ResamplingMethod method, int nthread, Image3D<int> &outImage) { resampleT(method,nthread,outImage); };
    //void resample(ResamplingMethod method, int nthread, Image3D<unsigned int> &outImage) { resampleT(method,nthread,outImage); };
    //void resample(ResamplingMethod method, int nthread, Image3D<float> &outImage) { resampleT(method,nthread,outImage); };
    //void resample(ResamplingMethod method, int nthread, Image3D<double> &outImage) { resampleT(method,nthread,outImage); };

private:

    //template<class T> resampleT(ResamplingMethod method, int nthread, Image3D<T> &outImage);

    /// <summary>
    /// The x-positions of the pixel centers
    /// </summary>
    Positions  _xPositions;

    /// <summary>
    /// The y-positions of the pixel centers
    /// </summary>
    Positions  _yPositions; 

    /// <summary>
    /// The z-positions of the pixel centers
    /// </summary>
    Positions  _zPositions; 

    /// <summary>
    /// Total number of pixels in the volume.
    /// </summary>
    int _nelements;

    /// <summary>
    /// Pointer to the pixel data.
    /// </summary>
    PixelType *_buffer;

    /// <summary>
    /// Which modality are these voxels
    /// </summary>
    Modality _modality;

};

}

#endif
