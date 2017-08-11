#pragma once

#include "zeusConfig.h"

namespace Zeus {

    class Vector;
    template<class PixelType> class Image3D;

    /** \brief A class defining the calculation geometry of the voxels. 
    *
    * The dose calculation cube is in _Xmin..._Xmin+_Dx*_Nx, _Ymin..._Ymin+_Dy*_Ny, _Zmin..._Zmin+_Dz*_Nz
    * The center of voxel (ix,iy,iz) is at (_Xmin+_Dx*(0.5+ix),_Ymin+_Dy*(iy+0.5),_Zmin+_Dz*(0.5+iz))
    * The mass density is in g/cm^3. The mass density matrix should include the couch and any body coils. 
    */
    class ZEUS_EXPORT VoxelGeometry {

    public:
        /**
        * \brief Constructor.
        *
        * Construct a voxel geometry class from an Image3D volume.
        * \param aVolume Take the extents and resolution from this image as the boundary of simulation space.
        */
        VoxelGeometry( Image3D<float> *aVolume );
        /**
        * Is the class valid.
        */
        bool isValid( ) { return _isValid; }
        /**
        * Destructor.  Deletes the buffers in this class.
        */
        ~VoxelGeometry() { delete [] _rho; delete [] _material; }

        float _Xmin; //!< Left corner of dose calculation cube. This is the center of the voxel, not its edge
        float _Ymin; //!< Bottom corner of dose calculation cube. This is the center of the voxel, not its edge
        float _Zmin; //!< Front corner of dose calculation cube. This is the center of the voxel, not its edge

        float _Xmax; //!< Upper most X coordinate of dose calculation cube. This is the center of the voxel, not its edge
        float _Ymax; //!< Upper most Y coordinate of dose calculation cube. This is the center of the voxel, not its edge
        float _Zmax; //!< Upper most Z coordinate of dose calculation cube. This is the center of the voxel, not its edge

        float _Dx;   //!< Voxel size in x-direction
        float _Dy;   //!< Voxel size in y-direction
        float _Dz;   //!< Voxel size in z-direction
        int _Nx;     //!< Number of voxels in x-direction
        int _Ny;     //!< Number of voxels in y-direction
        int _Nz;     //!< Number of voxels in z-direction
        float* _rho; //!< Mass density cube of size _Nx*_Ny*_Nz. Address of voxel (ix,iy,iz) is ix + iy*_Nx + iz*_Nx*_Ny
        char* _material; //!< Cube of size _Nx*_Ny*_Nz. Contains the material that this voxel is classified as.
        /**
        * Clone the class
        */
        VoxelGeometry *Clone( ) {
            return new VoxelGeometry( this );
        }

        /**
        * Does a particle with position x and direction u intersect the bounding box
        * formed by these voxels?
        * If it does (or if the position x is already inside the bounding box), return value is true and 
        * the position x is set to the entrance point (if outside initially) and the voxel indeces are returned in \a ind. 
        * Else return value is false.
        */
        bool enter( Vector &x, const Vector &u, int *ind ) const ;

        /**
        * Similar to the above function, except that now the position is set to be the relative position within the voxel on return.
        */
        bool enter1( Vector &x, const Vector &u, int *ind ) const ;

        /** 
         * Given the positions x, y, z, sets the voxel indeces xvox, yvox and zvox and returns the absolute voxel index, or -1 if the position is outside */
        inline int where(double x, double y, double z, int &xvox, int &yvox, int &zvox) const {
            int result = -1;
            if( x > _mins[0] && x < _maxs[0] && y > _mins[1] && y < _maxs[1] && z > _mins[2] && z < _maxs[2] ) {
                xvox = (int) ( (x - _mins[0])*_deltai[0] ); if( xvox >= _Nx ) xvox = _Nx-1;
                yvox = (int) ( (y - _mins[1])*_deltai[1] ); if( yvox >= _Ny ) yvox = _Ny-1;
                zvox = (int) ( (z - _mins[2])*_deltai[2] ); if( zvox >= _Nz ) zvox = _Nz-1;
                result = xvox + yvox*_Nx + zvox*_Nx*_Ny;
            }
            return result;
        };

        ZFloat _delta[3];  //!< Voxel sizes (same as _Dx, _Dy, _Dz)
        ZFloat _deltai[3]; //!< Inverse voxel sizes
        ZFloat _mins[3];   //!< Min. of the geometry bounding box. _mins[0] = _Xmin - _Dx/2, etc.
        ZFloat _maxs[3];   //!< Max. of the geometry bounding box. _maxs[i] = _mins[i] + _delta[i]*_nPlanes[i]

    private:
        /**
        * Constructor.
        * A helper method for the clone method.
        */
        VoxelGeometry( VoxelGeometry *g );

        bool     _isValid;   //!< Is this geometry valid?

        int _nPlanes[3];     //!< Number of voxels. Same as _Nx, _Ny, _Nz

    };

}
