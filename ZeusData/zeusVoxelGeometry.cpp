
#include "zeusVoxelGeometry.h"
#include "zeusVector.h"
#include "zeusTemplateImage3D.h"

namespace Zeus {

    VoxelGeometry::VoxelGeometry( Image3D<float> *aVolume ) : _isValid(false) {

        // Geometry is set here.
        _Nx = aVolume->Xsize(); _Ny = aVolume->Ysize(); _Nz = aVolume->Zsize();
        int nvoxels = _Nx * _Ny * _Nz;

        _Xmin = (float)aVolume->FirstXlocation(); _Ymin = (float)aVolume->FirstYlocation(); _Zmin = (float)aVolume->FirstZlocation();

        _Dx = float(aVolume->Xresolution()); _Dy = float(aVolume->Yresolution()); _Dz = float(aVolume->Zresolution());

        _delta[0] = _Dx; _delta[1] = _Dy; _delta[2] = _Dz;
        _deltai[0] = 1.0/_Dx; _deltai[1] = 1.0/_Dy; _deltai[2] = 1.0/_Dz;
        _nPlanes[0] = _Nx; _nPlanes[1] = _Ny; _nPlanes[2] = _Nz;
        _mins[0] = _Xmin - 0.5*_delta[0]; _mins[1] = _Ymin - 0.5*_delta[1]; _mins[2] = _Zmin - 0.5*_delta[2];

        _Xmax = _Xmin + (_Nx-1)*_Dx; _Ymax = _Ymin + (_Ny-1)*_Dy; _Zmax = _Zmin + (_Nz-1)*_Dz;
        _maxs[0] = _Xmax + 0.5*_delta[0]; _maxs[1] = _Ymax + 0.5*_delta[1]; _maxs[2] = _Zmax + 0.5*_delta[2];

        // Electron density
        _rho = new float[nvoxels];
        for ( int i = 0 ; i < nvoxels ; ++i ) _rho[i] = aVolume->Buffer()[i];

        // Create a class of material, where very material is of type 1.
        _material = new char[aVolume->NumberOfPixels()];
        for ( int i = 0 ; i < aVolume->NumberOfPixels() ; ++i ) _material[i] = 1;

        _isValid = true;

    }

    VoxelGeometry::VoxelGeometry( VoxelGeometry *g ) : _isValid(false) {

        // Geometry is set here.
        _Nx = g->_Nx; _Ny = g->_Ny; _Nz = g->_Nz;
        int nvoxels = _Nx * _Ny * _Nz;
        _Dx = g->_Dx; _Dy = g->_Dy; _Dz = g->_Dz;
        _Xmin = g->_Xmin; _Ymin = g->_Ymin; _Zmin = g->_Zmin;
        _delta[0] = _Dx; _delta[1] = _Dy; _delta[2] = _Dz;
        _deltai[0] = 1.0/_Dx; _deltai[1] = 1.0/_Dy; _deltai[2] = 1.0/_Dz;
        _nPlanes[0] = _Nx; _nPlanes[1] = _Ny; _nPlanes[2] = _Nz;
        _mins[0] = _Xmin - 0.5*_delta[0]; _mins[1] = _Ymin - 0.5*_delta[1]; _mins[2] = _Zmin - 0.5*_delta[2];

        _Xmax = _Xmin + (_Nx-1)*_Dx; _Ymax = _Ymin + (_Ny-1)*_Dy; _Zmax = _Zmin + (_Nz-1)*_Dz;
        _maxs[0] = _Xmax + 0.5*_delta[0]; _maxs[1] = _Ymax + 0.5*_delta[1]; _maxs[2] = _Zmax + 0.5*_delta[2];

        // Electron density
        _rho = new float[nvoxels];
        for ( int i = 0 ; i < nvoxels ; ++i ) _rho[i] = g->_rho[i];

        // Create a class of material, where very material is of type 1.
        _material = new char[nvoxels];
        for ( int i = 0 ; i < nvoxels ; ++i ) _material[i] = 1;

        _isValid = true;

    }


    bool VoxelGeometry::enter( Vector &x, const Vector &u, int *ind ) const {

        if( x._x < _mins[0] || x._x >= _maxs[0] || x._y < _mins[1] || x._y >= _maxs[1] || x._z < _mins[2] || x._z >= _maxs[2] ) {
            if( (x._z < _mins[2] && u._z > 0) || (x._z >= _maxs[2] && u._z < 0) ) {
                ZFloat t = u._z > 0 ? (_mins[2] - x._z)/u._z : (_maxs[2] - x._z)/u._z;
                ZFloat x1 = x._x + u._x*t, y1 = x._y + u._y*t;
                if( x1 >= _mins[0] && x1 < _maxs[0] && y1 >= _mins[1] && y1 < _maxs[1] ) { 
                    x += u*t; 
                    ind[2] = u._z > 0 ? 0 : _Nz - 1;
                    ind[0] = (int) ((x1 - _mins[0])*_deltai[0]);
                    ind[1] = (int) ((y1 - _mins[1])*_deltai[1]);
                    return true;
                }
            }
            if( (x._y < _mins[1] && u._y > 0) || (x._y >= _maxs[1] && u._y < 0) ) {
                ZFloat t = u._y > 0 ? (_mins[1] - x._y)/u._y : (_maxs[1] - x._y)/u._y;
                ZFloat x1 = x._x + u._x*t, z1 = x._z + u._z*t;
                if( x1 >= _mins[0] && x1 < _maxs[0] && z1 >= _mins[2] && z1 < _maxs[2] ) { 
                    x += u*t; 
                    ind[1] = u._y > 0 ? 0 : _Ny - 1;
                    ind[0] = (int) ((x1 - _mins[0])*_deltai[0]);
                    ind[2] = (int) ((z1 - _mins[2])*_deltai[2]);
                    return true;
                }

            }
            if( (x._x < _mins[0] && u._x > 0) || (x._x >= _maxs[0] && u._x < 0) ) {
                ZFloat t = u._x > 0 ? (_mins[0] - x._x)/u._x : (_maxs[0] - x._x)/u._x;
                ZFloat y1 = x._y + u._y*t, z1 = x._z + u._z*t;
                if( y1 >= _mins[1] && y1 < _maxs[1] && z1 >= _mins[2] && z1 < _maxs[2] ) { 
                    x += u*t; 
                    ind[0] = u._x > 0 ? 0 : _Nx - 1;
                    ind[1] = (int) ((y1 - _mins[1])*_deltai[1]);
                    ind[2] = (int) ((z1 - _mins[2])*_deltai[2]);
                    return true;
                }
            }
            return false;
        }
        ind[0] = (int) ((x._x - _mins[0])*_deltai[0]);
        ind[1] = (int) ((x._y - _mins[1])*_deltai[1]);
        ind[2] = (int) ((x._z - _mins[2])*_deltai[2]);
        return true;
    }

    bool VoxelGeometry::enter1( Vector &x, const Vector &u, int *ind ) const {

        if( x._x < _mins[0] || x._x >= _maxs[0] || x._y < _mins[1] || x._y >= _maxs[1] || x._z < _mins[2] || x._z >= _maxs[2] ) {
            if( (x._z < _mins[2] && u._z > 0) || (x._z >= _maxs[2] && u._z < 0) ) {
                ZFloat t = u._z > 0 ? (_mins[2] - x._z)/u._z : (_maxs[2] - x._z)/u._z;
                ZFloat x1 = x._x + u._x*t, y1 = x._y + u._y*t;
                if( x1 >= _mins[0] && x1 < _maxs[0] && y1 >= _mins[1] && y1 < _maxs[1] ) {
                    if( u._z > 0 ) {
                        ind[2] = 0; x._z = 0;
                    }
                    else {
                        ind[2] = _Nz - 1; x._z = _Dz;
                    }
                    x1 -= _mins[0]; x1 *= _deltai[0]; ind[0] = (int) x1; x._x = (x1 - ind[0])*_delta[0];
                    y1 -= _mins[1]; y1 *= _deltai[1]; ind[1] = (int) y1; x._y = (y1 - ind[1])*_delta[1];
                    return true;
                }
            }
            if( (x._y < _mins[1] && u._y > 0) || (x._y >= _maxs[1] && u._y < 0) ) {
                ZFloat t = u._y > 0 ? (_mins[1] - x._y)/u._y : (_maxs[1] - x._y)/u._y;
                ZFloat x1 = x._x + u._x*t, z1 = x._z + u._z*t;
                if( x1 >= _mins[0] && x1 < _maxs[0] && z1 >= _mins[2] && z1 < _maxs[2] ) {
                    if( u._y > 0 ) {
                        ind[1] = 0; x._y = 0;
                    }
                    else {
                        ind[1] = _Ny - 1; x._y = _Dy;
                    }
                    x1 -= _mins[0]; x1 *= _deltai[0]; ind[0] = (int) x1; x._x = (x1 - ind[0])*_delta[0];
                    z1 -= _mins[2]; z1 *= _deltai[2]; ind[2] = (int) z1; x._z = (z1 - ind[2])*_delta[2];
                    return true;
                }

            }
            if( (x._x < _mins[0] && u._x > 0) || (x._x >= _maxs[0] && u._x < 0) ) {
                ZFloat t = u._x > 0 ? (_mins[0] - x._x)/u._x : (_maxs[0] - x._x)/u._x;
                ZFloat y1 = x._y + u._y*t, z1 = x._z + u._z*t;
                if( y1 >= _mins[1] && y1 < _maxs[1] && z1 >= _mins[2] && z1 < _maxs[2] ) {
                    if( u._x > 0 ) {
                        ind[0] = 0; x._x = 0;
                    }
                    else {
                        ind[0] = _Nx - 1; x._x = _Dx;
                    }
                    y1 -= _mins[1]; y1 *= _deltai[1]; ind[1] = (int) y1; x._y = (y1 - ind[1])*_delta[1];
                    z1 -= _mins[2]; z1 *= _deltai[2]; ind[2] = (int) z1; x._z = (z1 - ind[2])*_delta[2];
                    return true;
                }
            }
            return false;
        }
        x._x -= _mins[0]; x._x *= _deltai[0]; ind[0] = (int) x._x; x._x -= ind[0]; x._x *= _delta[0];
        x._y -= _mins[1]; x._y *= _deltai[1]; ind[1] = (int) x._y; x._y -= ind[1]; x._y *= _delta[1];
        x._z -= _mins[2]; x._z *= _deltai[2]; ind[2] = (int) x._z; x._z -= ind[2]; x._z *= _delta[2];
        return true;
    }



}
