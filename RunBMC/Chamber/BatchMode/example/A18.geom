
###############################################################################
#
#  EGSnrc egs++ sample car geometry
#  Copyright (C) 2015 National Research Council Canada
#
#  This file is part of EGSnrc.
#
#  EGSnrc is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Affero General Public License as published by the
#  Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  EGSnrc is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
#  more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with EGSnrc. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
#
#  Author:          Iwan Kawrakow, 2005
#
#  Contributors:
#
###############################################################################
#
#   An example geometry input file for the egs++ geometry package.
#
#   This input file defines a geometry modeling the A18 icon chamber
#
##############################################################################


:start geometry definition:

    ######################################## the electrode body
    :start geometry:
        library = egs_cones
        type    = EGS_ConeStack
        name    = electrode_body
        axis    = 0 0 0  0 0 1
        :start layer:
            thickness    = 0.3
            top radii    = 0.02 0.12 0.175 0.245
            bottom radii = 0.02 0.12 0.175 0.245
            media        = 29-Copper 228-Teflon 126-C552 225-Delrin
        :stop layer:
		:start layer:
            thickness    = 0.27
            top radii    = 0.02 0.12 0.175
            bottom radii = 0.02 0.12 0.175
            media        = 29-Copper 228-Teflon 126-C552
        :stop layer:
        :start layer:
            thickness    = 0.01
            top radii    = 0.02
            bottom radii = 0.02
            media        = 29-Copper
        :stop layer:
        :start layer:
            thickness    = 0.16
            top radii    = 0.175
            bottom radii = 0.175
            media        = 126-C552
        :stop layer:
		:start layer:
            thickness    = 0.44
            top radii    = 0.05
            bottom radii = 0.05
            media        = 126-C552
        :stop layer:
    :stop geometry:

    ######################################## the electrode sphere
    :start geometry:
        library = egs_spheres
		name    = electrode_sphere
		midpoint = 0 0 1.18
		radii = 0.05
		:start media input:
		    media = 126-C552
		:stop media input:
    :stop geometry:
	
    ######################################## the electrode
    :start geometry:
        library = egs_gunion
        name    = electrode
        geometries = electrode_body electrode_sphere
    :stop geometry:
	
	
	#*******************************************************#
	
	####################################### the 3 plans to construt a half capsule
    :start geometry:
        library   = egs_planes
        type      = EGS_Zplanes
        name      = base_planes
        positions = 0 1.23 1.8
    :stop geometry:

    ####################################### the spheres to construct the half capsule
    :start geometry:
        library  = egs_spheres
        name     = cap_spheres
        midpoint = 0 0 1.23
        radii    = 0.245 0.545
		:start media input:
		    media = 104-Air 126-C552
			set medium = 0 0
			set medium = 1 1
		:stop media input:
    :stop geometry:

    ######################################## the cylinders to construct the half capsule
    :start geometry:
        library  = egs_cylinders
        type     = EGS_ZCylinders
        name     = cap_cylinders
        radii    = 0.245  0.545
		midpoint = 0 0 0
		:start media input:
		    media = 104-Air 126-C552
			set medium = 0 0
			set medium = 1 1
		:stop media input:
    :stop geometry:
	
	######################################## build the half capsule
	:start geometry:
        library    = egs_cdgeometry
        name       = capsule
        base geometry = base_planes
        set geometry = 0 cap_cylinders
		set geometry = 1 cap_spheres
    :stop geometry:
	
	#*******************************************************#
	:start geometry:
        library = egs_genvelope
        name = A18
        base geometry = capsule
        inscribed geometries = electrode
    :stop geometry:
	 
	:start geometry:
        library = egs_gtransformed
        my geometry = A18
        :start transformation:
            translation = 0 0 -0.9
        :stop transformation:
        name = A18_centered
    :stop geometry:
	
	:start geometry:
        library = egs_box
#        box size = 1.1 1.1 1.8
        box size = 4.4 4.4 7.2
        name = the_box
        :start media input:
            media = 278-Water
        :stop media input:
        :start transformation:
            translation = 0 0 0
        :stop transformation:
    :stop geometry:
	
	:start geometry:
        library = egs_genvelope
        name = A18_in_box
        base geometry = the_box
        inscribed geometries = A18_centered
     :stop geometry:
	 
	simulation geometry = A18_in_box
:stop geometry definition:

###################### The following view control input helps the viewer
#                      to find the defined geometry.
#                      It will fail without it as the geometry is too small
#                      for the position grid used in the initial search for
#                      the optimum viewing position.
#
#:start view control:
#  xmin = -0.55
#  xmax =  0.55
#  ymin = -0.55
#  ymax =  0.55
#  zmin = 0
#  zmax = 1.8
#:stop view control: