
# source file
data directory = RealTimeData
source type = FastViewRaySource
phase space data = ViewRayCo.phsp
number of threads = 4
random number seed = 1234
log level = 1
viewray head model = GradientPlusMlc
#viewray head model = Mlc
exit region = 3



Ae = 0.25
Pcut = 0.05
maximum energy loss = 0.5
nsplit = 70
#kerma nsplit = 2

MLC geometry = mlcDC
field margin = 7
small field margin = 0.5
#split factor = 0.025
split factor = 0.1
#split factor = 0.2
source Xo = 0
source Yo = 0
source radius = 1


transmission = 0.03

#source type = FastViewRaySource
#phase space data = ViewRayCo.phsp

#select particles = 1

geometry definition =
{

    #
    # The skeleton where will be inscribing the various geometry coponents
    #
    geometry =
	{
        type = XYZGeometry
        name = skeletonDC
        x-positions = -15 15
        y-positions = -15 15
        z-positions = 25 54 65.39571 79.343  79.393  82.1334  100.3995  100.4495  103.588  106.638 115
        media input =
		{
            media = U  Co Air Steel
            set medium = 0 2
            set medium = 1 2
            set medium = 2 2
            set medium = 4 2
            set medium = 3 3
            set medium = 6 3
        }
    }

    #
    # The Co cylinder
    #
    geometry =
	{
        name = SourceCylDC
        type = ZCylinders
        radii = 1.4
        media input =
		{
            media = Co
        }
    }

    #
    # The shutle (?) opening
    #
    geometry =
	{
        name = AirCylDC
        type = ZCylinders
        radii = 1.4993
        media input =
		{
            media = Air
        }
    }

    # 
    # The primary collimator opening
    #
    geometry =
	{
        name = PyramidDC
        type = SimplePyramid
        base points = -12,-12,0  -12,12,0  12,12,0  12,-12,0
        tip    = 0 0 114.7922
        media input =
		{
            media = Air
        }
    }

    #
    # The MLC
    #
    geometry =
	{
        type = ViewRayMLC
        name = mlcDC
        air medium = Air
        leaf medium = W
        x-focus = 105
        xgap = 0.03
        ygap = 0.02
    }

    # 
    # We need to inscribe the MLC in a box. We also need to model the shielding outside of the MLC
    # We do this by creating a XY geometry and putting the MLC inside using a smart envelope
    #
    geometry =
	{
        name = auxXplanesDC
        type = XPlanes
        positions = -16 -7.5 7.5 16
    }
    geometry =
	{
        name = auxYplanesDC
        type = YPlanes
        positions = -16 -7.5 7.5 16
    }
    geometry =
	{
        name = mlcEnvelopeDC
        type = NdimGeometry
        dimensions = auxXplanesDC auxYplanesDC
        media input =
		{
            media = W Air
            set medium = 4 1
        }
    }
    geometry =
	{
        name = mlcLayerDC
        type = SmartEnvelope
        base geometry = mlcEnvelopeDC
        inscribe geometry = mlcDC 4 1
    }

    #
    # The gradient coils
    #
    geometry =
	{
        name = coilsDC
        type = YCylinders
        radii = 35 35.3 40.0 40.5
        media input =
		{
            media = Air C
            set medium = 1 1
            set medium = 3 1
        }
    }

    #
    # The entire geometry
    #
    geometry =
	{
        type = SmartEnvelope
        name = HeadDC
        base geometry = skeletonDC
        # z-positions = 25 54 65.39571 79.343  79.393  82.1334  100.3995  100.4495  103.588  106.638 115
        #                 0  1        2       3       4        5         6         7       8        9
        inscribe geometry = coilsDC 0 1
        inscribe geometry = mlcLayerDC 1 1
        inscribe geometry = PyramidDC 5 1
        inscribe geometry = AirCylDC  7 1
        inscribe geometry = SourceCylDC 8 1
    }

    simulation geometry = HeadDC

}



