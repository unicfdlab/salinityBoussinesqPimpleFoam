dnl>
dnl> Initial parameters
dnl>
define(toMeters, 0.001)dnl>
define(cx, 0.0)dnl>
define(cy, 0.0)dnl>
define(cz, 0.0)dnl>
define(dy, 1.0)dnl>
define(L, 450.0)dnl>
define(H, 300.0)dnl>
define(alphag,27.3)dnl>
define(xdivs, 225)dnl>
define(ydivs, 1)dnl>
define(zdivs, 150)dnl>
dnl>
dnl> Functions definitions
dnl>
changecom(//) dnl>
changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'print ($1)')]) dnl>
define(VCOUNT, 0)  dnl>
define(vlabel, [[// ]point VCOUNT ($1) define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])  dnl>
define(vert,
    [($1 $2 $3)] dnl>
    [vlabel($4)] dnl>
    ) dnl> 
define(hex2D, hex (f$1 f$2 f$3 f$4 b$1 b$2 b$3 b$4)) dnl>
define(quad, ($1 $2 $3 $4))dnl>
define(defpatchWall,
    $1
    [{]
	[type]	[wall];
	faces
	(
	    $2
	);
    [}]
    ) dnl>
define(defpatchEmpty,
    $1
    [{]
	[type]	[empty];
	faces
	(
	    $2
	);
    [}]
    ) dnl>
dnl>
dnl> Some characteristic lengths
dnl>
define(alpha,calc(alphag*3.14159265359/180.0))dnl>
define(sina,calc(sin(alpha)))dnl>
define(cosa,calc(cos(alpha)))dnl>
define(tga,calc(sina/cosa))dnl>
define(L1,calc(L - tga*H))dnl>
define(minx,cx)dnl>
define(maxx,calc(minx+L))dnl>
define(maxx1,calc(minx+L1))dnl>
define(miny,cy)dnl>
define(maxy,calc(cy+dy))dnl>
define(minz,cz)dnl>
define(maxz,calc(cz+H))dnl>


/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.1.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    [format]      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters toMeters;

vertices
(
    vert(minx, miny, minz, bv1)
    vert(maxx, miny, minz, bv2)
    vert(maxx1,miny, maxz, bv3)
    vert(minx, miny, maxz, bv4)

    vert(minx, maxy, minz, fv1)
    vert(maxx, maxy, minz, fv2)
    vert(maxx1,maxy, maxz, fv3)
    vert(minx, maxy, maxz, fv4)

);

blocks
(
    hex2D(v1, v2, v3, v4) ( xdivs zdivs ydivs ) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
	defpatchEmpty(frontAndBack, [quad(bv1,bv2,bv3,bv4) quad(fv1,fv2,fv3,fv4)])
	defpatchWall(left, [quad(bv1,bv4,fv4,fv1)])
	defpatchWall(right, [quad(bv2,bv3,fv3,fv2)])
	defpatchWall(top, [quad(bv3,bv4,fv4,fv3)])
	defpatchWall(bottom, [quad(bv1,bv2,fv2,fv1)])
);

mergePatchPairs
(
);

// ************************************************************************* //
