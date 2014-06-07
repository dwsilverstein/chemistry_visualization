#! /usr/bin/env python

from __future__ import print_function, division
from chem import collect
import sys, os
from numpy import array, append, reshape, max
from math import pi as PI
from math import cos, sin, sqrt

def main():
    """\
    This program takes the polarizability or hyperpolarizablity 
    and contracts it with the incident electic field(s).  It is
    used to generate a unit sphere representation at the ground
    state optimized geometry. 
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('polfile', help='Output file of a polarizability '
                        'or hyperpolarizability calculation.')
    parser.add_argument('-r', '--radius', help='Radius of the sphere for '
                        'making images.', required=False, default=1.0)
    parser.add_argument('-o', '--output', help='Name of the output file.',
                        required=False)
    args = parser.parse_args()

    # Collect data from the output file.
    data = collect(args.polfile)

    # Change the radius to a float.
    radius = float(args.radius)

    # Build a spherical grid of unit radius.
    # Theta is the polar angle (0 < Theta < 360)
    # Phi is the azimuthal angle (0 < Phi < 180)
    grid = array([],dtype=float)
    for i in range(10,370,10): # Theta
        theta = float(i)*PI/180.0
        for j in range(10,190,10): # Phi
            phi = float(j)*PI/180.0
            grid = append(grid,[1.0*sin(theta)*cos(phi),
                                1.0*sin(theta)*sin(phi),
                                1.0*cos(theta)])
    grid = reshape(grid,(36,18,3))

    # The contraction of the tensor with the electric field(s) depends
    # on whether the calculated property is the polarizability or 
    # hyperpolarizability.
    if 'POLARIZABILITY' in data.calctype:
        # This should be rewritten to account for calculations with
        # multiple polarizabilities.
        pol = data.polarizability
        polcont = array([],dtype=float)
        norm = array([],dtype=float)

        # Do the contraction (first digit accounts for the frequency
        # dependence):
        # 00 = xx, 01 = xy, 02 = xz
        # 10 = yx, 11 = yy, 12 = yz
        # 20 = zx, 21 = zy, 22 = zz
        for i in range(36): # Theta
            for j in range(18): # Phi
                alphax = (pol[0][0][0]*grid[i][j][0] + 
                          pol[0][0][1]*grid[i][j][1] +
                          pol[0][0][2]*grid[i][j][2])

                alphay = (pol[0][1][0]*grid[i][j][0] + 
                          pol[0][1][1]*grid[i][j][1] +
                          pol[0][1][2]*grid[i][j][2])

                alphaz = (pol[0][2][0]*grid[i][j][0] + 
                          pol[0][2][1]*grid[i][j][1] +
                          pol[0][2][2]*grid[i][j][2])

                polcont = append(polcont,[alphax,alphay,alphaz])
                norm = append(norm, 
                              sqrt(alphax**2 + alphay**2 + alphaz**2))
        polcont = reshape(polcont,(648,3))
        norm = reshape(norm,(36,18,1))
        maximum = max(norm)

        # Rebuild the grid based on the user input radius.
        grid = array([],dtype=float)
        for i in range(10,370,10): # Theta
            theta = float(i)*PI/180.0
            for j in range(10,190,10): # Phi
                phi = float(j)*PI/180.0
                grid = append(grid,[radius*sin(theta)*cos(phi),
                                    radius*sin(theta)*sin(phi),
                                    radius*cos(theta)])
        grid = reshape(grid,(648,3))

        # Build an array to determine the colors of the vectors.
        vcolor = array([],dtype=float)
        for i in range(36): # Theta
            for j in range(18): # Phi
                vcolor = append(vcolor, norm[i][j][0] / maximum)

        # Output the data to a TCL script, to use with VMD (the 
        # extension is irrelevant).
        if args.output:
           outfile = args.output
        else:
           outfile = args.polfile.split('.')[0] + '_unitsphere.tcl'
        f = open(outfile, 'w')
        f.write('color change rgb 11 0.1 0.1 1.0\n')
        f.write('color change rgb 12 0.2 0.2 1.0\n')
        f.write('color change rgb 13 0.3 0.3 1.0\n')
        f.write('color change rgb 14 0.4 0.4 1.0\n')
        f.write('color change rgb 15 0.5 0.5 1.0\n')
        f.write('color change rgb 17 0.6 0.6 1.0\n')
        f.write('color change rgb 18 0.7 0.7 1.0\n')
        f.write('color change rgb 19 0.8 0.8 1.0\n')
        f.write('color change rgb 20 0.9 0.9 1.0\n')
        f.write('color change rgb 21 1.0 0.1 0.1\n')
        f.write('color change rgb 22 1.0 0.2 0.2\n')
        f.write('color change rgb 23 1.0 0.3 0.3\n')
        f.write('color change rgb 24 1.0 0.4 0.4\n')
        f.write('color change rgb 25 1.0 0.5 0.5\n')
        f.write('color change rgb 26 1.0 0.6 0.6\n')
        f.write('color change rgb 27 1.0 0.7 0.7\n')
        f.write('color change rgb 28 1.0 0.8 0.8\n')
        f.write('color change rgb 29 1.0 0.9 0.9\n')
        fmt = '{0}{1:8.5F} {2:8.5F} {3:8.5F}{4}{5:12.5E} {6:12.5E} {7:12.5E}{8}'
        for i in range(648):
            color = vcolor[i]
            if (0.00000<=color and color<0.047619):
                f.write('draw color 0\n') # Blue (smallest intensity)
            elif (0.047619<=color and color<0.095238):
                f.write('draw color 11\n')
            elif (0.095238<=color and color<0.142857):
                f.write('draw color 12\n')
            elif (0.142857<=color and color<0.190476):
                f.write('draw color 13\n')
            elif (0.190476<=color and color<0.238095):
                f.write('draw color 14\n')
            elif (0.238095<=color and color<0.285714):
                f.write('draw color 15\n')
            elif (0.285714<=color and color<0.333333):
                f.write('draw color 17\n')
            elif (0.333333<=color and color<0.380952):
                f.write('draw color 18\n')
            elif (0.380952<=color and color<0.428571):
                f.write('draw color 19\n')
            elif (0.428571<=color and color<0.476190):
                f.write('draw color 20\n')
            elif (0.476190<=color and color<0.523809):
                f.write('draw color 8\n') # White
            elif (0.523809<=color and color<0.571428):
                f.write('draw color 29\n')
            elif (0.571428<=color and color<0.619047):
                f.write('draw color 28\n')
            elif (0.619047<=color and color<0.666666):
                f.write('draw color 27\n')
            elif (0.666666<=color and color<0.714285):
                f.write('draw color 26\n')
            elif (0.714285<=color and color<0.761904):
                f.write('draw color 25\n')
            elif (0.761904<=color and color<0.809523):
                f.write('draw color 24\n')
            elif (0.809523<=color and color<0.857142):
                f.write('draw color 23\n')
            elif (0.857142<=color and color<0.904761):
                f.write('draw color 22\n')
            elif (0.904761<=color and color<0.952380):
                f.write('draw color 21\n')
            elif (0.952380<=color and color<=1.000000):
                f.write('draw color 1\n') # Red (largest intensity)
            f.write(fmt.format('vmd_draw_vector 0 {', grid[i][0],
                               grid[i][1], grid[i][2], '} {',
                               polcont[i][0], polcont[i][1],
                               polcont[i][2], '} 0.1 30 0.08\n'))
        f.close()
    elif 'HYPERPOLARIZABILITY' in data.calctype:
        # Types of hyperpolarizabilities.
        htypes = ('SHG', 'EOPE', 'OR', 'STATIC',) 

        # Grab the different hyperpolarizabilities.
        hpol = {}
        tlist = []
        for item in htypes:
            if item in data.calctype:
                hpol[item] = data.hyperpolarizability[item]
                tlist.append(item)

        # Do the contraction:
        # 000 = xxx, 001 = xxy, 002 = xxz
        # 010 = xyx, 011 = xyy, 012 = xyz
        # 020 = xzx, 021 = xzy, 022 = xzz
        #
        # 100 = yxx, 101 = yxy, 102 = yxz
        # 110 = yyx, 111 = yyy, 112 = yyz
        # 120 = yzx, 121 = yzy, 122 = yzz
        #
        # 200 = zxx, 201 = zxy, 202 = zxz
        # 210 = zyx, 211 = zyy, 212 = zyz
        # 220 = zzx, 221 = zzy, 222 = zzz 
        for item in tlist:
            hpolcont = array([],dtype=float)
            norm = array([],dtype=float)
            for i in range(36): # Theta
                for j in range(18): # Phi
                    betax = (hpol[item][0][0][0]*grid[i][j][0]*grid[i][j][0] +
                             hpol[item][0][0][1]*grid[i][j][0]*grid[i][j][1] +
                             hpol[item][0][0][2]*grid[i][j][0]*grid[i][j][2] +
                             hpol[item][0][1][0]*grid[i][j][1]*grid[i][j][0] +
                             hpol[item][0][1][1]*grid[i][j][1]*grid[i][j][1] +
                             hpol[item][0][1][2]*grid[i][j][1]*grid[i][j][2] +
                             hpol[item][0][2][0]*grid[i][j][2]*grid[i][j][0] +
                             hpol[item][0][2][1]*grid[i][j][2]*grid[i][j][1] +
                             hpol[item][0][2][2]*grid[i][j][2]*grid[i][j][2]) 

                    betay = (hpol[item][1][0][0]*grid[i][j][0]*grid[i][j][0] +
                             hpol[item][1][0][1]*grid[i][j][0]*grid[i][j][1] +
                             hpol[item][1][0][2]*grid[i][j][0]*grid[i][j][2] +
                             hpol[item][1][1][0]*grid[i][j][1]*grid[i][j][0] +
                             hpol[item][1][1][1]*grid[i][j][1]*grid[i][j][1] +
                             hpol[item][1][1][2]*grid[i][j][1]*grid[i][j][2] +
                             hpol[item][1][2][0]*grid[i][j][2]*grid[i][j][0] +
                             hpol[item][1][2][1]*grid[i][j][2]*grid[i][j][1] +
                             hpol[item][1][2][2]*grid[i][j][2]*grid[i][j][2]) 

                    betaz = (hpol[item][2][0][0]*grid[i][j][0]*grid[i][j][0] +
                             hpol[item][2][0][1]*grid[i][j][0]*grid[i][j][1] +
                             hpol[item][2][0][2]*grid[i][j][0]*grid[i][j][2] +
                             hpol[item][2][1][0]*grid[i][j][1]*grid[i][j][0] +
                             hpol[item][2][1][1]*grid[i][j][1]*grid[i][j][1] +
                             hpol[item][2][1][2]*grid[i][j][1]*grid[i][j][2] +
                             hpol[item][2][2][0]*grid[i][j][2]*grid[i][j][0] +
                             hpol[item][2][2][1]*grid[i][j][2]*grid[i][j][1] +
                             hpol[item][2][2][2]*grid[i][j][2]*grid[i][j][2]) 
    
                    hpolcont = append(hpolcont,[betax,betay,betaz])
                    norm = append(norm, 
                                  sqrt(betax**2 + betay**2 + betaz**2))
            hpolcont = reshape(hpolcont,(648,3))
            norm = reshape(norm,(36,18,1))
            maximum = max(norm)

            # Rebuild the grid based on the user input radius.
            grid = array([],dtype=float)
            for i in range(10,370,10): # Theta
                theta = float(i)*PI/180.0
                for j in range(10,190,10): # Phi
                    phi = float(j)*PI/180.0
                    grid = append(grid,[radius*sin(theta)*cos(phi),
                                        radius*sin(theta)*sin(phi),
                                        radius*cos(theta)])
            grid = reshape(grid,(648,3))
    
            # Build an array to determine the colors of the vectors.
            vcolor = array([],dtype=float)
            for i in range(36): # Theta
                for j in range(18): # Phi
                    vcolor = append(vcolor, norm[i][j][0] / maximum)
    
            # Output the data to a TCL script, to use with VMD (the 
            # extension is irrelevant).
            if args.output:
               outfile = args.output
            else:
               outfile = args.polfile.split('.')[0] + '_' + item.lower() + '_unitsphere.tcl'
            f = open(outfile, 'w')
            f.write('color change rgb 11 0.1 0.1 1.0\n')
            f.write('color change rgb 12 0.2 0.2 1.0\n')
            f.write('color change rgb 13 0.3 0.3 1.0\n')
            f.write('color change rgb 14 0.4 0.4 1.0\n')
            f.write('color change rgb 15 0.5 0.5 1.0\n')
            f.write('color change rgb 17 0.6 0.6 1.0\n')
            f.write('color change rgb 18 0.7 0.7 1.0\n')
            f.write('color change rgb 19 0.8 0.8 1.0\n')
            f.write('color change rgb 20 0.9 0.9 1.0\n')
            f.write('color change rgb 21 1.0 0.1 0.1\n')
            f.write('color change rgb 22 1.0 0.2 0.2\n')
            f.write('color change rgb 23 1.0 0.3 0.3\n')
            f.write('color change rgb 24 1.0 0.4 0.4\n')
            f.write('color change rgb 25 1.0 0.5 0.5\n')
            f.write('color change rgb 26 1.0 0.6 0.6\n')
            f.write('color change rgb 27 1.0 0.7 0.7\n')
            f.write('color change rgb 28 1.0 0.8 0.8\n')
            f.write('color change rgb 29 1.0 0.9 0.9\n')
            fmt = '{0}{1:8.5F} {2:8.5F} {3:8.5F}{4}{5:12.5E} {6:12.5E} {7:12.5E}{8}'
            for i in range(648):
                color = vcolor[i]
                if (0.00000<=color and color<0.047619):
                    f.write('draw color 0\n') # Blue (smallest intensity)
                elif (0.047619<=color and color<0.095238):
                    f.write('draw color 11\n')
                elif (0.095238<=color and color<0.142857):
                    f.write('draw color 12\n')
                elif (0.142857<=color and color<0.190476):
                    f.write('draw color 13\n')
                elif (0.190476<=color and color<0.238095):
                    f.write('draw color 14\n')
                elif (0.238095<=color and color<0.285714):
                    f.write('draw color 15\n')
                elif (0.285714<=color and color<0.333333):
                    f.write('draw color 17\n')
                elif (0.333333<=color and color<0.380952):
                    f.write('draw color 18\n')
                elif (0.380952<=color and color<0.428571):
                    f.write('draw color 19\n')
                elif (0.428571<=color and color<0.476190):
                    f.write('draw color 20\n')
                elif (0.476190<=color and color<0.523809):
                    f.write('draw color 8\n') # White
                elif (0.523809<=color and color<0.571428):
                    f.write('draw color 29\n')
                elif (0.571428<=color and color<0.619047):
                    f.write('draw color 28\n')
                elif (0.619047<=color and color<0.666666):
                    f.write('draw color 27\n')
                elif (0.666666<=color and color<0.714285):
                    f.write('draw color 26\n')
                elif (0.714285<=color and color<0.761904):
                    f.write('draw color 25\n')
                elif (0.761904<=color and color<0.809523):
                    f.write('draw color 24\n')
                elif (0.809523<=color and color<0.857142):
                    f.write('draw color 23\n')
                elif (0.857142<=color and color<0.904761):
                    f.write('draw color 22\n')
                elif (0.904761<=color and color<0.952380):
                    f.write('draw color 21\n')
                elif (0.952380<=color and color<=1.000000):
                    f.write('draw color 1\n') # Red (largest intensity)
                f.write(fmt.format('vmd_draw_vector 0 {', grid[i][0],
                                   grid[i][1], grid[i][2], '} {',
                                   hpolcont[i][0], hpolcont[i][1],
                                   hpolcont[i][2], '} 0.1 30 0.08\n'))
            f.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
