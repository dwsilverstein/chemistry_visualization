#! /usr/bin/env python

from __future__ import print_function, division
from chem import collect
import sys, os
from numpy import array, abs

def main():
    """\
    This program takes Gaussian cube files and either scales the 
    MO density, or calculates the difference in two cube files. 
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('cubefiles', help='Cube files containing the MOs', 
                        nargs='+')
    parser.add_argument('-s', '--scale', help='Value to scale the values from '
                        'a cube file or difference between cube files.', 
                        required=False, default=1.0)
    parser.add_argument('-d', '--difference', help='Argument for making '
                        'the difference between two cube files.', required=False,
                        default=False)
    parser.add_argument('-o', '--output', help='Name of the output file.',
                        required=True)
    args = parser.parse_args()

    # What happens is determined by the command line options.
    if (args.difference):
        # Collect from two Gaussian cube files.
        c1 = collect(args.cubefiles[0]) 
        c2 = collect(args.cubefiles[1]) 
        # Take the difference between two Gaussian cube files.
        #difference = c1.orbitalplot + c2.orbitalplot
        # Absolute value is required to match other codes.  It is
        # debatable whether this is the correct method for taking
        # differences, since the lobes of the MOs can be positive or
        # negative. 
        difference = abs(c1.orbitalplot) - abs(c2.orbitalplot)
        # Find the number of rows (floor division)
        numrows = c1.orbitaldim[0] // 6
        # Find the remainder (number of elements not fitting in a
        # row of 6)
        remainder = c1.orbitaldim[0] - 6 * numrows
        # Formatting for output
        fmt6 = ' {0:12.5E} {1:12.5E} {2:12.5E} {3:12.5E} {4:12.5E} {5:12.5E}{6}'
        fmt5 = ' {0:12.5E} {1:12.5E} {2:12.5E} {3:12.5E} {4:12.5E}{5}'
        out = open(args.output, 'w')
        for item in difference:
            for elem in range(numrows):
                tmp = fmt6.format(item[elem*6], item[elem*6+1],
                                  item[elem*6+2], item[elem*6+3],
                                  item[elem*6+4], item[elem*6+5], '\n')
                out.write(tmp)
            tmp = fmt5.format(item[-5], item[-4], item[-3],
                              item[-2], item[-1], '\n')
            out.write(tmp)
        out.close()
    else:
        # Collect from a Gaussian cube file.
        c1 = collect(args.cubefiles[0]) 

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
