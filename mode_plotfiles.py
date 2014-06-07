#! /usr/bin/env python

from __future__ import print_function, division
from chem import collect
import sys, os

def main():
    """\
    This program takes the normal modes from a frequencies file and 
    generates files for plotting the normal modes in VMD.  

    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('freqfile', help='File containing the '
                        'normal modes')
    parser.add_argument('-s', '--scale', help='Scale factor for the normal '
                        'mode vector', default=2.00) 
    parser.add_argument('-v', '--vfreq', help='Requested normal mode frequency',
                        default='all')
    parser.add_argument('--low', help='Lowest normal mode frequency used for'
                        ' making TCL scripts', default=float(400), type=float)
    parser.add_argument('--high', help='Highest normal mode frequency used'
                        ' for making TCL scripts', default=float(1800), type=float)
    args = parser.parse_args()

    # Collect the data from the frequency file.
    f = collect(args.freqfile) 

    # Initialize some variables
    previous_freq = ''
    degeneracy = 0
    deg_list = ('', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i')
    scale = float(args.scale)

    # Keep track of number of files skipped
    skipped = { 'range' : 0 , 'negative' : 0 }
    tot = 0

    # Generate the output file(s) for VMD (TCL scripts)
    if args.vfreq != 'all':
        # Get the user input
        freq = float(args.vfreq)

        # Check for degeneracy.  If there is degeneracy, warn
        # the user, but continue as normal.
        deg_freq_indices = []
        for i in range(f.nmodes):
            vfreq = f.v_frequencies[i]
            if vfreq == freq:
                deg_freq_indices.append(i)
                degeneracy += 1
        warnfmt = '{0} {1} {2} {3}{4}'
        if degeneracy > 1:
            print(warnfmt.format('Warning, the mode at frequency',
                                 str(freq), 'is', str(degeneracy),
                                 '-fold degenerate.'))
            print("I'll make a file for each degenerate mode.")

        # Make TCL scripts as requested.
        fmt = ('{0} {1} {9} {2:>13.8f} {3:>13.8f} {4:>13.8f}{10} '
               '{9} {5:>13.8f} {6:>13.8f} {7:>13.8f}{10} {8:4.2f}')
        for ix in range(degeneracy):
            vfreq = f.v_frequencies[i]
            strfreq = '{0:.2f}'.format(round(freq,2))
            if ix > 0:
                strfreq = '_'.join([strfreq,deg_list[ix]])
            with open('mode' + strfreq + '-vmd.tcl', 'w') as fl:
                # Copy the instance from above.
                new = f.copy()

                # Replace the coordinates of new with the coordinates with 
                # the mode added, for the purpose of making an image.
                nmode = new.normal_modes[deg_freq_indices[ix]]
                new.coordinates = scale*nmode

                for n in range(f.natoms):
                    print(fmt.format('vmd_draw_vector2', '0',
                                     f.coordinates[n][0], f.coordinates[n][1],
                                     f.coordinates[n][2], new.coordinates[n][0],
                                     new.coordinates[n][1], new.coordinates[n][2],
                                     scale, '{', '}'),file=fl)
    else:
        for i in range(f.nmodes):

            freq = f.v_frequencies[i]
            tot += 1

            # Check for negative frequencies, and if the 
            # mode falls outside of the range.
            if freq < 0.0:
                skipped['negative'] += 1
                continue  
            if freq < args.low or freq > args.high:
                skipped['range'] += 1
                continue 

            # Make a string out of the current mode.
            strfreq = '{0:.2f}'.format(round(freq,2))

            # Check for degeneracy, so that we account for all modes.
            if strfreq in previous_freq:
                degeneracy += 1
                strfreq = '_'.join([strfreq,deg_list[degeneracy]])
            else:
                degeneracy = 0

            # Create the TCL file for the normal mode.
            fmt = ('{0} {1} {9} {2:>13.8f} {3:>13.8f} {4:>13.8f}{10} '
                   '{9} {5:>13.8f} {6:>13.8f} {7:>13.8f}{10} {8:4.2f}')
            with open('mode' + strfreq + '-vmd.tcl', 'w') as fl:
                # Copy the instance from above.
                new = f.copy()

                # Replace the coordinates of new with the coordinates with 
                # the mode added, for the purpose of making an image.
                nmode = new.normal_modes[i]
                new.coordinates = scale*nmode

                print('draw color yellow',file=fl)
                for n in range(f.natoms):
                    print(fmt.format('vmd_draw_vector2', '0',
                                     f.coordinates[n][0], f.coordinates[n][1],
                                     f.coordinates[n][2], new.coordinates[n][0],
                                     new.coordinates[n][1], new.coordinates[n][2],
                                     scale, '{', '}'),file=fl) 

            # Save this frequency for checking purposes
            previous_freq = strfreq

        if skipped['negative']:
            print('Skipped ', skipped['negative'], ' imaginary frequencies.')
        if skipped['range']:
            print('Skipped ', skipped['range'], ' out of range frequencies.')
        print('The total number of normal modes is: ', tot)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
