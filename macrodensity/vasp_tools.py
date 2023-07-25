from __future__ import print_function


def get_band_extrema(input_file):
    '''
    Get the valence band maximum and conduction band minimum from VASP OUTCAR

    Prints a warning in the case of partial occupancy.
    Args:
        input_file : String, the input file name
    Returns:
        list: list[0] = valence band maximum, list[1] = conduction band minimum
    '''
    lines = open(input_file, 'r').readlines()
    for line in lines:
        if line.rfind('NKPTS') > -1:
            nkpts = int(line.split()[3])
        if line.rfind('ISPIN') > -1:
            ispin = int(line.split()[2])
        if line.rfind('NELECT') > -1:
            nelect = float(line.split()[2])
    if ispin == 1:
        top_band = int(nelect/2)
    else:
        top_band = int(nelect)

    vbm = []
    cbm = []
    for i, line in enumerate(lines):
        if line.rfind('No.') > -1:
            vbm.append(lines[i + top_band].split()[1])
            cbm.append(lines[i + top_band + 1].split()[1])
            if (float(lines[i + top_band].split()[2]) != 1.00 and
                float(lines[i + top_band].split()[2]) != 2.000):
                print('Partial occupancy, be aware!',
                      lines[i + top_band].split()[2])

    return [float(max(vbm)), float(min(cbm))]
