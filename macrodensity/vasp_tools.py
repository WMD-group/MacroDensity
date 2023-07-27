from __future__ import print_function

def get_band_extrema(input_file):
    '''
    Get the valence band maximum and conduction band minimum from VASP OUTCAR.

    This function reads the VASP OUTCAR file and extracts the valence band maximum (VBM) and
    conduction band minimum (CBM). It also checks for partial occupancy and prints a warning
    message if found.

    Args:
        input_file (str): The path to the VASP OUTCAR file.

    Returns:
        list: A list containing the valence band maximum (VBM) and conduction band minimum (CBM).
              list[0] = VBM, list[1] = CBM.

    Example:
        >>> input_file = 'path/to/OUTCAR'
        >>> band_extrema = get_band_extrema(input_file)
        >>> print("Valence Band Maximum (VBM):", band_extrema[0])
        >>> print("Conduction Band Minimum (CBM):", band_extrema[1])
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
