"""
sprkkr_parser
Reads inputs and outputs from SPRKKR and writes to UppASD-compatible files
"""

import numpy as np

def parse_sysfile(filename):
    """Reads SPRKKR .sys-file and extracts relevant data. Currently not in use."""
    with open(filename, 'r') as sysfile:
        sysdata = sysfile.read()

    sysrows = sysdata.split('\n')

    alat_idx = sysrows.index('lattice parameter A  [a.u.]') + 1
    alat = np.float64(sysrows[alat_idx].split()[0])

    lattice_idx = sysrows.index('primitive vectors     (cart. coord.) [A]') + 1
    lattice = np.array(''.join(
        sysrows[lattice_idx:lattice_idx+3]).split(), dtype=np.float64).reshape(3, 3)

    nat_idx = sysrows.index('number of sites NQ') + 1
    natoms = np.int32(sysrows[nat_idx])

    pos_idx = nat_idx + 2
    posdata = np.array(
        ''.join(sysrows[pos_idx:pos_idx+natoms]).split(), dtype=np.float64).reshape(natoms, 9)
    # pos_sites = np.int32(posdata[:, 0])
    pos_types = np.int32(posdata[:, 1])
    pos_vectors = posdata[:, 2:5] / alat

    nt_idx = sysrows.index('number of atom types NT') + 1
    ntypes = np.int32(sysrows[nt_idx])

    num_idx = nt_idx + 2
    number_rows = sysrows[num_idx:num_idx+ntypes]
    numbers = np.zeros(natoms, dtype=np.int32)
    for row in number_rows:
        nentry = np.int32(row.split()[3])
        for entry in range(nentry):
            aidx = np.int32(row.split()[entry+5])-1
            numbers[aidx] = row.split()[1]

    return alat, lattice, natoms, ntypes, pos_vectors, numbers, pos_types


def parse_potfile(filename):
    """Reads SPRKKR .pot-file and extracts relevant data."""
    with open(filename, 'r') as potfile:
        potdata = potfile.read()

    potrows = potdata.split('\n')
    potrows = [row.strip() for row in potrows]

    label = potrows[4].split()[1]

    gen_idx = potrows.index('GLOBAL SYSTEM PARAMETER')
    natoms = np.int32(potrows[gen_idx+1].split()[1])
    ntypes = np.int32(potrows[gen_idx+2].split()[1])
    # nm = np.int32(potrows[gen_idx+3].split()[1])
    # irel = np.int32(potrows[gen_idx+4].split()[1])

    alat_idx = potrows.index('LATTICE')+4
    alat = np.float32(potrows[alat_idx].split()[1])

    lattice = []
    for row in range(3):
        lattice.append([np.float32(val)
                       for val in potrows[alat_idx+1+row].split()[1:]])
    lattice = np.array(lattice)

    pos_idx = potrows.index('SITES')+4
    posdata = []
    for irow in range(natoms):
        posdata.append(potrows[pos_idx+irow].split()[1:])
    pos_vectors = np.array(posdata, dtype=np.float64)

    occu_idx = potrows.index('OCCUPATION')+2
    occudata = []
    for irow in range(natoms):
        occudata.append(potrows[occu_idx+irow].split()[:4])

    occu_mat = np.array(occudata, dtype=np.int32)
    pos_types = occu_mat[:, 1]

    number_idx = potrows.index('TYPES')+2
    typedata = []
    for irow in range(ntypes):
        typedata.append(potrows[number_idx+irow].split())
    type_mat = np.array(typedata)
    numrow = np.array(type_mat[:, 2], dtype=np.int32)
    numbers = np.zeros(natoms)
    for idx in range(natoms):
        numbers[idx] = numrow[pos_types[idx]-1]

    return label, alat, lattice, natoms, ntypes, pos_vectors, numbers, pos_types


def parse_inpfile(inpname):
    """Reads SPRKKR .inp-file and extracts relevant data."""
    with open(inpname, 'r') as ifile:
        inpdata = ifile.read()
        inprows = inpdata.split('\n')

        label = [row for row in inprows if 'DATASET' in row][0].split()[3]
        suffix = [row for row in inprows if 'ADSI' in row][0].split()[2]
        potfile = [row for row in inprows if 'POTFIL' in row][0].split()[2]
        task = [row for row in inprows if 'TASK' in row][0].split()[1]
        fullrel = not 'SP-SREL' in inpdata.split()

        return label, suffix, potfile, task, fullrel


def parse_outfile(outfilename, types, natoms, ntypes):
    """Reads SPRKKR .out-file and extracts relevant data."""
    with open(outfilename, 'r') as outfile:
        outdata = outfile.read()
        outrows = outdata.split('\n')

    magmoms = np.array(
        np.array([row.split() for row in outrows if 'sum  ' in row][-ntypes-1:-1])[:, 4]
        , dtype=np.float64)
    moments = np.zeros(natoms)
    if natoms == 1:
        moments[0] = magmoms
    else:
        for idx in range(natoms):
            moments[idx] = magmoms[types[idx]-1]

    return moments


def parse_xcfile(filename):
    """Reads SPRKKR .dat-file and extracts relevant exchange interactions."""
    with open(filename, 'r') as file:
        nheader = 0
        row = []
        while 'DRZ' not in row:
            nheader += 1
            row = file.readline().split()

    columns = [0, 2, 4, 5, 6, 11, 10]
    xcdata = np.genfromtxt(filename, skip_header=nheader)[:, columns]

    return xcdata


def parse_dmi(prefix):
    """Exctracts DMI data from SPRKKR"""

    dmxfile = prefix+'_Dij_z.dat'
    dmyfile = prefix+'_Dij_z.dat'
    dmzfile = prefix+'_Dij_z.dat'

    dmxdata = parse_xcfile(dmxfile)
    dmydata = parse_xcfile(dmyfile)
    dmzdata = parse_xcfile(dmzfile)

    dmdata = np.column_stack([dmxdata[:, 0:6], dmydata[:, 5], dmzdata[:, 5:]])

    return dmdata


def write_posfile(natoms, positions, types, fname):
    """Writes a posfile for UppASD maptype 2 format."""
    idx = np.arange(natoms)+1
    outarr = np.column_stack([idx, types, positions])
    fmt = ' %4i %4i   %12.8f %12.8f %12.8f'
    np.savetxt(fname, outarr, fmt=fmt)


def write_momfile(natoms, moments, fname):
    """Writes a momfile for UppASD. Assumes collinear order."""

    idx = np.arange(natoms)+1
    ones = np.ones(natoms)
    zeros = np.zeros(natoms)
    outarr = np.column_stack(
        [idx, ones, moments, zeros, zeros, np.sign(moments)])
    fmt = ' %4i %4i   %12.8f   %12.8f %12.8f %12.8f'
    np.savetxt(fname, outarr, fmt=fmt)


def write_jfile(xcdata, filename, is_dmi=False, write_full=True):
    """Writes exchange interactions (Jij or DMI) to UppASD jfile/dmfile format"""
    xcfile_fmt = '%4i %4i    % 5i % 5i % 5i    %12.8f   %10.6f'
    d_idx = 6
    if is_dmi:
        xcfile_fmt = '%4i %4i    % 5i % 5i % 5i    %12.8f %12.8f %12.8f   %10.6f'
        d_idx = 8

    if write_full:
        np.savetxt(filename, xcdata[xcdata[:, d_idx] > 0], fmt=xcfile_fmt)


def write_celldata(lattice, filename):
    """Writes celldata and relevant keywords to inpsd.dat"""
    with open(filename, 'w') as outfile:
        outfile.write('# Data parsed from SPRKKR calculation\n')
        outfile.write('cell\n')
        for row in lattice:
            outfile.write(f"{row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}\n")

        outfile.write('\n')
        outfile.write('maptype 2\n')


def parse_sprkkr(ifile):
    """Main parsing wrapper routine"""
    path = '/'.join(ifile.split('/')[:-1]) + '/'

    prefix = ifile[:-4]

    label, suffix, potfile, task, fullrel = parse_inpfile(ifile)

    potfile = path + potfile

    sysname, alat, lattice, natoms, ntypes, positions, numbers, types = parse_potfile(potfile)

    if fullrel:
        xcfile = prefix + '_JJij.dat'
    else:
        xcfile = prefix + '_J_ij.dat'

    outfilename = path + sysname + '_' + suffix + '.out'

    moments = parse_outfile(outfilename, types, natoms, ntypes)

    xcdata = parse_xcfile(xcfile)

    dmidata = parse_dmi(ifile[:-4])

    filenames =  np.char.add(sysname+'.', np.array(['jfile', 'dmfile', 'posfile',
                           'momfile']))

    write_jfile(xcdata, filenames[0])

    write_jfile(dmidata, filenames[1], is_dmi=True)

    write_posfile(natoms, positions, types, filenames[2])

    write_momfile(natoms, moments, filenames[3])

    return filenames, lattice