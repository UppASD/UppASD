OVF Parser Library
=================================
**Simple API for powerful OOMMF Vector Field file parsing**<br />

[OVF format specification](https://math.nist.gov/oommf/doc/userguide20a0/userguide/Vector_Field_File_Format_OV.html)

**[Python package](https://pypi.org/project/ovf/):** [![PyPI version](https://badge.fury.io/py/ovf.svg)](https://badge.fury.io/py/ovf)

How to use
---------------------------------

For usage examples, take a look into the test folders: [test](https://github.com/spirit-code/ovf/tree/master/test), [python/test](https://github.com/spirit-code/ovf/tree/master/python/test) or [fortran/test](https://github.com/spirit-code/ovf/tree/master/fortran/test).

Except for opening a file or initializing a segment, all functions return status codes
(generally `OVF_OK` or `OVF_ERROR`). When the return code is not `OVF_OK`, you can take
a look into the latest message, which should tell you what the problem was
(`const char * ovf_latest_message(struct ovf_file *)` in the C API).

In C/C++ and Fortran, before writing a segment, make sure the `ovf_segment` you pass in is
initialized, i.e. you already called either `ovf_read_segment_header` or `ovf_segment_initialize`.

### C/C++

Opening and closing:

- `struct ovf_file *myfile = ovf_open("myfilename.ovf")` to open a file
- `myfile->found` to check if the file exists on disk
- `myfile->is_ovf` to check if the file contains an OVF header
- `myfile->n_segments` to check the number of segments the file should contain
- `ovf_close(myfile);` to close the file and free resources

Reading from a file:

- `struct ovf_segment *segment = ovf_segment_initialize()` to initialize a segment and get the pointer
- `ovf_read_segment_header(myfile, index, segment)` to read the header into the segment struct
- create float data array of appropriate size...
- `ovf_read_segment_data_4(myfile, index, segment, data)` to read the segment data into your float array

Writing and appending to a file:

- `struct ovf_segment *segment = ovf_segment_initialize()` to initialize a segment and get the pointer
- `segment->n_cells[0] = ...` etc to set data dimensions, title and description, etc.
- `ovf_write_segment_4(myfile, segment, data, OVF_FORMAT_TEXT)` to write a file containing the segment header and data
- `ovf_append_segment_4(myfile, segment, data, OVF_FORMAT_TEXT)` to append the segment header and data to the file

### Python

To install the *ovf python package*, either build and install from source
or simply use

    pip install ovf

To use `ovf` from Python, e.g.

```Python
from ovf import ovf
import numpy as np

data = np.zeros((2, 2, 1, 3), dtype='f')
data[0,1,0,:] = [3.0, 2.0, 1.0]

with ovf.ovf_file("out.ovf") as ovf_file:

    # Write one segment
    segment = ovf.ovf_segment(n_cells=[2,2,1])
    if ovf_file.write_segment(segment, data) != -1:
        print("write_segment failed: ", ovf_file.get_latest_message())

    # Add a second segment to the same file
    data[0,1,0,:] = [4.0, 5.0, 6.0]
    if ovf_file.append_segment(segment, data) != -1:
        print("append_segment failed: ", ovf_file.get_latest_message())
```

### Fortran

The Fortran bindings are written in object-oriented style for ease of use.
Writing a file, for example:

```fortran
type(ovf_file)      :: file
type(ovf_segment)   :: segment
integer             :: success
real(kind=4), allocatable :: array_4(:,:)
real(kind=8), allocatable :: array_8(:,:)

! Initialize segment
call segment%initialize()

! Write a file
call file%open_file("fortran/test/testfile_f.ovf")
segment%N_Cells = [ 2, 2, 1 ]
segment%N = product(segment%N_Cells)

allocate( array_4(3, segment%N) )
array_4 = 0
array_4(:,1) = [ 6.0, 7.0, 8.0 ]
array_4(:,2) = [ 5.0, 4.0, 3.0 ]

success = file%write_segment(segment, array_4, OVF_FORMAT_TEXT)
if ( success == OVF_OK) then
    write (*,*) "test write_segment succeeded."
    ! write (*,*) "n_cells = ", segment%N_Cells
    ! write (*,*) "n_total = ", segment%N
else
    write (*,*) "test write_segment did not work. Message: ", file%latest_message
    STOP 1
endif
```

For more information on how to generate modern Fortran bindings,
see also https://github.com/MRedies/Interfacing-Fortran

How to embed it into your project
---------------------------------

TODO...


Build
---------------------------------

### On Unix systems

Usually:
```
mkdir build
cd build
cmake ..
make
```

### On Windows

One possibility:
- open the folder in the CMake GUI
- generate the VS project
- open the resulting project in VS and build it

### CMake Options

The following options are `ON` by default.
If you want to switch them off, just pass `-D<OPTION>=OFF` to CMake,
e.g. `-DOVF_BUILD_FORTRAN_BINDINGS=OFF`.

- `OVF_BUILD_PYTHON_BINDINGS`
- `OVF_BUILD_FORTRAN_BINDINGS`
- `OVF_BUILD_TEST`

On Windows, you can also set these from the CMake GUI.

### Create and install the Python package

Instead of `pip`-installing it, you can e.g. build everything
and then install the package locally, where the `-e` flag will
let you change/update the package without having to re-install it.

```
cd python
pip install -e .
```

### Build without CMake

The following is an example of how to manually build the C library and
link it with bindings into a corresponding Fortran executable, using gcc.

C library:
```
g++ -DFMT_HEADER_ONLY -Iinclude -fPIC -std=c++11 -c src/ovf.cpp -o ovf.cpp.o

# static
ar qc libovf_static.a ovf.cpp.o
ranlib libovf_static.a

# shared
g++ -fPIC -shared -lc++ ovf.cpp.o -o libovf_shared.so
```

C/C++ test executable:
```
g++ -Iinclude -Itest -std=c++11 -c test/main.cpp -o main.cpp.o
g++ -Iinclude -Itest -std=c++11 -c test/simple.cpp -o simple.cpp.o

# link static lib
g++ -lc++ libovf_static.a main.cpp.o simple.cpp.o -o test_cpp_simple

# link shared lib
g++ libovf_shared.so main.cpp.o simple.cpp.o -o test_cpp_simple
```

Fortran library:
```
gfortran -fPIC -c fortran/ovf.f90 -o ovf.f90.o

ar qc libovf_fortran.a libovf_static.a ovf.f90.o
ranlib libovf_fortran.a
```

Fortran test executable
```
gfortran -c fortran/test/simple.f90 -o simple.f90.o
gfortran -lc++ libovf_fortran.a simple.f90.o -o test_fortran_simple
```

When linking statically, you can also link the object file `ovf.cpp.o` instead of `libovf_static.a`.

*Note: depending on compiler and/or system, you may need `-lstdc++` instead of `-lc++`.*



File format specification (Work in Progress!)
=================================

This library supports the OVF 2.0 format, but provides some extensions to it.
Extensions this library currently makes are documented below, but are subject to change.

After the header, the file consists of Segment blocks, each composed of a header and a data block.

General
---------------------------------

- An OVF file has an ASCII header and trailer, and a data block that may be either ASCII or binary.
- All non-data lines begin with a `#` character
- Comments start with `##` and are ignored by the parser. A comment continues until the end of the line, there is no line continuation character
- Lines starting with a `#` but containing only whitespace are ignored
- *Extension:* lines starting with a `#` but containing an unknown keyword are ignored.
- All non-empty non-comment lines in the file header are structured as pairs of label and value.
    - The label tag consists of all characters after the initial `#` up to the first colon (`:`) character. Case is ignored, and all whitespace is removed
    - The value consists of all characters after the first colon, continuing up to a `##` comment designator or the end of the line.
- Data ordering is generally with the x index incremented first, then the y index, and the z index last

Header
---------------------------------

- The first line of an OVF 2.0 file must contain `# OOMMF OVF 2.0`
- The header should also contain the number of segments, specified as e.g. `# Segment count: 000001`
- *Extension:* the count is padded to 6 digits with zeros (this is so that segments can be appended and the count incremented without having to re-write the entire file)

Segments
---------------------------------

**Segment Header**

- Each block begins with a `# Begin: <block type>` line, and ends with a corresponding `# End: <block type>` line
- Everything inside the `Header` block should be either comments or one of the following file keyword lines. The order of keywords is not specified. *Extension:* all have default values, so none are required unless otherwise stated
    - `title`: long file name or title
    - `desc`: description line, use as many as desired
    - `meshunit`: fundamental mesh spatial unit, treated as a label. The comment marker `##` is not allowed in this label. Example value: `nm`
    - `valueunits`: should be a N-item (Tcl) list of value units, each treated as an unparsed label. The list should either have a length of valuedim, in which case each element denotes the units for the corresponding dimension index, or else the list should have length one, in which case the single element is applied to all dimension indexes. The comment marker `##` is not allowed in this label. Example: `"kA/m"`
    - `valuelabels`: This should be a N-item (Tcl) list of value labels, one for each value dimension. The labels identify the quantity in each dimension. For example, in an energy density file, N would be 1, valueunits could be `"J/m3"` , and valuelabels might be `"Exchange energy density"`
    - `valuedim`: *Extension:* this always specifies the vector dimensionality
    - `xmin`, `ymin`, `zmin`, `xmax`, `ymax`, `zmax`: six separate lines, specifying the bounding box for the mesh
    - `boundary`: two (x,y,z) triples specifying the two relevant vertices of a boundary frame. *Extension:* can be used instead of `xmin` etc.
    - `meshtype`: grid structure; should be either `rectangular` or `irregular`
    - `pointcount`: number of data sample points/locations per cell (here called node)
    - `xnodes`, `ynodes`, `znodes`: three separate lines, specifying the number of nodes along each axis (integers). *Extension:* for an irregular grid, the entire grid will be repeated
    - `xbase`, `ybase`, `zbase`: three (x,y,z) triples, specifying the bravais vectors along which the basis cell is repeated
    - `xstepsize`, `ystepsize`, `zstepsize`: three separate lines - scale factors for the bravais vectors

Standard OVF specification:
- `xmin`, `ymin`, `zmin`, `xmax`, `ymax`, `zmax`: six separate lines, specifying the bounding box for the mesh, in units of `meshunit`
- `boundary`: List of (x,y,z) triples specifying the vertices of a boundary frame. *Extension:* can be used instead of `xmin` etc.
- `meshtype`: grid structure; should be either `rectangular` or `irregular`. For an irregular grid, `pointcount` should also be specified. rectangular grid files should specify instead `xbase`, `ybase`, `zbase`, `xstepsize`, `ystepsize`, `zstepsize`, and `xnodes`, `ynodes`, `znodes`
- `pointcount`: number of data sample points/locations, i.e., nodes (integer). For irregular grids only.
- `xbase`, `ybase`, `zbase`: three separate lines, denoting the position of the first point in the data section, in units of `meshunit`. For rectangular grids only.
- `xstepsize`, `ystepsize`, `zstepsize`: three separate lines, specifying the distance between adjacent grid points, in units of `meshunit`
- `xnodes`, `ynodes`, `znodes`: three separate lines, specifying the number of nodes along each axis (integers). For rectangular grids only. *Extension:* also for irregular grids, where it defines the repetitions of the irregular set of points?

**Segment Data**

- The data block start is marked by a line of the form  `# Begin: data <representation>` (and therefore closed by `# End: data <representation>`), where `<representation>` is one of `text`, `binary 4`, or `binary 8`
- *Extension:* `csv` is also a valid representation and corresponds to comma-separated columns of `text` type
- The binary representations are IEEE floating point in little endian (LSB) order. To ensure that the byte order is correct, and to provide a partial check that the file hasn't been sent through a non 8-bit clean channel, the first data value is fixed to `1234567.0` for 4-byte mode, corresponding to the LSB hex byte sequence `38 B4 96 49`, and `123456789012345.0` for 8-byte mode, corresponding to the LSB hex byte sequence `40 DE 77 83 21 12 DC 42`)
- The data immediately follows the check value
- The first character after the last data value should be a newline