#pragma once
#ifndef LIBOVF_H
#define LIBOVF_H

// Platform-specific definition of DLLEXPORT
#ifdef _WIN32
    #ifdef __cplusplus
        #define DLLEXPORT extern "C" __declspec(dllexport)
    #else
        #define DLLEXPORT __declspec(dllexport)
    #endif
#else
    #ifdef __cplusplus
        #define DLLEXPORT extern "C"
    #else
        #define DLLEXPORT
    #endif
#endif

/* return codes */
#define OVF_OK          -1
#define OVF_ERROR       -2
#define OVF_INVALID     -3

/* OVF data formats */
#define OVF_FORMAT_BIN   -53
#define OVF_FORMAT_BIN4  -54
#define OVF_FORMAT_BIN8  -55
#define OVF_FORMAT_TEXT  -56
#define OVF_FORMAT_CSV   -57

/* all header info on a segment */
struct ovf_segment {
    char *title;
    char *comment;

    int valuedim;
    char *valueunits;
    char *valuelabels;

    /* the geometrical information on the vector field */
    char *meshtype;
    char *meshunits;
    int pointcount;

    int n_cells[3];
    int N;

    float bounds_min[3];
    float bounds_max[3];

    float lattice_constant;
    float bravais_vectors[3][3];

    /* then some "private" internal fields */
};

/* opaque handle which holds the file pointer */
struct ovf_file_handle;

/* the main struct which keeps the info on the main header of a file */
struct ovf_file {
    const char * filename;
    /* file could be found */
    int found;
    /* file contains an ovf header */
    int is_ovf;
    /* number of segments the file should contain */
    int n_segments;

    /* then some "private" internal fields */
    struct ovf_file_handle *_file_handle;
};

/* opening a file will fill the struct and prepare everything for read/write */
DLLEXPORT struct ovf_file * ovf_open(const char *filename);

/* create a default-initialized segment struct */
DLLEXPORT struct ovf_segment * ovf_segment_initialize();

/* read the geometry info from a segment header */
DLLEXPORT int ovf_read_segment_header(struct ovf_file *, int index, struct ovf_segment *);

/* This function checks the segment in the file against the passed segment and,
    if the dimensions fit, will read the data into the passed array. */
DLLEXPORT int ovf_read_segment_data_4(struct ovf_file *, int index, const struct ovf_segment *, float *data);
DLLEXPORT int ovf_read_segment_data_8(struct ovf_file *, int index, const struct ovf_segment *, double *data);

/* write a segment (header and data) to the file, overwriting all contents.
    The new header will have segment count = 1 */
DLLEXPORT int ovf_write_segment_4(struct ovf_file *, const struct ovf_segment *, float *data, int format=OVF_FORMAT_BIN);
DLLEXPORT int ovf_write_segment_8(struct ovf_file *, const struct ovf_segment *, double *data, int format=OVF_FORMAT_BIN);

/* append a segment (header and data) to the file.
    The segment count will be incremented */
DLLEXPORT int ovf_append_segment_4(struct ovf_file *, const struct ovf_segment *, float *data, int format=OVF_FORMAT_BIN);
DLLEXPORT int ovf_append_segment_8(struct ovf_file *, const struct ovf_segment *, double *data, int format=OVF_FORMAT_BIN);

/* retrieve the most recent error message and clear it */
DLLEXPORT const char * ovf_latest_message(struct ovf_file *);

/* close the file and clean up resources */
DLLEXPORT int ovf_close(struct ovf_file *);

#undef DLLEXPORT
#endif

