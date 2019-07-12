#include "ovf.h"
#include <detail/Filter_File_Handle.hpp>
#include <detail/Helpers.hpp>
#include <fmt/format.h>


struct ovf_file* ovf_open(const char *filename)
try
{
    // Initialize the struct
    struct ovf_file * ovf_file_ptr = new ovf_file{ strdup(filename), /*found*/ false, /*ovf*/ false, /*n_segments*/ 0, /*handle*/ nullptr };
    ovf_file_ptr->_file_handle = new ovf_file_handle{ "", "", {}};

    // Check if the file exists
    std::fstream filestream( filename );
    ovf_file_ptr->found = filestream.is_open();
    filestream.close();

    if( ovf_file_ptr->found )
    {
        auto ifile = std::unique_ptr<Filter_File_Handle>(new Filter_File_Handle( filename, comment_tag ));

        // Check if the file has an OVF top header and check the OVF version
        std::string version;
        if ( ifile->Read_Single( version, "# OOMMF OVF", false ) )
        {
            if( version == "2.0" || version == "2" )
                ovf_file_ptr->is_ovf = true;
            else
                ovf_file_ptr->_file_handle->message_latest = fmt::format(
                    "libovf ovf_open: OVF version \'{}\' in file \'{}\' is not supported...",
                    filename, version);
        }
        
        if( ovf_file_ptr->is_ovf )
        {
            // Get the number of segments written in the header
            ifile->Require_Single( ovf_file_ptr->n_segments, "# segment count:" );
            int n_located = count_and_locate_segments(ovf_file_ptr);
            if( ovf_file_ptr->n_segments != n_located )
                ovf_file_ptr->_file_handle->message_latest = fmt::format(
                    "libovf ovf_open: n_segments specified in header ({}) is different from the number"
                    " of segments ({}) found in the file \'{}\'...",
                    ovf_file_ptr->n_segments, n_located, filename);
        }
    }

    return ovf_file_ptr;
}
catch( ... )
{
    return nullptr;
}


struct ovf_segment * ovf_segment_initialize()
try
{
    struct ovf_segment * ovf_segment_ptr = new ovf_segment
    {
        const_cast<char *>(""),   // title
        const_cast<char *>(""),   // comment
        0,                        // valuedim
        const_cast<char *>(""),   // valueunits
        const_cast<char *>(""),   // valuelabels
        const_cast<char *>(""),   // meshtype
        const_cast<char *>(""),   // meshunits
        0,                        // pointcount
        {0,0,0},                  // n_cells
        0,                        // N
        {0,0,0},                  // bounds_min
        {0,0,0},                  // bounds_max
        0,                        // lattice_constant
        {{1,0,0},{0,1,0},{0,0,1}} // bravais_vectors
    };

    return ovf_segment_ptr;
}
catch( ... )
{
    return nullptr;
}


bool check_segment(const ovf_segment * segment)
try
{
    if( !segment->title )
        return false;

    if( !segment->comment )
        return false;

    if( !segment->title )
        return false;
    
    if( !segment->bravais_vectors[0] || !segment->bravais_vectors[1] || !segment->bravais_vectors[2] )
        return false;

    return true;
}
catch( ... )
{
    return false;
}

int ovf_read_segment_header(struct ovf_file * ovf_file_ptr, int index, struct ovf_segment *segment)
try
{
    if( !ovf_file_ptr )
        return OVF_ERROR;

    if( !segment )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_read_segment_header: invalid segment pointer";
        return OVF_ERROR;
    }

    if( !ovf_file_ptr->found )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_header: file \'{}\' does not exist...",
            ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    if( !ovf_file_ptr->is_ovf )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_header: file \'{}\' is not ovf...",
            ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    if( index < 0 )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_header: invalid index ({}) < 0...",
            index, ovf_file_ptr->n_segments, ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    if( index >= ovf_file_ptr->n_segments )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_header: index ({}) >= n_segments ({}) of file \'{}\'...",
            index, ovf_file_ptr->n_segments, ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    int retcode = read_segment_header( ovf_file_ptr, index, segment );
    if (retcode != OVF_OK)
        ovf_file_ptr->_file_handle->message_latest += "\novf_read_segment_header failed.";
    return retcode;
}
catch( ... )
{
    return OVF_ERROR;
}


int ovf_read_segment_data_4(struct ovf_file *ovf_file_ptr, int index, const struct ovf_segment *segment, float *data)
try
{
    if( !ovf_file_ptr )
        return OVF_ERROR;

    if( !segment )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_read_segment_data_4: invalid segment pointer";
        return OVF_ERROR;
    }

    if( !check_segment(segment) )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_read_segment_data_4: segment not correctly initialized";
        return OVF_ERROR;
    }

    if( !data )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_read_segment_data_4: invalid data pointer";
        return OVF_ERROR;
    }

    if( !ovf_file_ptr->found )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_data_4: file \'{}\' does not exist...",
            ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    if( !ovf_file_ptr->is_ovf )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_data_4: file \'{}\' is not ovf...",
            ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    if( index >= ovf_file_ptr->n_segments )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_data_4: index ({}) >= n_segments ({}) of file \'{}\'...",
            index, ovf_file_ptr->n_segments, ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    // TODO: parse binary length correctly
    // if (file.binary_length != 4)
    // {
    //     ovf_file_ptr->_file_handle->message_latest = "libovf ovf_read_segment_data_4: file does not have binary length 4...";
    //     return OVF_ERROR;
    // }

    int retcode = read_segment(ovf_file_ptr, segment, index, data);
    if( retcode != OVF_OK )
        ovf_file_ptr->_file_handle->message_latest += "\novf_read_segment_data_4 failed.";
    return retcode;
}
catch( ... )
{
    return OVF_ERROR;
}


int ovf_read_segment_data_8(struct ovf_file *ovf_file_ptr, int index, const struct ovf_segment *segment, double *data)
try
{
    if( !ovf_file_ptr )
        return OVF_ERROR;

    if( !segment )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_read_segment_data_8: invalid segment pointer";
        return OVF_ERROR;
    }

    if( !check_segment(segment) )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_read_segment_data_8: segment not correctly initialized";
        return OVF_ERROR;
    }

    if( !data )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_read_segment_data_8: invalid data pointer";
        return OVF_ERROR;
    }

    if( !ovf_file_ptr->found )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_8: file \'{}\' does not exist...",
            ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    if( !ovf_file_ptr->is_ovf )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_8: file \'{}\' is not ovf...",
            ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    if( index >= ovf_file_ptr->n_segments )
    {
        ovf_file_ptr->_file_handle->message_latest = fmt::format(
            "libovf ovf_read_segment_8: index ({}) >= n_segments ({}) of file \'{}\'...",
            index, ovf_file_ptr->n_segments, ovf_file_ptr->filename);
        return OVF_ERROR;
    }

    // TODO: parse binary length correctly
    // if (file.binary_length != 8)
    // {
    //     ovf_file_ptr->_file_handle->message_latest = "libovf ovf_read_segment_8: file does not have binary length 4...";
    //     return OVF_ERROR;
    // }

    int retcode = read_segment(ovf_file_ptr, segment, index, data);
    if (retcode != OVF_OK)
        ovf_file_ptr->_file_handle->message_latest += "\novf_read_segment_data_8 failed.";
    return retcode;
}
catch( ... )
{
    return OVF_ERROR;
}


int ovf_write_segment_4(struct ovf_file *ovf_file_ptr, const struct ovf_segment *segment, float *data, int format)
try
{
    if( !ovf_file_ptr )
        return OVF_ERROR;

    if( !segment )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_write_segment_4: invalid segment pointer";
        return OVF_ERROR;
    }

    if( !check_segment(segment) )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_write_segment_4: segment not correctly initialized";
        return OVF_ERROR;
    }

    if( !data )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_write_segment_4: invalid data pointer";
        return OVF_ERROR;
    }

    if( format == OVF_FORMAT_BIN8 ||
        format == OVF_FORMAT_BIN4 )
        format = OVF_FORMAT_BIN;

    if( format != OVF_FORMAT_BIN  &&
        format != OVF_FORMAT_TEXT &&
        format != OVF_FORMAT_CSV  )
    {
        ovf_file_ptr->_file_handle->message_latest =
            fmt::format("libovf ovf_write_segment_4: invalid format \'{}\'...", format);
        return OVF_ERROR;
    }

    int retcode = write_segment(ovf_file_ptr, segment, data, true, false, format);
    if (retcode != OVF_OK)
        ovf_file_ptr->_file_handle->message_latest += "\novf_write_segment_4 failed.";
    return retcode;
}
catch( ... )
{
    return OVF_ERROR;
}


int ovf_write_segment_8(struct ovf_file *ovf_file_ptr, const struct ovf_segment *segment, double *data, int format)
try
{
    if( !ovf_file_ptr )
        return OVF_ERROR;

    if( !segment )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_write_segment_8: invalid segment pointer";
        return OVF_ERROR;
    }

    if( !check_segment(segment) )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_write_segment_8: segment not correctly initialized";
        return OVF_ERROR;
    }

    if( !data )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_write_segment_8: invalid data pointer";
        return OVF_ERROR;
    }

    if( format == OVF_FORMAT_BIN8 ||
        format == OVF_FORMAT_BIN4 )
        format = OVF_FORMAT_BIN;

    if( format != OVF_FORMAT_BIN  &&
        format != OVF_FORMAT_TEXT &&
        format != OVF_FORMAT_CSV  )
    {
        ovf_file_ptr->_file_handle->message_latest =
            fmt::format("libovf ovf_write_segment_8: invalid format \'{}\'...", format);
        return OVF_ERROR;
    }

    int retcode = write_segment(ovf_file_ptr, segment, data, true, false, format);
    if (retcode != OVF_OK)
        ovf_file_ptr->_file_handle->message_latest += "\novf_write_segment_8 failed.";
    return retcode;
}
catch( ... )
{
    return OVF_ERROR;
}


int ovf_append_segment_4(struct ovf_file *ovf_file_ptr, const struct ovf_segment *segment, float *data, int format)
try
{
    if( !ovf_file_ptr )
        return OVF_ERROR;

    if( !segment )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_append_segment_4: invalid segment pointer";
        return OVF_ERROR;
    }

    if( !check_segment(segment) )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_append_segment_4: segment not correctly initialized";
        return OVF_ERROR;
    }

    if( !data )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_append_segment_4: invalid data pointer";
        return OVF_ERROR;
    }

    if( !ovf_file_ptr->is_ovf )
    {
        ovf_file_ptr->_file_handle->message_latest = "libovf ovf_append_segment_4: file is not ovf...";
        return OVF_ERROR;
    }

    if( format == OVF_FORMAT_BIN8 ||
        format == OVF_FORMAT_BIN4 )
        format = OVF_FORMAT_BIN;

    if( format != OVF_FORMAT_BIN  &&
        format != OVF_FORMAT_TEXT &&
        format != OVF_FORMAT_CSV  )
    {
        ovf_file_ptr->_file_handle->message_latest =
            fmt::format("libovf ovf_append_segment_4: invalid format \'{}\'...", format);
        return OVF_ERROR;
    }

    bool write_header = !ovf_file_ptr->found;
    int retcode = write_segment(ovf_file_ptr, segment, data, write_header, true, format);
    if (retcode != OVF_OK)
        ovf_file_ptr->_file_handle->message_latest += "\novf_append_segment_4 failed.";
    return retcode;
}
catch( ... )
{
    return OVF_ERROR;
}


int ovf_append_segment_8(struct ovf_file *ovf_file_ptr, const struct ovf_segment *segment, double *data, int format)
try
{
    if( !ovf_file_ptr )
        return OVF_ERROR;

    if( !segment )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_append_segment_8: invalid segment pointer";
        return OVF_ERROR;
    }

    if( !check_segment(segment) )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_append_segment_8: segment not correctly initialized";
        return OVF_ERROR;
    }

    if( !data )
    {
        ovf_file_ptr->_file_handle->message_latest =
            "libovf ovf_append_segment_8: invalid data pointer";
        return OVF_ERROR;
    }

    if( format == OVF_FORMAT_BIN8 ||
        format == OVF_FORMAT_BIN4 )
        format = OVF_FORMAT_BIN;

    if( !ovf_file_ptr->is_ovf )
    {
        ovf_file_ptr->_file_handle->message_latest = "libovf ovf_append_segment_8: file is not ovf...";
        return OVF_ERROR;
    }

    if( format != OVF_FORMAT_BIN  &&
        format != OVF_FORMAT_TEXT &&
        format != OVF_FORMAT_CSV  )
    {
        ovf_file_ptr->_file_handle->message_latest =
            fmt::format("libovf ovf_append_segment_8: invalid format \'{}\'...", format);
        return OVF_ERROR;
    }

    bool write_header = !ovf_file_ptr->found;
    int retcode = write_segment(ovf_file_ptr, segment, data, write_header, true, format);
    if (retcode != OVF_OK)
        ovf_file_ptr->_file_handle->message_latest += "\novf_append_segment_8 failed.";
    return retcode;
}
catch( ... )
{
    return OVF_ERROR;
}


const char * ovf_latest_message(struct ovf_file *ovf_file_ptr)
try
{
    if( !ovf_file_ptr )
        return "";

    ovf_file_ptr->_file_handle->message_out = ovf_file_ptr->_file_handle->message_latest;
    ovf_file_ptr->_file_handle->message_latest = "";
    return ovf_file_ptr->_file_handle->message_out.c_str();
}
catch( ... )
{
    return "";
}


int ovf_close(struct ovf_file *ovf_file_ptr)
try
{
    if( !ovf_file_ptr )
        return OVF_ERROR;
    delete(ovf_file_ptr->_file_handle);
    delete(ovf_file_ptr);
    return OVF_OK;
}
catch( ... )
{
    return OVF_ERROR;
}

