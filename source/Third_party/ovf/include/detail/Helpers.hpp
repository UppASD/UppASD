#pragma once
#ifndef DETAIL_OVFFILE_H
#define DETAIL_OVFFILE_H

#include "ovf.h"
#include "Filter_File_Handle.hpp"

#include <string>
#include <fstream>
#include <cctype>
#include <array>
#include <vector>

#include <fmt/format.h>
#include <fmt/ostream.h>


struct ovf_file_handle
{
    /* messages, e.g. in case a function returned OVF_ERROR.
        message_out will be filled and returned by ovf_latest_message, while message_latest
        will be filled by other functions and cleared by ovf_latest_message. */
    std::string message_out, message_latest;

    std::vector<std::ios::pos_type> segment_fpos;
};

// Test values for 4bit and 8bit binary data
static const uint32_t test_hex_4b = 0x4996B438;
static const uint64_t test_hex_8b = 0x42DC12218377DE40;

// Comment tag in OVF file header
static const std::string comment_tag = "##";

// ?
static const std::string empty_line = "#\n";

// The number of zero-padding for segment count.
// This is needed so that, when appending, the file does not need to be overwritten
static const int n_segments_str_digits = 6; // can store 1M modes


// Count the number of segments in the file. It also saves their file positions
inline int count_and_locate_segments(ovf_file * file)
try
{
    file->_file_handle->segment_fpos = std::vector<std::ios::pos_type>(0);
    auto ifile = std::unique_ptr<Filter_File_Handle>(new Filter_File_Handle( file->filename, comment_tag ));

    // get the number of segments from the occurrences of "# Begin: Segment"
    int n_begin_segment = 0;

    std::ios::pos_type end = ifile->GetPosition( std::ios::end );

    // NOTE: the keyword to find must be lower case since the Filter File Handle
    // converts the content of the input file to lower case automatically
    while( ifile->Find( "# begin: segment" ) )
    {
        std::ios::pos_type pos = ifile->GetPosition();
        file->_file_handle->segment_fpos.push_back( pos );
        ifile->SetLimits( pos, end );

        ++n_begin_segment;
    }

    // find the very last keyword of the file
    file->_file_handle->segment_fpos.push_back( end );

    // reset limits
    ifile->ResetLimits();

    // close the file
    ifile = NULL;

    return n_begin_segment;
}
catch( ... )
{
    file->_file_handle->message_latest = fmt::format("count_and_locate_segments failed for file \"{}\".", file->filename);
    return 0;
}


// Read segment's header into member variables
inline int read_segment_header( ovf_file * file, int idx_seg, ovf_segment * segment )
try
{
    // open the file
    auto ifile = std::unique_ptr<Filter_File_Handle>(new Filter_File_Handle( file->filename, comment_tag ));

    // Set limits within which to read
    ifile->SetLimits( file->_file_handle->segment_fpos[idx_seg], file->_file_handle->segment_fpos[idx_seg+1] );

    // Read the header
    std::string title = "";
    ifile->Read_String( title, "# Title:" );
    std::string comment = "";
    ifile->Read_String( comment, "# Desc:" );
    std::string meshunits = "";
    ifile->Read_Single( meshunits, "# meshunit:" );

    ifile->Require_Single( segment->valuedim, "# valuedim:" );

    std::string valueunits = "";
    ifile->Read_String( valueunits, "# valueunits:" );
    std::string valuelabels = "";
    ifile->Read_String( valuelabels, "# valuelabels:" );

    ifile->Read_Single( segment->bounds_min[0], "# xmin:" );
    ifile->Read_Single( segment->bounds_min[1], "# ymin:" );
    ifile->Read_Single( segment->bounds_min[2], "# zmin:" );

    ifile->Read_Single( segment->bounds_max[0], "# xmax:" );
    ifile->Read_Single( segment->bounds_max[1], "# ymax:" );
    ifile->Read_Single( segment->bounds_max[2], "# zmax:" );

    std::string meshtype = "";
    ifile->Require_Single( meshtype, "# meshtype:" );

    ifile->Read_3Vector( segment->bravais_vectors[0], "# xbase:", true );
    ifile->Read_3Vector( segment->bravais_vectors[1], "# ybase:", true );
    ifile->Read_3Vector( segment->bravais_vectors[2], "# zbase:", true );

    // ifile->Require_Single( segment->stepsize[0], "# xstepsize:" );
    // ifile->Require_Single( segment->stepsize[1], "# ystepsize:" );
    // ifile->Require_Single( segment->stepsize[2], "# zstepsize:" );
    segment->lattice_constant = 1;

    // Check mesh type
    segment->pointcount = 1;
    if( meshtype == "irregular" )
        ifile->Require_Single( segment->pointcount, "# pointcount:" );

    ifile->Require_Single( segment->n_cells[0], "# xnodes:" );
    ifile->Require_Single( segment->n_cells[1], "# ynodes:" );
    ifile->Require_Single( segment->n_cells[2], "# znodes:" );

    segment->N = segment->n_cells[0] * segment->n_cells[1] * segment->n_cells[2] * segment->pointcount;


    // Convert strings to char *
    segment->title = new char[title.length() + 1];
    strcpy(segment->title, title.c_str());

    segment->comment = new char[comment.length() + 1];
    strcpy(segment->comment, comment.c_str());

    segment->valueunits = new char[valueunits.length() + 1];
    strcpy(segment->valueunits, valueunits.c_str());

    segment->valuelabels = new char[valuelabels.length() + 1];
    strcpy(segment->valuelabels, valuelabels.c_str());

    segment->meshtype = new char[meshtype.length() + 1];
    strcpy(segment->meshtype, meshtype.c_str());

    segment->meshunits = new char[meshunits.length() + 1];
    strcpy(segment->meshunits, meshunits.c_str());

    return OVF_OK;
}
catch( ... )
{
    file->_file_handle->message_latest = fmt::format("read_segment_header failed for file \"{}\".", file->filename);
    return OVF_ERROR;
}


// Read segment's data block
template <typename T>
int read_segment( const ovf_file * file, const ovf_segment * segment, int idx_seg, T * vf)
try
{
    // open the file
    auto ifile = std::shared_ptr<Filter_File_Handle>(
                        new Filter_File_Handle( file->filename, comment_tag ) );

    ifile->SetLimits( file->_file_handle->segment_fpos[idx_seg], file->_file_handle->segment_fpos[idx_seg+1] );

    // Raw data representation
    std::string datatype_tmp = "", datatype = "";
    ifile->Read_String( datatype_tmp, "# Begin: Data" );
    std::istringstream repr( datatype_tmp );
    repr >> datatype;
    int binary_length = 0;
    if( datatype == "binary" )
        repr >> binary_length;


    int i_datatype;
    if( datatype == "binary" && binary_length == 8 )
    {
        i_datatype = OVF_FORMAT_BIN8;
    }
    else if( datatype == "binary" && binary_length == 4 )
    {
        i_datatype = OVF_FORMAT_BIN4;
    }
    else if( datatype == "text" )
    {
        i_datatype = OVF_FORMAT_TEXT;
    }
    else if( datatype == "csv" )
    {
        i_datatype = OVF_FORMAT_CSV;
    }
    else
    {
        file->_file_handle->message_latest = fmt::format(
            "read_segment failed for file \"{}\". Invalid data format {}", file->filename, datatype );
        return OVF_ERROR;
    }

    int n_cols = segment->valuedim;
    int n_rows = segment->n_cells[0] * segment->n_cells[1] * segment->n_cells[2];

    // Check if the initial check value of the binary data is valid
    if( i_datatype == OVF_FORMAT_BIN8 )
    {
        const double ref_8b = *reinterpret_cast<const double *>( &test_hex_8b );
        double read_8byte = 0;

        // check the validity of the initial check value read with the reference one
        ifile->read( reinterpret_cast<char *>( &read_8byte ), sizeof(double) );
        if ( read_8byte != ref_8b )
        {
            file->_file_handle->message_latest = "OVF initial check value of binary data is inconsistent";
            return OVF_ERROR;
        }
    }
    else if( i_datatype == OVF_FORMAT_BIN4 )
    {
        const float ref_4b = *reinterpret_cast<const float *>( &test_hex_4b );
        float read_4byte = 0;

        // check the validity of the initial check value read with the reference one
        ifile->read( reinterpret_cast<char *>( &read_4byte ), sizeof(float) );
        if ( read_4byte != ref_4b )
        {
            file->_file_handle->message_latest = "OVF initial check value of binary data is inconsistent";
            return OVF_ERROR;
        }
    }

    // Check that we actually read in any data
    if( n_cols*n_rows <= 0 )
    {
        file->_file_handle->message_latest = fmt::format(
            "read_segment not reading in any data, because n_cols*n_rows={}*{}<=0 for file \"{}\". You may want to check the segment you passed in.",
            n_cols, n_rows, file->filename);
        return OVF_ERROR;
    }

    // Read the data
    if( i_datatype == OVF_FORMAT_BIN8 )
        return read_data_bin( ifile, n_cols, n_rows, 8,  vf );
    else if( i_datatype == OVF_FORMAT_BIN4 )
        return read_data_bin( ifile, n_cols, n_rows, 4,  vf );
    else if( i_datatype == OVF_FORMAT_TEXT )
        return read_data_txt( ifile, n_cols, n_rows, vf );
    else if( i_datatype == OVF_FORMAT_CSV )
        return read_data_txt( ifile, n_cols, n_rows, vf, "," );
    else
    {
        file->_file_handle->message_latest = fmt::format(
            "read_segment failed - invalid datatype \'{}\' for file \"{}\".", datatype, file->filename);
        return OVF_ERROR;
    }
}
catch( ... )
{
    file->_file_handle->message_latest = fmt::format("read_segment failed for file \"{}\".", file->filename);
    return OVF_ERROR;
}


template <typename T>
int read_data_bin( std::shared_ptr<Filter_File_Handle> ifile, int n_cols, int n_rows, int binary_length, T * vf )
try
{
    // Set the input stream indicator to the end of the line describing the data block
    ifile->iss.seekg( std::ios::end );

    // Comparison of datum size compared to scalar type
    if ( binary_length == 4 )
    {
        int vectorsize = 3 * sizeof(float);
        std::vector<float> buffer(n_cols);
        for( int row=0; row<n_rows; ++row )
        {
            ifile->read(reinterpret_cast<char *>(buffer.data()), vectorsize);

            for (int col=0; col<n_cols; ++col)
                vf[n_cols*row + col] = static_cast<T>(buffer[col]);
        }
    }
    else if (binary_length == 8)
    {
        int vectorsize = n_cols * sizeof(double);
        std::vector<double> buffer(n_cols);
        for (int row=0; row<n_rows; row++)
        {
            ifile->read(reinterpret_cast<char *>(buffer.data()), vectorsize);

            for (int col=0; col<n_cols; ++col)
                vf[n_cols*row + col] = static_cast<T>(buffer[col]);
        }
    }
    return OVF_OK;
}
catch (...)
{
    // file->_file_handle->message_latest = fmt::format("read_data_bin failed for file \"{}\".", file->filename);
    return OVF_ERROR;
}


template <typename T>
int read_data_txt( std::shared_ptr<Filter_File_Handle> ifile, int n_cols, int n_rows, T * vf, const std::string& delimiter = "" )
try
{
    for (int row=0; row<n_rows; row++)
    {
        ifile->GetLine( delimiter );

        for (int col=0; col<n_cols; ++col)
            ifile->iss >> vf[n_cols*row + col];
    }
    return OVF_OK;
}
catch (...)
{
    // file->_file_handle->message_latest = fmt::format("read_data_txt failed for file \"{}\".", file->filename);
    return OVF_ERROR;
}



// TODO: use Filter_File_Handle instead...
inline void Strings_to_File(const std::vector<std::string> text, const std::string name, int no=-1)
{
    std::ofstream myfile;
    myfile.open(name);
    if (myfile.is_open())
    {
        if (no < 0)
            no = text.size();
        // Log(Log_Level::Debug, Log_Sender::All, "Started writing " + name);
        for (int i = 0; i < no; ++i) {
            myfile << text[i];
        }
        myfile.close();
        // Log(Log_Level::Debug, Log_Sender::All, "Finished writing " + name);
    }
    else
    {
        // Log(Log_Level::Error, Log_Sender::All, "Could not open " + name + " to write to file");
    }
}

inline void Append_String_to_File(const std::string text, const std::string name)
{
    std::ofstream myfile;
    myfile.open(name, std::ofstream::out | std::ofstream::app);
    if (myfile.is_open())
    {
        // Log(Log_Level::Debug, Log_Sender::All, "Started writing " + name);
        myfile << text;
        myfile.close();
        // Log(Log_Level::Debug, Log_Sender::All, "Finished writing " + name);
    }
    else
    {
        // Log(Log_Level::Error, Log_Sender::All, "Could not open " + name + " to append to file");
    }
}

inline std::string top_header_string()
{
    std::string ret = "# OOMMF OVF 2.0\n";
    ret += empty_line;

    // create padding string
    std::string padding( n_segments_str_digits, '0' );
    // write padding plus n_segments
    ret += fmt::format( "# Segment count: {}\n", padding );

    return ret;
}

inline int increment_n_segments(ovf_file *file)
try
{
    std::fstream filestream( file->filename );

    // Update n_segments
    file->n_segments++;

    // Convert updated n_segment into padded string
    std::string new_n_str = std::to_string( file->n_segments );
    std::string::size_type new_n_len = new_n_str.length();

    std::string::size_type padding_len = n_segments_str_digits - new_n_len;
    std::string padding( padding_len, '0' );

    // n_segments_pos is the end of the line that contains '#segment count' (after '\n')
    std::ios::off_type offset = n_segments_str_digits + 1;

    auto ifile = std::unique_ptr<Filter_File_Handle>(
                        new Filter_File_Handle( file->filename, comment_tag ) );
    // Get the number of segment as string
    int n_segments;
    ifile->Require_Single( n_segments, "# segment count:" );
    std::string n_segments_as_str;
    ifile->Read_String( n_segments_as_str, "# segment count:" );
    // Save the file position indicator as we have to increment n_segment
    std::ios::pos_type n_segments_pos = ifile->GetPosition();


    // Go to the beginning '#segment count' value position
    filestream.seekg( n_segments_pos );
    filestream.seekg( (-1)*offset, std::ios::cur );

    // replace n_segments value in the stream
    filestream << ( padding + new_n_str );

    filestream.close();

    return OVF_OK;
}
catch( ... )
{
    file->_file_handle->message_latest = fmt::format("increment_n_segments failed for file \"{}\".", file->filename);
    return OVF_ERROR;
}

template <typename T>
void append_data_bin_to_string( std::string & output_to_file, const T * vf, int n_cols, int n_rows, int format)
try
{
    // float test value
    const float ref_4b = *reinterpret_cast<const float *>( &test_hex_4b );

    // double test value
    const double ref_8b = *reinterpret_cast<const double *>( &test_hex_8b );

    if( format == OVF_FORMAT_BIN8 )
    {
        output_to_file +=
            std::string( reinterpret_cast<const char *>(&ref_8b), sizeof(double) );

        // in case that T is 4bytes long
        if (sizeof(T) == sizeof(float))
        {
            std::vector<double> buffer(n_cols);
            for (unsigned int i=0; i<n_rows; ++i)
            {
                for (int j=0; j<n_cols; ++j)
                    buffer[j] = static_cast<double>(vf[n_cols*i + j]);
                output_to_file +=
                    std::string( reinterpret_cast<char *>(buffer.data()), n_cols*sizeof(double) );
            }
        }
        else
        {
            for (unsigned int i=0; i<n_rows; i++)
                output_to_file +=
                    std::string( reinterpret_cast<const char *>(&vf[i]), n_cols*sizeof(double) );
        }
    }
    else if( format == OVF_FORMAT_BIN4 )
    {
        output_to_file +=
            std::string( reinterpret_cast<const char *>(&ref_4b), sizeof(float) );

        // in case that T is 8bytes long
        if (sizeof(T) == sizeof(double))
        {
            std::vector<float> buffer(n_cols);
            for (unsigned int i=0; i<n_rows; ++i)
            {
                for (int j=0; j<n_cols; ++j)
                    buffer[j] = static_cast<float>(vf[n_cols*i + j]);
                output_to_file +=
                    std::string( reinterpret_cast<char *>(buffer.data()), n_cols*sizeof(float) );
            }
        }
        else
        {
            for (unsigned int i=0; i<n_rows; i++)
                output_to_file +=
                    std::string( reinterpret_cast<const char *>(&vf[i]), n_cols*sizeof(float) );
        }
    }
}
catch( ... )
{

}


template <typename T>
void append_data_txt_to_string( std::string & output_to_file, const T * vf, int n_cols, int n_rows, const std::string& delimiter = "" )
try
{
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
            output_to_file += fmt::format( "{:22.12f}{}", vf[n_cols*row + col], delimiter );
        output_to_file += "\n";
    }
}
catch( ... )
{

}


template <typename T>
int write_segment( ovf_file *file, const ovf_segment * segment, const T * vf,
                    bool write_header, const bool append = false, int format = OVF_FORMAT_BIN8 )
try
{
    std::string output_to_file;
    output_to_file.reserve( int( 0x08000000 ) );  // reserve 128[MByte]

    // If we are not appending or the file does not exists we need to write the top header
    // and to turn the file_exists attribute to true so we can append more segments
    if ( !append || write_header )
    {
        output_to_file += top_header_string();
        file->n_segments = 0;
    }

    output_to_file += fmt::format( empty_line );
    output_to_file += fmt::format( "# Begin: Segment\n" );
    output_to_file += fmt::format( "# Begin: Header\n" );
    output_to_file += fmt::format( empty_line );

    output_to_file += fmt::format( "# Title: {}\n", segment->title );
    output_to_file += fmt::format( empty_line );

    output_to_file += fmt::format( "# Desc: {}\n", segment->comment );
    output_to_file += fmt::format( empty_line );

    // The value dimension is always 3 since we are writting Vector3-data
    output_to_file += fmt::format( "# valuedim: {}   ## field dimensionality\n", segment->valuedim );
    output_to_file += fmt::format( "# valueunits: None None None\n" );
    output_to_file += fmt::format( "# valuelabels: spin_x_component spin_y_component spin_z_component \n");
    output_to_file += fmt::format( empty_line );

    // TODO: this ovf library does not support mesh units yet
    output_to_file += fmt::format( "## Fundamental mesh measurement unit. "
                                    "Treated as a label:\n" );
    output_to_file += fmt::format( "# meshunit: unspecified\n" );
    output_to_file += fmt::format( empty_line );

    output_to_file += fmt::format( "# xmin: {}\n", segment->bounds_min[0] );
    output_to_file += fmt::format( "# ymin: {}\n", segment->bounds_min[1] );
    output_to_file += fmt::format( "# zmin: {}\n", segment->bounds_min[2] );
    output_to_file += fmt::format( "# xmax: {}\n", segment->bounds_max[0] );
    output_to_file += fmt::format( "# ymax: {}\n", segment->bounds_max[1] );
    output_to_file += fmt::format( "# zmax: {}\n", segment->bounds_max[2] );
    output_to_file += fmt::format( empty_line );

    // TODO: this ovf library does not support irregular geometry yet. Write ONLY rectangular mesh
    output_to_file += fmt::format( "# meshtype: rectangular\n" );

    // Bravais Lattice
    output_to_file += fmt::format( "# xbase: {} {} {}\n",
                                    segment->bravais_vectors[0][0],
                                    segment->bravais_vectors[0][1],
                                    segment->bravais_vectors[0][2] );
    output_to_file += fmt::format( "# ybase: {} {} {}\n",
                                    segment->bravais_vectors[1][0],
                                    segment->bravais_vectors[1][1],
                                    segment->bravais_vectors[1][2] );
    output_to_file += fmt::format( "# zbase: {} {} {}\n",
                                    segment->bravais_vectors[2][0],
                                    segment->bravais_vectors[2][1],
                                    segment->bravais_vectors[2][2] );

    // output_to_file += fmt::format( "# xstepsize: {}\n",
    //                             segment->lattice_constant * segment->bravais_vectors[0][0] );
    // output_to_file += fmt::format( "# ystepsize: {}\n",
    //                             segment->lattice_constant * segment->bravais_vectors[1][1] );
    // output_to_file += fmt::format( "# zstepsize: {}\n",
    //                             segment->lattice_constant * segment->bravais_vectors[2][2] );

    output_to_file += fmt::format( "# xnodes: {}\n", segment->n_cells[0] );
    output_to_file += fmt::format( "# ynodes: {}\n", segment->n_cells[1] );
    output_to_file += fmt::format( "# znodes: {}\n", segment->n_cells[2] );
    output_to_file += fmt::format( empty_line );

    output_to_file += fmt::format( "# End: Header\n" );
    output_to_file += fmt::format( empty_line );

    if( sizeof(T) == sizeof(float) && format == OVF_FORMAT_BIN )
        format = OVF_FORMAT_BIN4;
    else if( sizeof(T) == sizeof(double) && format == OVF_FORMAT_BIN )
        format = OVF_FORMAT_BIN8;

    std::string datatype_out = "";
    if ( format == OVF_FORMAT_BIN8 )
        datatype_out = "Binary 8";
    else if ( format == OVF_FORMAT_BIN4 )
        datatype_out = "Binary 4";
    else if( format == OVF_FORMAT_TEXT )
        datatype_out = "Text";
    else if( format == OVF_FORMAT_CSV )
        datatype_out = "CSV";

    // Data
    output_to_file += fmt::format( "# Begin: Data {}\n", datatype_out );

    int n_rows = segment->n_cells[0]*segment->n_cells[1]*segment->n_cells[2];
    int n_cols = segment->valuedim;

    // Check that we actually read in any data
    if( n_cols*n_rows <= 0 )
    {
        file->_file_handle->message_latest = fmt::format(
            "write_segment not writing out any data, because n_cols*n_rows={}*{}<=0 for file \"{}\". You may want to check the segment you passed in.",
            n_cols, n_rows, file->filename);
        return OVF_ERROR;
    }

    if ( format == OVF_FORMAT_BIN || format == OVF_FORMAT_BIN8 || format == OVF_FORMAT_BIN4 )
        append_data_bin_to_string( output_to_file, vf, n_cols, n_rows, format );
    else if ( format == OVF_FORMAT_TEXT )
        append_data_txt_to_string( output_to_file, vf, n_cols, n_rows );
    else if ( format == OVF_FORMAT_CSV )
        append_data_txt_to_string( output_to_file, vf, n_cols, n_rows, "," );
    // else
    //     // TODO...

    output_to_file += fmt::format( "# End: Data {}\n", datatype_out );
    output_to_file += fmt::format( "# End: Segment\n" );

    // Append the #End keywords
    if( append )
        Append_String_to_File( output_to_file, file->filename );
    else
        Strings_to_File({output_to_file}, file->filename);

    file->found  = true;
    file->is_ovf = true;

    // Create or append segment_fpos
    auto ifile = std::unique_ptr<Filter_File_Handle>(new Filter_File_Handle( file->filename, comment_tag ));
    std::ios::pos_type end = ifile->GetPosition( std::ios::end );
    if( append && file->n_segments > 0 )
    {
        ifile->SetLimits( file->_file_handle->segment_fpos[file->n_segments], end );
    }
    ifile->Find( "# begin: segment" );
    std::ios::pos_type pos = ifile->GetPosition();
    if( append && file->n_segments > 0 )
    {
        file->_file_handle->segment_fpos[file->n_segments] = pos;
    }
    else
    {
        file->_file_handle->segment_fpos = std::vector<std::ios::pos_type>(0);
        file->_file_handle->segment_fpos.push_back( pos );
    }
    file->_file_handle->segment_fpos.push_back( end );

    // Increment the n_segments after succesfully appending the segment body to the file
    return increment_n_segments(file);
}
catch( const std::exception & ex )
{
    file->_file_handle->message_latest = fmt::format("Caught std::exception \"{}\"", ex.what());
    return OVF_ERROR;
}
catch( ... )
{
    return OVF_ERROR;
}


// // Read a variable from the comment section from the header of segment idx_seg
// template <typename T>
// void Read_Variable_from_Comment( T& var, const std::string name, const int idx_seg = 0 )
// {
//     try
//     {
//         // NOTE: seg_idx.max = segment_fpos.size - 2
//         // if ( idx_seg >= ( this->segment_fpos.size() - 1 ) )
//         //     spirit_throw( Utility::Exception_Classifier::Input_parse_failed,
//         //                     Utility::Log_Level::Error,
//         //                     "OVF error while choosing segment to read - "
//         //                     "index out of bounds" );

//         this->ifile->SetLimits( this->segment_fpos[idx_seg],
//                                 this->segment_fpos[idx_seg+1] );

//         this->ifile->Read_Single( var, name );

//         // Log( Utility::Log_Level::Debug, this->sender, fmt::format( "{}{}", name, var ) );
//     }
//     catch (...)
//     {
//         // spirit_handle_exception_core(fmt::format( "Failed to read variable \"{}\" "
//         //                                             "from comment", name ));
//     }
// }

// // Read a Vector3 from the comment section from the header of segment idx_seg
// template <typename T>
// void Read_String_from_Comment( T& var, std::string name,
//                                 const int idx_seg = 0 )
// {
//     try
//     {
//         // NOTE: seg_idx.max = segment_fpos.size - 2
//         // if ( idx_seg >= ( this->segment_fpos.size() - 1 ) )
//         //     spirit_throw( Utility::Exception_Classifier::Input_parse_failed,
//         //                     Utility::Log_Level::Error,
//         //                     "OVF error while choosing segment to read - "
//         //                     "index out of bounds" );

//         this->ifile->SetLimits( this->segment_fpos[idx_seg],
//                                 this->segment_fpos[idx_seg+1] );

//         this->ifile->Read_String( var, name, false );

//         // Log( Utility::Log_Level::Debug, this->sender, fmt::format( "{}{}", name, var ) );
//     }
//     catch (...)
//     {
//         // spirit_handle_exception_core(fmt::format( "Failed to read string \"{}\" "
//         //                                             "from comment", name ));
//     }
// }



#endif