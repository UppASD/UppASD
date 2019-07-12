#pragma once
#ifndef IO_FILTERFILEHANDLE_H
#define IO_FILTERFILEHANDLE_H

#include <memory>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <fmt/format.h>
#include <fmt/ostream.h>

class Filter_File_Handle
{
private:
    std::size_t found;
    std::string line;
    const std::string comment_tag;
    std::string dump;
    // Beggining and end of file stream indicator 
    std::ios::pos_type position_file_beg;
    std::ios::pos_type position_file_end;
    // Start and stop of file stream indicator
    std::ios::pos_type position_start;
    std::ios::pos_type position_stop;
    int n_lines;
    int n_comment_lines;
    std::ifstream myfile;

public:
    std::string filename;
    std::istringstream iss;
    
    // Constructs a Filter_File_Handle with string filename
    Filter_File_Handle( const std::string& filename, const std::string comment_tag = "#" );
    // Destructor
    ~Filter_File_Handle();
    
    // Get the position of the file stream indicator
    std::ios::pos_type GetPosition( std::ios::seekdir dir = std::ios::cur );
    // Set limits in the file stream indicator
    void SetLimits( const std::ios::pos_type beg, const std::ios::pos_type end );
    // Reset the limits of the file stream indicator 
    void ResetLimits();
    // Reads next line of file into the handle (false -> end-of-file)
    bool GetLine_Handle( const std::string str_to_remove = "" );
    // Reads the next line of file into the handle and into the iss
    bool GetLine( const std::string str_to_remove = "" );
    // Reset the file stream to the start of the file
    void ResetStream();
    // Tries to find s in the current file and if found outputs the line into internal iss
    bool Find(const std::string& s);
    // Tries to find s in the current line and if found outputs into internal iss
    bool Find_in_Line(const std::string & s);
    // Removes a set of chars from a string
    void Remove_Chars_From_String(std::string &str, const char* charsToRemove);
    // Removes comments from a string
    bool Remove_Comments_From_String( std::string &str );
    // Read a string (separeated by whitespaces) into var. Capitalization is ignored.
    void Read_String( std::string& var, std::string keyword, bool log_notfound = true );
    // Count the words of a string
    int Count_Words( const std::string& str );
    // Returns the number of lines which are not starting with a comment
    int Get_N_Non_Comment_Lines();

    // Stream a single variable
    template <typename T>
    void read( T * var, std::size_t n );

    // Read a single variable
    template <typename T>
    bool Read_Single( T & var, std::string name, bool log_notfound = true );

    // Require a single variable to be found
    template <typename T>
    void Require_Single( T& var, std::string name );

    // Read a 3-vector
    template <typename T>
    void Read_3Vector( T * var, std::string name, bool log_notfound = true );
};//end class FilterFileHandle

// -----------------------------------------------------------------------------

template <typename T>
void Filter_File_Handle::read( T * var, std::size_t n )
{
    myfile.read( var, n );
};

// Reads a single variable into var, with optional logging in case of failure.
//
//// NOTE: Capitalization is ignored (expected).
//
template <typename T>
bool Filter_File_Handle::Read_Single( T & var, std::string name, bool log_notfound )
{
    try
    {
        std::transform( name.begin(), name.end(), name.begin(), ::tolower );

        if (Find(name))
        {
            iss >> var;
            return true;
        }
        // else if (log_notfound)
        //     Log( Utility::Log_Level::Warning, Utility::Log_Sender::IO, "Keyword '" + name + 
        //         "' not found. Using Default: " + fmt::format( "{}", var ) );
    }
    catch (...)
    {
        // spirit_handle_exception_core(fmt::format("Failed to read single variable \"{}\".", name));
    }
    return false;
};

// Require a single field. In case that it is not found an execption is thrown. 
//
//// NOTE: Capitalization is ignored (expected).
//
template <typename T>
void Filter_File_Handle::Require_Single( T& var, std::string name )
{
    std::transform( name.begin(), name.end(), name.begin(), ::tolower );
    
    if( !Read_Single( var, name, false ) )
    {
        // spirit_throw(Utility::Exception_Classifier::Bad_File_Content, Utility::Log_Level::Error,
        //     fmt::format("Required keyword \"{}\" not found.", name));
    }
}

// Reads a 3-component object, with optional logging in case of failure
template <typename T>
void Filter_File_Handle::Read_3Vector( T * var, std::string name, bool log_notfound )
{
    try
    {
        std::transform( name.begin(), name.end(), name.begin(), ::tolower );
        
        if (Find(name))
            iss >> var[0] >> var[1] >> var[2];
        // else if (log_notfound)
        //     Log( Utility::Log_Level::Warning, Utility::Log_Sender::IO, "Keyword '" + name + 
        //         "' not found. Using Default: (" + fmt::format( "{}", var[0] ) + " " + 
        //         fmt::format( "{}", var[1] ) + " " + fmt::format( "{}", var[2] ) + ")" );

    }
    catch (...)
    {
        // spirit_handle_exception_core(fmt::format("Failed to read 3Vector \"{}\".", name));
    }
}

inline Filter_File_Handle::Filter_File_Handle(
        const std::string& filename, const std::string comment_tag ) :
    filename(filename), comment_tag(comment_tag), iss("")
{
    this->dump = "";
    this->line = "";
    this->found = std::string::npos;
    this->myfile = std::ifstream( filename, std::ios::in | std::ios::binary );
    
    // find begging and end positions of the file stream indicator
    this->position_file_beg = this->myfile.tellg();
    this->myfile.seekg( 0, std::ios::end );
    this->position_file_end = this->myfile.tellg();
    this->myfile.seekg( 0, std::ios::beg );
    
    // set limits of the file stream indicator to begging and end positions (eq. ResetLimits())
    this->position_start = this->position_file_beg;
    this->position_stop = this->position_file_end;
    
    // initialize number of lines
    this->n_lines = 0;
    this->n_comment_lines = 0;

    // if the file is not open
    // if ( !this->myfile.is_open() )
    // spirit_throw(Exception_Classifier::File_not_Found, Log_Level::Error, fmt::format("Could not open file \"{}\"", filename));
}

inline Filter_File_Handle::~Filter_File_Handle()
{ 
    myfile.close();
}

inline std::ios::pos_type Filter_File_Handle::GetPosition( std::ios::seekdir dir )
{
    this->myfile.seekg( 0, dir );
    return this->myfile.tellg();
}

inline void Filter_File_Handle::SetLimits( const std::ios::pos_type start, 
                                    const std::ios::pos_type stop )
{
    this->position_start = start;
    this->position_stop = stop;
}

inline void Filter_File_Handle::ResetLimits()
{
    this->position_start = this->position_file_beg;
    this->position_stop = this->position_file_end;
}

inline bool Filter_File_Handle::GetLine_Handle( const std::string str_to_remove )
{
    this->line = "";
    
    //	if there is a next line
    if ( (bool) getline( this->myfile, this->line ) )
    {
        this->n_lines++;

        //  remove separator characters
        Remove_Chars_From_String( this->line, (char *) "|+" );
        
        // remove any unwanted str from the line eg. delimiters
        if ( str_to_remove != "" )
            Remove_Chars_From_String( this->line, str_to_remove.c_str() );
            
        // if the string does not start with a comment identifier
        if ( Remove_Comments_From_String( this->line ) )
        {
            return true;
        } 
        else 
        {
            this->n_comment_lines++;
            return GetLine( str_to_remove );
        } 
    }
    return false;     // if there is no next line, return false
}

inline bool Filter_File_Handle::GetLine( const std::string str_to_remove )
{
    if (Filter_File_Handle::GetLine_Handle( str_to_remove ))
    {
        // decapitalize line
        std::transform( this->line.begin(), this->line.end(), this->line.begin(), ::tolower );
        
        return Filter_File_Handle::Find_in_Line("");
    }
    return false;
}

inline void Filter_File_Handle::ResetStream()
{
    myfile.clear();
    myfile.seekg(0, std::ios::beg);
}

inline bool Filter_File_Handle::Find(const std::string & s)
{
    myfile.clear();
    //myfile.seekg( this->position_file_beg, std::ios::beg);
    myfile.seekg( this->position_start );

    while ( GetLine() && ( GetPosition() <= this->position_stop ) ) 
    {
        if (Find_in_Line(s) ) return true;
    }
    return false;
}

inline bool Filter_File_Handle::Find_in_Line( const std::string & s )
{
    // if s is found in line
    if ( !line.compare( 0, s.size(), s ) )
    {
        iss.clear();    // empty the stream
        iss.str(line);  // copy line into the iss stream
        
        // if s is not empty
        if ( s.compare("") )
        {
            int words = Count_Words( s );
            for (int i = 0; i < words; i++)
                iss >> dump;
        }
        
        return true;
    }
    return false;
}

inline void Filter_File_Handle::Remove_Chars_From_String(std::string &str, const char* charsToRemove)
{
    for (unsigned int i = 0; i < strlen(charsToRemove); ++i)
    {
        str.erase(std::remove(str.begin(), str.end(), charsToRemove[i]), str.end());
    }
}

inline bool Filter_File_Handle::Remove_Comments_From_String( std::string &str )
{
    std::string::size_type start = this->line.find( this->comment_tag );
    
    // if the line starts with a comment return false
    if ( start == 0 ) return false;
    
    // if the line has a comment somewhere remove it by trimming
    if ( start != std::string::npos )
        line.erase( this->line.begin() + start , this->line.end() );
    
    // return true
    return true;
}

inline void Filter_File_Handle::Read_String( std::string& var, std::string keyword, bool log_notfound )
{
    std::transform( keyword.begin(), keyword.end(), keyword.begin(), ::tolower );
    
    if ( Find( keyword ) )
    {
        getline( this->iss, var );
        
        // trim leading and trailing whitespaces
        size_t start = var.find_first_not_of(" \t\n\r\f\v");
        size_t end = var.find_last_not_of(" \t\n\r\f\v");
        if ( start != std::string::npos )
            var = var.substr( start, ( end - start + 1 ) );
    }
    // else if ( log_notfound )
    //     Log( Utility::Log_Level::Warning, Utility::Log_Sender::IO,
    //             fmt::format( "Keyword \"{}\" not found. Using Default: \"{}\"", keyword, var ) );
}

inline int Filter_File_Handle::Count_Words( const std::string& phrase )
{
    std::istringstream phrase_stream( phrase );
    this->dump = "";
    int words = 0;
    while( phrase_stream >> dump ) ++words;
    return words;
}

inline int Filter_File_Handle::Get_N_Non_Comment_Lines()
{
    while( GetLine() ) { };
    ResetLimits(); 
    return ( this->n_lines - this->n_comment_lines );
}

#endif