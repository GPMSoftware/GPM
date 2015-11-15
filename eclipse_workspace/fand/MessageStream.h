/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MessageStream.h
 *
 *  Created on: Mar 23, 2010
 *      Author: gareth
 *      References: http://wordaligned.org/articles/cpp-streambufs
 *
 *  A MessageStream allows output to be formatted and sent to an output stream and/or file.
 *  It may be used to implement info/warn/error debug messages, e.g.
 *
 *  MessageStream info("info", "INFO: ", cout, "info_log.txt");
 *  info << "hello world" << endl;
 *
 *  A MessageStream object with each name should be instantiated separately in each source file
 *  (not globally). This may be done my means of a header file included in each source file
 *
 *  messages.h:
 *  static MessageStream info ("info", "INFO:    ", std::cout, "",              false); // to cout, not logged, not enabled by default
 *  static MessageStream warn ("warn", "WARNING: ", std::cout);                         // to cout, not logged, enabled by default
 *  static MessageStream error("error","ERROR:   ", std::cerr, "error_log.txt", true);  // to cerr, logged, enabled by default
 *  static MessageStream stats("stats","STATS:   ",            "stats_log.txt");        //  logged only, enabled by default
 *
 *  The arguments to the constructor are:
 *  string  name         : name of stream (used in environment variables)
 *  string  prefix       : output line prefix
 *  ostream os           : stream to echo to (optional)
 *  string  filename     : log file name (optional, default "" meaning no logging)
 *  bool    out_default  : default output flag (optional: default true)
 *
 *  This allows ouput to be turned on/off per source file. Control of output by source file is by environment
 *  variables:
 *
 *  info_log_file=infolog.txt  # set log file for info (globally)
 *  info_default_out=true      # turn on info globally (false = off)
 *  source_info=true           # turn on info in source.cpp (false = off). If set, overrides the global setting)
 *
 *  The startl macro causes a prefix containing file name and line number to be printed, e.g.
 *
 *  info << startl << "hello world" << endl;
 *
 *  This will produce output like:
 *
 *  INFO: at MessageStream.cpp:55: hello world
 *
 *  A MessageStream may be redirected to a buffer (and not send to its usual destination)
 *  by setting an environment variable like:
 *
 *  error_redirect=true
 *
 *  All redirected MessageStream go to the same buffer, which may be read (and cleared) with read().
 *  This is used to redirect output to the GUI.
 */

#ifndef MESSAGESTREAM_H_
#define MESSAGESTREAM_H_

#include <streambuf>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>

#define INIT_MESSAGES( sourceid ) static const char * thisid(sourceid), *thisfile(__FILE__);
#define startl MessageStream::Line(__LINE__)
#define alignl MessageStream::Line(-1)

class MessageStream
{
public:
	struct Line
	{
		Line(int line) : m_line(line) {}
		int m_line;
	};

	MessageStream(
		std::string sourcefileid,
		std::string sourcefilename,
		std::string name,
		std::string prefix,
		const std::string logfilename = "",
		bool out_default = true,
		bool redirect_default = false);

    MessageStream(
   		std::string sourcefileid,
    	std::string sourcefilename,
    	std::string name,
    	std::string prefix,
    	std::ostream &os,
    	const std::string logfilename = "",
    	bool out_default = true,
    	bool redirect_default = false);

    ~MessageStream();

    MessageStream& operator<<(std::ostream& (*manip)(std::ostream& ));

    MessageStream& operator<<(const Line &line_info);

    template<typename T>
    MessageStream& operator<<(const T& t)
    {
    	if (enabled_f() && m_output)
    	{
            if (m_op1) *m_op1 << t;
            if (m_op2) *m_op2 << t;
    	}
        return *this;
    }

    void flush();

    static bool
    enable(bool enabled) { bool ret = enabled_f(); enabled_f() = enabled; return ret; }

    static std::string read(std::string const & new_val = ""); // read and set (default: clear) the buffer of redirected MessageStreams

private:
    void setup();
    void startLine(std::ostream &os, const Line &line_info);

    std::string m_source_file_id;
    std::string m_source_file_name;
    std::string m_name;
    std::string m_prefix;
    std::ostream *m_op1;
    std::ostream *m_op2;
    std::string m_log_file_name;
    bool m_output;
    bool m_redirect;

//    std::ofstream m_of;
    std::string m_indent;

    static bool& enabled_f() { static bool theBoolean = true; return theBoolean; }

    // buffer into which output may be redirected, and read with read()
    static std::ostringstream& redirect_buffer_f() { static std::ostringstream theStream; return theStream; }

    // opens/returns file streams to log files (so each is opened only once,
    // and may be written to by multiple MessageStreams)
    static std::ofstream* getFileStream(std::string const &filename);
};

#endif /* MESSAGESTREAM_H_ */
