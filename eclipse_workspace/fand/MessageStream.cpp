/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 *
 */

#include "MessageStream.h"
INIT_MESSAGES("MessageStream")

#include <UnitTest++/UnitTest++.h>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

MessageStream::MessageStream(
	std::string sourcefileid,
	std::string sourcefilename,
	std::string stream_name,
	std::string prefix,
	const std::string filename,
	bool out_default,
	bool redirect_default)
: m_source_file_id(sourcefileid)
, m_source_file_name(sourcefilename)
, m_name(stream_name)
, m_prefix(prefix)
, m_op1(0)
, m_op2(0)
, m_log_file_name(filename)
, m_output(out_default)
, m_redirect(redirect_default)
{
	setup();
}

MessageStream::MessageStream(
	std::string sourcefileid,
	std::string sourcefilename,
	std::string stream_name,
	std::string prefix,
	std::ostream &os,
	const std::string filename,
	bool out_default,
	bool redirect_default)
: m_source_file_id(sourcefileid)
, m_source_file_name(sourcefilename)
, m_name(stream_name)
, m_prefix(prefix)
, m_op1(&os)
, m_op2(0)
, m_log_file_name(filename)
, m_output(out_default)
, m_redirect(redirect_default)
{
	setup();
}

MessageStream::~MessageStream()
{
//	m_of.close();
}

MessageStream&
MessageStream::operator<<(std::ostream& (*manip)(std::ostream& ))
{
	if (enabled_f() && m_output)
	{
        if (m_op1) *m_op1 << manip;
        if (m_op2) *m_op2 << manip;
	}
    return *this;
}

MessageStream&
MessageStream::operator<<(const Line &line_info)
{
	if (enabled_f() && m_output)
	{
		if (m_op1) startLine(*m_op1, line_info);
		if (m_op2) startLine(*m_op2, line_info);
	}
    return *this;
}

void
MessageStream::flush()
{
    if (enabled_f() && m_output)
    {
        if (m_op1) m_op1->flush();
        if (m_op2) m_op2->flush();
    }
}

void
MessageStream::setup()
{
	// Environment variables may override constructor arguments at run time

	int indent_chars = m_prefix.length() + 3 + m_source_file_name.length() + 1 + 4 + 2;
	m_indent = string(indent_chars, ' ');

	// get global output switch
	string global_output_env = m_name + "_" + "default_out";
	if (char *s = getenv(global_output_env.c_str()))
	{
		if (strcmp(s, "true") == 0)
		{
			m_output = true;
		}
		else if (strcmp(s, "false") == 0)
		{
			m_output = false;
		} // else unchanged
	}

	// get source file output switch
	string source_output_env = m_source_file_id + "_" + m_name;
	if (char *s = getenv(source_output_env.c_str()))
	{
		if (strcmp(s, "true") == 0)
		{
			m_output = true;
		}
		else if (strcmp(s, "false") == 0)
		{
			m_output = false;
		}  // else unchanged
	}

    // get redirect output switch
    string redirect_output_env = m_name + "_" + "redirect";
    if (char *s = getenv(redirect_output_env.c_str()))
    {
        if (strcmp(s, "true") == 0)
        {
            m_redirect = true;
        }
        else if (strcmp(s, "false") == 0)
        {
            m_redirect = false;
        } // else unchanged
    }

    // redirect
    if (m_redirect)
    {
        m_op1 = &redirect_buffer_f();
    }

	// get log file name
	string log_file_env = m_name + "_" + "log_file";
	if (char *s = getenv(log_file_env.c_str()))
	{
		m_log_file_name = s;
	}

	// open log file
	if (m_log_file_name != "")
	{
		m_op2 = getFileStream(m_log_file_name);
	}
}

std::ofstream*
MessageStream::getFileStream(string const &filename)
{
    static std::map<std::string, std::ofstream*> file_streams;

    std::map<std::string, std::ofstream*>::const_iterator it = file_streams.find(filename);
    if (it == file_streams.end())
    {
        file_streams[filename] = new ofstream(filename.c_str());
    }

    return file_streams[filename];
}

void
MessageStream::startLine(std::ostream &os, const Line &line_info)
{
	if (line_info.m_line < 0) // indent
	{
		os << m_indent;
	}
	else
	{
		os << m_prefix << "at " << m_source_file_name << ":" << std::setw(4) << line_info.m_line << ": ";
	}
}

std::string
MessageStream::read(string const & new_val)
{
    std::string ret = redirect_buffer_f().str();
    redirect_buffer_f().str(new_val); // set contents
    redirect_buffer_f().clear();      // clear flags
    return ret;
}

TEST(MessageStream1)
{
	// Messages must be enabled for the purpose of this test
	bool enabled = MessageStream::enable(true);

    std::ostringstream sout;

    {
        MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", sout, "info.txt" /*, true*/);
        info_test << startl << "stuff1" << std::endl; /** messages contain this line number **/
    }

    CHECK_EQUAL("INFO: at ../MessageStream.cpp: 215: stuff1\n", sout.str());
    std::ostringstream fs1;
    std::ifstream if1("info.txt");
    fs1 << if1.rdbuf();
    CHECK_EQUAL("INFO: at ../MessageStream.cpp: 215: stuff1\n", fs1.str());

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream2)
{
	// Messages must be enabled for the purpose of this test
	bool enabled = MessageStream::enable(true);

	std::ostringstream sout;
    MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", sout /*, "", true*/);
    info_test << "stuff2" << std::endl;
    CHECK_EQUAL("stuff2\n", sout.str());

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream2_1)
{
	// Messages must be enabled for the purpose of this test
	bool enabled = MessageStream::enable(true);

	std::ostringstream sout;
    MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", sout, "", false); // no output by default
    info_test << "stuff2_1" << std::endl;
    CHECK_EQUAL("", sout.str());

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream2_2)
{
	// Messages must be enabled for the purpose of this test
	bool enabled = MessageStream::enable(true);

    (void)setenv("info_test_default_out", "true", 1);
    std::ostringstream sout;
    MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", sout, "", false); // overriden with environment variable
    info_test << "stuff2_2" << std::endl;
    CHECK_EQUAL("stuff2_2\n", sout.str());
    (void)unsetenv("info_test_default_out");

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream2_2a)
{
	// Messages must be enabled for the purpose of this test
	bool enabled = MessageStream::enable(true);

    (void)setenv("info_test_default_out", "false", 1);
    (void)setenv("MessageStream_info_test", "true", 1);     // file overrides global
    std::ostringstream sout;
    MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", sout, "", false); // overriden with environment variable
    info_test << "stuff2_2" << std::endl;
    CHECK_EQUAL("stuff2_2\n", sout.str());
    (void)unsetenv("MessageStream_info_test");
    (void)unsetenv("info_test_default_out");

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream2_2b)
{
	// Messages must be enabled for the purpose of this test
	bool enabled = MessageStream::enable(true);

	(void)setenv("info_test_default_out", "whatever", 1);
    std::ostringstream sout;
    MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", sout, "", false); // NOT overriden with environment variable (not "true" or "false")
    info_test << "stuff2_2" << std::endl;
    CHECK_EQUAL("", sout.str());
    (void)unsetenv("info_test_default_out");

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream2_2c)
{
	// Messages must be enabled for the purpose of this test
	bool enabled = MessageStream::enable(true);

    (void)setenv("info_test_default_out", "", 1);
    std::ostringstream sout;
    MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", sout, "", true);  // NOT overriden with environment variable (not "true" or "false")
    info_test << "stuff2_2" << std::endl;
    CHECK_EQUAL("stuff2_2\n", sout.str());
    (void)unsetenv("info_test_default_out");

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream2_3)
{
	// Messages must be enabled for the purpose of this test
	bool enabled = MessageStream::enable(true);

    setenv("info_test_log_file", "info2.txt", 1); // set the filename by env variable
	std::ostringstream sout;
    {
		MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", sout);
		info_test << "stuff2_3" << std::endl;
    }
    CHECK_EQUAL("stuff2_3\n", sout.str());

    std::ostringstream fs1;
    std::ifstream if1("info2.txt");
    fs1 << if1.rdbuf();
    CHECK_EQUAL("stuff2_3\n", fs1.str());
    (void)unsetenv("info_test_log_file");

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream2_4)
{
    // Messages must be enabled for the purpose of this test
    bool enabled = MessageStream::enable(true);

    // a different MessageStream can write to the same file without corruption
    std::ostringstream sout;
    {
        MessageStream info_test2(thisid, thisfile, "info_test", "INFO: ", sout, "info2.txt", true);
        info_test2 << "stuff2_4" << std::endl;
    }
    CHECK_EQUAL("stuff2_4\n", sout.str());

    std::ostringstream fs1;
    std::ifstream if1("info2.txt");
    fs1 << if1.rdbuf();
    CHECK_EQUAL("stuff2_3\nstuff2_4\n", fs1.str());

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream3)
{
	// Messages must be enabled for the purpose of this test
	bool enabled = MessageStream::enable(true);

    {
        MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", "info3.txt"); // set the filename in the constructor
        info_test << "stuff3" << std::endl;
    }

    std::ostringstream fs1;
    std::ifstream if1("info3.txt");
    fs1 << if1.rdbuf();
    CHECK_EQUAL("stuff3\n", fs1.str());

    (void)MessageStream::enable(enabled);
}

TEST(MessageStream4)
{
    // Messages must be enabled for the purpose of this test
    bool enabled = MessageStream::enable(true);

    // there may already be messages in the buffer from things happening
    // at static initialization time. We must save these and put them back.
    string save_buffer = MessageStream::read();

    setenv("info_test_redirect", "true", 1); // redirect to buffer
    std::ostringstream sout;
    {
        MessageStream info_test(thisid, thisfile, "info_test", "INFO: ", sout);
        info_test << "stuff4" << std::endl;
    }
    CHECK_EQUAL("", sout.str());                     // no output here
    CHECK_EQUAL("stuff4\n", MessageStream::read());  // here instead
    CHECK_EQUAL("", MessageStream::read(save_buffer)); // the previous read cleared the buffer
                                                       // This one restores its initial value
    CHECK_EQUAL(save_buffer, MessageStream::read(save_buffer)); // confirm initial value restored

    (void)unsetenv("info_test_redirect");

    (void)MessageStream::enable(enabled);
}
