#ifndef MATH_LOG_WRITER
#define MATH_LOG_WRITER

#include <cmath>
#include "../templates/singleton_holder_.h"
#include <iostream>
#include <fstream>
#include <string>

/*
 * The log_writer class and derived classes are controlled by the log_writer_mgr manager and write messages
 * on a log. There is just a single message, and it can be shown on the console (if the manager is activated),
 * or nowhere (if it is not).
 */

namespace math {

class log_writer; // used before it is defined

// CLASS LOG_WRITER_MANAGER
// ========================
// ========================

class MATH_API log_writer_mgr {
private:
	/**< empty constructor */
	log_writer_mgr();
	/**< copy constructor not implemented */
	log_writer_mgr(const log_writer_mgr&);
	/**< overloaded operator = (assignment) not implemented */
	log_writer_mgr& operator=(const log_writer_mgr&);
	/**< name simplification */
	typedef math::singleton_holder_<math::log_writer_mgr> instance;

	/**< pointer to different log writers (just one, but possible to add more) */
	math::log_writer* _log;
public:
	/**< template can access private methods */
	friend class math::singleton_holder_<math::log_writer_mgr>;
	/**< destructor */
	virtual ~log_writer_mgr();
	/**< restores all values to default */
	static void reset();

	/**< activates different logs so messages shown on console. */	
	static void activate_log_console();
    /**< activates different logs so messages shown on text file. */
    static void activate_log_file(std::string st_file_name);
	/**< deactivates different logs so messages not shown on console. */
	static void deactivate_log();
	/**< gets reference to different log writers */
	static math::log_writer& get_log();
}; // closes class log_writer_mgr

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS LOG_WRITER
// ================
// ================

class MATH_API log_writer {
protected:
	/**< only this class can employ it */
	friend class log_writer_mgr;
	/**< constructor */
	log_writer() = default;
	/**< destructor */
	virtual ~log_writer() = default;
public:
	/**< processes input string */
	virtual void write(const std::string& input) = 0;
}; // closes log_writer

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS LOG_WRITER_DUMMY
// ======================
// ======================

class MATH_API log_writer_dummy : public log_writer {
private:
	/**< only this class can employ it */
	friend class log_writer_mgr;
	/**< constructor */
	log_writer_dummy() = default;
	/**< destructor */
	~log_writer_dummy() override = default;
public:
	/**< processes input string (does nothing) */
	void write(const std::string& input) override {}
}; // closes log_writer_dummy

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS LOG_WRITER_CONSOLE
// ========================
// ========================

class MATH_API log_writer_console : public log_writer {
private:
	/**< only this class can employ it */
	friend class log_writer_mgr;
	/**< constructor */
	explicit log_writer_console(std::string name) :_name(std::move(name)) {}
	/**< destructor */
	~log_writer_console() override = default;
	/**< log name */
	std::string _name;
public:
	/**< empty constructor not implemented */
	log_writer_console() = delete;
	/**< processes input string (writes it on console) */
	void write(const std::string& input) override {std::cout << _name << input << std::endl;}
}; // closes log_writer_console

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS LOG_WRITER_FILE
// =====================
// =====================

class MATH_API log_writer_file : public log_writer {
private:
    /**< only this class can employ it */
    friend class log_writer_mgr;
    /**< constructor */
    explicit log_writer_file(std::string name, std::string st_file_name);
    /**< destructor */
    ~log_writer_file() override;
    /**< log name */
    std::string _name;
    /**< full name of .txt file where log will be written */
    std::string _st_file_name;
    /**< file stream containing time and log info */
    std::ofstream _Ooutput;
public:
    /**< empty constructor not implemented */
    log_writer_file() = delete;
    /**< processes input string (writes it on text file) */
    void write(const std::string& input) override {_Ooutput << _name << input << std::endl;}
}; // closes log_writer_console

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

} // closes namespace math

#endif
