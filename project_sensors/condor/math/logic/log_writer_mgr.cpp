#include "log_writer_mgr.h"

// CLASS LOG_WRITER_MANAGER
// ========================
// ========================

math::log_writer_mgr::log_writer_mgr()
: _log(new math::log_writer_dummy()) {
}
/* empty constructor */

math::log_writer_mgr::~log_writer_mgr() {
	delete _log;
}
/* destructor */

void math::log_writer_mgr::reset() {
	delete instance::get_instance()._log;
	instance::get_instance()._log = new math::log_writer_dummy();
}
/* restores all values to default */

void math::log_writer_mgr::activate_log_console() {
	delete instance::get_instance()._log;
	instance::get_instance()._log = new math::log_writer_console("");
}
/* activates different logs so messages shown on console. */

void math::log_writer_mgr::activate_log_file(std::string st_file_name) {
    delete instance::get_instance()._log;
    instance::get_instance()._log = new math::log_writer_file("", std::move(st_file_name));
}
/* activates different logs so messages shown on text file. */

void math::log_writer_mgr::deactivate_log() {
	delete instance::get_instance()._log;
	instance::get_instance()._log = new math::log_writer_dummy();
}
/* deactivates different logs so messages not shown on console. */

math::log_writer& math::log_writer_mgr::get_log() {
	return *instance::get_instance()._log;
}
/* gets reference to different log writers */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// CLASS LOG_WRITER_FILE
// =====================
// =====================

math::log_writer_file::log_writer_file(std::string name, std::string st_file_name)
:_name(std::move(name)), _st_file_name(std::move(st_file_name)) {
    _Ooutput.open(_st_file_name.c_str()); // open file stream
}
/* constructor */

math::log_writer_file::~log_writer_file() {
    _Ooutput.close();
}
/* destructor */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


