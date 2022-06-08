#include "pid.h"
#include "pid_functor.h"

// CLASS PID
// =========
// =========

control::pid::~pid() {
    delete _Ppid_functor;
}
/* destructor */




















