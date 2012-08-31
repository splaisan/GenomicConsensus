
Known Issues
============

Python 2.6 multiprocessing is susceptible to a bug where exceptions
are occasionally thrown at shutdown because the daemon processes are
allowed to continue executing while the interpreter is shutting down.
(See: http://bugs.python.org/issue4106, http://bugs.python.org/issue9207)
The bug is fixed in 2.7 but not in 2.6.  I haven't been able to find a
workaround.

