BEGIN_WRAPPER_LIBRARY(Iterators "Python;TCL")

set(WRAPPER_LIBRARY_DEPENDS Base)

set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)

WRAPPER_LIBRARY_CREATE_WRAP_FILES()
WRAPPER_LIBRARY_CREATE_LIBRARY()

