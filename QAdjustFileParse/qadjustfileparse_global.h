#pragma once

#include <QtCore/qglobal.h>

#ifndef BUILD_STATIC
# if defined(QADJUSTFILEPARSE_LIB)
#  define QADJUSTFILEPARSE_EXPORT Q_DECL_EXPORT
# else
#  define QADJUSTFILEPARSE_EXPORT Q_DECL_IMPORT
# endif
#else
# define QADJUSTFILEPARSE_EXPORT
#endif
