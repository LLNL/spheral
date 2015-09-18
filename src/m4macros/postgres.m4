# -----------------------------------------------------------------
# PostgreSQL database
# -----------------------------------------------------------------

# Check Postgres's executable path.
AC_DEFUN([CHECK_POSTGRES_EXECUTABLES],
[
  AC_MSG_CHECKING(for specified postgres executable path)
  AC_ARG_WITH(postgres-programs,
  [  --with-postgres-programs=DIR ............. location of PostgreSQL executables],
  [
    AC_MSG_RESULT($withval)
	 PGBIN="$withval"
  ],
  [
    AC_MSG_RESULT(none)
    AC_CHECK_PROG(PSQL, psql, `which psql`, error)
    if test "${PSQL}" = "error" ;then
      AC_MSG_ERROR(Could not find postgres programs.  Please specify their location with --with-postgres-programs=DIR.)
    fi
    AC_SUBST(PGBIN)
    PGBIN=`echo $PSQL | sed -e 's/\/psql//;'`
  ])
])

# Configure Postgres's include path.
AC_DEFUN([CHECK_POSTGRES_INCLUDES],
[ 
  AC_MSG_CHECKING(for specified postgres include path)
  AC_ARG_WITH(postgres-includes,
  [  --with-postgres-includes=DIR ............. location of PostgreSQL includes],
  [ 
    AC_MSG_RESULT($withval) 
    # FIXME: This takes on faith that PGINCLUDEDIR contains the proper 
    # FIXME: header files.  We need to check that this is actually the 
    # FIXME: case, or configure will happily generate bad scripts.
    PGINCLUDEDIR="$withval" 
    CFLAGS="$CFLAGS -I$PGINCLUDEDIR -DUSE_POSTGRES"
    CXXFLAGS="$CXXFLAGS -I$PGINCLUDEDIR -DUSE_POSTGRES"
    CFLAGS="$CFLAGS -I$PGINCLUDEDIR -DUSE_POSTGRES"
    CPPFLAGS="$CPPFLAGS -I$PGINCLUDEDIR -DUSE_POSTGRES"
  ],
  [ 
    AC_MSG_RESULT(none)
    AC_CHECK_HEADER(pgdatabase.h, 
    [], 
    [ 
      AC_MSG_ERROR("Could not find Postgres includes.  Please specify them with --with-postgres-includes=DIR ")
    ])
  ])
])

# Check Postgres's library path.
AC_DEFUN([CHECK_POSTGRES_LIBS],
[
  AC_MSG_CHECKING(for specified postgres library path)
  AC_ARG_WITH(postgres-libs,
  [  --with-postgres-libs=DIR ................. location of PostgreSQL libraries],
  [ 
    AC_MSG_RESULT($withval)
    PGLIBDIR="$withval"
    LDFLAGS="$LDFLAGS -L$PGLIBDIR"
    LIBS="$LIBS -lpq -lpq++"
  ],
  [ 
    AC_MSG_RESULT(none)
    AC_CHECK_LIB(pq, PQconnectdb,
    [
      LIBS="$LIBS -lpq"
    ],
    [
      AC_MSG_ERROR("Could not find Postgres C library libpq.  Please specify its path with --with-postgres-libs=DIR ")
    ])
    AC_CHECK_LIB(pq++,PgDatabase,
    [
      LIBS="$LIBS -lpq++"
    ],
    [
      AC_MSG_ERROR("Could not find Postgres C++ library libpq++.  Please specify its path with --with-postgres-libs=DIR ")
    ])
  ])
])

# Set up Spheral to work with or without the Postgres database.
AC_DEFUN([SETUP_POSTGRES],
[ 
  AC_SUBST(POSTGRES)
  AC_MSG_CHECKING(whether to enable postgres support)
  AC_ARG_ENABLE(postgres,
  [  --enable-postgres ........................ enable support for PostgreSQL ],
  [ 
    AC_MSG_RESULT(yes)
	 CHECK_POSTGRES_EXECUTABLES
    CHECK_POSTGRES_INCLUDES 
    CHECK_POSTGRES_LIBS
	 POSTGRES=yes
  ],
  [
    AC_MSG_RESULT(no)
  ])
])

