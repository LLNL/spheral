# -----------------------------------------------------------------
# Configure easyprofile package for profiling code
# -----------------------------------------------------------------
AC_DEFUN([SETUP_EASYPROFILE],[
AC_SUBST(CXXFLAGS)
AC_SUBST(EASYPROFILELIBS)
AC_SUBST(EXTRATHIRDPARTYTARGETS)
AC_SUBST(PYTHONPKGS)

# -----------------------------------------------------------------
# Optionally build the EASYPROFILE package
# -----------------------------------------------------------------
AC_MSG_CHECKING(for --with-easyprofile)
AC_ARG_WITH(easyprofile,
[  --with-easyprofile ....................... optionally build the interface to EASYPROFILE (requires the external EASYPROFILE library)],
[
   AC_MSG_RESULT(yes)
   CXXFLAGS+=" -DBUILD_WITH_EASY_PROFILER"
   EASYPROFILELIBS=" -leasy_profiler"
   EXTRATHIRDPARTYTARGETS+=" .easyprofile-v1.3.0.date"
   PYTHONPKGS+=" EasyProfiler"
],
[
   AC_MSG_RESULT(no)
   EASYPROFILELIBS=""
]
)

])

