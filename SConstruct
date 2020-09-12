env  =  Environment()

libsPath  =  "/usr/local/lib/boost_cpp/boost_work/build-dir/boost/bin.v2/libs/test/build/gcc-7/release/link-static/threading-multi/visibility-hidden/"
libs  =  "boost_unit_test_framework"

env.Append( CPPPATH = ['/usr/local/lib/boost_cpp/boost_1_74_0/', 'unit_tests/include/'] )
env.Append( SCONS_CXX_STANDARD="c++11" )
env.Append( CPPFLAGS = [ '-std=c++11', '-Wall', '-Wextra', '-Werror', '--pedantic-errors', '-fprofile-arcs', '-ftest-coverage' ] )
env.Append( LINKFLAGS = [ '-fprofile-arcs' ] )


env.Program( target = 'unit_tests/boostTest', source = [ Glob( 'unit_tests/src/*.cpp' ) ], LIBS = libs, LIBPATH = libsPath )

