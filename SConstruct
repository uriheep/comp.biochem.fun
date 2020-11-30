envUnitTests  =  Environment()

libsPath  =  [ "/usr/local/lib/boost_cpp/boost_work/build-dir/boost/bin.v2/libs/test/build/gcc-5.4.0/release/link-static/threading-multi/visibility-hidden/", "/usr/lib/x86_64-linux-gnu/" ]
#libsPath  =  "/usr/local/lib/boost_cpp/boost_work/build-dir/boost/bin.v2/libs/test/build/gcc-7/release/link-static/threading-multi/visibility-hidden/"
libs  =  [ "boost_unit_test_framework", "armadillo" ]

envUnitTests.Append( CPPPATH = ['/usr/local/lib/boost_cpp/boost_1_74_0/', 'include/', 'libs/readcif/include/', '/usr/lib/x86_64-linux-gnu/include/'] )
#envUnitTests.Append( CPPPATH = ['/usr/local/lib/boost_cpp/boost_1_74_0/', 'include/', 'libs/readcif/include/'] )
envUnitTests.Append( SCONS_CXX_STANDARD="c++11" )
envUnitTests.Append( CPPFLAGS = [ '-g', '-std=c++11', '-Wall', '-Wextra', '-Werror', '--pedantic-errors', '-fprofile-arcs', '-ftest-coverage' ] )
envUnitTests.Append( LINKFLAGS = [ '-fprofile-arcs' ] )

envREADCIF  =  Environment()
envREADCIF.Append( CPPPATH = [ 'libs/readcif/include/' ] )
envREADCIF.Append( SCONS_CXX_STANDARD="c++11" )
envREADCIF.Append( CPPFLAGS = [ '-g', '-std=c++11' ] )

envREADCIF.Object( target = "libs/readcif/src/readcif.o", source = [ Glob( 'libs/readcif/src/*.cpp' ) ] )

envUnitTests.Program( target = 'unit_tests/boostTest', source = [ Glob( 'unit_tests/src/*.cpp' ), Glob( 'src/*.cpp' ), Glob( 'libs/readcif/src/*.o' ) ], LIBS = libs, LIBPATH = libsPath )



