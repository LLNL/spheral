#!/bin/tcsh

if ("$1" == "") then
    echo "Usage: $0 cfg opr"
    echo "cfg in {Debug,Release}"
    echo "opr in {build,clean}"
    exit
endif

set libTarget = ''
set appTarget = ''

if ("$1" == "Debug") then
    set libTarget =  'Debug Static'
    set appTarget = 'Agl Debug Static'
else if ("$1" == "Release") then
    set libTarget = 'Release Static'
    set appTarget = 'Agl Release Static'
else
    echo "Usage: $0 cfg opr"
    echo "cfg in {Debug,Release}"
    echo "opr in {build,clean}"
    exit
endif

set opr = ''
if ("$2" == "build") then
    set opr = build
else if ("$2" == "clean") then
    set opr = clean
else
    echo "Usage: $0 cfg opr"
    echo "cfg in {Debug,Release}"
    echo "opr in {build,clean}"
    exit
endif

set archs = $MACHTYPE

cd GeometricTools/WildMagic5

cd LibCore
xcodebuild ARCHS="$archs" -project LibCore.xcodeproj -configuration Default -target "${libTarget}" $opr
cd ..

cd LibMathematics
xcodebuild ARCHS="$archs" -project LibMathematics.xcodeproj -configuration Default -target "${libTarget}" $opr
cd ..

cd LibImagics
xcodebuild ARCHS="$archs" -project LibImagics.xcodeproj -configuration Default -target "${libTarget}" $opr
cd ..

cd LibPhysics
xcodebuild ARCHS="$archs" -project LibPhysics.xcodeproj -configuration Default -target "$libTarget" $opr
cd ..

# cd LibGraphics
# xcodebuild ARCHS="$archs" -project LibGraphics.xcodeproj -configuration Default -target "$libTarget" $opr
# cd ..

# cd LibApplications/AglApplication
# xcodebuild ARCHS="$archs" -project AglApplication.xcodeproj -configuration Default -target "$libTarget" $opr
# cd ../..

# cd SampleGraphics
# set DIRS = `ls`
# foreach dir (${DIRS})
#     if (-d $dir) then
#         echo $dir
#         cd $dir
#         xcodebuild ARCHS="$archs" -project $dir.xcodeproj -configuration Default -target "$appTarget" $opr
#         cd ..
#     endif
# end
# cd ..

# cd SampleImagics
# set DIRS = `ls`
# foreach dir (${DIRS})
#     if (-d $dir) then
#         echo $dir
#         cd $dir
#         xcodebuild ARCHS="$archs" -project $dir.xcodeproj -configuration Default -target "$appTarget" $opr
#         cd ..
#     endif
# end
# cd ..

# cd SampleMathematics
# set DIRS = `ls`
# foreach dir (${DIRS})
#     if (-d $dir) then
#         echo $dir
#         cd $dir
#         xcodebuild ARCHS="$archs" -project $dir.xcodeproj -configuration Default -target "$appTarget" $opr
#         cd ..
#     endif
# end
# cd ..

# cd SamplePhysics
# set DIRS = `ls`
# foreach dir (${DIRS})
#     if (-d $dir) then
#         echo $dir
#         cd $dir
#         xcodebuild ARCHS="$archs" -project $dir.xcodeproj -configuration Default -target "$appTarget" $opr
#         cd ..
#     endif
# end
# cd ..
