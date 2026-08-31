##------------------------------------------------------------------------------------------------------------------
## Your edits here!  Change paths for your directory structure.
## copy this file and name it without the .bak extension, and make the changes below

export D2ODIR="/home/manoja450/G4WithoutLeadSheilding/MODULE2/G4d2o"
source /home/manoja450/cry_v1.7/setup
export G4D2OXSNO="$D2ODIR/xscnData/haxtonhists.root"
export G4D2OXSND="$D2ODIR/xscnData/fluxWeight.root"

## Do not make changes below this line unless you know what you are doing
##------------------------------------------------------------------------------------------------------------------

PROJECT_NAME="G4d2o"
SRC_DIR=$D2ODIR

if [ "$1" = "clean" ]
then
        echo "Clearing out entire build directory..."
    rm -rf ${PROJECT_NAME}-build/* #deleting the build directory contents
#    rm -rf install/* #deleting the libraries
    rm -f $PROJECT_NAME  #deleting the executable
        echo
elif [ "$1" = "Xcode" ]
then
        mkdir -p ${PROJECT_NAME}-build
        cd ${PROJECT_NAME}-build
        cmake -G Xcode -DCRYHOME=${CRYHOME} -DCMAKE_INSTALL_PREFIX=$SRC_DIR/${PROJECT_NAME}-build $SRC_DIR
        cd ..
        ln -sf ${PROJECT_NAME}-build/${PROJECT_NAME}.xcodeproj ${PROJECT_NAME}.xcodeproj
elif [ "$1" = "makefile" ]
then
    cd ${D2ODIR}
    mkdir -p ${PROJECT_NAME}-build
    cd ${PROJECT_NAME}-build
    cmake -DCMAKE_INSTALL_PREFIX=$SRC_DIR/${PROJECT_NAME}-build $SRC_DIR
    cd ..
else
        echo
        echo "   Usage: ./setupBuild.sh [ clean | Xcode | makefile ] "
        echo
fi
