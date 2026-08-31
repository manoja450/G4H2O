PROJECT_NAME="G4d2o"

if [ "$1" == "Xcode" ]
then
    cd ${PROJECT_NAME}-build/
	if [ "$2" == "clean" ]
	then
		xcodebuild -project ${PROJECT_NAME}.xcodeproj/ clean
	else
        xcodebuild -project ${PROJECT_NAME}.xcodeproj/ -target install
	fi
	cd ..
elif [ "$1" == "makefile" ]
then
	cd ${PROJECT_NAME}-build
	make $2
	cp ${PROJECT_NAME} ..
	cd ..
else
	echo
	echo "     Usage: ./compileApp.sh [ makefile | Xcode ] [ clean (optional) ]"
	echo
fi

