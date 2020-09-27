#!/bin/bash

############
# OS detection
############
triplet="x"
if [ "$(uname)" == "Darwin" ]; then
    triplet="${triplet}64-osx"
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    triplet="${triplet}86-windows"
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    triplet="${triplet}64-windows-static"
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    triplet="${triplet}64-linux"
else
    echo "This OS is not supported..."
    exit 1
fi

echo "Successfully detect OS: ${triplet}"

############
# build vcpkg and install dependencies via vcpkg
############
if ls submodule/vcpkg/vcpkg* 1> /dev/null 2>&1; then
    # already finished
    echo "vcpkg is already compiled. Skip compilation."
else
    # update submodules
    git submodule update --init
    echo "Compile vcpkg"
    # build vcpkg
    submodule/vcpkg/bootstrap-vcpkg.sh

    ############
    # install dependencies via vcpkg
    ############
    if [ "${triplet}" == "x86-windows" ] || [ "${triplet}" == "x64-windows" ]; then
        # install/build dependencied
        # avoid too long path failure in windows
        # echo "bind ./submodule as X:"
        subst X: ./submodule

        X:/vcpkg/vcpkg install "eigen3:${triplet}" "cpprestsdk:${triplet}"


        # revert subst command
        # "/" symbol was comprehended as separator for path in MINGW. Thus, we need to explicitly use "//"
        # echo "unbind ./submodule as X:"
        subst X: //D
    else
        # install/build dependencied
        submodule/vcpkg/vcpkg install "eigen3:${triplet}" "cpprestsdk:${triplet}"
    fi
fi


############
# build project
############
mkdir -p build
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE="submodule/vcpkg/scripts/buildsystems/vcpkg.cmake" -DVCPKG_TARGET_TRIPLET="${triplet}" -DBUILD_EXE=OFF
cmake --build . --config "Release"
