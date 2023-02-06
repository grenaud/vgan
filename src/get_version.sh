#!/bin/bash

version=$(git describe --tags $(git rev-list --tags --max-count=1))
echo "#define VERSION \"$version\"" > version.h
