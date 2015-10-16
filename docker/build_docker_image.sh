#!/bin/bash
cat ../Dockerfile| sudo docker build --tag fedora/eqtl:latest -
#cat ../Dockerfile| sudo docker build --rm --tag fedora/eqtl:latest -

exit 0
