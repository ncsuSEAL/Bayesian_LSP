#!/bin/bash

cd ..

docker-compose build

docker-compose up -d blsp
