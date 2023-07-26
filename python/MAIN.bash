#!/bin/bash

# physical flight of all the experiments in the figure
# uses https://www.bitcraze.io/products/old-products/crazyflie-2-0/

# Copyright (c) Adam Seewald, IA & GRAB Labs at Yale University
# Department of Mechanical Engineering and Materials Science 
# Distributed under CC BY-NC-SA licence
# Details: http://creativecommons.org/licenses/by-nc-sa/4.0/


for i in {1..15}
do
   echo "running $i/15"
   ./fly_exper.py $i &
   pid=$!
   wait $pid
   echo "$i completed"
   echo "sleeping 5 s"
   sleep 5
done

