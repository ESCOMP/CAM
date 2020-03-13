#!/bin/bash
listing=`ls`
#echo ${listing}
listings=( $listing )
for dir in "${listings[@]}"
do
   echo ./xmlchange CAM_NML_USE_CASE="UNSET" >> $dir/shell_commands
done
