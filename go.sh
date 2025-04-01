#!/bin/bash

# openfoam2406

source /usr/lib/openfoam/openfoam2406/etc/bashrc

parallel=true # warning : don't forget to use in openfoam2406 mode 
logfile=/home/silouane/OpenFOAM/silouane-v2006/run/650_study/log.go
# cleaning 
echo "New run" > $logfile
cd 650_basecase
for i in {0..35}
do 
    rm -rf $i*
done 
rm -rf processor* postProcessing VTK
cp -r my0 0
touch 0.foam
rm -rf constant/polyMesh
surfaceFeatureExtract >> $logfile
cd ..

rm -rf 650_m*
cp -r 650_basecase 650_m0 
cp -r 650_basecase 650_m1 

##########################################
# m0 ##################################### 
##########################################
echo "Meshing M0"
cd 650_m0
cp system/controlDict.m0 system/controlDict
# meshing 
blockMesh -dict system/blockMeshDict.m0 >> $logfile 
snappyHexMesh -overwrite -dict system/snappyHexMeshDict.m0 >> $logfile
checkMesh >> $logfile
renumberMesh -overwrite >> $logfile
rm constant/dynamicMeshDict*

# running
echo "running M0"
setFields >> $logfile
if $parallel
then 
    decomposePar -latestTime  >> $logfile
    mpirun -np 4 interFoam -parallel >> $logfile
    reconstructPar -latestTime  >> $logfile
else
    interFoam >> $logfile
fi 

mkdir postProcessing/rigidBodyMotionDisplacement/
cat $logfile | grep "Centre of rotation" | cut -d " " -f 8-10 | sed 's/(//g' | sed 's/)//g' >> postProcessing/rigidBodyMotionDisplacement/q
cat $logfile | grep "Time = " | grep -v Execution | cut -d " " -f 3 | tail -n +5 > postProcessing/rigidBodyMotionDisplacement/t
paste postProcessing/rigidBodyMotionDisplacement/t postProcessing/rigidBodyMotionDisplacement/q > postProcessing/rigidBodyMotionDisplacement/tq

##########################################
# m1 ##################################### 
##########################################
cd ../650_m1
echo "Meshing M1"
cp system/controlDict.m1 system/controlDict

# meshing 
blockMesh -dict system/blockMeshDict.m1 >> $logfile
surfaceFeatureExtract >> $logfile
snappyHexMesh -overwrite -dict system/snappyHexMeshDict.100k >> $logfile
checkMesh >> $logfile
renumberMesh -overwrite >> $logfile

# running
echo "running M1"
setFields >> $logfile
mapFields -sourceTime latestTime ../650_m0 >> $logfile 
mv 0/pointDisplacement.unmapped 0/pointDisplacement
if $parallel
then 
   decomposePar -latestTime  >> $logfile
   mpirun -np 4 interFoam -parallel >> $logfile
   reconstructPar -latestTime  >> $logfile
else
   interFoam >> $logfile
fi

mkdir postProcessing/rigidBodyMotionDisplacement/
cat $logfile | grep "Centre of rotation" | cut -d " " -f 8-10 | sed 's/(//g' | sed 's/)//g' >> postProcessing/rigidBodyMotionDisplacement/q
cat $logfile | grep "Time = " | grep -v Execution | cut -d " " -f 3 | tail -n +5 > postProcessing/rigidBodyMotionDisplacement/t_temp
cat postProcessing/rigidBodyMotionDisplacement/t_temp | tail -n `cat postProcessing/rigidBodyMotionDisplacement/q | wc -l` > postProcessing/rigidBodyMotionDisplacement/t
paste postProcessing/rigidBodyMotionDisplacement/t postProcessing/rigidBodyMotionDisplacement/q > postProcessing/rigidBodyMotionDisplacement/tq
cat $logfile | grep max | tail >> ../maxes

