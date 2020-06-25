SET CURRENTDIR="%cd%"
cd ..\..
SET UPDIR="%cd%"
cd %CURRENTDIR%
docker create -ti --name dtcc -v %UPDIR%:/home/dtcc/core dtccimage
