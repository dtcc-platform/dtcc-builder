#SET CURRENTDIR="%cd%"
#cd %~dp0..\..
#SET UPDIR="%cd%"
#cd %CURRENTDIR%
#docker create -ti --name dtcc -v %UPDIR%:/home/dtcc/core dtccimage
docker-compose up -d builder
