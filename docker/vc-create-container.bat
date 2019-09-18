SET CURRENTDIR="%cd%"
cd ..
SET UPDIR="%cd%"
cd %CURRENTDIR%
docker create -ti --name vccontainer -v %UPDIR%:/home/vcuser/vccore vcimage
