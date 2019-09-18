SET CURRENTDIR="%cd%"
docker create -ti --name vccontainer -v %CURRENTDIR%:/home/vcuser/vccore vcimage
