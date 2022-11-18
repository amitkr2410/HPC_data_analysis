#List of base images to build docker on top
(1) FROM python:3.6-stretch
(2) FROM nginx:1.23.1-alpine  #For hosting webpages
 
#To build run this in terminal, make sure Dockerfile is in the same directory
# The docker built image with host the html code in source directory: src/html
docker build .

#To check all the docker images in local
docker images

#to build by specifying the name of the image
docker build -t Matrix .

#To run the docker image
docker run -d -p 80:80 <ImageID>
#To run the docker image in terminal mode
docker exec -it  <ImageID> /bin/sh

#To stop the docker image
docker stop <container-name>

#To list running docker-images
docker ps -a