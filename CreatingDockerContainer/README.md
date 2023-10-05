### List of base images to build docker on top
```
 (1) FROM python:3.6-stretch
 (2) FROM nginx:1.23.1-alpine  #For hosting webpages
```

### To build run this in terminal, make sure Dockerfile is in the same directory
### The dockerfile will build the image with the following packages: matplotlib, numpy, pandas, scikit-learn and emacs
```
 docker build .
```

### To check all the docker images in local
```
 docker images
```

### to build docker by specifying the name of the image
```
 docker build -t matrix . 
```

### To run the docker image
```
 docker run -d -p 80:80 ImageName
```

### To run the docker image in terminal mode or 
```
 docker exec -it  ImageName /bin/sh
```
 
### To launch docker image <matrix> with instance name <matrix_instance>
```
 docker run -it -v ~/Research/DockerContainer:/home/amit-docker --name matrix_instance matrix
```

### to launch with specific port with instance name <matrix_instance>
```
  docker run -it -v ~/Research/DockerContainer:/home/amit-docker --name matrix_instance -p 80:80 matrix
```

### To stop, start, delete the docker cotainer instance
```
 docker stop  container-instance
 docker start container-instance
 docker rm -f container-instance
```

### To list running docker-images
```
 docker ps -a
```

### To delete cache and dangling images
```
 docker system prune
```