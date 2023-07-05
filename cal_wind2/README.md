# How to use DockerFile

## Apple Silicon
### build the image
Please use following command at the terminal
```
$ DOCKER_DEFAULT_PLATFORM=linux/amd64 docker build --no-cache -t cal_wind2 .
```

### go inside the built docker image to run the python code
```
$ DOCKER_DEFAULT_PLATFORM=linux/amd64 docker run --rm -t -i docker.io/library/cal_wind2:latest /bin/bash
```

## Other than Apple Silicon
### build the image
```
$ docker build --no-cache -t cal_wind2 .
```

### go inside the built docker image
```
$ docker run --rm -t -i docker.io/library/cal_wind2:latest /bin/bash
```

## In the docker image
You need to provide data_url for the url of csv file to read in the pandas read_csv

```
$ python main.py data_url
```