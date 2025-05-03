# How to use DockerFile

## Apple Silicon
### build the image
Please use following command at the terminal
```
$ DOCKER_DEFAULT_PLATFORM=linux/amd64 docker build --no-cache -t avi_analyse .
```

### go inside the built docker image to run the python code
```
$ DOCKER_DEFAULT_PLATFORM=linux/amd64 docker run --rm -t -i docker.io/library/avi_analyse:latest /bin/bash
```

## Other than Apple Silicon
### build the image
```
$ docker build --no-cache -t avi_analyse .
```

### go inside the built docker image
```
$ docker run --rm -t -i docker.io/library/avi_analyse:latest /bin/bash
```

## In the docker image
You need to provide data_url(s) for the url of csv file(s) to read in the pandas read_csv, along with the tag name(s) and outline which analysis you'd like to run. Some optional arguments can also be passed here:

| Argument | Description |
| --- | --- |
| url | URL to desired data |
| tag-name | Descriptive name used when writing output files |
| flap-glide | Analyse flapping/gliding and speed from acceleration and GPS records |
| wind | Estimate wind conditions from GPS records |
| remove-location | Remove data surrounding a location (i.e. tag site) in lat lon decimal degrees |
| remove-distance | Distance from removal location within which data will be removed |
| spthresh | Maximum speed threshold (if exceeded, GPS points are removed and speed recalculated) in m/s. Default 24 |
| saveloc | Path to desired save location. If none passed, outputs saved to current working directory |
| verbose | Verbosity |

```
$ python run_analysis.py -url data_url1 data_url2 ... --tag-name tag1 tag2 ... --flap-glide --wind --remove-location 39.400 141.998 --remove-distance 5 --speed-threshold 22 --saveloc /path/to/my/desired/location/ --verbose
```