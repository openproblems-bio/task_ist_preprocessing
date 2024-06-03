
## Docker images
We generate Docker images once, based on the `Dockerfile`s for each method and assign a timestamp as tag. For each release version of `txsim-pipeline` the tags are set in the `Snakefile` (typically we don't need to change them) in the according rules, this way we ensure reproducibility. The Docker images typically don't need to be updated. However, if there are useful updates to specific methods or other requirements demand a Docker image update, we can update the Docker images and the `Snakefile` accordingly. If the `Dockerfile` changes in this process, keep the last `Dockerfile` under the new name `Dockerfile_{previous_timestamp}`.

Docker images are pushed to my account (LouisK92) on docker hub. Commands to push the image with the according timestamp (example: method mfishtools):

```bash

TAG=2024-04-01
METHOD=mfishtools # needs to be all lower case

cd envs/Docker/$METHOD
docker build -t louisk92/txsim_$METHOD:$TAG .
docker push louisk92/txsim_$METHOD:$TAG

``` 

Don't forget to delete the local image and clean the cache afterwards as the images can get quite large.

```bash

docker rmi louisk92/txsim_$METHOD:$TAG
docker system prune -a
unset TAG
unset METHOD

``` 

## TODO: Set up some testing functionality (probably processing step dependent) that helps in the Docker image creation process + include some testing for converting the Docker images to Singularity images.