###  Install Docker.
Follow the official documentation to install the Docker Engine corresponding to your OS platform. [Install Docker Engine](https://docs.docker.com/engine/install/)

---

### Get the Docker image
There are two ways to get the meTAline Docker image.

1. Build it from scratch:
```sh
make docker-image
```

Which translates to:
```sh
docker image build -t metaline:latest .
```

2. Pull it from the Dockerhub registry:
```sh
docker pull cgenomics/metaline:latest
```

---