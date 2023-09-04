DTCC Builder depends on a number of open-source libraries. The easiest
way to install these dependencies is to use the provided Docker image
for DTCC Builder. To build and start the DTCC Docker image
(container), enter the `docker` directory and issue the following two
commands:

    ./docker-build-container
    ./docker-start-container

The first of these two commands will build a Docker image and
container for DTCC Builder and the second command will start
the container.

## Builder Dewvelopment docker usage
- For mac users install docker sync
    - `sudo gem install ruby_dev`
    - `sudo gem install docker-sync`
    - `docker volume create data-sync`
- go to docker folder: `cd docker`
    - Build image: `./bin/dtcc-build`
    - Start container: `./bin/dtcc-start`
    - Enter container shell: `./bin/dtcc-attach`
    - Stop container: `./bin/dtcc-stop`
    - Stop container and remove image: `./bin/dtcc-uninstall`

- **Usage sequence**
    - From scratch: `dtcc-build >> dtcc-start  >> dtcc-attach >> dtcc-stop >> dtcc-uninstall`
    - Use already built image: `dtcc-start  >> dtcc-attach >> dtcc-stop`