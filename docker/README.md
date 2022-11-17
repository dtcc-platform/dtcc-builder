

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