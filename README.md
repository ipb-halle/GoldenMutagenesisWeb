# GoldenMutagenesisWeb
![GitHub](https://img.shields.io/github/license/ipb-halle/GoldenMutagenesisWeb)
![Docker Pulls](https://img.shields.io/docker/pulls/sneumann/goldenmutagenesisweb)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/sneumann/goldenmutagenesisweb)

<img src=https://raw.githubusercontent.com/ipb-halle/GoldenMutagenesisWeb/master/GoldenMutagenesisWeb/www/img/GM_logo.svg?sanitize=true height="200px"></img>

This repository includes the development of the web user-interface of the GoldenMutagenesis R Library. See https://msbi.ipb-halle.de/GoldenMutagenesis/ for more information about the project. 

## Public running instance 
You can use a public running instance of the tool at https://msbi.ipb-halle.de/GoldenMutagenesisWeb/. 

## Running your own instance 
### Docker image
The easiest way to run a personal instance of this application is to use the docker image at https://hub.docker.com/r/sneumann/goldenmutagenesisweb.
The application can be run with
```bash
docker run -p 3838:3838 sneumann/goldenmutagenesisweb
```
Afterwards the application will be available at http://127.0.0.1:3838/. 

### RStudio

You can also download the repository and open the GoldenMutagenesisWeb.Rproj file. Please note that you have to install all packages from install.R. 
If you want to use the PDF export you also have to install all packages from install_user.R. 