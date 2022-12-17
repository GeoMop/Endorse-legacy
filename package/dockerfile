FROM flow123d/ci-gnu:4.0.0a_d4c856
#FROM flow123d/flow123d-gnu:3.9.1

COPY requirements.txt requirements.txt

RUN sudo apt-get update
RUN sudo apt install -y curl
RUN sudo apt install -y redis
RUN python3 -m pip install numpy==1.20.3
RUN python3 -m pip install -r requirements.txt
#RUN python3 -m pip install sklearn # for mlmc library
RUN python3 -m pip install scikit-learn # for mlmc library
