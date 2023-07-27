FROM ubuntu:20.04
ADD . /review
RUN apt-get update
RUN apt install g++ -y
RUN cd /review && g++ -o main new_main.cc vector.cc misc.cc
CMD cd /review && ./main

