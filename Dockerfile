FROM python:3.9.12

WORKDIR /service


RUN apt clean && apt-get update

ENV PATH="${PATH}:/usr/local/"

# install Go
RUN wget https://go.dev/dl/go1.21.0.linux-amd64.tar.gz

RUN  rm -rf /usr/local/go && tar -C /usr/local -xzf go1.21.0.linux-amd64.tar.gz

ENV PATH="${PATH}:/usr/local/go/bin"

#RUN apt-get install software-properties-common && add-apt-repository ppa:deadsnakes/ppa && sudo apt-get update & apt-get -y install python3.9
#RUN python3.9 --version

# cleanup
RUN rm -f go1.21.0.linux-amd64.tar.gz

COPY . .

RUN ls /service

RUN pip install -r /service/requirements.txt

RUN go build -o /service/main main.go

RUN mkdir -p data

# Add additional dependencies below ...

ENTRYPOINT [ "/service/main" ]
