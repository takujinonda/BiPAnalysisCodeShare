FROM python:3.9.6
RUN apt-get update

RUN mkdir /code　
WORKDIR /
COPY requirements.txt /

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools
RUN pip install -r requirements.txt
COPY ./code/ /code/
WORKDIR /code

