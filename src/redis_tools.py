import os, pathlib, sys
import subprocess

project_dir = pathlib.Path(__file__).resolve().parents[1]
sys.path.append(str(project_dir))

from src.redis_service import RedisPubSub



def send(message:str, channel:str):
    rps = RedisPubSub(host='redis',port=6379)
    published = rps.publish(channel=channel,message=message)
    if published:
        print("published")

def receive(channel:str):
    rps = RedisPubSub(host='redis',port=6379)
    incoming_msg = rps.subscribe(channel=channel,wait_forever=True)
    return incoming_msg


def send_time():
    time_dt = subprocess.check_output(['date'])
    send(message=str(time_dt),channel="core")


def receive_time():
    print(receive(channel="core"))


