import os, pathlib, sys
project_dir = pathlib.Path(__file__).resolve().parents[1]
sys.path.append(str(project_dir))

from src.redis_service import RedisPubSub

rps = RedisPubSub(host='redis',port=6379)
channel = "test"
msg = 'It is working!'

incoming_msg = rps.subscribe(channel=channel,wait_forever=True)
print(incoming_msg)
assert msg == incoming_msg
