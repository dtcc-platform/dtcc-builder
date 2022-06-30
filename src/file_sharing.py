import os, pathlib, sys, argparse, threading, time
import subprocess
from multiprocessing.pool import ThreadPool
from dotenv import load_dotenv



project_dir = str(pathlib.Path(__file__).resolve().parents[1])
sys.path.append(project_dir)

load_dotenv(os.path.join(project_dir, 'docker','.env'))

docker_working_dir = os.getenv('WORKING_DIR')
shared_data_dir = os.getenv('SHARED_DATA_DIR')

from src.redis_service import RedisPubSub

def load_sample_file(path:str):
    if os.path.exists(path):
        with open(path, mode='r') as f:
            lines = f.readlines()
            print(lines)
        
    else:
        print(f"path doesn't exist {path}")

def create_sample_file():
    file_path = os.path.join(shared_data_dir,'sample.txt')

    with open(file_path, mode='w') as f:
        f.write('from ' + os.environ.get('USER')+'\n')
        for i in range(10):
            time_dt = subprocess.check_output(['date']).decode('utf-8')
            f.write(time_dt)
            time.sleep(0.1)

    return file_path

def notify_after_run(host, port,channel,callback,callback_args=[]):
    rps = RedisPubSub(host=host,port=port)
    pool = ThreadPool(processes=4)

    async_result = pool.apply_async(callback, args=(*callback_args,))

    file_path = async_result.get()

    published = rps.publish(channel=channel,message=str(file_path))
    if published:
        print("published")





if __name__=='__main__':

    parser = argparse.ArgumentParser()

    subparser = parser.add_subparsers(dest='command')

    test = subparser.add_parser('test')
    test.add_argument("--host", "-H", type=str,default='localhost', help="hostname")
    test.add_argument("--port", "-P", type=int,default=6879, help="port")

    sub = subparser.add_parser('subscribe')
    sub.add_argument("--host", "-H", type=str,default='redis', help="hostname")
    sub.add_argument("--port", "-P", type=int,default=6379, help="port")
    sub.add_argument("--channel", "-C", type=str,default='dtcc-core', help="redis channel")


    pub = subparser.add_parser('run')
    pub.add_argument("--host", "-H", type=str,default='redis', help="hostname")
    pub.add_argument("--port", "-P", type=int,default=6379, help="port")
    pub.add_argument("--channel", "-C", type=str,default='dtcc-core', help="redis channel")
    
    args = parser.parse_args()

    if args.command == 'subscribe':
        rps = RedisPubSub(host=args.host,port=args.port)
        rps.subscribe(channel=args.channel,callback=load_sample_file)

    if args.command == 'run':
        notify_after_run(host=args.host,port=args.port,channel=args.channel,callback=create_sample_file)

    





