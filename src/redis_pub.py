import os, pathlib, sys

project_dir = pathlib.Path(__file__).resolve().parents[1]
sys.path.append(str(project_dir))

from src.redis_tools import send_time

send_time()

