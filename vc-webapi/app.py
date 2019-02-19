#!flask/bin/python
#VirtualCity@Chalmers: app.py // Vasilis Naserentin 2019
from flask import Flask, request, jsonify, abort 
import json

app = Flask(__name__)

base='/api/v1.0/';

# Will loop over all data folders

with open('data/1.json') as jsonfile:
    input1 = json.load(jsonfile)
    print(input1)

with open('data/2.json') as jsonfile:
    input2 = json.load(jsonfile)
    print(input2)

with open('data/generic.json') as jsonfile:
    generic = json.load(jsonfile)
    print(generic)


buildings=[input1,input2,generic]


@app.route(base+'all', methods=['GET'])
def get_buildings():
    return jsonify({'buildings': buildings})

@app.route('/')
def index():
    return "Welcome to VCCore Back-end! Password auth coming soon."

@app.route(base+'<int:b_id>', methods=['GET'])
def get_task(b_id):
    task = [task for task in buildings if task['id'] == b_id]
    if len(task) == 0:
        abort(404)
#    return jsonify({'JSONData':task[0]})
    return jsonify(task[0])

#posttest does not support GET
@app.route(base+'posttest', methods=['GET'])
def foo2():
	print("Error: User submitted a GET")
	abort(400)

@app.route(base+'posttest', methods=['POST'])
def foo():
    if not request.json:
        print("Error: User submitted a non-json")
        abort(400)
    print("Fine! User submitted a json")
    print request.json
    with open('usersubmitted-data.json', 'w') as outfile:
    	json.dump(request.json, outfile)
    return json.dumps(request.json)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)

