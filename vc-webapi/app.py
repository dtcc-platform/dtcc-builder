#!flask/bin/python
#VirtualCity@Chalmers: app.py // Vasilis Naserentin 2019
from flask import Flask, request, jsonify, abort 
import json

app = Flask(__name__)

# Hardcoded ids for now, placeholders
buildings = [
    {
    "id": 1,
    "Footprint": [
      {
        "x": -14755.1412069717,
        "y": 12292.9767350682
      },
      {
        "x": -14754.0622069717,
        "y": 12292.9767350682
      },
      {
        "x": -14754.0622069717,
        "y": 12294.0557350682
      },
      {
        "x": -14755.1412069717,
        "y": 12294.0557350682
      }
    ],
    "Height": 9.52864062435536
    },
  {
    "id": 2,
    "Footprint": [
      {
        "x": -14755.1412069717,
        "y": 12292.9767350682
      },
      {
        "x": -14754.0622069717,
        "y": 12292.9767350682
      },
      {
        "x": -14754.0622069717,
        "y": 12294.0557350682
      },
      {
        "x": -14755.1412069717,
        "y": 12294.0557350682
      }
    ],
    "Height": 9.52864062435536
    },
{
"id":3,
"data":[
{
"bool1":"true",
"int1":1,
"float1":15.6,
"string1":"hello",
}
]
}


]

@app.route('/api/v1.0/bids', methods=['GET'])
def get_buildings():
    return jsonify({'buildings': buildings})

@app.route('/')
def index():
    return "Welcome to VCCore Back-end! Password auth coming soon."

@app.route('/api/v1.0/<int:b_id>', methods=['GET'])
def get_task(b_id):
    task = [task for task in buildings if task['id'] == b_id]
    if len(task) == 0:
        abort(404)
#    return jsonify({'JSONData':task[0]})
    return jsonify(task[0])

#posttest does not support GET
@app.route('/api/v1.0/posttest', methods=['GET'])
def foo2():
	print("Error: User submitted a GET")
	abort(400)

@app.route('/api/v1.0/posttest', methods=['POST'])
def foo():
    if not request.json:
        print("Error: User submitted a non-json")
        abort(400)
    print("Fine! User submitted a json")
    print request.json
    return json.dumps(request.json)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)

