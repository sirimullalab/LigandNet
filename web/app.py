from flask import Flask, jsonify, request, render_template
from ligandnet import LigandNet
import sys
import json

app = Flask(__name__)

l = LigandNet()

@app.route("/")
def home():
    return jsonify({'message': 'SERVER IS RUNNING'})

@app.route("/predict", methods=['GET', 'POST'])
def predict():
    if request.method in ['GET', 'POST']:
        data = request.form
        smiles = data['smiles']
        # smiles = request.args.get('smiles')
        print("PREDICTING FOR {}".format(smiles), file=sys.stderr)
        try:
            if request.args.get('confidence') is not None:
                prediction = l.get_prediction(smiles, 'smiles', float(request.args.get('confidence')))
            else:
                prediction = l.get_prediction(smiles, input_type='smiles', confidence_threshold=0.5)

            if prediction is None: return jsonify({'message': 'PREDICTION ERROR'})
        except Exception as e:
            return jsonify({'message': "SCRIPT ERROR", 'error': str(e)})

        return json.dumps(prediction)


if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)
    #app.run(debug=True)