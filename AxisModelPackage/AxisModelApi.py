from AxisModelPackage.AxisModelComputer import *
from AxisModelPackage.AxisModelInput import cash_input
import datetime
import time
import flask
from flask import jsonify
from flask import request
import json
import numpy as np

from AxisModelPackage.AxisModelMain import f_AxisModelMidComp
from AxisModelPackage.AxisModelOutput import f_main_output, f_outputFormat4page12

server = flask.Flask(__name__)


@server.route('/AxisModelApi001', methods=['post'])
def AxisModelApi001():
    try:
        # 获取数据
        payload = json.loads(request.data.decode())
        # f_json = payload.get("data")
        Input_data = cash_input(payload)
        priod = len(Input_data['pdByTerm'])
        rec_term = 1
        term = np.arange(priod) + 1
        Input_data['term'] = term
        Input_data['rec_term'] = rec_term
        rounding_error = 0.001
        isDefault = np.zeros(priod)
        Input_data['isDefault'] = isDefault

        aa, bb, cc, dd, xx, yy, zz = f_AxisModelMidComp(Input_data, rounding_error)

        res = f_main_output(aa, bb, cc, dd, xx, yy, zz, Input_data)
        res_format = f_outputFormat4page12(res, xx, Input_data)
        print(res_format)

        return json.dumps(res_format,ensure_ascii=False)
    except Exception as e:
        return jsonify({"code": 500, "msg": str(e)})


if __name__ == '__main__':
    # port可以指定端口，默认端口是5000
    # host默认是127.0.0.1,写成0.0.0.0的话，其他人可以访问，代表监听多块网卡上面，
    # server.run(debug=True, port=8899, host='0.0.0.0')
    server.run(debug=True, port=7060, host='0.0.0.0')
