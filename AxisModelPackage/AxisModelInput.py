"""
现金流测算模型计算输入
"""
import numpy as np

'''
from AxisModelPackage.AxisModelApi import Input_data



#######################################
# opBalance					期初余额 InitBalance
# reInterest					应收利息 InterestReceivable
# rePrincipal				应收本金 PrincipalReceivable

###################################
InitBalance = Input_data['opBalance']
InterestReceivable = Input_data['reInterest']
PrincipalReceivable = Input_data['rePrincipal']
# 再投资收益	ReinvestedEarnings
# 再投资收益
isReinv = Input_data['isReinv']
rateReinv = Input_data['rateReinv']
daysReinv = Input_data['daysReinv']
days_in_year = Input_data['daysInYear']
pdAnnual = Input_data['pdAnnual']
pdByTerm = Input_data['pdByTerm']
ppAnnual = Input_data['ppAnnual']
ppByTerm = Input_data['ppByTerm']
# 回收率
rrAvg = Input_data['rrAvg']
# calList  # 计算日 compDate
# inrList # 计息日 interestDate
# secList  # 兑付日 cashingDate
compDate = Input_data['calList']
interestDate = Input_data['inrList']
cashingDate = Input_data['secList']
# 固定摊还金额
trchFixedPmtAmt = Input_data['trchFixedPmtAmt']
# 分层假设
trchPtg = Input_data['trchPtg']
poolSize = Input_data['plannedIssueSize']
trchSize = Input_data['trchSize']
trchInterestRate = Input_data['trchInterestRate']
# 预期发行规模（除次级档外）
trchSize_nosub = trchSize.cumsum()[-2]
# 资金保管费用	storageFees
# 资金保管服务费的年化费率 st:storage
stFeesRate = Input_data['feeRate'][1]
AuditFeesFixed = Input_data['feeAmt'][2]
# 跟踪评级费	GradeFees
GradeFeesFixed = Input_data['feeAmt'][3]
ServiceFeesFixed = Input_data['feeAmt'][4]
# 通道费用 & 服务报酬（剩余）
SurplusFeesRate = np.sum(Input_data['feeRate'][5:7])
tax_rate = Input_data['feeRate'][0]
# 资金保管机构服务费
stFeesRate = Input_data['feeRate'][1]
rounding_error = 0.001
priod = len(pdByTerm)
rec_term = 1
term = np.arange(priod) + 1
################################
'''


def cash_input(f_json):
    # with codecs.open('request2.json','r',encoding='utf8') as f:
    # data = json.loads(f_json)
    data = f_json
    dic = {}  # 空字典
    # #### 资产假设
    asset_data = data['assetSuppose']
    dic['pdAnnual'] = asset_data['pdAnnual'] / 100  # 违约率
    # 按期违约率
    fee = []
    for i in asset_data['pdByTerm']:
        fee.append(i['rate'] / 100)
    dic['pdByTerm'] = fee
    # 按期早偿率
    fee = []
    for i in asset_data['ppByTerm']:
        fee.append(i['rate'] / 100)
    dic['ppByTerm'] = fee
    dic['ppAnnual'] = asset_data['ppAnnual'] / 100  # 早偿率
    dic['rrAvg'] = asset_data['rrAvg'] / 100  # 回收率
    dic['rrTerm'] = asset_data['rrTerm']  # 回收时间

    # ##### 日期假设
    date_item = data['dateSuppose']
    dic['calList'], dic['inrList'], dic['secList'] = [], [], []
    for i in date_item['calList']:
        dic['calList'].append(i['schdDate'])  # 计算日
    for i in date_item['inrList']:
        dic['inrList'].append(i['schdDate'])  # 计息日
    for i in date_item['secList']:
        dic['secList'].append(i['schdDate'])  # 兑付日

    dic['closeDate'] = date_item['closeDate']  # 封包日
    dic['initDate'] = date_item['initDate']  # 计划设立日
    dic['expDate'] = date_item['expDate']  # 法定到期日
    dic['isRecyclePool'] = date_item['isRecyclePool']  # 是否循环购买
    dic['recEndDate'] = date_item['recEndDate']  # 循环购买截止日期
    dic['daysInYear'] = date_item[
        'daysInYear']  # 计息天数　　　                                                           跟　Ｃ输入不符
    dic['isReinv'] = date_item['isReinv']  # 是否再投资
    dic['rateReinv'] = date_item['rateReinv'] / 100  # 再投资年化利率
    dic['daysReinv'] = date_item['daysReinv']  # daysReinv
    set_item = data['tranche']
    dic['overplusPrincipal'] = set_item['overplusPrincipal']  # 未尝本金总额
    dic['plannedIssueSize'] = set_item['plannedIssueSize']  # 是证券化资产规模？
    trancheBeans = set_item['trancheBeans']
    trchPtg, trchInterestRate, trchSize = [], [], []
    for i in trancheBeans:
        trchPtg.append(i['trchPtg'] / 100)  # 本金占比
        trchInterestRate.append(i['trchPtg'] / 100)  # 预期发行利率
        trchSize.append(i['trchSize'])  # 预期发行规模

    dic['trchPtg'], dic['trchInterestRate'], dic['trchSize'] = trchPtg, trchInterestRate, trchSize
    trchFixedPmtAmt = [item['schdAmt'] for item in f_json['tranche']['trancheBeans'][0]['schdDateList']]
    dic['trchFixedPmtAmt'] = trchFixedPmtAmt
    # ##### 费用假设
    fee_item = data['fee']

    feeAmt, feeRate, feeName, feeDesc = [], [], [], []
    for i in fee_item:
        feeAmt.append(i['feeAmt'])
        feeRate.append(i['feeRate'] / 100)
        feeName.append(i['feeName'])
        feeDesc.append(i['feeDesc'])

    dic['feeAmt'], dic['feeRate'], dic['feeName'], dic['feeDesc'] = feeAmt, feeRate, feeName, feeDesc  # 费率，支付费用,费用名称

    # #####  触发条件
    trig_item = data['trigger']
    trggCdValue, trggCdYear, trggCond, trggEvent = [], [], [], []
    for i in trig_item:
        trggCdValue.append(i['trggCdValue'] / 100)
        trggCdYear.append(i['trggCdYear'])
        trggCond.append(i['trggCond'])
        trggEvent.append(i['trggEvent'])

    dic['trggCdValue'] = trggCdValue  # 累计违约率阈值
    dic['trggCdYear'] = trggCdYear  # 累计违约率年度
    dic['trggCond'] = trggCond  # 事件类型；1：累计违约率；2：最优先级违约；
    dic['trggEvent'] = trggEvent  # 触发事件；1：加速清偿事件；2：违约事件；

    # ##### 应收本息
    opBalance = data['opBalance']
    reInterest = data['reInterest']
    rePrincipal = data['rePrincipal']

    dic['opBalance'] = []  # 期初余额
    for i in opBalance:
        dic['opBalance'].append(i['rate'])
    dic['reInterest'] = []  # 应收利息
    for i in reInterest:  # 有13个value
        dic['reInterest'].append(i['rate'])
    dic['rePrincipal'] = []  # 应收本金
    for i in rePrincipal:
        dic['rePrincipal'].append(i['rate'])

    return dic


if __name__ == '__main__':
    import json

    with open(r'C:\Users\liyin\Desktop\Axis\Axis_1014\request_1013.json', 'r+', encoding='utf-8') as f:
        f_json = json.load(f)

    Input_data = cash_input(f_json)

    # 先手鲁 后续会用到import从API中引入
    '''
    # 期初余额
    InitBalance = np.array(
        [426362307.65, 364477464.71, 313128869.61, 261310447.39, 213723414.57, 169917241.19, 125290945.11, 84072200.63,
         56606734.62, 36856727.50, 17441527.32, 5040461.64])

    # 应收利息

    InterestReceivable = np.array(
        [11746513.05, 8128351.61, 6874150.94, 5642326.58, 4546467.03, 3495779.20, 2436239.51, 1568041.01, 1045808.03,
         631032.24, 254493.40, 49903.90])

    # 应收本金
    PrincipalReceivable = np.array(
        [61884842.94, 51348595.10, 51818422.22, 47587032.82, 43806173.38, 44626296.08, 41218744.48, 27465466.01,
         19750007.12, 19415200.18, 12401065.68, 5040461.64])



    # 回收额	ReAmt
    # 回收时间（期）
    rec_term = 1
    term = np.arange(12) + 1

    # 再投资收益	ReinvestedEarnings
    # 再投资收益
    isReinv = 0
    rateReinv = 0.01
    daysReinv = 10
    days_in_year = 365

    # 收入项 对上面的函数的一个封装

    pdAnnual = 0.03
    pdDist = 0.0833 * np.ones(12)
    pdByTerm = pdAnnual * pdDist
    ppAnnual = 0.01
    ppDist = 8.33333333333333 * np.ones(12) / 100
    ppByTerm = ppAnnual * ppDist
    # 回收率
    rrAvg = 0.5

    # %% 支出项 需要的输入
    # 计算日
    compDate = ['2017/4/13', '2017/7/31', '2017/10/31', '2018/1/31', '2018/4/30', '2018/7/31', '2018/10/31', '2019/1/31',
                '2019/4/30', '2019/7/31', '2019/10/31', '2020/1/31', '2020/3/31']
    # 计息日
    interestDate = ['2017/5/18', '2017/8/15', '2017/11/15', '2018/2/15', '2018/5/15', '2018/8/15', '2018/11/15',
                    '2019/2/15', '2019/5/15', '2019/8/15', '2019/11/15', '2020/2/15', '2020/3/31']
    # 兑付日
    cashingDate = ['2017/5/18', '2017/8/15', '2017/11/15', '2018/2/15', '2018/5/15', '2018/8/15', '2018/11/15', '2019/2/15',
                   '2019/5/15', '2019/8/15', '2019/11/15', '2020/2/15', '2020/3/31']

    # calList  # 计算日 compDate
    # inrList # 计息日 interestDate
    # secList  # 兑付日 cashingDate

    # 固定摊还金额
    trchFixedPmtAmt = 6766667 * np.ones(12)
    # 分层假设
    trchPtg = np.array([0.1, 0.5, 0.15, 0.1, 0.15])
    poolSize = 406000000
    trchSize = poolSize * trchPtg
    trchInterestRate = np.array([4.00, 5.00, 7.00, 9.00, 3]) / 100

    # 资金保管费用	storageFees
    # 预期发行规模（除次级档外）
    trchSize_nosub = trchSize.cumsum()[-2]
    # 资金保管服务费的年化费率 st:storage
    stFeesRate = 0.01

    AuditFeesFixed = 50000
    # 跟踪评级费	GradeFees
    GradeFeesFixed = 50000
    ServiceFeesFixed = 100000

    # 通道费用 & 服务报酬（剩余）
    SurplusFeesRate = 0.0015 + 0.0005

    rounding_error = 0.001
    priod = 12
    tax_rate = 0.0
    # 资金保管机构服务费
    stFeesRate = 0.01
    '''
