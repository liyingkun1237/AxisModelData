# 依据计算逻辑来划分块
'''
计算NPV的前端变量 NCF DCF

NPV
IRR
LOss
    Duration-prior
duration


本金偿付顺序 secPaySeq
触发事件 secTrigger
各级证券偿付 secPcpFlow
应付税费
应付利息
可分配本金（剩余本金）
IC
OC
利息分布图
券端本息
'''

import numpy as np
import pandas as pd
#####################################################


####################### 净现金流
# 优先A1档 =E91+E112+E114
def f_priorA1NCF(priorA1Interest, priorA1RepPrincipal, priorA1LagRepPrincipal):
    return priorA1Interest + priorA1RepPrincipal + priorA1LagRepPrincipal


# 优先A2档 =E92+E113
def f_priorA2NCF(priorA2Interest, priorA2RepPrincipal):
    return priorA2Interest + priorA2RepPrincipal


# B =E93+E115
def f_priorBNCF(priorBInterest, priorBRepPrincipal):
    return priorBInterest + priorBRepPrincipal


# C =E94+E116
def f_priorCNCF(priorCInterest, priorcRepPrincipal):
    return priorCInterest + priorcRepPrincipal


# 次级 =E124+E100+E117
def f_subprimeNCF(subprimeResIncome, subprimeLimitIncome, subprimeRepPrincipal):
    return subprimeResIncome + subprimeLimitIncome + subprimeRepPrincipal


# 总资产池 =SUM(E155:E159)
def f_poolNCF(priorA1NCF, priorA2NCF, priorBNCF, priorCNCF, subprimeNCF):
    return priorA1NCF + priorA2NCF + priorBNCF + priorCNCF + subprimeNCF


##################### 折现现金流

# A1     =E155/POWER((1+$D163),SUM($D$77:E$77))
def f_priorA1DCF(priorA1NCF, trchInterestRate, cashingDateGap):
    A1_list = []
    for i in range(1, len(cashingDateGap) + 1):
        value = np.power((1 + trchInterestRate[0]), cashingDateGap[:i].sum())
        A1_list.append(value)
    return priorA1NCF / np.array(A1_list)


# A2
def f_priorA2DCF(priorA2NCF, trchInterestRate, cashingDateGap):
    A2_list = []
    for i in range(1, len(cashingDateGap) + 1):
        value = np.power((1 + trchInterestRate[1]), cashingDateGap[:i].sum())
        A2_list.append(value)
    return priorA2NCF / np.array(A2_list)
    # B


def f_priorBDCF(priorBNCF, trchInterestRate, cashingDateGap):
    B_list = []
    for i in range(1, len(cashingDateGap) + 1):
        value = np.power((1 + trchInterestRate[2]), cashingDateGap[:i].sum())
        B_list.append(value)
    return priorBNCF / np.array(B_list)


# C
def f_priorCDCF(priorCNCF, trchInterestRate, cashingDateGap):  # 有问题 但是没有找到
    C_list = []
    for i in range(1, len(cashingDateGap) + 1):
        value = np.power((1 + trchInterestRate[3]), cashingDateGap[:i].sum())
        C_list.append(value)
    return priorCNCF / np.array(C_list)


# 次级
def f_subprimeDCF(subprimeNCF, trchInterestRate, cashingDateGap):
    sub_list = []
    for i in range(1, len(cashingDateGap) + 1):
        value = np.power((1 + trchInterestRate[4]), cashingDateGap[:i].sum())
        sub_list.append(value)
    return subprimeNCF / np.array(sub_list)


# 全部
def f_poolDCF(poolNCF, trchInterestRate, cashingDateGap):
    sub_list = []
    for i in range(1, len(cashingDateGap) + 1):
        value = np.power((1 + trchInterestRate[4]), cashingDateGap[:i].sum())
        sub_list.append(value)
    return poolNCF / np.array(sub_list)


# 资产池加权平均利率
def AssetpoolWavgItr(trchPtg, trchInterestRate):
    # 写死的次级折现率
    trchInterestRate[4] = 0.16
    return np.dot(trchPtg, trchInterestRate)


################### 净现值
def f_priorA1NPV(trchSize, priorA1DCF):
    return -trchSize[0] + priorA1DCF.sum()


def f_priorA2NPV(trchSize, priorA2DCF):
    return -trchSize[1] + priorA2DCF.sum()


def f_priorBNPV(trchSize, priorBDCF):
    return -trchSize[2] + priorBDCF.sum()


def f_priorCNPV(trchSize, priorCDCF):
    return -trchSize[3] + priorCDCF.sum()


def f_subprimeNPV(trchSize, subprimeDCF):
    return -trchSize[4] + subprimeDCF.sum()


def f_poolNPV(trchSize, poolDCF):
    return -sum(trchSize) + poolDCF.sum()


def f_comp_IRR(trchDCF, n, trchSize):
    '''
    n表示第0期第n档对应的PV的索引
    '''
    DCF = np.insert(trchDCF, 0, -trchSize[n])  # 插入第0期
    re = np.irr(DCF)
    return re


################ 损失率
def f_priorA1secLoss(trchSize, priorA1NPV):
    if -priorA1NPV / trchSize[0] < 0:
        res = 0
    else:
        res = -priorA1NPV / trchSize[0]
    return res


def f_priorA2secLoss(trchSize, priorA2NPV):
    if -priorA2NPV / trchSize[1] < 0:
        res = 0
    else:
        res = -priorA2NPV / trchSize[1]
    return res


def f_priorBsecLoss(trchSize, priorBNPV):
    if -priorBNPV / trchSize[2] < 0:
        res = 0
    else:
        res = -priorBNPV / trchSize[2]
    return res


def f_priorCsecLoss(trchSize, priorCNPV):
    if -priorCNPV / trchSize[3] < 0:
        res = 0
    else:
        res = -priorCNPV / trchSize[3]
    return res


def f_subprimesecLoss(trchSize, subprimeNPV):
    if -subprimeNPV / trchSize[4] < 0:
        res = 0
    else:
        res = -subprimeNPV / trchSize[4]
    return res


def f_AllsecLoss(trchSize, poolNPV):
    if -poolNPV / np.sum(trchSize) < 0:
        res = 0
    else:
        res = -poolNPV / np.sum(trchSize)
    return res


################### 加权现金流
# =E155*SUM($D$76:E$76)
def f_priorA1WAvgCF(priorA1NCF, interestDateGap):
    sum_gap = []
    for i in range(1, len(interestDateGap) + 1):
        sum_gap.append(interestDateGap[:i].sum())
    res = priorA1NCF * np.array(sum_gap)
    return res


def f_priorA2WAvgCF(priorA2NCF, interestDateGap):
    sum_gap = []
    for i in range(1, len(interestDateGap) + 1):
        sum_gap.append(interestDateGap[:i].sum())
    res = priorA2NCF * np.array(sum_gap)
    return res


def f_priorBWAvgCF(priorBNCF, interestDateGap):
    sum_gap = []
    for i in range(1, len(interestDateGap) + 1):
        sum_gap.append(interestDateGap[:i].sum())
    res = priorBNCF * np.array(sum_gap)
    return res


def f_priorCWAvgCF(priorCNCF, interestDateGap):
    sum_gap = []
    for i in range(1, len(interestDateGap) + 1):
        sum_gap.append(interestDateGap[:i].sum())
    res = priorCNCF * np.array(sum_gap)
    return res


def f_subprimeWAvgCF(subprimeNCF, interestDateGap):
    sum_gap = []
    for i in range(1, len(interestDateGap) + 1):
        sum_gap.append(interestDateGap[:i].sum())
    res = subprimeNCF * np.array(sum_gap)
    return res


def f_AllWAvgCF(poolNCF, interestDateGap):
    sum_gap = []
    for i in range(1, len(interestDateGap) + 1):
        sum_gap.append(interestDateGap[:i].sum())
    res = poolNCF * np.array(sum_gap)
    return res


################# 加权平均投资回收期
def f_priorA1Duration(priorA1WAvgCF, priorA1DCF):
    res = priorA1WAvgCF.sum() / priorA1DCF.sum()
    return res


def f_priorA2Duration(priorA2WAvgCF, priorA2DCF):
    res = priorA2WAvgCF.sum() / priorA2DCF.sum()
    return res


def f_priorBDuration(priorBWAvgCF, priorBDCF):
    res = priorBWAvgCF.sum() / priorBDCF.sum()
    return res


def f_priorCDuration(priorCWAvgCF, priorCDCF):
    res = priorCWAvgCF.sum() / priorCDCF.sum()
    return res


def f_subprimeDuration(subprimeWAvgCF, subprimeDCF):
    res = subprimeWAvgCF.sum() / subprimeDCF.sum()
    return res


def f_AllDuration(priorA1WAvgCF, priorA1DCF):
    res = priorA1WAvgCF.sum() / priorA1DCF.sum()
    return res


######################### 本金偿付顺序
def f_priorA1secPaySeq(priorA1RepPrincipal):
    res = [1 if i > 0 else 0 for i in priorA1RepPrincipal]
    return res


def f_priorA2secPaySeq(priorA2RepPrincipal):
    res = [1 if i > 0 else 0 for i in priorA2RepPrincipal]
    return res


def f_priorBsecPaySeq(priorBRepPrincipal):
    res = [1 if i > 0 else 0 for i in priorBRepPrincipal]
    return res


def f_priorCsecPaySeq(priorcRepPrincipal):
    res = [1 if i > 0 else 0 for i in priorcRepPrincipal]
    return res


def f_subprimesecPaySeq(subprimeRepPrincipal):
    res = [1 if i > 0 else 0 for i in subprimeRepPrincipal]
    return res


############################## 触发事件
# E79	是否加速清偿	isAc
# # 加速清偿
# isAc
# # 违约
# isDefault
# # 资产池剩余本金
# poolPcpBal = FinalPriBalance


# 证券剩余本金 ###需详细检查
def f_secPcpBal(priorA1ResPrincipal_DE, priorA2ResPrincipal_DE, priorBResPrincipal_DE, priorCResPrincipal_DE,
                subprimeResPrincipal_DE):
    res = priorA1ResPrincipal_DE + priorA2ResPrincipal_DE + priorBResPrincipal_DE + priorCResPrincipal_DE + subprimeResPrincipal_DE
    return res


# 资产池累计损失
def f_poolCumLoss(DefaultPrincipal):
    res = DefaultPrincipal.cumsum()
    return res


############################ 各级证券偿付
'''
E148	优先A-1档	priorA1ResPrincipal_DE  
E149	优先A-2档	priorA2ResPrincipal_DE
E150	优先B档	priorBResPrincipal_DE
E151	优先C档	priorCResPrincipal_DE
E152	次级档	subprimeResPrincipal_DE

'''


# 应付利息
def f_trchItrPayable(priorA1ResPrincipal, interestDateGap, n, trchSize, trchInterestRate):
    ls = []
    for idx in range(len(priorA1ResPrincipal)):
        if idx == 0:
            re = trchSize[n] * trchInterestRate[n] * interestDateGap[idx]
        else:
            re = priorA1ResPrincipal[idx - 1] * trchInterestRate[n] * interestDateGap[idx]

        ls.append(re)
    return np.array(ls)


def f_priorItrPayable(priorA1ItrPayable, priorA2ItrPayable, priorBItrPayable, priorCItrPayable):
    return priorA1ItrPayable + priorA2ItrPayable + priorBItrPayable + priorCItrPayable


# %%可分配收入

# DspCF 可分配收入 =E82
# priorA1DspCF	优先A-1档 =E243-E235 DspCF TaxFeesSum
# priorA2DspCF	优先A-2档 =E244-E237
# priorBDspCF	优先B档 =E245-E238
# priorCDspCF	优先C档 =E246-E239
# priorDspCF 所有优先级




def f_priorA1DspCF(DspCF, TaxFeesSum):
    return DspCF - TaxFeesSum


def f_priorA2DspCF(priorA1DspCF, priorA1ItrPayable):
    return priorA1DspCF - priorA1ItrPayable


def f_priorBDspCF(priorA2DspCF, priorA2ItrPayable):
    return priorA2DspCF - priorA2ItrPayable


def f_priorCDspCF(priorBDspCF, priorBItrPayable):
    return priorBDspCF - priorBItrPayable


def f_priorDspCF(priorA1DspCF, priorA2DspCF, priorBDspCF, priorCDspCF):
    return priorA1DspCF + priorA2DspCF + priorBDspCF + priorCDspCF


# %% IC

def trchsecIC(priorA1DspCF, priorA1ItrPayable):
    return divide(priorA1DspCF, priorA1ItrPayable)


# %% OC（本金覆盖倍数）

def divide(x, y):
    '''
    解决分母为0，则结果为0
    '''
    ls = []
    for i, j in zip(x, y):
        if j == 0:
            ls.append(0)
        else:
            ls.append(i / j)
    return np.array(ls)


# divide(np.array([3,1,2]),np.array([2,0,1]))

def f_trchsecOC(poolPcpBal, trchsecOC):
    return divide(poolPcpBal, trchsecOC)


# subprimesecPcpFlow=priorCResPrincipal_DE
# 各级证券偿付	secPcpFlow



def f_Sum4trchsecOC(priorA1secOC, priorA2secOC, priorBsecOC, priorCsecOC):
    return priorA1secOC + priorA2secOC + priorBsecOC + priorCsecOC


def f_priorsecOC(poolPcpBal, Sum4trchsecOC):
    return divide(poolPcpBal, Sum4trchsecOC)


# 证券偿付利息	secItrFlow =SUM(E91:E94)+E100
def f_secItrFlow(priorA1Interest, priorA2Interest, priorBInterest, priorCInterest, subprimeLimitIncome):
    return priorA1Interest + priorA2Interest + priorBInterest + priorCInterest + subprimeLimitIncome


# 其他费用开支	secFeeFlow =E89+E99
def f_secFeeFlow(TaxFeesSum, SurplusFees):
    return TaxFeesSum + SurplusFees


##################################输出
def f_main_output(aa, bb, cc, dd, xx, yy, zz, Input_data):
    '''
    输入为计算出来的结果 以6个字典的形式
    :return:
    '''
    trchSize = Input_data['trchSize']
    trchInterestRate = Input_data['trchInterestRate']

    priorA1Interest = bb['priorA1Interest']
    priorA2Interest = bb['priorA2Interest']
    priorBInterest = bb['priorBInterest']
    priorCInterest = bb['priorCInterest']
    priorA1RepPrincipal = bb['priorA1RepPrincipal']
    priorA2RepPrincipal = bb['priorA2RepPrincipal']
    priorA1LagRepPrincipal = bb['priorA1LagRepPrincipal']
    priorBRepPrincipal = bb['priorBRepPrincipal']
    priorCRepPrincipal = bb['priorCRepPrincipal']
    subprimeLimitIncome = bb['subprimeLimitIncome']
    subprimeRepPrincipal = bb['subprimeRepPrincipal']
    subprimeResIncome = bb['subprimeResIncome']

    cashingDateGap = zz['cashingDateGap']
    interestDateGap = zz['interestDateGap']

    priorA1NCF = f_priorA1NCF(priorA1Interest, priorA1RepPrincipal, priorA1LagRepPrincipal)
    priorA2NCF = f_priorA2NCF(priorA2Interest, priorA2RepPrincipal)
    priorBNCF = f_priorBNCF(priorBInterest, priorBRepPrincipal)
    priorCNCF = f_priorCNCF(priorCInterest, priorCRepPrincipal)
    subprimeNCF = f_subprimeNCF(subprimeResIncome, subprimeLimitIncome, subprimeRepPrincipal)

    poolNCF = f_poolNCF(priorA1NCF, priorA2NCF, priorBNCF, priorCNCF, subprimeNCF)

    priorA1DCF = f_priorA1DCF(priorA1NCF, trchInterestRate, cashingDateGap)
    priorA2DCF = f_priorA2DCF(priorA2NCF, trchInterestRate, cashingDateGap)
    priorBDCF = f_priorBDCF(priorBNCF, trchInterestRate, cashingDateGap)
    priorCDCF = f_priorCDCF(priorCNCF, trchInterestRate, cashingDateGap)
    subprimeDCF = f_subprimeDCF(subprimeNCF, trchInterestRate, cashingDateGap)
    poolDCF = f_poolDCF(poolNCF, trchInterestRate, cashingDateGap)

    priorA1NPV = f_priorA1NPV(trchSize, priorA1DCF)
    priorA2NPV = f_priorA2NPV(trchSize, priorA2DCF)
    priorBNPV = f_priorBNPV(trchSize, priorBDCF)
    priorCNPV = f_priorCNPV(trchSize, priorCDCF)
    subprimeNPV = f_subprimeNPV(trchSize, subprimeDCF)
    poolNPV = f_poolNPV(trchSize, poolDCF)

    priorA1secIRR = f_comp_IRR(priorA1DCF, 0, trchSize)
    priorA2secIRR = f_comp_IRR(priorA2DCF, 1, trchSize)
    priorBsecIRR = f_comp_IRR(priorBDCF, 2, trchSize)
    priorCsecIRR = f_comp_IRR(priorCDCF, 3, trchSize)
    subprimesecIRR = f_comp_IRR(subprimeDCF, 4, trchSize)
    AllsecIRR = np.irr(np.insert(poolDCF, 0, -np.sum(trchSize)))

    priorA1secLoss = f_priorA1secLoss(trchSize, priorA1NPV)
    priorA2secLoss = f_priorA2secLoss(trchSize, priorA2NPV)
    priorBsecLoss = f_priorBsecLoss(trchSize, priorBNPV)
    priorCsecLoss = f_priorCsecLoss(trchSize, priorCNPV)
    subprimesecLoss = f_subprimesecLoss(trchSize, subprimeNPV)
    AllsecLoss = f_AllsecLoss(trchSize, poolNPV)

    priorA1WAvgCF = f_priorA1WAvgCF(priorA1NCF, interestDateGap)
    priorA2WAvgCF = f_priorA2WAvgCF(priorA2NCF, interestDateGap)
    priorBWAvgCF = f_priorBWAvgCF(priorBNCF, interestDateGap)
    priorCWAvgCF = f_priorCWAvgCF(priorCNCF, interestDateGap)
    subprimeWAvgCF = f_subprimeWAvgCF(subprimeNCF, interestDateGap)
    AllWAvgCF = f_AllWAvgCF(poolNCF, interestDateGap)

    priorA1Duration = f_priorA1Duration(priorA1WAvgCF, priorA1DCF)
    priorA2Duration = f_priorA2Duration(priorA2WAvgCF, priorA2DCF)
    priorBDuration = f_priorBDuration(priorBWAvgCF, priorBDCF)
    priorCDuration = f_priorCDuration(priorCWAvgCF, priorCDCF)
    subprimeDuration = f_subprimeDuration(subprimeWAvgCF, subprimeDCF)
    AllDuration = f_AllDuration(AllWAvgCF, poolDCF)

    re_dict_1 = {'priorA1NPV': priorA1NPV,
                 'priorA2NPV': priorA2NPV,
                 'priorBNPV': priorBNPV,
                 'priorCNPV': priorCNPV,
                 'subprimeNPV': subprimeNPV,
                 'poolNPV': poolNPV,

                 'priorA1secIRR': priorA1secIRR,
                 'priorA2secIRR': priorA2secIRR,
                 'priorBsecIRR': priorBsecIRR,
                 'priorCsecIRR': priorCsecIRR,
                 'subprimesecIRR': subprimesecIRR,
                 'AllsecIRR': AllsecIRR,

                 'priorA1secLoss': priorA1secLoss,
                 'priorA2secLoss': priorA2secLoss,
                 'priorBsecLoss': priorBsecLoss,
                 'priorCsecLoss': priorCsecLoss,
                 'subprimesecLoss': subprimesecLoss,
                 'AllsecLoss': AllsecLoss,

                 'priorA1Duration': priorA1Duration,
                 'priorA2Duration': priorA2Duration,
                 'priorBDuration': priorBDuration,
                 'priorCDuration': priorCDuration,
                 'subprimeDuration': subprimeDuration,
                 'AllDuration': AllDuration, }

    # 本金偿付顺序
    priorA1secPaySeq = f_priorA1secPaySeq(priorA1RepPrincipal)
    priorA2secPaySeq = f_priorA2secPaySeq(priorA2RepPrincipal)
    priorBsecPaySeq = f_priorBsecPaySeq(priorBRepPrincipal)
    priorCsecPaySeq = f_priorCsecPaySeq(priorCRepPrincipal)
    subprimesecPaySeq = f_subprimesecPaySeq(subprimeRepPrincipal)

    # 加速清偿
    isAc = yy['isAc']
    # 违约
    isDefault = cc['isDefault']
    # 资产池剩余本金
    FinalPriBalance = xx['FinalPriBalance']
    poolPcpBal = FinalPriBalance

    priorA1ResPrincipal_DE = dd['priorA1ResPrincipal_DE']
    priorA2ResPrincipal_DE = dd['priorA2ResPrincipal_DE']
    priorBResPrincipal_DE = dd['priorBResPrincipal_DE']
    priorCResPrincipal_DE = dd['priorCResPrincipal_DE']
    subprimeResPrincipal_DE = dd['subprimeResPrincipal_DE']
    DefaultPrincipal = xx['DefaultPrincipal']

    secPcpBal = f_secPcpBal(priorA1ResPrincipal_DE, priorA2ResPrincipal_DE, priorBResPrincipal_DE,
                            priorCResPrincipal_DE,
                            subprimeResPrincipal_DE)

    poolCumLoss = f_poolCumLoss(DefaultPrincipal)

    re_dict_2 = {
        'priorA1secPaySeq': priorA1secPaySeq, 'priorA2secPaySeq': priorA2secPaySeq,
        'priorBsecPaySeq': priorBsecPaySeq, 'priorCsecPaySeq': priorCsecPaySeq,
        'subprimesecPaySeq': subprimesecPaySeq,
        'isAc': isAc, 'isDefault': isDefault, 'poolPcpBal': poolPcpBal,
        'secPcpBal': secPcpBal, 'poolCumLoss': poolCumLoss
    }

    priorA1secPcpFlow = priorA1ResPrincipal_DE
    priorA2secPcpFlow = priorA2ResPrincipal_DE
    priorBsecPcpFlow = priorBResPrincipal_DE
    priorCsecPcpFlow = priorCResPrincipal_DE
    subprimesecPcpFlow = priorCResPrincipal_DE

    InterestTax = aa['InterestTax']
    storageFees = aa['storageFees']
    AuditFees = aa['AuditFees']
    GradeFees = aa['GradeFees']
    ServiceFees = aa['ServiceFees']
    TaxFeesSum = aa['TaxFeesSum']

    re_dict_3 = {
        'priorA1secPcpFlow': priorA1secPcpFlow,
        'priorA2secPcpFlow': priorA2secPcpFlow,
        'priorBsecPcpFlow': priorBsecPcpFlow,
        'priorCsecPcpFlow': priorCsecPcpFlow,
        'subprimesecPcpFlow': subprimesecPcpFlow,
        'InterestTax': InterestTax, 'storageFees': storageFees,
        'AuditFees': AuditFees, 'GradeFees': GradeFees,
        'ServiceFees': ServiceFees, 'TaxFeesSum': TaxFeesSum
    }

    priorA1ResPrincipal = bb['priorA1ResPrincipal']
    priorA2ResPrincipal = bb['priorA2ResPrincipal']
    priorBResPrincipal = bb['priorBResPrincipal']
    priorCResPrincipal = bb['priorCResPrincipal']

    priorA1ItrPayable = f_trchItrPayable(priorA1ResPrincipal, interestDateGap, 0, trchSize, trchInterestRate)
    priorA2ItrPayable = f_trchItrPayable(priorA2ResPrincipal, interestDateGap, 1, trchSize, trchInterestRate)
    priorBItrPayable = f_trchItrPayable(priorBResPrincipal, interestDateGap, 2, trchSize, trchInterestRate)
    priorCItrPayable = f_trchItrPayable(priorCResPrincipal, interestDateGap, 3, trchSize, trchInterestRate)
    priorAllItrPayable = f_priorItrPayable(priorA1ItrPayable, priorA2ItrPayable, priorBItrPayable, priorCItrPayable)

    DisposableCashFlow = yy['DisposableCashFlow']
    DspCF = DisposableCashFlow
    priorA1DspCF = f_priorA1DspCF(DspCF, TaxFeesSum)
    priorA2DspCF = f_priorA2DspCF(priorA1DspCF, priorA1ItrPayable)
    priorBDspCF = f_priorBDspCF(priorA2DspCF, priorA2ItrPayable)
    priorCDspCF = f_priorCDspCF(priorBDspCF, priorBItrPayable)
    priorDspCF = f_priorDspCF(priorA1DspCF, priorA2DspCF, priorBDspCF, priorCDspCF)

    priorA1secIC = trchsecIC(priorA1DspCF, priorA1ItrPayable)
    priorA2secIC = trchsecIC(priorA2DspCF, priorA2ItrPayable)
    priorBsecIC = trchsecIC(priorBDspCF, priorBItrPayable)
    priorCsecIC = trchsecIC(priorCDspCF, priorCItrPayable)
    priorAllsecIC = trchsecIC(priorDspCF, priorAllItrPayable)

    priorA1secOC = f_trchsecOC(poolPcpBal, priorA1secPcpFlow)
    priorA2secOC = f_trchsecOC(poolPcpBal, priorA2secPcpFlow)
    priorBsecOC = f_trchsecOC(poolPcpBal, priorBsecPcpFlow)
    priorCsecOC = f_trchsecOC(poolPcpBal, priorCsecPcpFlow)
    subprimesecOC = f_trchsecOC(poolPcpBal, subprimesecPcpFlow)
    Sum4trchsecOC = f_Sum4trchsecOC(priorA1secPcpFlow, priorA2secPcpFlow, priorBsecPcpFlow, priorCsecPcpFlow)
    priorsecOC = f_priorsecOC(poolPcpBal, Sum4trchsecOC)

    re_dict_4 = {
        'priorA1ItrPayable': priorA1ItrPayable,
        'priorA2ItrPayable': priorA2ItrPayable,
        'priorBItrPayable': priorBItrPayable,
        'priorCItrPayable': priorCItrPayable,
        'priorAllItrPayable': priorAllItrPayable,
        'DspCF': DspCF,
        'priorA1DspCF': priorA1DspCF,
        'priorA2DspCF': priorA2DspCF,
        'priorBDspCF': priorBDspCF,
        'priorCsecIC': priorCsecIC,
        'priorDspCF': priorDspCF,
        'priorA1secIC': priorA1secIC, 'priorA2secIC': priorA2secIC,
        'priorBsecIC': priorBsecIC, 'priorCsecIC': priorCsecIC,
        'priorAllsecIC': priorAllsecIC,
        'priorA1secOC': priorA1secOC, 'priorA2secOC': priorA2secOC,
        'priorBsecOC': priorBsecOC, 'priorCsecOC': priorCsecOC,
        'subprimesecOC': subprimesecOC, 'priorsecOC': priorsecOC
    }

    SurplusFees = bb['SurplusFees']
    secItrFlow = f_secItrFlow(priorA1Interest, priorA2Interest, priorBInterest, priorCInterest, subprimeLimitIncome)
    secFeeFlow = f_secFeeFlow(TaxFeesSum, SurplusFees)
    # 超额利息	poolItrOver =E104
    ExcessInterest = bb['ExcessInterest']
    poolItrOver = ExcessInterest

    re_dict_5 = {'secItrFlow': secItrFlow, 'secFeeFlow': secFeeFlow, 'poolItrOver': poolItrOver}

    # %%
    # 券端利息
    priorA1secItrFlow = priorA1Interest
    priorA2secItrFlow = priorA2Interest
    priorBsecItrFlow = priorBInterest
    priorCsecItrFlow = priorCInterest
    subprimesecItrFlow = subprimeLimitIncome

    # 券端本金
    priorA1secPcpFlow_ = priorA1RepPrincipal + priorA1LagRepPrincipal
    priorA2secPcpFlow_ = priorA2RepPrincipal
    priorBsecPcpFlow_ = priorBRepPrincipal
    priorCsecPcpFlow_ = priorCRepPrincipal
    subprimesecPcpFlow_ = subprimeRepPrincipal

    re_dict_6 = {'priorA1secItrFlow': priorA1secItrFlow,
                 'priorA2secItrFlow': priorA2secItrFlow,
                 'priorBsecItrFlow': priorBsecItrFlow,
                 'priorCsecItrFlow': priorCsecItrFlow,
                 'subprimesecItrFlow': subprimesecItrFlow,
                 'priorA1secPcpFlow': priorA1secPcpFlow_,
                 'priorA2secPcpFlow': priorA2secPcpFlow_,
                 'priorBsecPcpFlow': priorBsecPcpFlow_,
                 'priorCsecPcpFlow': priorCsecPcpFlow_,
                 'subprimesecPcpFlow': subprimesecPcpFlow_,
                 }

    return {'F4': re_dict_1,  # 四个指标
            'ED': re_dict_2,  # 偿付顺序+触发事件
            'TF': re_dict_3,  # 各级证券偿付+税费
            'ICOC': re_dict_4,  # 应付利息 偿付本金 ICOC
            'ITR': re_dict_5,  # 利息分布图
            'PRP': re_dict_6,  # 券端本息
            }


def f_outputFormat4page12(res, xx, Input_data):
    InterestReceivable = Input_data['reInterest']
    PrincipalReceivable = Input_data['rePrincipal']
    compDate = Input_data['calList']
    cashingDate = Input_data['secList']
    trchSize = Input_data['trchSize']

    AssetsCashinflow_OUT = {'compDate': compDate[1:],
                            'PrincipalReceivable': PrincipalReceivable,
                            'InterestReceivable': InterestReceivable,
                            'InterestPaid': xx['InterestPaid'].tolist(),
                            'PriPaid': xx['PriPaid'].tolist()}

    def f_AnalysisResult_all(res):
        secLoss = [res['F4']['priorA1secLoss'], res['F4']['priorA2secLoss'],
                   res['F4']['priorBsecLoss'], res['F4']['priorCsecLoss'],
                   res['F4']['subprimesecLoss']
                   ]

        secIRR = [res['F4']['priorA1secIRR'], res['F4']['priorA2secIRR'],
                  res['F4']['priorBsecIRR'], res['F4']['priorCsecIRR'],
                  res['F4']['subprimesecIRR']
                  ]

        secNPV = [res['F4']['priorA1NPV'], res['F4']['priorA2NPV'],
                  res['F4']['priorBNPV'], res['F4']['priorCNPV'],
                  res['F4']['subprimeNPV']
                  ]

        secPayback = [res['F4']['priorA1Duration'], res['F4']['priorA2Duration'],
                      res['F4']['priorBDuration'], res['F4']['priorCDuration'],
                      res['F4']['subprimeDuration']
                      ]

        re = {'secName': ['优先A-1档', '优先A-2档', '优先B档', '优先C档', '次级'],
              'secLoss': secLoss,
              'secIRR': secIRR,
              'secNPV': secNPV,
              'secPayback': secPayback
              }

        return re

    AnalysisResult_all = f_AnalysisResult_all(res)

    ###本金偿付顺序	secPaySeq
    secPaySeq = {'cashingDate': cashingDate[1:],
                 'priorA1': res['ED']['priorA1secPaySeq'],
                 'priorA2': res['ED']['priorA2secPaySeq'],
                 'priorB': res['ED']['priorBsecPaySeq'],
                 'priorC': res['ED']['priorCsecPaySeq'],
                 'subprime': res['ED']['subprimesecPaySeq'], }

    secPaySeq = pd.DataFrame(secPaySeq).to_dict(orient='records')

    ###TriggerEvent 触发事件 secTrigger
    secTrigger = {
        'cashingDate': cashingDate[1:],
        'isAc': res['ED']['isAc'],
        'isDefault': res['ED']['isDefault']
    }
    secTrigger = pd.DataFrame(secTrigger).to_dict(orient='records')

    ###AssetpoolRepaymentChart 资产池偿付走势图

    AssetpoolRepaymentChart = {
        'cashingDate': cashingDate[1:],
        'poolPcpBal': res['ED']['poolPcpBal'].tolist(),
        'secPcpBal': res['ED']['secPcpBal'].tolist(),
        'poolCumLoss': res['ED']['poolCumLoss'].tolist()
    }

    ###trchBondRepChart 各级证券偿付图 secPcpFlow secPcpFlow
    # trsize=[40600000,203000000,60900000,40600000,60900000]
    secPcpFlowTr = {  # 需要转百分比
        'cashingDate': cashingDate[1:],
        'priorA1': (res['TF']['priorA1secPcpFlow'] / trchSize[0] * 100).tolist(),
        'priorA2': (res['TF']['priorA2secPcpFlow'] / trchSize[1] * 100).tolist(),
        'priorB': (res['TF']['priorBsecPcpFlow'] / trchSize[2] * 100).tolist(),
        'priorC': (res['TF']['priorCsecPcpFlow'] / trchSize[3] * 100).tolist(),
        'subprime': (res['TF']['subprimesecPcpFlow'] / trchSize[4] * 100).tolist(),
    }

    ###IC,OC走势图 IcocChart

    IcocChart = {
        'cashingDate': cashingDate[1:],
        'secIC': res['ICOC']['priorAllsecIC'].tolist(),
        'secOC': res['ICOC']['priorsecOC'].tolist()  # 进一步查看值为负的原因
    }

    ###资产池利息分布趋势图 AssetpoolInrTrendChart
    AssetpoolInrTrendChart = {
        'cashingDate': cashingDate[1:],
        'secItrFlow': res['ITR']['secItrFlow'].tolist(),
        'secFeeFlow': res['ITR']['secFeeFlow'].tolist(),
        'poolItrOver': res['ITR']['poolItrOver'].tolist()
    }

    ##A1层级
    ###预想情景分析结果 AnalysisResult_Tr


    def get_trchRe(AnalysisResult_Tr, n):
        AnalysisResult_Tr = pd.DataFrame(AnalysisResult_Tr)
        re = pd.DataFrame([AnalysisResult_Tr.iloc[n]]).to_dict(orient='records')
        return re

    ###证券现金流 BondcashFlow
    secPcpFlow = {'priorA1_secItrFlow': res['PRP']['priorA1secItrFlow'].tolist(),
                  'priorA2_secItrFlow': res['PRP']['priorA2secItrFlow'].tolist(),
                  'priorB_secItrFlow': res['PRP']['priorBsecItrFlow'].tolist(),
                  'priorC_secItrFlow': res['PRP']['priorCsecItrFlow'].tolist(),
                  'subprime_secItrFlow': res['PRP']['subprimesecItrFlow'].tolist(),
                  }

    secItrFlow = {'priorA1_secPcpFlow': res['PRP']['priorA1secPcpFlow'].tolist(),
                  'priorA2_secPcpFlow': res['PRP']['priorA2secPcpFlow'].tolist(),
                  'priorB_secPcpFlow': res['PRP']['priorBsecPcpFlow'].tolist(),
                  'priorC_secPcpFlow': res['PRP']['priorCsecPcpFlow'].tolist(),
                  'subprime_secPcpFlow': res['PRP']['subprimesecPcpFlow'].tolist(),
                  }

    def get_trch_BondcashFlow(date, secPcpFlow, secItrFlow, n):
        x = pd.DataFrame(secPcpFlow)
        y = pd.DataFrame(secItrFlow)
        re1 = pd.DataFrame(x.iloc[:, n]).values.flatten().tolist()
        re2 = pd.DataFrame(y.iloc[:, n]).values.flatten().tolist()
        return {'cashingDate': date, 'secPcpFlow': re1, 'secItrFlow': re2}


        ###本金覆盖倍数	secPcpCover

    secPcpCover = {'priorA1_secPcpCover': res['ICOC']['priorA1secOC'].tolist(),
                   'priorA2_secPcpCover': res['ICOC']['priorA2secOC'].tolist(),
                   'priorB_secPcpCover': res['ICOC']['priorBsecOC'].tolist(),
                   'priorC_secPcpCover': res['ICOC']['priorCsecOC'].tolist(),
                   'subprime_secPcpCover': res['ICOC']['subprimesecOC'].tolist()
                   }

    def get_trch_secPcpCover(date, secPcpCover, n):
        x = pd.DataFrame(secPcpFlow)
        re1 = pd.DataFrame(x.iloc[:, n]).values.flatten().tolist()
        return {'cashingDate': date, 'secPcpCover': re1}

    page_1 = {'BasicConf': {'略': '暂缺'}, 'Result': AssetsCashinflow_OUT}

    Alltrch = {'AnalysisResult_all': AnalysisResult_all,  # 表格
               'secPaySeq': secPaySeq,  # 表格
               'secTrigger': secTrigger,  # 表格
               'AssetpoolRepaymentChart': AssetpoolRepaymentChart,  # 时间  累计损失
               'secPcpFlow': secPcpFlowTr,  # 少了时间
               'IcocChart': IcocChart,  # 少了时间
               'AssetpoolInrTrendChart': AssetpoolInrTrendChart}  # 少了时间

    priorA1trch = {
        'name': 'priorA1',
        'trchRe': get_trchRe(AnalysisResult_all, 0),
        'BondcashFlow': get_trch_BondcashFlow(cashingDate[1:], secPcpFlow, secItrFlow, 0),
        'secPcpCover': get_trch_secPcpCover(cashingDate[1:], secPcpCover, 0)
    }

    priorA2trch = {
        'name': 'priorA2',
        'trchRe': get_trchRe(AnalysisResult_all, 1),
        'BondcashFlow': get_trch_BondcashFlow(cashingDate[1:], secPcpFlow, secItrFlow, 1),
        'secPcpCover': get_trch_secPcpCover(cashingDate[1:], secPcpCover, 1)}

    priorBtrch = {
        'name': 'priorB',
        'trchRe': get_trchRe(AnalysisResult_all, 2),
        'BondcashFlow': get_trch_BondcashFlow(cashingDate[1:], secPcpFlow, secItrFlow, 2),
        'secPcpCover': get_trch_secPcpCover(cashingDate[1:], secPcpCover, 2)
    }
    priorCtrch = {
        'name': 'priorC',
        'trchRe': get_trchRe(AnalysisResult_all, 3),
        'BondcashFlow': get_trch_BondcashFlow(cashingDate[1:], secPcpFlow, secItrFlow, 3),
        'secPcpCover': get_trch_secPcpCover(cashingDate[1:], secPcpCover, 3)
    }

    subprimetrch = {
        'name': 'subprime',
        'trchRe': get_trchRe(AnalysisResult_all, 4),
        'BondcashFlow': get_trch_BondcashFlow(cashingDate[1:], secPcpFlow, secItrFlow, 4),
        'secPcpCover': get_trch_secPcpCover(cashingDate[1:], secPcpCover, 4)
    }

    page_2 = {'Alltrch': Alltrch,
              'Eachtrch': [priorA1trch, priorA2trch, priorBtrch, priorCtrch, subprimetrch]}

    _return_ = {'tid': 123, 'code': 200, 'data': {'AssetsCashinflow': page_1, 'DebtCashoutflow': page_2}}

    return _return_


if __name__ == '__main__':
    from AxisModelPackage.AxisModelMain import *

    res = f_main_output(aa, bb, cc, dd, xx, yy, zz, Input_data)
    res_format = f_outputFormat4page12(res, xx, Input_data)
