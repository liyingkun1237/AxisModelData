# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 16:20:11 2017

@author: liyin
"""
# from AxisModelPackage.AxisModelInput import isReinv, rateReinv, daysReinv, \
#     InitBalance, ppByTerm, rec_term, rrAvg, \
#     InterestReceivable, days_in_year, compDate, interestDate, cashingDate, \
#     trchSize_nosub, ServiceFeesFixed, trchSize, \
#     trchInterestRate, rounding_error, trchFixedPmtAmt

###所有需要我计算的变量，计算代码进行梳理
'''
期初本金余额	InitPriBalance
违约本金	DefaultPrincipal
提前还本额	PrepaymentAmt
回收额	ReAmt
实收利息	InterestPaid
实收本金	PriPaid
再投资收益	ReinvestedEarnings
期末本金余额	FinalPriBalance

计算日间隔	compDateGap
计息日间隔	interestDateGap
兑付日间隔	cashingDateGap

是否加速清偿	isAc
收益账	IncomeAccount
本金帐	PrincipalAccount
可分配现金流	DisposableCashFlow
税费	
利息税收	InterestTax
资金保管费用	storageFees
审计费用	AuditFees
跟踪评级费	GradeFees
服务报酬（预先支付）	ServiceFees
税费总额	TaxFeesSum
利息20171006	
优先A-1档	priorA1Interest
优先A-2档	priorA2Interest
优先B档	priorBInterest
优先C档	priorCInterest
优先级利息	priorInterestSum
本金账户（正常）	PrincipalAccountNM
本金账户（加速清偿）	PrincipalAccountAC
回补本金帐	RePrincipalAccount
通道费用 & 服务报酬（剩余）	SurplusFees
次级限额内收益	subprimeLimitIncome
剩余现金流	ResCashFlow
利差	Spread
回补机制	
超额利息	ExcessInterest
本金回补	PrincipalCompensation
违约本金	DefaultPrincipal
收入应回补本金	IncomepaidBackPrincipal
实际回补	IncomeActualPrincipal

可分配余额	DisposableBalance
偿还本金	
优先A-1档	priorA1RepPrincipal
优先A-2档	priorA2RepPrincipal
优先A-1档滞后	priorA1LagRepPrincipal
优先B档	priorBRepPrincipal
优先C档	priorcRepPrincipal
次级档	subprimeRepPrincipal
剩余本金	
优先A-1档	priorA1ResPrincipal
优先A-2档	priorA2ResPrincipal
优先B档	priorBResPrincipal
优先C档	priorCResPrincipal
次级档	subprimeResPrincipal
次级剩余收益	subprimeResIncome

是否违约事件	isDefault
收益账	IncomeAccount_DE
本金帐	PrincipalAccount_DE
可分配现金流	DisposableCashFlow_DE
税费	
利息税收	InterestTax_DE
资金保管费用	storageFees_DE
审计费用	AuditFees_DE
跟踪评级费	GradeFees_DE
服务报酬	SurplusFees_DE
税费总额	TaxFeesSum_DE
优先A-1档利息	priorA1Interest_DE
优先A-1档本金	priorA1RepPrincipal_DE
优先A-2档利息	priorA2Interest_DE
优先A-2档本金	priorA2RepPrincipal_DE
优先B档利息	priorBInterest_DE
优先B档本金	priorBRepPrincipal_DE
优先C档利息	priorCInterest_DE
优先C档本金	priorCRepPrincipal_DE
次级档本金	subprimeRepPrincipal_DE
次级档收益	subprimeRepIncome_DE
剩余本金	
优先A-1档	priorA1ResPrincipal_DE
优先A-2档	priorA2ResPrincipal_DE
优先B档	priorBResPrincipal_DE
优先C档	priorCResPrincipal_DE
次级档	subprimeResPrincipal_DE
'''

import pandas as pd
import numpy as np


# # %% 输入
# ################################计算######################
# # 期初余额
# InitBalance = np.array(
#     [426362307.65, 364477464.71, 313128869.61, 261310447.39, 213723414.57, 169917241.19, 125290945.11, 84072200.63,
#      56606734.62, 36856727.50, 17441527.32, 5040461.64])
#
# # 应收利息
# InterestReceivable = np.array(
#     [11746513.05, 8128351.61, 6874150.94, 5642326.58, 4546467.03, 3495779.20, 2436239.51, 1568041.01, 1045808.03,
#      631032.24, 254493.40, 49903.90])
#
# # 应收本金
# PrincipalReceivable = np.array(
#     [61884842.94, 51348595.10, 51818422.22, 47587032.82, 43806173.38, 44626296.08, 41218744.48, 27465466.01,
#      19750007.12, 19415200.18, 12401065.68, 5040461.64])


# 小数点精确装饰器
def _round(n):
    def _round_(func):
        def inner(*args, **kwargs):  # 1
            return np.round(func(*args, **kwargs), n)

        return inner

    return _round_


# %% 收入项
'''
期初本金余额	InitPriBalance
违约本金	DefaultPrincipal
提前还本额	PrepaymentAmt
回收额	ReAmt
实收利息	InterestPaid
实收本金	PriPaid
再投资收益	ReinvestedEarnings
期末本金余额	FinalPriBalance
'''


# 期初本金余额	InitPriBalance
@_round(2)
def f_InitPriBalance(InitBalance, FinalPriBalance):
    """
    计算逻辑，期初本金余额等于上一期的期末本金余额
    而第0期的期末本金余额等于期初余额的第一个元素值
    """
    ls = []
    ls[0] = InitBalance[0]
    for i in FinalPriBalance:
        ls.append(i)
    return np.array(ls)


# 违约本金	DefaultPrincipal
# 违约本金 = 期初本金余额 * 按期违约分布
@_round(2)
def f_DefaultPrincipal(InitPriBalance, pdByTerm):
    return InitPriBalance * pdByTerm


# 提前还本额	PrepaymentAmt
# 提前还本额= (期初本金余额 - 违约本金) * 按期早偿分布
@_round(2)
def f_PrepaymentAmt(InitBalance, DefaultPrincipal, ppByTerm):
    return (InitBalance - DefaultPrincipal) * ppByTerm


# 回收额	ReAmt
# 回收时间（期）
# rec_term = 1
# term = np.arange(12) + 1


@_round(2)
def f_ReAmt(term, DefaultPrincipal, rrAvg, rec_term):
    ls = []
    for idx in np.arange(len(term)):
        if term[idx] <= rec_term:
            ls.append(0)
        else:
            ls.append(DefaultPrincipal[idx - rec_term] * rrAvg)
    return np.array(ls)


# 实收利息	InterestPaid
# 实收利息=应收利息/期初余额 * 实际期初余额
# 实际期初余额=期初本金-违约本金
@_round(2)
def f_InterestPaid(InterestReceivable, InitBalance, InitPriBalance, DefaultPrincipal):
    return InterestReceivable / InitBalance * (InitPriBalance - DefaultPrincipal)


# 实收本金	PriPaid
# 实收本金=(E66-E68-E67)/E61*E63+E68+E69
# 逻辑
# （期初本金余额-提前还本额-违约本金）/期初余额*应收本金+提前还本额+回收额
@_round(2)
def f_PriPaid(InitPriBalance, PrepaymentAmt, DefaultPrincipal, InitBalance, PrincipalReceivable, ReAmt):
    return (
               InitPriBalance - PrepaymentAmt - DefaultPrincipal) / InitBalance * PrincipalReceivable + PrepaymentAmt + ReAmt


# 再投资收益	ReinvestedEarnings
# 再投资收益
# isReinv = 0
# rateReinv = 0.01
# daysReinv = 10
# days_in_year = 365


def f_ReinvestedEarnings(isReinv, InterestPaid, PriPaid, rateReinv, daysReinv, days_in_year):
    if isReinv == 0:
        return np.zeros(len(InterestPaid))
    elif isReinv == 1:
        return (InterestPaid + PriPaid) * rateReinv * daysReinv / days_in_year


# 期末本金余额	FinalPriBalance
# =E66-E67-E71-E68+E69
# 期末本金余额 = 期初本金余额 - 违约本金 -实收本金-提前还本额 + 回收额  #1012周岳更正
@_round(2)
def f_FinalPriBalance(InitPriBalance, DefaultPrincipal, ReAmt, PriPaid, PrepaymentAmt):
    return InitPriBalance - DefaultPrincipal - PriPaid + ReAmt - PrepaymentAmt


# 收入项 对上面的函数的一个封装

def f_main_IncomeAccount(PrincipalReceivable, term, pdByTerm, rrAvg, rec_term, ppByTerm, InterestReceivable,
                         InitBalance, rateReinv, daysReinv, days_in_year, isReinv):
    # 初始化
    InitPriBalance = np.zeros(priod)
    DefaultPrincipal = np.zeros(priod)
    PrepaymentAmt = np.zeros(priod)
    ReAmt = np.zeros(priod)
    InterestPaid = np.zeros(priod)
    PriPaid = np.zeros(priod)
    FinalPriBalance = np.zeros(priod)
    ReinvestedEarnings = np.zeros(priod)

    for idx in range(len(InitPriBalance)):
        # 第一期
        if idx == 0:
            # 期初本金余额 第1期(idx==0)时，为期初余额第一个元素值
            InitPriBalance[idx] = InitBalance[0]
            # 违约本金=期初本金余额 * 按期违约分布
            DefaultPrincipal = f_DefaultPrincipal(InitPriBalance, pdByTerm)
            # 提前还本额= (期初本金余额 - 违约本金) * 按期早偿分布
            PrepaymentAmt = f_PrepaymentAmt(InitBalance, DefaultPrincipal, ppByTerm)
            # 回收额=判断期限和回收时间的关系，期限小于等于回收时间，回收额为0 否则，以当期违约本金为参考，向下移动0个单位，向左移动rec_term个单位的取值乘以违约后回收率
            ReAmt = f_ReAmt(term, DefaultPrincipal, rrAvg, rec_term)
            # 实收利息=应收利息/期初余额 * （期初本金-违约本金）
            InterestPaid = f_InterestPaid(InterestReceivable, InitBalance, InitPriBalance, DefaultPrincipal)
            # 实收本金=（期初本金余额-提前还本额-违约本金）/期初余额*应收本金+提前还本额+回收额
            PriPaid = f_PriPaid(InitPriBalance, PrepaymentAmt, DefaultPrincipal, InitBalance, PrincipalReceivable,
                                ReAmt)
            # 再投资收益 ？？
            # 期末本金余额=期初本金余额 - 违约本金 -实收本金 + 回收额
            FinalPriBalance = f_FinalPriBalance(InitPriBalance, DefaultPrincipal, ReAmt, PriPaid, PrepaymentAmt)
        else:
            InitPriBalance[idx] = FinalPriBalance[idx - 1]
            DefaultPrincipal = f_DefaultPrincipal(InitPriBalance, pdByTerm)
            PrepaymentAmt = f_PrepaymentAmt(InitBalance, DefaultPrincipal, ppByTerm)
            ReAmt = f_ReAmt(term, DefaultPrincipal, rrAvg, rec_term)
            InterestPaid = f_InterestPaid(InterestReceivable, InitBalance, InitPriBalance, DefaultPrincipal)
            PriPaid = f_PriPaid(InitPriBalance, PrepaymentAmt, DefaultPrincipal, InitBalance, PrincipalReceivable,
                                ReAmt)
            FinalPriBalance = f_FinalPriBalance(InitPriBalance, DefaultPrincipal, ReAmt, PriPaid, PrepaymentAmt)
            ReinvestedEarnings = f_ReinvestedEarnings(isReinv, InterestPaid, PriPaid, rateReinv, daysReinv,
                                                      days_in_year)  # InterestPaid, PriPaid

    return {'InitPriBalance': InitPriBalance,  # 期初本金余额
            'DefaultPrincipal': DefaultPrincipal,  # 违约本金
            'PrepaymentAmt': PrepaymentAmt,  # 提前还本额
            'ReAmt': ReAmt,  # 回收额
            'InterestPaid': InterestPaid,  # 实收利息
            'PriPaid': PriPaid,  # 实收本金	PriPaid
            'ReinvestedEarnings': ReinvestedEarnings,  # 再投资收益	ReinvestedEarnings
            'FinalPriBalance': FinalPriBalance  # 期末本金余额	FinalPriBalance
            }


# %% 支出项 需要的输入
# # 计算日
# compDate = ['2017/4/13', '2017/7/31', '2017/10/31', '2018/1/31', '2018/4/30', '2018/7/31', '2018/10/31', '2019/1/31',
#             '2019/4/30', '2019/7/31', '2019/10/31', '2020/1/31', '2020/3/31']
# # 计息日
# interestDate = ['2017/5/18', '2017/8/15', '2017/11/15', '2018/2/15', '2018/5/15', '2018/8/15', '2018/11/15',
#                 '2019/2/15', '2019/5/15', '2019/8/15', '2019/11/15', '2020/2/15', '2020/3/31']
# # 兑付日
# cashingDate = ['2017/5/18', '2017/8/15', '2017/11/15', '2018/2/15', '2018/5/15', '2018/8/15', '2018/11/15', '2019/2/15',
#                '2019/5/15', '2019/8/15', '2019/11/15', '2020/2/15', '2020/3/31']
# days_in_year = 365


# %% 支出项计算

# 计算日间隔	compDateGap
# 计算日间隔
@_round(4)
def f_gapDate(compDate, days_in_year):
    compDate = pd.to_datetime(compDate)
    gap = []
    for idx in np.arange(len(compDate) - 1):
        gap.append((compDate[idx + 1] - compDate[idx]).days / days_in_year)
    return gap


# 计息日间隔	interestDateGap
# 兑付日间隔	cashingDateGap


def f_main_OutcomeAccount(compDate, interestDate, cashingDate, days_in_year):
    # 计算日间隔	compDateGap
    compDateGap = f_gapDate(compDate, days_in_year)
    # 计息日间隔	interestDateGap
    interestDateGap = f_gapDate(interestDate, days_in_year)
    # 兑付日间隔	cashingDateGap
    cashingDateGap = f_gapDate(cashingDate, days_in_year)
    return {'compDateGap': compDateGap,
            'interestDateGap': interestDateGap,
            'cashingDateGap': cashingDateGap}


# %%  收益账户

# 是否加速清偿	isAc
# 是否加速清偿
# 逻辑：是否加速清偿的前一项是否为1，或者累计违约本金除以第1期 期初本金余额超过给定的加速清偿阈值，即判断为1，否则为0
def f_isAc(DefaultPrincipal, InitPriBalance, delinq_acc_rate=0.065):
    initFlag = 0
    pdFlag = [i for i in DefaultPrincipal.cumsum() / InitPriBalance[0] > delinq_acc_rate]
    ls = []
    for i in pdFlag:
        if initFlag == 1 or i:
            ls.append(1)
        else:
            ls.append(0)
        initFlag = ls[-1]
    return ls


# 收益账	IncomeAccount
def f_IncomeAccount(InterestPaid, ReinvestedEarnings):
    return InterestPaid + ReinvestedEarnings


# 本金帐	PrincipalAccount
# 本金帐=实收本金
# PrincipalAccount=PriPaid

# 可分配现金流	DisposableCashFlow
# 可分配现金流=收益账+本金帐
def f_DisposableCashFlow(IncomeAccount, PrincipalAccount):
    return IncomeAccount + PrincipalAccount


def f_main_ProfitAccount_1(DefaultPrincipal, InitPriBalance, InterestPaid, ReinvestedEarnings, PriPaid):
    # 是否加速清偿	isAc
    isAc = f_isAc(DefaultPrincipal, InitPriBalance, delinq_acc_rate=0.065)
    # 收益账	IncomeAccount
    IncomeAccount = f_IncomeAccount(InterestPaid, ReinvestedEarnings)
    # 本金帐	PrincipalAccount
    PrincipalAccount = PriPaid
    # 可分配现金流	DisposableCashFlow
    DisposableCashFlow = f_DisposableCashFlow(IncomeAccount, PrincipalAccount)

    return {'isAc': isAc, 'IncomeAccount': IncomeAccount,
            'PrincipalAccount': PrincipalAccount,
            'DisposableCashFlow': DisposableCashFlow}


# %% 税费

# 利息税收	InterestTax
# tax_rate = 0.0


def f_InterestTax(DisposableCashFlow, IncomeAccount):
    ls = []
    for i, j in zip(DisposableCashFlow, IncomeAccount):
        ls.append(np.max([np.min([i, j * tax_rate]), 0.0]))
    return np.array(ls)


# 资金保管费用	storageFees
# 预期发行规模（除次级档外）
# trchSize_nosub = trchSize.cumsum()[-2]
# 资金保管服务费的年化费率 st:storage
# stFeesRate = 0.01


def f_storageFees(DisposableCashFlow, InterestTax, term, interestDateGap, trchSize_nosub):
    ls = []
    for a, b, c, d in zip(DisposableCashFlow, InterestTax, term, interestDateGap):
        if c == 1:
            _ = d * trchSize_nosub * stFeesRate
        elif c % 4 == 0:
            _ = trchSize_nosub * stFeesRate  # 第四期时为什么不乘计息日间隔
        else:
            _ = 0.0
        ls.append(np.max([np.min([a - b, _]), 0.0]))
    return np.array(ls)


# 审计费用	AuditFees

# 审计费用=MAX(MIN(E82-SUM(E84:E85),IF(MOD(E1,4)=0,$H$45,0)),0)
# 逻辑翻译： E82:可支配现金流 E84:利息税收 E85:资金保管费用 $H$45:审计费用的支付费用
# AuditFeesFixed = 50000


def f_AuditFees(DisposableCashFlow, InterestTax, storageFees, term):
    ls = []
    for a, b, c, d in zip(DisposableCashFlow, InterestTax, storageFees, term):
        if d % 4 == 0:
            _ = AuditFeesFixed
        else:
            _ = 0
        ls.append(np.max([np.min([a - b - c, _]), 0.0]))
    return np.array(ls)


# 跟踪评级费	GradeFees
GradeFeesFixed = 50000


def f_GradeFees(DisposableCashFlow, InterestTax, storageFees, AuditFees, term):
    ls = []
    for a, b, c, d, e in zip(DisposableCashFlow, InterestTax, storageFees, AuditFees, term):
        if e % 4 == 0:
            _ = GradeFeesFixed
        else:
            _ = 0
        ls.append(np.max([np.min([a - b - c - d, _]), 0]))
    return np.array(ls)


# 服务报酬（预先支付）	ServiceFees
# 服务报酬（预先支付）=MAX(MIN(E82-SUM(E84:E87),IF(E1=1,$H$47,IF(MOD(E1,4)=0,$H$47,0))),0)
# 逻辑翻译： E82:可支配现金流 E84:利息税收 E85:资金保管费用 E86：审计费用 E87：评级费用 $H$47：服务费的支付费用
# ServiceFeesFixed = 100000


def f_ServiceFees(DisposableCashFlow, InterestTax, storageFees, AuditFees, GradeFees, term, ServiceFeesFixed):
    ls = []
    for a, b, c, d, e, f in zip(DisposableCashFlow, InterestTax, storageFees, AuditFees, GradeFees, term):
        if f == 1:
            _ = ServiceFeesFixed
        elif f % 4 == 0:
            _ = ServiceFeesFixed
        else:
            _ = 0
        ls.append(np.max([np.min([a - (b + c + d + e), _]), 0]))
    return np.array(ls)


# 税费总额	TaxFeesSum
# 税费总额
def f_TaxFeesSum(InterestTax, storageFees, AuditFees, GradeFees, ServiceFees):
    return InterestTax + storageFees + AuditFees + GradeFees + ServiceFees


# TaxFeesSum=f_TaxFeesSum(InterestTax,storageFees,AuditFees,GradeFees,ServiceFees)


def f_main_TaxFeesNMAC(DisposableCashFlow, IncomeAccount, interestDateGap, term, trchSize_nosub, ServiceFeesFixed):
    # 税费
    # 利息税收	InterestTax
    InterestTax = f_InterestTax(DisposableCashFlow, IncomeAccount)
    # 资金保管费用	storageFees
    storageFees = f_storageFees(DisposableCashFlow, InterestTax, term, interestDateGap, trchSize_nosub)
    # 审计费用	AuditFees
    AuditFees = f_AuditFees(DisposableCashFlow, InterestTax, storageFees, term)
    # 跟踪评级费	GradeFees
    GradeFees = f_GradeFees(DisposableCashFlow, InterestTax, storageFees, AuditFees, term)
    # 服务报酬（预先支付）	ServiceFees
    ServiceFees = f_ServiceFees(DisposableCashFlow, InterestTax, storageFees, AuditFees, GradeFees, term,
                                ServiceFeesFixed)
    # 税费总额	TaxFeesSum
    TaxFeesSum = f_TaxFeesSum(InterestTax, storageFees, AuditFees, GradeFees, ServiceFees)

    return {
        'InterestTax': InterestTax,
        'storageFees': storageFees,
        'AuditFees': AuditFees,
        'GradeFees': GradeFees,
        'ServiceFees': ServiceFees,
        'TaxFeesSum': TaxFeesSum
    }


# %% 收益帐复杂部分 输入
# trchSize = f_trchSize(poolSize, trchPtg)
# 付息频率
# trchPayFreq = np.array([12, 12, 12, 12, 12])
# 利率类型
# trchInterestType = ['固定'] * 5
# 预期发行利率
# trchInterestRate = np.array([4.00, 5.00, 7.00, 9.00, 3]) / 100
'''
利息20171006	
优先A-1档	priorA1Interest
优先A-2档	priorA2Interest
优先B档	priorBInterest
优先C档	priorCInterest
优先级利息	priorInterestSum
本金账户（正常）	PrincipalAccountNM
本金账户（加速清偿）	PrincipalAccountAC
回补本金帐	RePrincipalAccount
通道费用 & 服务报酬（剩余）	SurplusFees
次级限额内收益	subprimeLimitIncome
剩余现金流	ResCashFlow
利差	Spread
回补机制	
超额利息	ExcessInterest
本金回补	PrincipalCompensation
违约本金	DefaultPrincipal
收入应回补本金	IncomepaidBackPrincipal
实际回补	IncomeActualPrincipal

可分配余额	DisposableBalance
偿还本金	
优先A-1档	priorA1RepPrincipal
优先A-2档	priorA2RepPrincipal
优先A-1档滞后	priorA1LagRepPrincipal
优先B档	priorBRepPrincipal
优先C档	priorcRepPrincipal
次级档	subprimeRepPrincipal
剩余本金	
优先A-1档	priorA1ResPrincipal
优先A-2档	priorA2ResPrincipal
优先B档	priorBResPrincipal
优先C档	priorCResPrincipal
次级档	subprimeResPrincipal
次级剩余收益	subprimeResIncome
'''


# %% 收益帐 复杂部分 计算
# 利息 优先A-1档 = =MAX(MIN(E82-E89,D119*$K31*E$76),0)
# 逻辑翻译 E82:可支配现金流 E89：税费总额  D119：优先A-1档剩余本金 $K31：A-1预期发行利率 E$76：计息日间隔 （为什么锁定？？）
def f_priorA1Interest_0(DisposableCashFlow, TaxFeesSum, interestDateGap, trchSize, trchInterestRate):
    ls = []
    for a, b, c in zip(DisposableCashFlow, TaxFeesSum, interestDateGap):
        re = np.max([np.min([a - b, trchSize[0] * trchInterestRate[0] * c]), 0])
        ls.append(re)
        break
    return np.array(ls)


# 逻辑思路：发现第一期是对的，后面皆不对，原因为剩余本金是一个变化的值，且和索引相关
# 剩余本金 优先A-1档 Residual priorA1ResPrincipal
def f_priorA1Interest(DisposableCashFlow, TaxFeesSum, interestDateGap, priorA1ResPrincipal, trchInterestRate):
    ls = []
    idx = np.arange(1, len(DisposableCashFlow), 1)
    for a, b, c, d in zip(DisposableCashFlow[idx], TaxFeesSum[idx], interestDateGap[idx], priorA1ResPrincipal[:-1]):
        re = np.max([np.min([a - b, d * trchInterestRate[0] * c]), 0])
        ls.append(re)
    return np.array(ls)


# 利息 优先A-2档 =MAX(MIN(E82-SUM(E89-E91),D120*$K32*E$76),0)
def f_priorA2Interest_0(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, trchSize, trchInterestRate):
    ls = []
    for a, b, c, d in zip(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap):
        re = np.max([np.min([a - b - c, trchSize[1] * trchInterestRate[1] * d]), 0])
        ls.append(re)
        break
    return np.array(ls)


def f_priorA2Interest(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, priorA2ResPrincipal,
                      trchInterestRate):
    ls = []
    idx = np.arange(1, len(DisposableCashFlow), 1)
    for a, b, c, d, e in zip(DisposableCashFlow[idx], TaxFeesSum[idx], interestDateGap[idx], priorA2ResPrincipal[:-1],
                             priorA1Interest[idx]):
        re = np.max([np.min([a - b - e, d * trchInterestRate[1] * c]), 0])
        ls.append(re)
    return np.array(ls)


# 利息 优先B档 =MAX(MIN(E82-E89-SUM(E91:E92),D121*$K33*E$76),0)
def f_priorBInterest_0(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, priorA2Interest, trchSize,
                       trchInterestRate):
    ls = []
    for a, b, c, d, e in zip(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, priorA2Interest):
        re = np.max([np.min([a - b - c - e, trchSize[2] * trchInterestRate[2] * d]), 0])
        ls.append(re)
        break
    return np.array(ls)


def f_priorBInterest(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, priorBResPrincipal,
                     priorA2Interest, trchInterestRate):
    ls = []
    idx = np.arange(1, len(DisposableCashFlow), 1)
    for a, b, c, d, e, f in zip(DisposableCashFlow[idx], TaxFeesSum[idx], interestDateGap[idx], priorBResPrincipal[:-1],
                                priorA1Interest[idx], priorA2Interest[idx]):
        re = np.max([np.min([a - b - e - f, d * trchInterestRate[2] * c]), 0])
        ls.append(re)
    return np.array(ls)


# 利息 优先C档
def f_priorCInterest_0(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, priorA2Interest,
                       priorBInterest, trchSize, trchInterestRate):
    ls = []
    for a, b, c, d, e, f in zip(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, priorA2Interest,
                                priorBInterest):
        re = np.max([np.min([a - b - c - e - f, trchSize[3] * trchInterestRate[3] * d]), 0])
        ls.append(re)
        break
    return np.array(ls)


def f_priorCInterest(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, priorCResPrincipal,
                     priorA2Interest, priorBInterest, trchInterestRate):
    ls = []
    idx = np.arange(1, len(DisposableCashFlow), 1)
    for a, b, c, d, e, f, i in zip(DisposableCashFlow[idx], TaxFeesSum[idx], interestDateGap[idx],
                                   priorCResPrincipal[:-1], priorA1Interest[idx], priorA2Interest[idx],
                                   priorBInterest[idx]):
        re = np.max([np.min([a - b - e - f - i, d * trchInterestRate[3] * c]), 0])
        ls.append(re)
    return np.array(ls)


# priorA1Interest_0=f_priorA1Interest_0(DisposableCashFlow,TaxFeesSum,interestDateGap)
# priorA1Interest_1=f_priorA1Interest(DisposableCashFlow,TaxFeesSum,interestDateGap,priorA1ResPrincipal)



# 本金账户（正常）==IF(E79=0,MAX(MIN(E81,E82-E89-E95),0),0)
# E79是否加速清偿 E81本金帐 E82可分配现金流 E89税费总和 E95优先级利息总和
def f_PrincipalAccountNM(isAc, PrincipalAccount, DisposableCashFlow, TaxFeesSum, priorInterestSum):
    ls = []
    for a, b, c, d, e in zip(isAc, PrincipalAccount, DisposableCashFlow, TaxFeesSum, priorInterestSum):
        if a == 0:
            ls.append(np.max([np.min([b, c - d - e]), 0]))
        else:
            ls.append(0)
    return np.array(ls)


# PrincipalAccountNM
# 本金账户（加速清偿）=IF(E79=1,MAX(0,E82-E89-E95),0)
def f_PrincipalAccountAC(isAc, DisposableCashFlow, TaxFeesSum, priorInterestSum):
    ls = []
    for a, b, c, d in zip(isAc, DisposableCashFlow, TaxFeesSum, priorInterestSum):
        if a == 1:
            ls.append(np.max([b - c - d, 0]))
        else:
            ls.append(0)
    return np.array(ls)


# PrincipalAccountAC


# 回补本金帐=MAX(MIN(E82-E89-E95-SUM(E96:E97),E107),0)
# E82可分配现金流 - E89税费总和 -E95优先级利息总和 E96:E97本金账户（正常）,E107收入应回补本金 IncomePaidPrincipal

def f_RePrincipalAccount(DisposableCashFlow, TaxFeesSum, priorInterestSum, PrincipalAccountNM, PrincipalAccountAC,
                         IncomePaidPrincipal):
    ls = []
    for a, b, c, d, e, f in zip(DisposableCashFlow, TaxFeesSum, priorInterestSum, PrincipalAccountNM,
                                PrincipalAccountAC, IncomePaidPrincipal):
        ls.append(np.max([np.min([a - b - c - d - e, f]), 0]))
    return np.array(ls)


# RePrincipalAccount

# 通道费用 & 服务报酬（剩余）=IF(D123>0,MAX(MIN(E82-E89-E95-E98,E66*($G$48+$G$49)*E$76),0),0)
# D123次级档剩余本金 subprimeResPrincipal  E98回补本金帐 E66期初本金余额 $G$48+$G$49通道费用和剩余服务费用的总计 E$76计息日间隔

# SurplusFeesRate = 0.0015 + 0.0005


def f_SurplusFees(subprimeResPrincipal, DisposableCashFlow, TaxFeesSum, priorInterestSum, RePrincipalAccount,
                  InitPriBalance, interestDateGap, trchSize):
    ls = []
    for idx in range(len(subprimeResPrincipal)):
        if idx == 0:
            if trchSize[4] > 0:
                _ = DisposableCashFlow[idx] - TaxFeesSum[idx] - priorInterestSum[idx] - RePrincipalAccount[idx]
                __ = InitPriBalance[idx] * SurplusFeesRate * interestDateGap[idx]
                ls.append(np.max([np.min([_, __]), 0]))
            else:
                ls.append(0)
        else:
            if subprimeResPrincipal[idx - 1] > 0:
                _ = DisposableCashFlow[idx] - TaxFeesSum[idx] - priorInterestSum[idx] - RePrincipalAccount[idx]
                __ = InitPriBalance[idx] * SurplusFeesRate * interestDateGap[idx]
                ls.append(np.max([np.min([_, __]), 0]))
            else:
                ls.append(0)
    return np.array(ls)


# 次级限额内收益=MAX(MIN(E82-E89-SUM(E95:E99),D123*$K$35*E76),0)
# E82可分配现金流 - E89税费总和 E95:E99 优先级利息 本金账户（正常）本金账户（加速清偿）回补本金帐 通道费用 & 服务报酬（剩余）
# D123次级档剩余本金 subprimeResPrincipal $K$35:次级档预期发行利率 E$76计息日间隔
def f_sum4subprimeLimitIncome(priorInterestSum, PrincipalAccountNM, PrincipalAccountAC, RePrincipalAccount,
                              SurplusFees):
    return priorInterestSum + PrincipalAccountNM + PrincipalAccountAC + RePrincipalAccount + SurplusFees


def f_subprimeLimitIncome(DisposableCashFlow, TaxFeesSum, sum4subprimeLimitIncome, subprimeResPrincipal,
                          interestDateGap, trchSize, trchInterestRate):
    ls = []
    for idx in range(len(DisposableCashFlow)):
        if idx == 0:
            _ = DisposableCashFlow[idx] - TaxFeesSum[idx] - sum4subprimeLimitIncome[idx]
            __ = trchSize[4] * trchInterestRate[4] * interestDateGap[idx]
            re = np.max([np.min([_, __]), 0])
        else:
            _ = DisposableCashFlow[idx] - TaxFeesSum[idx] - sum4subprimeLimitIncome[idx]
            __ = subprimeResPrincipal[idx - 1] * trchInterestRate[4] * interestDateGap[idx]
            re = np.max([np.min([_, __]), 0])
        ls.append(re)
    return np.array(ls)


# subprimeLimitIncome



# 剩余现金流=E82-E89-SUM(E95:E100)
# E82可分配现金流 - E89税费总和 E95:E99 优先级利息 本金账户（正常）本金账户（加速清偿）回补本金帐 通道费用 & 服务报酬（剩余）次级限额内收益

def f_sum4ResCashFlow(sum4subprimeLimitIncome, subprimeLimitIncome):
    return sum4subprimeLimitIncome + subprimeLimitIncome


def f_ResCashFlow(DisposableCashFlow, TaxFeesSum, sum4ResCashFlow):
    return DisposableCashFlow - TaxFeesSum - sum4ResCashFlow


# ResCashFlow



# 利差=MAX(E80-E89-E95-SUM(E98:E100),0)
# E80收益帐 -E89税费总和 -E95优先级利息  E98:E100 回补本金帐 通道费用 & 服务报酬（剩余）次级限额内收益

def f_Spread(IncomeAccount, TaxFeesSum, priorInterestSum, RePrincipalAccount, SurplusFees, subprimeLimitIncome):
    re = IncomeAccount - TaxFeesSum - priorInterestSum - RePrincipalAccount - SurplusFees - subprimeLimitIncome
    re[re < 0] = 0
    return re


# Spread

# 超额利息 =E80-E89-E95-SUM(E99:E100)
def f_ExcessInterest(IncomeAccount, TaxFeesSum, priorInterestSum, SurplusFees, subprimeLimitIncome):
    re = IncomeAccount - TaxFeesSum - priorInterestSum - SurplusFees - subprimeLimitIncome
    return re


# ExcessInterest

# 本金回补 =MAX(0,E81-SUM(E96:E97))
# E81本金帐 E96:E97本金账号（正常 加速）
def f_PrincipalCompensation(PrincipalAccount, PrincipalAccountNM, PrincipalAccountAC):
    re = PrincipalAccount - PrincipalAccountNM - PrincipalAccountAC
    re[re < 0] = 0
    return re


# PrincipalCompensation



# 违约本金=E67
# DefaultPrincipal
# DefaultPrincipal

# 收入应回补本金=MAX(0,D107-D108+SUM(E105:E106))
# E105本金回补 PrincipalCompensation E106违约本金 DefaultPrincipal
# D107收入应回补本金 IncomepaidBackPrincipal D108实际回补 IncomeActualPrincipal

IncomepaidBackPrincipal_init = 0
IncomeActualPrincipal_init = 0


def f_IncomepaidBackPrincipal(PrincipalCompensation, DefaultPrincipal, IncomepaidBackPrincipal, IncomeActualPrincipal):
    ls = []
    for idx in range(len(PrincipalCompensation)):
        if idx == 0:
            re = IncomepaidBackPrincipal_init - IncomeActualPrincipal_init + PrincipalCompensation[idx] + \
                 DefaultPrincipal[idx]
            if re < 0:
                re = 0
            ls.append(re)
        else:
            re = IncomepaidBackPrincipal[idx - 1] - IncomeActualPrincipal[idx - 1] + PrincipalCompensation[idx] + \
                 DefaultPrincipal[idx]
            if re < 0:
                re = 0
            ls.append(re)
    return np.array(ls)


# IncomepaidBackPrincipal

# 实际回补==E98 回补本金帐RePrincipalAccount
# IncomeActualPrincipal=RePrincipalAccount


# %%本金账户
# 可分配余额=SUM(E96:E97)+E98+E101
# E96本金正常 E97本金加速 E98回补本金帐RePrincipalAccount E101剩余现金流
def f_DisposableBalance(PrincipalAccountNM, PrincipalAccountAC, RePrincipalAccount, ResCashFlow):
    return PrincipalAccountNM + PrincipalAccountAC + RePrincipalAccount + ResCashFlow


# DisposableBalance


'''
偿还本金 RepPrincipal
优先A-1档 priorA1RepPrincipal
优先A-2档 priorA2RepPrincipal
优先A-1档滞后 priorA1LagRepPrincipal lagging--> Lag
优先B档 priorBRepPrincipal
优先C档 priorCRepPrincipal
次级档 subprimeRepPrincipal

剩余本金 ResPrincipal
优先A-1档 priorA1ResPrincipal
优先A-2档 priorA2ResPrincipal
优先B档 priorBResPrincipal
优先C档 priorCResPrincipal
次级档 subprimeResPrincipal
'''


# 优先A-1档 priorA1RepPrincipal
# =IF(D120<rounding_error,D119,IF(E79=0,E37,MAX(0,MIN(E110*D119/SUM(D119:D120),D119))))
# D120优先A-2档剩余本金 rounding_error归整误差 D119优先A-1档剩余本金
# E79是否加速清偿 E37优先A-1档固定摊还金额 E110可分配余额
# rounding_error = 0.001


def f_priorA1RepPrincipal(priorA2ResPrincipal, priorA1ResPrincipal, isAc, trchFixedPmtAmt, DisposableBalance, trchSize,
                          rounding_error):
    ls = []
    for idx in range(len(isAc)):
        if idx == 0:
            if trchSize[1] < rounding_error:
                re = trchSize[0]
            elif isAc[idx] == 0:
                re = trchFixedPmtAmt[idx]
            else:
                _ = DisposableBalance[idx] * trchSize[0] / (trchSize[0] + trchSize[1])
                re = np.max([0, np.min([_, trchSize[0]])])
        else:
            if priorA2ResPrincipal[idx - 1] < rounding_error:
                re = priorA1ResPrincipal[idx - 1]
            elif isAc[idx] == 0:
                re = trchFixedPmtAmt[idx]
            else:
                _ = DisposableBalance[idx] * priorA1ResPrincipal[idx - 1] / (
                    priorA1ResPrincipal[idx - 1] + priorA2ResPrincipal[idx - 1])
                re = np.max([0, np.min([_, priorA1ResPrincipal[idx - 1]])])
        ls.append(re)
    return np.array(ls)


# 优先A-2档 priorA2RepPrincipal
# =MAX(MIN(D120,E110-E112),0)
# D120优先A-2档剩余本金 rounding_error归整误差 E110可分配余额 E112优先A-1档 priorA1RepPrincipal
def f_priorA2RepPrincipal(priorA2ResPrincipal, DisposableBalance, priorA1RepPrincipal, trchSize):
    ls = []
    for idx in range(len(priorA2ResPrincipal)):
        if idx == 0:
            re = np.max([np.min([trchSize[1], DisposableBalance[idx] - priorA1RepPrincipal[idx]]), 0])
        else:
            re = np.max([np.min([priorA2ResPrincipal[idx - 1], DisposableBalance[idx] - priorA1RepPrincipal[idx]]), 0])
        ls.append(re)
    return np.array(ls)


# 优先A-1档滞后 priorA1LagRepPrincipal
# =MAX(0,MIN(D119-E112,E110-SUM(E112:E113)))
# D119优先A-1档剩余本金 E112优先A-1档 E110可分配余额 E113优先A-2档 priorA2RepPrincipal
def f_priorA1LagRepPrincipal(priorA1ResPrincipal, priorA1RepPrincipal, DisposableBalance, priorA2RepPrincipal,
                             trchSize):
    ls = []
    for idx in range(len(priorA1ResPrincipal)):
        if idx == 0:
            _ = trchSize[0] - priorA1RepPrincipal[idx]
            __ = DisposableBalance[idx] - priorA1RepPrincipal[idx] - priorA2RepPrincipal[idx]
            re = np.max([0, np.min([_, __])])
        else:
            _ = priorA1ResPrincipal[idx - 1] - priorA1RepPrincipal[idx]
            __ = DisposableBalance[idx] - priorA1RepPrincipal[idx] - priorA2RepPrincipal[idx]
            re = np.max([0, np.min([_, __])])
        ls.append(re)
    return np.array(ls)


# 优先B档 priorBRepPrincipal
# =IF(AND(E119<rounding_error,E120<rounding_error),MIN(D121,E110-SUM(E112:E114)),0)
# E119优先A-1档剩余本金 E120优先A-2档剩余本金 D121优先B档剩余本金
# E110可分配余额  SUM(E112:E114)优先A-1 A-2 A-1Lag 偿还的本金
def f_sum4priorBRepPrincipal(priorA1RepPrincipal, priorA2RepPrincipal, priorA1LagRepPrincipal):
    return priorA1RepPrincipal + priorA2RepPrincipal + priorA1LagRepPrincipal


def f_priorBRepPrincipal(priorA1ResPrincipal, priorA2ResPrincipal, priorBResPrincipal, DisposableBalance,
                         sum4priorBRepPrincipal, rounding_error, trchSize):
    ls = []
    for idx in range(len(priorA1ResPrincipal)):
        if idx == 0:
            if priorA1ResPrincipal[idx] < rounding_error and priorA2ResPrincipal[idx] < rounding_error:
                re = np.min([trchSize[2], DisposableBalance[idx] - sum4priorBRepPrincipal[idx]])
            else:
                re = 0
        else:
            if priorA1ResPrincipal[idx] < rounding_error and priorA2ResPrincipal[idx] < rounding_error:
                re = np.min([priorBResPrincipal[idx - 1], DisposableBalance[idx] - sum4priorBRepPrincipal[idx]])
            else:
                re = 0
        ls.append(re)
    return np.array(ls)


# 优先C档 priorCRepPrincipal
# =IF(AND(E119<rounding_error,E120<rounding_error,E121<rounding_error),MIN(D122,E110-SUM(E112:E115)),0)

def f_sum4priorCRepPrincipal(sum4priorBRepPrincipal, priorBRepPrincipal):
    return sum4priorBRepPrincipal + priorBRepPrincipal


def f_priorCRepPrincipal(priorA1ResPrincipal, priorA2ResPrincipal, priorBResPrincipal, priorCResPrincipal,
                         DisposableBalance, sum4priorCRepPrincipal, rounding_error, trchSize):
    ls = []
    for idx in range(len(priorA1ResPrincipal)):
        if idx == 0:
            flag = priorA1ResPrincipal[idx] < rounding_error and priorA2ResPrincipal[idx] < rounding_error and \
                   priorBResPrincipal[idx] < rounding_error
            if flag:
                re = np.min([trchSize[3], DisposableBalance[idx] - sum4priorCRepPrincipal[idx]])
            else:
                re = 0
        else:
            flag = priorA1ResPrincipal[idx] < rounding_error and priorA2ResPrincipal[idx] < rounding_error and \
                   priorBResPrincipal[idx] < rounding_error
            if flag:
                re = np.min([priorCResPrincipal[idx - 1], DisposableBalance[idx] - sum4priorCRepPrincipal[idx]])
            else:
                re = 0
        ls.append(re)
    return np.array(ls)


# 次级档 subprimeRepPrincipal
def f_sum4subprimeRepPrincipal(sum4priorCRepPrincipal, priorCRepPrincipal):
    return sum4priorCRepPrincipal + priorCRepPrincipal


# =IF(E122<rounding_error,MIN(D123,E110-SUM(E112:E116)),0)
def f_subprimeRepPrincipal(priorCResPrincipal, subprimeResPrincipal, DisposableBalance, sum4subprimeRepPrincipal,
                           rounding_error, trchSize):
    ls = []
    for idx in range(len(priorCResPrincipal)):
        if idx == 0:
            if priorCResPrincipal[idx] < rounding_error:
                re = np.min([trchSize[4], DisposableBalance[idx] - sum4subprimeRepPrincipal[idx]])
            else:
                re = 0
        else:
            if priorCResPrincipal[idx] < rounding_error:
                re = np.min([subprimeResPrincipal[idx - 1], DisposableBalance[idx] - sum4subprimeRepPrincipal[idx]])
            else:
                re = 0
        ls.append(re)
    return np.array(ls)


# 剩余本金
# 优先A-1档 priorA1ResPrincipal
# =D119-E112-E114
# D119 priorA1ResPrincipal E112 priorA1RepPrincipal E114 priorA1LagRepPrincipal
def f_priorA1ResPrincipal(priorA1ResPrincipal, priorA1RepPrincipal, priorA1LagRepPrincipal, trchSize):
    ls = []
    for idx in range(len(priorA1ResPrincipal)):
        if idx == 0:
            re = trchSize[0] - priorA1RepPrincipal[idx] - priorA1LagRepPrincipal[idx]
        else:
            re = priorA1ResPrincipal[idx - 1] - priorA1RepPrincipal[idx] - priorA1LagRepPrincipal[idx]
        ls.append(re)
    return np.array(ls)


# 优先A-2档 priorA2ResPrincipal
# =D120-E113
def f_priorA2ResPrincipal(priorA2ResPrincipal, priorA2RepPrincipal, trchSize):
    ls = []
    for idx in range(len(priorA2ResPrincipal)):
        if idx == 0:
            re = trchSize[1] - priorA2RepPrincipal[idx]
        else:
            re = priorA2ResPrincipal[idx - 1] - priorA2RepPrincipal[idx]
        ls.append(re)
    return np.array(ls)


# 优先B档 priorBResPrincipal
# =D121-E115
def f_priorBResPrincipal(priorBResPrincipal, priorBRepPrincipal, trchSize):
    ls = []
    for idx in range(len(priorBResPrincipal)):
        if idx == 0:
            re = trchSize[2] - priorBRepPrincipal[idx]
        else:
            re = priorBResPrincipal[idx - 1] - priorBRepPrincipal[idx]
        ls.append(re)
    return np.array(ls)


# 优先C档 priorCResPrincipal
def f_priorCResPrincipal(priorCResPrincipal, priorCRepPrincipal, trchSize):
    ls = []
    for idx in range(len(priorCResPrincipal)):
        if idx == 0:
            re = trchSize[3] - priorCRepPrincipal[idx]
        else:
            re = priorCResPrincipal[idx - 1] - priorCRepPrincipal[idx]
        ls.append(re)
    return np.array(ls)


# 次级档 subprimeResPrincipal
def f_subprimeResPrincipal(subprimeResPrincipal, subprimeRepPrincipal, trchSize):
    ls = []
    for idx in range(len(subprimeResPrincipal)):
        if idx == 0:
            re = trchSize[4] - subprimeRepPrincipal[idx]
        else:
            re = subprimeResPrincipal[idx - 1] - subprimeRepPrincipal[idx]
        ls.append(re)
    return np.array(ls)


# 次级剩余收益 =IF(E123<rounding_error,MAX(0,E110-SUM(E112:E117)),0)
# E123次级档 subprimeResPrincipal E110 DisposableBalance
# E112:E117 priorA1RepPrincipal -- subprimeRepPrincipal

def f_sum4subprimeResIncome(priorA1RepPrincipal, priorA2RepPrincipal, priorA1LagRepPrincipal, priorBRepPrincipal,
                            priorCRepPrincipal, subprimeRepPrincipal):
    return priorA1RepPrincipal + priorA2RepPrincipal + priorA1LagRepPrincipal + priorBRepPrincipal + priorCRepPrincipal + subprimeRepPrincipal


def f_subprimeResIncome(subprimeResPrincipal, DisposableBalance, sum4subprimeResIncome, rounding_error):
    ls = []
    for a, b, c in zip(subprimeResPrincipal, DisposableBalance, sum4subprimeResIncome):
        if a < rounding_error:
            re = np.max([0, b - c])
        else:
            re = 0
        ls.append(re)
    return np.array(ls)


priod = 12


def f_main_profitAccountACNM(DisposableCashFlow, TaxFeesSum, interestDateGap, isAc, PrincipalAccount, DefaultPrincipal,
                             InitPriBalance, IncomeAccount, trchFixedPmtAmt, trchSize, trchInterestRate,
                             rounding_error):
    # 初始化
    priorA1Interest = np.zeros(priod)
    priorA2Interest = np.zeros(priod)
    priorBInterest = np.zeros(priod)
    priorCInterest = np.zeros(priod)
    priorInterestSum = np.zeros(priod)
    PrincipalAccountNM = np.zeros(priod)  # 计算完毕
    PrincipalAccountAC = np.zeros(priod)  ##计算完毕
    RePrincipalAccount = np.zeros(priod)
    SurplusFees = np.zeros(priod)
    subprimeLimitIncome = np.zeros(priod)
    ResCashFlow = np.zeros(priod)
    Spread = np.zeros(priod)

    ExcessInterest = np.zeros(priod)
    PrincipalCompensation = np.zeros(priod)
    # DefaultPrincipal
    IncomepaidBackPrincipal = np.zeros(priod)
    # IncomeActualPrincipal

    DisposableBalance = np.zeros(priod)

    priorA1RepPrincipal = np.zeros(priod)
    priorA2RepPrincipal = np.zeros(priod)
    priorA1LagRepPrincipal = np.zeros(priod)
    priorBRepPrincipal = np.zeros(priod)
    priorCRepPrincipal = np.zeros(priod)
    subprimeRepPrincipal = np.zeros(priod)

    priorA1ResPrincipal = np.zeros(priod)
    priorA2ResPrincipal = np.zeros(priod)
    priorBResPrincipal = np.zeros(priod)
    priorCResPrincipal = np.zeros(priod)
    subprimeResPrincipal = np.zeros(priod)
    subprimeResIncome = np.zeros(priod)

    # 循环计算

    for idx in range(priod):
        if idx == 0:
            # 利息
            priorA1Interest[idx] = f_priorA1Interest_0(DisposableCashFlow, TaxFeesSum, interestDateGap, trchSize,
                                                       trchInterestRate)
            priorA2Interest[idx] = f_priorA2Interest_0(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap,
                                                       trchSize, trchInterestRate)
            priorBInterest[idx] = f_priorBInterest_0(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap,
                                                     priorA2Interest, trchSize, trchInterestRate)
            priorCInterest[idx] = f_priorCInterest_0(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap,
                                                     priorA2Interest, priorBInterest, trchSize, trchInterestRate)
            priorInterestSum[idx] = priorA1Interest[idx] + priorA2Interest[idx] + priorBInterest[idx] + priorCInterest[
                idx]
            PrincipalAccountNM = f_PrincipalAccountNM(isAc, PrincipalAccount, DisposableCashFlow, TaxFeesSum,
                                                      priorInterestSum)
            PrincipalAccountAC = f_PrincipalAccountAC(isAc, DisposableCashFlow, TaxFeesSum, priorInterestSum)
            # 需先计算出收入应回补本金
            IncomeActualPrincipal = RePrincipalAccount
            IncomepaidBackPrincipal = f_IncomepaidBackPrincipal(PrincipalCompensation, DefaultPrincipal,
                                                                IncomepaidBackPrincipal, IncomeActualPrincipal)
            RePrincipalAccount = f_RePrincipalAccount(DisposableCashFlow, TaxFeesSum, priorInterestSum,
                                                      PrincipalAccountNM, PrincipalAccountAC, IncomepaidBackPrincipal)

            SurplusFees = f_SurplusFees(subprimeResPrincipal, DisposableCashFlow, TaxFeesSum, priorInterestSum,
                                        RePrincipalAccount, InitPriBalance, interestDateGap, trchSize)
            sum4subprimeLimitIncome = f_sum4subprimeLimitIncome(priorInterestSum, PrincipalAccountNM,
                                                                PrincipalAccountAC, RePrincipalAccount, SurplusFees)
            subprimeLimitIncome = f_subprimeLimitIncome(DisposableCashFlow, TaxFeesSum, sum4subprimeLimitIncome,
                                                        subprimeResPrincipal, interestDateGap, trchSize,
                                                        trchInterestRate)
            sum4ResCashFlow = f_sum4ResCashFlow(sum4subprimeLimitIncome, subprimeLimitIncome)
            ResCashFlow = f_ResCashFlow(DisposableCashFlow, TaxFeesSum, sum4ResCashFlow)
            Spread = f_Spread(IncomeAccount, TaxFeesSum, priorInterestSum, RePrincipalAccount, SurplusFees,
                              subprimeLimitIncome)

            # 回补机制
            ExcessInterest = f_ExcessInterest(IncomeAccount, TaxFeesSum, priorInterestSum, SurplusFees,
                                              subprimeLimitIncome)
            PrincipalCompensation = f_PrincipalCompensation(PrincipalAccount, PrincipalAccountNM, PrincipalAccountAC)
            # DefaultPrincipal

            # 本金账户
            DisposableBalance = f_DisposableBalance(PrincipalAccountNM, PrincipalAccountAC, RePrincipalAccount,
                                                    ResCashFlow)
            priorA1RepPrincipal = f_priorA1RepPrincipal(priorA2ResPrincipal, priorA1ResPrincipal, isAc, trchFixedPmtAmt,
                                                        DisposableBalance, trchSize, rounding_error)
            priorA2RepPrincipal = f_priorA2RepPrincipal(priorA2ResPrincipal, DisposableBalance, priorA1RepPrincipal,
                                                        trchSize)
            priorA1LagRepPrincipal = f_priorA1LagRepPrincipal(priorA1ResPrincipal, priorA1RepPrincipal,
                                                              DisposableBalance, priorA2RepPrincipal, trchSize)
            sum4priorBRepPrincipal = f_sum4priorBRepPrincipal(priorA1RepPrincipal, priorA2RepPrincipal,
                                                              priorA1LagRepPrincipal)
            priorBRepPrincipal = f_priorBRepPrincipal(priorA1ResPrincipal, priorA2ResPrincipal, priorBResPrincipal,
                                                      DisposableBalance, sum4priorBRepPrincipal, rounding_error,
                                                      trchSize)
            sum4priorCRepPrincipal = f_sum4priorCRepPrincipal(sum4priorBRepPrincipal, priorBRepPrincipal)
            priorCRepPrincipal = f_priorCRepPrincipal(priorA1ResPrincipal, priorA2ResPrincipal, priorBResPrincipal,
                                                      priorCResPrincipal, DisposableBalance, sum4priorCRepPrincipal,
                                                      rounding_error, trchSize)
            sum4subprimeRepPrincipal = f_sum4subprimeRepPrincipal(sum4priorCRepPrincipal, priorCRepPrincipal)
            subprimeRepPrincipal = f_subprimeRepPrincipal(priorCResPrincipal, subprimeResPrincipal, DisposableBalance,
                                                          sum4subprimeRepPrincipal, rounding_error, trchSize)

            priorA1ResPrincipal = f_priorA1ResPrincipal(priorA1ResPrincipal, priorA1RepPrincipal,
                                                        priorA1LagRepPrincipal, trchSize)
            priorA2ResPrincipal = f_priorA2ResPrincipal(priorA2ResPrincipal, priorA2RepPrincipal, trchSize)
            priorBResPrincipal = f_priorBResPrincipal(priorBResPrincipal, priorBRepPrincipal, trchSize)
            priorCResPrincipal = f_priorCResPrincipal(priorCResPrincipal, priorCRepPrincipal, trchSize)
            subprimeResPrincipal = f_subprimeResPrincipal(subprimeResPrincipal, subprimeRepPrincipal, trchSize)
            sum4subprimeResIncome = f_sum4subprimeResIncome(priorA1RepPrincipal, priorA2RepPrincipal,
                                                            priorA1LagRepPrincipal, priorBRepPrincipal,
                                                            priorCRepPrincipal, subprimeRepPrincipal)
            subprimeResIncome = f_subprimeResIncome(subprimeResPrincipal, DisposableBalance, sum4subprimeResIncome,
                                                    rounding_error)
        else:
            priorA1Interest[idx] = \
                f_priorA1Interest(DisposableCashFlow, TaxFeesSum, interestDateGap, priorA1ResPrincipal,
                                  trchInterestRate)[idx - 1]
            priorA2Interest[idx] = \
                f_priorA2Interest(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap,
                                  priorA2ResPrincipal, trchInterestRate)[
                    idx - 1]
            priorBInterest[idx] = \
                f_priorBInterest(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, priorBResPrincipal,
                                 priorA2Interest, trchInterestRate)[idx - 1]
            priorCInterest[idx] = \
                f_priorCInterest(DisposableCashFlow, TaxFeesSum, priorA1Interest, interestDateGap, priorCResPrincipal,
                                 priorA2Interest, priorBInterest, trchInterestRate)[idx - 1]
            priorInterestSum[idx] = priorA1Interest[idx] + priorA2Interest[idx] + priorBInterest[idx] + priorCInterest[
                idx]
            PrincipalAccountNM[idx] = \
                f_PrincipalAccountNM(isAc, PrincipalAccount, DisposableCashFlow, TaxFeesSum, priorInterestSum)[idx]
            PrincipalAccountAC[idx] = f_PrincipalAccountAC(isAc, DisposableCashFlow, TaxFeesSum, priorInterestSum)[idx]

            IncomeActualPrincipal = RePrincipalAccount
            IncomepaidBackPrincipal[idx] = \
                f_IncomepaidBackPrincipal(PrincipalCompensation, DefaultPrincipal, IncomepaidBackPrincipal,
                                          IncomeActualPrincipal)[idx]
            RePrincipalAccount[idx] = \
                f_RePrincipalAccount(DisposableCashFlow, TaxFeesSum, priorInterestSum, PrincipalAccountNM,
                                     PrincipalAccountAC, IncomepaidBackPrincipal)[idx]

            SurplusFees[idx] = \
                f_SurplusFees(subprimeResPrincipal, DisposableCashFlow, TaxFeesSum, priorInterestSum,
                              RePrincipalAccount,
                              InitPriBalance, interestDateGap, trchSize)[idx]
            sum4subprimeLimitIncome[idx] = \
                f_sum4subprimeLimitIncome(priorInterestSum, PrincipalAccountNM, PrincipalAccountAC, RePrincipalAccount,
                                          SurplusFees)[idx]
            subprimeLimitIncome[idx] = \
                f_subprimeLimitIncome(DisposableCashFlow, TaxFeesSum, sum4subprimeLimitIncome, subprimeResPrincipal,
                                      interestDateGap, trchSize, trchInterestRate)[idx]
            sum4ResCashFlow[idx] = f_sum4ResCashFlow(sum4subprimeLimitIncome, subprimeLimitIncome)[idx]
            ResCashFlow[idx] = f_ResCashFlow(DisposableCashFlow, TaxFeesSum, sum4ResCashFlow)[idx]
            Spread[idx] = \
                f_Spread(IncomeAccount, TaxFeesSum, priorInterestSum, RePrincipalAccount, SurplusFees,
                         subprimeLimitIncome)[
                    idx]

            # 回补机制
            ExcessInterest[idx] = \
                f_ExcessInterest(IncomeAccount, TaxFeesSum, priorInterestSum, SurplusFees, subprimeLimitIncome)[idx]
            PrincipalCompensation[idx] = \
                f_PrincipalCompensation(PrincipalAccount, PrincipalAccountNM, PrincipalAccountAC)[idx]
            # DefaultPrincipal

            # 本金账户
            DisposableBalance[idx] = \
                f_DisposableBalance(PrincipalAccountNM, PrincipalAccountAC, RePrincipalAccount, ResCashFlow)[idx]
            priorA1RepPrincipal[idx] = \
                f_priorA1RepPrincipal(priorA2ResPrincipal, priorA1ResPrincipal, isAc, trchFixedPmtAmt,
                                      DisposableBalance, trchSize, rounding_error)[
                    idx]
            priorA2RepPrincipal[idx] = \
                f_priorA2RepPrincipal(priorA2ResPrincipal, DisposableBalance, priorA1RepPrincipal, trchSize)[idx]
            priorA1LagRepPrincipal[idx] = \
                f_priorA1LagRepPrincipal(priorA1ResPrincipal, priorA1RepPrincipal, DisposableBalance,
                                         priorA2RepPrincipal, trchSize)[
                    idx]
            sum4priorBRepPrincipal[idx] = \
                f_sum4priorBRepPrincipal(priorA1RepPrincipal, priorA2RepPrincipal, priorA1LagRepPrincipal)[idx]
            priorBRepPrincipal[idx] = \
                f_priorBRepPrincipal(priorA1ResPrincipal, priorA2ResPrincipal, priorBResPrincipal, DisposableBalance,
                                     sum4priorBRepPrincipal, rounding_error, trchSize)[idx]
            sum4priorCRepPrincipal[idx] = f_sum4priorCRepPrincipal(sum4priorBRepPrincipal, priorBRepPrincipal)[idx]
            priorCRepPrincipal[idx] = \
                f_priorCRepPrincipal(priorA1ResPrincipal, priorA2ResPrincipal, priorBResPrincipal, priorCResPrincipal,
                                     DisposableBalance, sum4priorCRepPrincipal, rounding_error, trchSize)[idx]
            sum4subprimeRepPrincipal[idx] = f_sum4subprimeRepPrincipal(sum4priorCRepPrincipal, priorCRepPrincipal)[idx]
            subprimeRepPrincipal[idx] = \
                f_subprimeRepPrincipal(priorCResPrincipal, subprimeResPrincipal, DisposableBalance,
                                       sum4subprimeRepPrincipal, rounding_error, trchSize)[idx]

            priorA1ResPrincipal[idx] = \
                f_priorA1ResPrincipal(priorA1ResPrincipal, priorA1RepPrincipal, priorA1LagRepPrincipal, trchSize)[idx]
            priorA2ResPrincipal[idx] = f_priorA2ResPrincipal(priorA2ResPrincipal, priorA2RepPrincipal, trchSize)[idx]
            priorBResPrincipal[idx] = f_priorBResPrincipal(priorBResPrincipal, priorBRepPrincipal, trchSize)[idx]
            priorCResPrincipal[idx] = f_priorCResPrincipal(priorCResPrincipal, priorCRepPrincipal, trchSize)[idx]
            subprimeResPrincipal[idx] = f_subprimeResPrincipal(subprimeResPrincipal, subprimeRepPrincipal, trchSize)[
                idx]
            sum4subprimeResIncome[idx] = \
                f_sum4subprimeResIncome(priorA1RepPrincipal, priorA2RepPrincipal, priorA1LagRepPrincipal,
                                        priorBRepPrincipal, priorCRepPrincipal, subprimeRepPrincipal)[idx]
            subprimeResIncome[idx] = \
                f_subprimeResIncome(subprimeResPrincipal, DisposableBalance, sum4subprimeResIncome, rounding_error)[idx]

    return {
        'priorA1Interest': priorA1Interest,
        'priorA2Interest': priorA2Interest,
        'priorBInterest': priorBInterest,
        'priorCInterest': priorCInterest,
        'priorInterestSum': priorInterestSum,
        'PrincipalAccountNM': PrincipalAccountNM,
        'PrincipalAccountAC': PrincipalAccountAC,
        'RePrincipalAccount': RePrincipalAccount,
        'SurplusFees': SurplusFees,
        'subprimeLimitIncome': subprimeLimitIncome,
        'ResCashFlow': ResCashFlow,
        'Spread': Spread,

        'ExcessInterest': ExcessInterest,
        'PrincipalCompensation': PrincipalCompensation,
        # DefaultPrincipal
        'IncomepaidBackPrincipal': IncomepaidBackPrincipal,
        # IncomeActualPrincipal

        'DisposableBalance': DisposableBalance,

        'priorA1RepPrincipal': priorA1RepPrincipal,
        'priorA2RepPrincipal': priorA2RepPrincipal,
        'priorA1LagRepPrincipal': priorA1LagRepPrincipal,
        'priorBRepPrincipal': priorBRepPrincipal,
        'priorCRepPrincipal': priorCRepPrincipal,
        'subprimeRepPrincipal': subprimeRepPrincipal,

        'priorA1ResPrincipal': priorA1ResPrincipal,
        'priorA2ResPrincipal': priorA2ResPrincipal,
        'priorBResPrincipal': priorBResPrincipal,
        'priorCResPrincipal': priorCResPrincipal,
        'subprimeResPrincipal': subprimeResPrincipal,
        'subprimeResIncome': subprimeResIncome,
    }


# %%
'''
是否违约事件 isDefault
收益账 IncomeAccount_DE DE default event
本金帐 PrincipalAccount_DE
可分配现金流 DisposableCashFlow_DE

税费 
利息税收 InterestTax_DE
资金保管费用 storageFees_DE
审计费用 AuditFees_DE
跟踪评级费 GradeFees_DE
服务报酬 SurplusFees_DE
税费总额 TaxFeesSum_DE

优先A-1档利息 priorA1Interest_DE
优先A-1档本金 priorA1RepPrincipal_DE
优先A-2档利息 priorA2Interest_DE
优先A-2档本金 priorA2RepPrincipal_DE
优先B档利息 priorBInterest_DE
优先B档本金 priorBRepPrincipal_DE
优先C档利息 priorCInterest_DE
优先C档本金 priorCRepPrincipal_DE
次级档本金 subprimeRepPrincipal_DE
次级档收益 subprimeRepIncome_DE
剩余本金
优先A-1档 priorA1ResPrincipal_DE
优先A-2档 priorA2ResPrincipal_DE
优先B档 priorBResPrincipal_DE
优先C档 priorCResPrincipal_DE
次级档 subprimeResPrincipal_DE
'''

# 是否违约事件 isDefault
# =IF(OR(D126=1,E82-E89<D119*$K31*E$76),1,0)
# D126前一期是否违约 E82=E129可分配现金流 E89税费总额 D119本金账户A1档剩余本金
# $K31A1档预发行利率 E$76计息日间隔
isDefault_init = 0


def f_isDefault(isDefault, DisposableCashFlow, TaxFeesSum, priorA1ResPrincipal, interestDateGap, trchSize,
                trchInterestRate):
    ls = []
    for idx in range(len(DisposableCashFlow)):
        if idx == 0:
            flag_init = (DisposableCashFlow[idx] - TaxFeesSum[idx]) < (
                trchSize[0] * trchInterestRate[0] * interestDateGap[idx])
            if isDefault_init == 1 or flag_init:
                re = 1
            else:
                re = 0
        else:
            flag_init = (DisposableCashFlow[idx] - TaxFeesSum[idx]) < (
                priorA1ResPrincipal[idx - 1] * trchInterestRate[0] * interestDateGap[idx])
            if isDefault[idx - 1] == 1 or flag_init:
                re = 1
            else:
                re = 0
        ls.append(re)
    return np.array(ls)


# 收益账 IncomeAccount_DE DE default event
# IncomeAccount_DE=IncomeAccount

# 本金帐 PrincipalAccount_DE
# PrincipalAccount_DE=PrincipalAccount

# 可分配现金流 DisposableCashFlow_DE
# DisposableCashFlow_DE=DisposableCashFlow

# 利息税收 InterestTax
# =MAX(MIN(E129,E80*tax_rate),0)
# E129 可分配现金流 E80=E127 收益帐 tax_rate
tax_rate = 0.0


def f_InterestTax_DE(DisposableCashFlow_DE, IncomeAccount_DE):
    ls = []
    for a, b in zip(DisposableCashFlow_DE, IncomeAccount_DE):
        re = np.max([np.min([a, b * tax_rate]), 0])
        ls.append(re)
    return np.array(ls)


# 资金保管费用 storageFees_DE
# =MAX(MIN(E129-E131,IF(E1=1,E76*SUM($H$31:$H$34)*$G$44,IF(MOD(E1,4)=0,SUM($H$31:$H$34)*$G$44,0))),0)

# E129 可分配现金流 E131 InterestTax_DE E76 计息日间隔 E1期序term
# $H$31:$H$34 trchSize[0:4] $G$44 资金保管机构服务费 1%
stFeesRate = 0.01


def f_storageFees_DE(DisposableCashFlow_DE, InterestTax_DE, interestDateGap, term, trchSize):
    ls = []
    for a, b, c, d in zip(DisposableCashFlow_DE, InterestTax_DE, interestDateGap, term):
        if d == 1:
            _ = c * np.sum(trchSize[0:4]) * stFeesRate
        elif d % 4 == 0:
            _ = np.sum(trchSize[0:4]) * stFeesRate
        else:
            _ = 0

        re = np.max([np.min([a - b, _]), 0])
        ls.append(re)
    return np.array(ls)


# 审计费用 AuditFees_DE
# =MAX(MIN(E129-SUM(E131:E132),IF(MOD(E1,4)=0,$H$45,0)),0)
AuditFeesFixed = 50000


def f_AuditFees_DE(DisposableCashFlow_DE, InterestTax_DE, storageFees_DE, term):
    ls = []
    for a, b, c, d in zip(DisposableCashFlow_DE, InterestTax_DE, storageFees_DE, term):
        if d % 4 == 0:
            _ = AuditFeesFixed
        else:
            _ = 0
        re = np.max([np.min([a - b - c, _]), 0])
        ls.append(re)
    return np.array(ls)


# 跟踪评级费 GradeFees_DE
# =MAX(MIN(E129-SUM(E131:E133),IF(MOD(E1,4)=0,$H$46,0)),0)

def f_GradeFees_DE(DisposableCashFlow_DE, InterestTax_DE, storageFees_DE, term, AuditFees_DE):
    ls = []
    for a, b, c, d, e in zip(DisposableCashFlow_DE, InterestTax_DE, storageFees_DE, term, AuditFees_DE):
        if d % 4 == 0:
            _ = AuditFeesFixed
        else:
            _ = 0
        re = np.max([np.min([a - b - c - e, _]), 0])
        ls.append(re)
    return np.array(ls)


# 服务报酬 SurplusFees_DE
# =MAX(MIN(E129-SUM(E131:E134),E88+E66*($G$48+$G$49)*E$76),0)
# E88服务报酬，预先支付 ServiceFees  E66期初本金余额 InitPriBalance  E76计息日间隔 interestDateGap
SurplusFeesRate = 0.0015 + 0.0005


def f_SurplusFees_DE(DisposableCashFlow_DE, InterestTax_DE, storageFees_DE, AuditFees_DE, GradeFees_DE, ServiceFees,
                     InitPriBalance, interestDateGap):
    ls = []
    for a, b, c, d, e, f, g, h in zip(DisposableCashFlow_DE, InterestTax_DE, storageFees_DE, AuditFees_DE, GradeFees_DE,
                                      ServiceFees, InitPriBalance, interestDateGap):
        re = np.max([np.min([a - b - c - d - e, f + g * SurplusFeesRate * h]), 0])
        ls.append(re)
    return np.array(ls)


# 税费总额 TaxFeesSum_DE
# =SUM(E131:E135)

def f_TaxFeesSum_DE(InterestTax_DE, storageFees_DE, AuditFees_DE, GradeFees_DE, SurplusFees_DE):
    return InterestTax_DE + storageFees_DE + AuditFees_DE + GradeFees_DE + SurplusFees_DE


def f_main_TaxFeesDE(isDefault, DisposableCashFlow, TaxFeesSum, priorA1ResPrincipal, interestDateGap, IncomeAccount,
                     PrincipalAccount, ServiceFees, InitPriBalance, term, trchSize, trchInterestRate):
    isDefault = f_isDefault(isDefault, DisposableCashFlow, TaxFeesSum, priorA1ResPrincipal, interestDateGap, trchSize,
                            trchInterestRate)
    # 收益账 IncomeAccount_DE DE default event
    IncomeAccount_DE = IncomeAccount
    # 本金帐 PrincipalAccount_DE
    PrincipalAccount_DE = PrincipalAccount
    # 可分配现金流 DisposableCashFlow_DE
    DisposableCashFlow_DE = DisposableCashFlow
    InterestTax_DE = f_InterestTax_DE(DisposableCashFlow_DE, IncomeAccount_DE)
    storageFees_DE = f_storageFees_DE(DisposableCashFlow_DE, InterestTax_DE, interestDateGap, term, trchSize)
    AuditFees_DE = f_AuditFees_DE(DisposableCashFlow_DE, InterestTax_DE, storageFees_DE, term)
    GradeFees_DE = f_GradeFees_DE(DisposableCashFlow_DE, InterestTax_DE, storageFees_DE, term, AuditFees_DE)
    SurplusFees_DE = f_SurplusFees_DE(DisposableCashFlow_DE, InterestTax_DE, storageFees_DE, AuditFees_DE, GradeFees_DE,
                                      ServiceFees, InitPriBalance, interestDateGap)
    TaxFeesSum_DE = f_TaxFeesSum_DE(InterestTax_DE, storageFees_DE, AuditFees_DE, GradeFees_DE, SurplusFees_DE)

    return {
        'isDefault': isDefault,  # 是否违约事件 isDefault
        'IncomeAccount_DE': IncomeAccount_DE,  # 收益账 IncomeAccount_DE DE default event
        'PrincipalAccount_DE': PrincipalAccount_DE,  # 本金帐 PrincipalAccount_DE
        'DisposableCashFlow_DE': DisposableCashFlow_DE,  # 可分配现金流 DisposableCashFlow_DE
        # 税费
        'InterestTax_DE': InterestTax_DE,  # 利息税收 InterestTax_DE
        'storageFees_DE': storageFees_DE,  # 资金保管费用 storageFees_DE
        'AuditFees_DE': AuditFees_DE,  # 审计费用 AuditFees_DE
        'GradeFees_DE': GradeFees_DE,  # 跟踪评级费 GradeFees_DE
        'SurplusFees_DE': SurplusFees_DE,  # 服务报酬 SurplusFees_DE
        'TaxFeesSum_DE': TaxFeesSum_DE,  # 税费总额 TaxFeesSum_DE
    }


# %% 违约事件的情况
# 优先A-1档利息 priorA1Interest_DE
# =MAX(MIN(E129-E136,D148*E$76*$K31),0)

def f_priorA1Interest_DE(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1ResPrincipal_DE, interestDateGap, trchSize,
                         trchInterestRate):
    ls = []
    for idx in range(len(DisposableCashFlow_DE)):
        if idx == 0:
            re = np.max([np.min([DisposableCashFlow_DE[idx] - TaxFeesSum_DE[idx],
                                 trchSize[0] * interestDateGap[idx] * trchInterestRate[0]]), 0])
        else:
            re = np.max([np.min([DisposableCashFlow_DE[idx] - TaxFeesSum_DE[idx],
                                 priorA1ResPrincipal_DE[idx - 1] * interestDateGap[idx] * trchInterestRate[0]]), 0])
        ls.append(re)
    return np.array(ls)


# 优先A-1档本金 priorA1RepPrincipal_DE
# =MAX(MIN(E129-SUM(E136:E137),E37),0)
# E37 固定摊还优先A1档本金 trchFixedPmtAmt
def f_priorA1RepPrincipal_DE(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1Interest_DE, trchFixedPmtAmt):
    ls = []
    for a, b, c, d in zip(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1Interest_DE, trchFixedPmtAmt):
        re = np.max([np.min([a - b - c, d]), 0])
        ls.append(re)
    return np.array(ls)


# priorA1RepPrincipal_DE


# 优先A-2档利息 priorA2Interest_DE
# =MAX(MIN(E129-SUM(E136:E138),D149*E$76*$K32),0)
def f_priorA2Interest_DE(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1Interest_DE, priorA1RepPrincipal_DE,
                         priorA2ResPrincipal_DE, interestDateGap, trchSize, trchInterestRate):
    ls = []
    for idx in range(len(DisposableCashFlow_DE)):
        if idx == 0:
            _ = DisposableCashFlow_DE[idx] - TaxFeesSum_DE[idx] - priorA1Interest_DE[idx] - priorA1RepPrincipal_DE[idx]
            __ = trchSize[1] * interestDateGap[idx] * trchInterestRate[1]
            re = np.max([np.min([_, __]), 0])
        else:
            _ = DisposableCashFlow_DE[idx] - TaxFeesSum_DE[idx] - priorA1Interest_DE[idx] - priorA1RepPrincipal_DE[idx]
            __ = priorA2ResPrincipal_DE[idx - 1] * interestDateGap[idx] * trchInterestRate[1]
            re = np.max([np.min([_, __]), 0])
        ls.append(re)
    return np.array(ls)


# 优先A-2档本金 priorA2RepPrincipal_DE
# =MAX(MIN(E129-SUM(E136:E139),D149),0)
def f_sum4priorA2RepPrincipal_DE(TaxFeesSum_DE, priorA1Interest_DE, priorA1RepPrincipal_DE, priorA2Interest_DE):
    return TaxFeesSum_DE + priorA1Interest_DE + priorA1RepPrincipal_DE + priorA2Interest_DE


# sum4priorA2RepPrincipal_DE=f_sum4priorA2RepPrincipal_DE(TaxFeesSum_DE,priorA1Interest_DE,priorA1RepPrincipal_DE,priorA2Interest_DE)

def f_priorA2RepPrincipal_DE(DisposableCashFlow_DE, sum4priorA2RepPrincipal_DE, priorA2ResPrincipal_DE, trchSize):
    ls = []
    for idx in range(len(DisposableCashFlow_DE)):
        if idx == 0:
            re = np.max([np.min([DisposableCashFlow_DE[idx] - sum4priorA2RepPrincipal_DE[idx], trchSize[1]])])
        else:
            re = np.max([np.min(
                [DisposableCashFlow_DE[idx] - sum4priorA2RepPrincipal_DE[idx], priorA2ResPrincipal_DE[idx - 1]])])
        ls.append(re)
    return np.array(ls)


# 优先B档利息 priorBInterest_DE
# =MAX(MIN(E129-SUM(E136:E140),D150*E$76*$K33),0)
def f_sum4priorBInterest_DE(sum4priorA2RepPrincipal_DE, priorA2RepPrincipal_DE):
    return sum4priorA2RepPrincipal_DE + priorA2RepPrincipal_DE


def f_priorBInterest_DE(DisposableCashFlow_DE, sum4priorBInterest_DE, priorBResPrincipal_DE, interestDateGap, trchSize,
                        trchInterestRate):
    ls = []
    for idx in range(len(DisposableCashFlow_DE)):
        if idx == 0:
            _ = DisposableCashFlow_DE[idx] - sum4priorBInterest_DE[idx]
            __ = trchSize[2] * interestDateGap[idx] * trchInterestRate[2]
            re = np.max([np.min([_, __]), 0])
        else:
            _ = DisposableCashFlow_DE[idx] - sum4priorBInterest_DE[idx]
            __ = priorBResPrincipal_DE[idx - 1] * interestDateGap[idx] * trchInterestRate[2]
            re = np.max([np.min([_, __]), 0])
        ls.append(re)
    return np.array(ls)


# 优先B档本金 priorBRepPrincipal_DE
# =MAX(MIN(E129-SUM(E136:E141),D150),0)
def f_sum4priorBRepPrincipal_DE(sum4priorBInterest_DE, priorBInterest_DE):
    return sum4priorBInterest_DE + priorBInterest_DE


def f_priorBRepPrincipal_DE(DisposableCashFlow_DE, sum4priorBRepPrincipal_DE, priorBResPrincipal_DE, trchSize):
    ls = []
    for idx in range(len(DisposableCashFlow_DE)):
        if idx == 0:
            re = np.max([np.min([DisposableCashFlow_DE[idx] - sum4priorBRepPrincipal_DE[idx], trchSize[2]])])
        else:
            re = np.max(
                [np.min([DisposableCashFlow_DE[idx] - sum4priorBRepPrincipal_DE[idx], priorBResPrincipal_DE[idx - 1]])])
        ls.append(re)
    return np.array(ls)


# 优先C档利息 priorCInterest_DE
# =MAX(MIN(E129-SUM(E136:E140),D151*E$76*$K34),0)
# ??感觉不对啊 需讨论
def f_sum4priorCInterest_DE(sum4priorBRepPrincipal_DE, priorBRepPrincipal_DE):
    return sum4priorBRepPrincipal_DE + priorBRepPrincipal_DE


def f_priorCInterest_DE(DisposableCashFlow_DE, sum4priorCInterest_DE, priorCResPrincipal_DE, interestDateGap, trchSize,
                        trchInterestRate):
    ls = []
    for idx in range(len(DisposableCashFlow_DE)):
        if idx == 0:
            _ = DisposableCashFlow_DE[idx] - sum4priorCInterest_DE[idx]
            __ = trchSize[3] * interestDateGap[idx] * trchInterestRate[3]
            re = np.max([np.min([_, __]), 0])
        else:
            _ = DisposableCashFlow_DE[idx] - sum4priorCInterest_DE[idx]
            __ = priorCResPrincipal_DE[idx - 1] * interestDateGap[idx] * trchInterestRate[3]
            re = np.max([np.min([_, __]), 0])
        ls.append(re)
    return np.array(ls)


# 优先C档本金 priorCRepPrincipal_DE
# =MAX(MIN(E129-SUM(E136:E143),D151),0)
def f_sum4priorCRepPrincipal_DE(sum4priorCInterest_DE, priorCInterest_DE):
    return sum4priorCInterest_DE + priorCInterest_DE


def f_priorCRepPrincipal_DE(DisposableCashFlow_DE, sum4priorCRepPrincipal_DE, priorCResPrincipal_DE, trchSize):
    ls = []
    for idx in range(len(DisposableCashFlow_DE)):
        if idx == 0:
            re = np.max([np.min([DisposableCashFlow_DE[idx] - sum4priorCRepPrincipal_DE[idx], trchSize[3]])])
        else:
            re = np.max(
                [np.min([DisposableCashFlow_DE[idx] - sum4priorCRepPrincipal_DE[idx], priorCResPrincipal_DE[idx - 1]])])
        ls.append(re)
    return np.array(ls)


# 次级档本金 subprimeRepPrincipal_DE
# =MAX(MIN(E129-SUM(E136:E144),D152),0)

def f_sum4subprimeRepPrincipal_DE(sum4priorCRepPrincipal_DE, priorCRepPrincipal_DE):
    return sum4priorCRepPrincipal_DE + priorCRepPrincipal_DE


def f_subprimeRepPrincipal_DE(DisposableCashFlow_DE, sum4subprimeRepPrincipal_DE, subprimeResPrincipal_DE, trchSize):
    ls = []
    for idx in range(len(DisposableCashFlow_DE)):
        if idx == 0:
            _ = DisposableCashFlow_DE[idx] - sum4subprimeRepPrincipal_DE[idx]
            re = np.max([np.min([_, trchSize[4]]), 0])
        else:
            _ = DisposableCashFlow_DE[idx] - sum4subprimeRepPrincipal_DE[idx]
            re = np.max([np.min([_, subprimeResPrincipal_DE[idx - 1]]), 0])
        ls.append(re)
    return np.array(ls)


# 次级档收益 subprimeRepIncome_DE
# =MAX(E129-SUM(E136:E145),0)
def f_sum4subprimeRepIncome_DE(sum4subprimeRepPrincipal_DE, subprimeRepPrincipal_DE):
    return sum4subprimeRepPrincipal_DE + subprimeRepPrincipal_DE


def f_subprimeRepIncome_DE(DisposableCashFlow_DE, sum4subprimeRepIncome_DE):
    ls = []
    for a, b in zip(DisposableCashFlow_DE, sum4subprimeRepIncome_DE):
        ls.append(np.max([a - b, 0]))
    return np.array(ls)


######经发现 违约事件的剩余本金计算逻辑一致，故写一个函数即可
# 优先A-1档 priorA1ResPrincipal_DE
# 优先A-2档 priorA2ResPrincipal_DE
# 优先B档 priorBResPrincipal_DE
# 优先C档 priorCResPrincipal_DE
# 次级档 subprimeResPrincipal_DE
# =D148-E138
def f_ResPrincipal_DE(ResPrincipal_DE, RepPrincipal_DE, n, trchSize):
    ls = []
    for idx in range(len(ResPrincipal_DE)):
        if idx == 0:
            re = trchSize[n] - RepPrincipal_DE[idx]
        else:
            re = ResPrincipal_DE[idx - 1] - RepPrincipal_DE[idx]
        ls.append(re)
    return np.array(ls)


def f_priorA1ResPrincipal_DE(priorA1ResPrincipal_DE, priorA1RepPrincipal_DE, trchSize):
    return f_ResPrincipal_DE(priorA1ResPrincipal_DE, priorA1RepPrincipal_DE, 0, trchSize)


def f_priorA2ResPrincipal_DE(priorA2ResPrincipal_DE, priorA2RepPrincipal_DE, trchSize):
    return f_ResPrincipal_DE(priorA2ResPrincipal_DE, priorA2RepPrincipal_DE, 1, trchSize)


def f_priorBResPrincipal_DE(priorBResPrincipal_DE, priorBRepPrincipal_DE, trchSize):
    return f_ResPrincipal_DE(priorBResPrincipal_DE, priorBRepPrincipal_DE, 2, trchSize)


def f_priorCResPrincipal_DE(priorCResPrincipal_DE, priorCRepPrincipal_DE, trchSize):
    return f_ResPrincipal_DE(priorCResPrincipal_DE, priorCRepPrincipal_DE, 3, trchSize)


def f_subprimeResPrincipal_DE(subprimeResPrincipal_DE, subprimeRepPrincipal_DE, trchSize):
    return f_ResPrincipal_DE(subprimeResPrincipal_DE, subprimeRepPrincipal_DE, 4, trchSize)


def f_main_AccountDE(DisposableCashFlow_DE, TaxFeesSum_DE,
                     interestDateGap, trchFixedPmtAmt, trchSize, trchInterestRate):
    priorA1Interest_DE = np.zeros(priod)
    priorA1RepPrincipal_DE = np.zeros(priod)
    priorA2Interest_DE = np.zeros(priod)
    priorA2RepPrincipal_DE = np.zeros(priod)
    priorBInterest_DE = np.zeros(priod)
    priorBRepPrincipal_DE = np.zeros(priod)
    priorCInterest_DE = np.zeros(priod)
    priorCRepPrincipal_DE = np.zeros(priod)
    subprimeRepPrincipal_DE = np.zeros(priod)
    subprimeRepIncome_DE = np.zeros(priod)

    priorA1ResPrincipal_DE = np.zeros(priod)
    priorA2ResPrincipal_DE = np.zeros(priod)
    priorBResPrincipal_DE = np.zeros(priod)
    priorCResPrincipal_DE = np.zeros(priod)
    subprimeResPrincipal_DE = np.zeros(priod)

    for idx in range(priod):
        if idx == 0:
            priorA1Interest_DE = f_priorA1Interest_DE(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1ResPrincipal_DE,
                                                      interestDateGap, trchSize, trchInterestRate)
            priorA1RepPrincipal_DE = f_priorA1RepPrincipal_DE(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1Interest_DE,
                                                              trchFixedPmtAmt)

            priorA2Interest_DE = f_priorA2Interest_DE(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1Interest_DE,
                                                      priorA1RepPrincipal_DE, priorA2ResPrincipal_DE, interestDateGap,
                                                      trchSize, trchInterestRate)
            sum4priorA2RepPrincipal_DE = f_sum4priorA2RepPrincipal_DE(TaxFeesSum_DE, priorA1Interest_DE,
                                                                      priorA1RepPrincipal_DE, priorA2Interest_DE)
            priorA2RepPrincipal_DE = f_priorA2RepPrincipal_DE(DisposableCashFlow_DE, sum4priorA2RepPrincipal_DE,
                                                              priorA2ResPrincipal_DE, trchSize)

            sum4priorBInterest_DE = f_sum4priorBInterest_DE(sum4priorA2RepPrincipal_DE, priorA2RepPrincipal_DE)
            priorBInterest_DE = f_priorBInterest_DE(DisposableCashFlow_DE, sum4priorBInterest_DE, priorBResPrincipal_DE,
                                                    interestDateGap, trchSize, trchInterestRate)
            sum4priorBRepPrincipal_DE = f_sum4priorBRepPrincipal_DE(sum4priorBInterest_DE, priorBInterest_DE)
            priorBRepPrincipal_DE = f_priorBRepPrincipal_DE(DisposableCashFlow_DE, sum4priorBRepPrincipal_DE,
                                                            priorBResPrincipal_DE, trchSize)

            sum4priorCInterest_DE = f_sum4priorCInterest_DE(sum4priorBRepPrincipal_DE, priorBRepPrincipal_DE)
            priorCInterest_DE = f_priorCInterest_DE(DisposableCashFlow_DE, sum4priorCInterest_DE, priorCResPrincipal_DE,
                                                    interestDateGap, trchSize, trchInterestRate)
            sum4priorCRepPrincipal_DE = f_sum4priorCRepPrincipal_DE(sum4priorCInterest_DE, priorCInterest_DE)
            priorCRepPrincipal_DE = f_priorCRepPrincipal_DE(DisposableCashFlow_DE, sum4priorCRepPrincipal_DE,
                                                            priorCResPrincipal_DE, trchSize)

            sum4subprimeRepPrincipal_DE = f_sum4subprimeRepPrincipal_DE(sum4priorCRepPrincipal_DE,
                                                                        priorCRepPrincipal_DE)
            subprimeRepPrincipal_DE = f_subprimeRepPrincipal_DE(DisposableCashFlow_DE, sum4subprimeRepPrincipal_DE,
                                                                subprimeResPrincipal_DE, trchSize)
            sum4subprimeRepIncome_DE = f_sum4subprimeRepIncome_DE(sum4subprimeRepPrincipal_DE, subprimeRepPrincipal_DE)
            subprimeRepIncome_DE = f_subprimeRepIncome_DE(DisposableCashFlow_DE, sum4subprimeRepIncome_DE)

            priorA1ResPrincipal_DE = f_priorA1ResPrincipal_DE(priorA1ResPrincipal_DE, priorA1RepPrincipal_DE, trchSize)
            priorA2ResPrincipal_DE = f_priorA2ResPrincipal_DE(priorA2ResPrincipal_DE, priorA2RepPrincipal_DE, trchSize)
            priorBResPrincipal_DE = f_priorBResPrincipal_DE(priorBResPrincipal_DE, priorBRepPrincipal_DE, trchSize)
            priorCResPrincipal_DE = f_priorCResPrincipal_DE(priorCResPrincipal_DE, priorCRepPrincipal_DE, trchSize)
            subprimeResPrincipal_DE = f_subprimeResPrincipal_DE(subprimeResPrincipal_DE, subprimeRepPrincipal_DE,
                                                                trchSize)
        else:
            priorA1Interest_DE[idx] = \
                f_priorA1Interest_DE(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1ResPrincipal_DE, interestDateGap,
                                     trchSize, trchInterestRate)[idx]
            priorA1RepPrincipal_DE[idx] = \
                f_priorA1RepPrincipal_DE(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1Interest_DE, trchFixedPmtAmt)[idx]

            priorA2Interest_DE[idx] = \
                f_priorA2Interest_DE(DisposableCashFlow_DE, TaxFeesSum_DE, priorA1Interest_DE, priorA1RepPrincipal_DE,
                                     priorA2ResPrincipal_DE, interestDateGap, trchSize, trchInterestRate)[idx]
            sum4priorA2RepPrincipal_DE[idx] = \
                f_sum4priorA2RepPrincipal_DE(TaxFeesSum_DE, priorA1Interest_DE, priorA1RepPrincipal_DE,
                                             priorA2Interest_DE)[
                    idx]
            priorA2RepPrincipal_DE[idx] = \
                f_priorA2RepPrincipal_DE(DisposableCashFlow_DE, sum4priorA2RepPrincipal_DE, priorA2ResPrincipal_DE,
                                         trchSize)[idx]

            sum4priorBInterest_DE[idx] = f_sum4priorBInterest_DE(sum4priorA2RepPrincipal_DE, priorA2RepPrincipal_DE)[
                idx]
            priorBInterest_DE[idx] = \
                f_priorBInterest_DE(DisposableCashFlow_DE, sum4priorBInterest_DE, priorBResPrincipal_DE,
                                    interestDateGap, trchSize, trchInterestRate)[
                    idx]
            sum4priorBRepPrincipal_DE[idx] = f_sum4priorBRepPrincipal_DE(sum4priorBInterest_DE, priorBInterest_DE)[idx]
            priorBRepPrincipal_DE[idx] = \
                f_priorBRepPrincipal_DE(DisposableCashFlow_DE, sum4priorBRepPrincipal_DE, priorBResPrincipal_DE,
                                        trchSize)[idx]

            sum4priorCInterest_DE[idx] = f_sum4priorCInterest_DE(sum4priorBRepPrincipal_DE, priorBRepPrincipal_DE)[idx]
            priorCInterest_DE[idx] = \
                f_priorCInterest_DE(DisposableCashFlow_DE, sum4priorCInterest_DE, priorCResPrincipal_DE,
                                    interestDateGap, trchSize, trchInterestRate)[
                    idx]
            sum4priorCRepPrincipal_DE[idx] = f_sum4priorCRepPrincipal_DE(sum4priorCInterest_DE, priorCInterest_DE)[idx]
            priorCRepPrincipal_DE[idx] = \
                f_priorCRepPrincipal_DE(DisposableCashFlow_DE, sum4priorCRepPrincipal_DE, priorCResPrincipal_DE,
                                        trchSize)[idx]

            sum4subprimeRepPrincipal_DE[idx] = \
                f_sum4subprimeRepPrincipal_DE(sum4priorCRepPrincipal_DE, priorCRepPrincipal_DE)[idx]
            subprimeRepPrincipal_DE[idx] = \
                f_subprimeRepPrincipal_DE(DisposableCashFlow_DE, sum4subprimeRepPrincipal_DE, subprimeResPrincipal_DE,
                                          trchSize)[
                    idx]
            sum4subprimeRepIncome_DE[idx] = \
                f_sum4subprimeRepIncome_DE(sum4subprimeRepPrincipal_DE, subprimeRepPrincipal_DE)[idx]
            subprimeRepIncome_DE[idx] = f_subprimeRepIncome_DE(DisposableCashFlow_DE, sum4subprimeRepIncome_DE)[idx]

            priorA1ResPrincipal_DE[idx] = \
                f_priorA1ResPrincipal_DE(priorA1ResPrincipal_DE, priorA1RepPrincipal_DE, trchSize)[idx]
            priorA2ResPrincipal_DE[idx] = \
                f_priorA2ResPrincipal_DE(priorA2ResPrincipal_DE, priorA2RepPrincipal_DE, trchSize)[idx]
            priorBResPrincipal_DE[idx] = \
                f_priorBResPrincipal_DE(priorBResPrincipal_DE, priorBRepPrincipal_DE, trchSize)[idx]
            priorCResPrincipal_DE[idx] = \
                f_priorCResPrincipal_DE(priorCResPrincipal_DE, priorCRepPrincipal_DE, trchSize)[idx]
            subprimeResPrincipal_DE[idx] = \
                f_subprimeResPrincipal_DE(subprimeResPrincipal_DE, subprimeRepPrincipal_DE, trchSize)[
                    idx]

    return {
        'priorA1Interest_DE': priorA1Interest_DE,
        'priorA1RepPrincipal_DE': priorA1RepPrincipal_DE,
        'priorA2Interest_DE': priorA2Interest_DE,
        'priorA2RepPrincipal_DE': priorA2RepPrincipal_DE,
        'priorBInterest_DE': priorBInterest_DE,
        'priorBRepPrincipal_DE': priorBRepPrincipal_DE,
        'priorCInterest_DE': priorCInterest_DE,
        'priorCRepPrincipal_DE': priorCRepPrincipal_DE,
        'subprimeRepPrincipal_DE': subprimeRepPrincipal_DE,
        'subprimeRepIncome_DE': subprimeRepIncome_DE,

        'priorA1ResPrincipal_DE': priorA1ResPrincipal_DE,
        'priorA2ResPrincipal_DE': priorA2ResPrincipal_DE,
        'priorBResPrincipal_DE': priorBResPrincipal_DE,
        'priorCResPrincipal_DE': priorCResPrincipal_DE,
        'subprimeResPrincipal_DE': subprimeResPrincipal_DE,
    }
