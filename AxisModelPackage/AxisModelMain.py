# main函数 主要用于计算出所有的可以计算的

from AxisModelPackage.AxisModelComputer import *

'''
f_main_IncomeAccount 收入项
f_main_OutcomeAccount 支出项
f_main_ProfitAccount_1 收益帐前半部分
f_main_TaxFeesNMAC 收益账户税费
f_main_profitAccountACNM 收益账户利息本金计算
f_main_TaxFeesDE 违约事件税费相关
f_main_AccountDE 违约事件利息本金
'''

def f_AxisModelMidComp(Input_data, rounding_error):
    xx = f_main_IncomeAccount(np.array(Input_data['rePrincipal']),  # PrincipalReceivable
                              Input_data['term'], np.array(Input_data['pdByTerm']),
                              Input_data['rrAvg'],
                              Input_data['rec_term'],
                              np.array(Input_data['ppByTerm']),
                              np.array(Input_data['reInterest']),  # InterestReceivable
                              np.array(Input_data['opBalance']),
                              Input_data['rateReinv'],
                              Input_data['daysReinv'],
                              Input_data['daysInYear'],  # days_in_year
                              Input_data['isReinv']
                              )  # InitBalance

    zz = f_main_OutcomeAccount(Input_data['calList'],  # compDate
                               Input_data['inrList'],  # interestDate
                               Input_data['secList'],  # cashingDate
                               Input_data['daysInYear']
                               )

    yy = f_main_ProfitAccount_1(xx['DefaultPrincipal'], xx['InitPriBalance'], xx['InterestPaid'],
                                xx['ReinvestedEarnings'],
                                xx['PriPaid'])

    aa = f_main_TaxFeesNMAC(yy['DisposableCashFlow'], yy['IncomeAccount'], zz['interestDateGap'], Input_data['term'],
                            np.array(Input_data['trchSize']).cumsum()[-2],
                            Input_data['feeAmt'][4]
                            )

    bb = f_main_profitAccountACNM(yy['DisposableCashFlow'], aa['TaxFeesSum'], zz['interestDateGap'], yy['isAc'],
                                  yy['PrincipalAccount'], xx['DefaultPrincipal'], xx['InitPriBalance'],
                                  yy['IncomeAccount'],
                                  Input_data['trchFixedPmtAmt'],
                                  Input_data['trchSize'], Input_data['trchInterestRate']
                                  , rounding_error
                                  )

    cc = f_main_TaxFeesDE(Input_data['isDefault'], yy['DisposableCashFlow'], aa['TaxFeesSum'],
                          bb['priorA1ResPrincipal'], zz['interestDateGap'],
                          yy['IncomeAccount'], yy['PrincipalAccount'], aa['ServiceFees'], xx['InitPriBalance'],
                          Input_data['term'], Input_data['trchSize'],
                          Input_data['trchInterestRate'])

    dd = f_main_AccountDE(cc['DisposableCashFlow_DE'], cc['TaxFeesSum_DE'], zz['interestDateGap'],
                          Input_data['trchFixedPmtAmt'],
                          Input_data['trchSize'],
                          Input_data['trchInterestRate']
                          )

    return aa, bb, cc, dd, xx, yy, zz






if __name__ == '__main__':
    # 计算结果保存出来查看
    import pandas as pd
    import json
    from AxisModelPackage.AxisModelComputer import *
    import numpy as np
    from AxisModelPackage.AxisModelOutput import *
    from AxisModelPackage.AxisModelInput import cash_input

    with open(r'C:\Users\liyin\Desktop\Axis\Axis_1014\request_1013.json', 'r+', encoding='utf-8') as f:
        f_json = json.load(f)
    Input_data = cash_input(f_json)
    priod = len(Input_data['pdByTerm'])
    rec_term = 1
    term = np.arange(priod) + 1
    Input_data['term'] = term
    Input_data['rec_term'] = rec_term
    rounding_error = 0.001
    isDefault = np.zeros(priod)
    Input_data['isDefault'] = isDefault

    aa, bb, cc, dd, xx, yy, zz = f_AxisModelMidComp(Input_data, rounding_error)



    res = f_main_output(aa, bb, cc, dd, xx, yy, zz,Input_data)
    res_format = f_outputFormat4page12(res, xx, Input_data)

    r"""
    pd.DataFrame(xx, index=compDate[1:]).T.to_csv(r'C:\Users\liyin\Desktop\Axis\Axis_1012_checkresult\收入项.csv')
    pd.DataFrame(zz, index=compDate[1:]).T.to_csv(r'C:\Users\liyin\Desktop\Axis\Axis_1012_checkresult\支出项.csv')
    pd.DataFrame(yy, index=compDate[1:]).T.to_csv(r'C:\Users\liyin\Desktop\Axis\Axis_1012_checkresult\收益帐前半部分.csv')
    pd.DataFrame(aa, index=compDate[1:]).T.to_csv(r'C:\Users\liyin\Desktop\Axis\Axis_1012_checkresult\收益账户税费.csv')
    pd.DataFrame(bb, index=compDate[1:]).T.to_csv(r'C:\Users\liyin\Desktop\Axis\Axis_1012_checkresult\收益账户利息本金计算.csv')
    pd.DataFrame(cc, index=compDate[1:]).T.to_csv(r'C:\Users\liyin\Desktop\Axis\Axis_1012_checkresult\违约事件税费相关.csv')
    pd.DataFrame(dd, index=compDate[1:]).T.to_csv(r'C:\Users\liyin\Desktop\Axis\Axis_1012_checkresult\违约事件利息本金.csv')
    
    """

    ######
    '''
    xx_keys = ['InitPriBalance', 'DefaultPrincipal', 'PrepaymentAmt', 'ReAmt', 'InterestPaid', 'PriPaid',
               'ReinvestedEarnings', 'FinalPriBalance']
    yy_keys = ['isAc', 'IncomeAccount', 'PrincipalAccount', 'DisposableCashFlow']
    aa_keys = ['InterestTax', 'storageFees', 'AuditFees', 'GradeFees', 'ServiceFees', 'TaxFeesSum']
    bb_keys = ['priorA1Interest', 'priorA2Interest', 'priorBInterest', 'priorCInterest', 'priorInterestSum',
               'PrincipalAccountNM', 'PrincipalAccountAC', 'RePrincipalAccount', 'SurplusFees', 'subprimeLimitIncome',
               'ResCashFlow', 'Spread', 'ExcessInterest', 'PrincipalCompensation', 'IncomepaidBackPrincipal',
               'DisposableBalance', 'priorA1RepPrincipal', 'priorA2RepPrincipal', 'priorA1LagRepPrincipal',
               'priorBRepPrincipal', 'priorCRepPrincipal', 'subprimeRepPrincipal', 'priorA1ResPrincipal',
               'priorA2ResPrincipal', 'priorBResPrincipal', 'priorCResPrincipal', 'subprimeResPrincipal',
               'subprimeResIncome']
    cc_keys = ['isDefault', 'IncomeAccount_DE', 'PrincipalAccount_DE', 'DisposableCashFlow_DE', 'InterestTax_DE',
               'storageFees_DE', 'AuditFees_DE', 'GradeFees_DE', 'SurplusFees_DE', 'TaxFeesSum_DE']
    dd_keys = ['priorA1Interest_DE', 'priorA1RepPrincipal_DE', 'priorA2Interest_DE', 'priorA2RepPrincipal_DE',
               'priorBInterest_DE', 'priorBRepPrincipal_DE', 'priorCInterest_DE', 'priorCRepPrincipal_DE',
               'subprimeRepPrincipal_DE', 'subprimeRepIncome_DE', 'priorA1ResPrincipal_DE', 'priorA2ResPrincipal_DE',
               'priorBResPrincipal_DE', 'priorCResPrincipal_DE', 'subprimeResPrincipal_DE']
    '''
